#' banc-l2 — Acquire and reroot per-neuron L2 skeletons via `pcg_skel`.
#'
#' Skips neurons already on disk. Honours `banc.test.ids` if set globally.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, CAVE L2 cache
#'
#' @section Writes:
#'   - `<banc.l2swc.save.path>/<root_id>.swc`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.

#######################################
### Acquire and reroot L2 skeletons ###
#######################################
source("banc/banc-startup.R")
redo <- FALSE

# Direct us to the BANC dataset
bancr::choose_banc()

# Read IDs — use current-version root IDs (canonical for the active BANC pipeline)
banc.ids <- banctable_query_cached() %>%
  banc_filter_neurons() %>%
  dplyr::filter(!is.na(root_id))
banc.root.ids <- na.omit(unique(banc.ids$root_id))
message(sprintf("Target neurons (root_%s): %d", banc.version, length(banc.root.ids)))

# Read written skels
dir.create(banc.l2swc.save.path, showWarnings = FALSE, recursive = TRUE)
banc.done.list <- list.files(banc.l2swc.save.path, pattern = "\\.swc", full.names = TRUE)
banc.l2.done <- gsub("\\.swc","",basename(banc.done.list))
banc.l2.todo <- setdiff(banc.root.ids, banc.l2.done)

# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.l2.todo <- intersect(banc.l2.todo, banc.test.ids)
message(length(banc.l2.todo), " L2 skeletons to fetch, fetched: ", length(banc.l2.done))

if (length(banc.l2.todo) == 0) {
  message("All L2 skeletons up to date. Nothing to do.")
} else {

# Batch for parallel processing
multiplier <- 10
banc.l2.todo <- sample(banc.l2.todo)
upper <- ifelse((numCores*multiplier)<length(banc.l2.todo),numCores*multiplier,length(banc.l2.todo))
batches <- split(banc.l2.todo, round(seq(from = 1, to = upper, length.out = length(banc.l2.todo))))

# Get L2 skeletons in for loop
by.query <- foreach::foreach(batch = seq_along(batches),
                              .combine = 'c',
                              .errorhandling='pass') %do% {
  
  # Get neuron IDs                              
  neuron.ids <- batches[[batch]]
  
  # Get L2 skeletons
  banc.l2.skels <- banc_read_l2skel(neuron.ids, OmitFailures = TRUE)
  failed.ids <- setdiff(neuron.ids, names(banc.l2.skels))
  if (length(failed.ids)) {
    message(sprintf("  CAVE fetch failed for %d/%d neurons: %s",
                    length(failed.ids), length(neuron.ids),
                    paste(head(failed.ids, 5), collapse = ", ")))
  }
  banc.l2.skels[,"id"] <- names(banc.l2.skels)
  neuron.ids <- names(banc.l2.skels)
  if(!length(banc.l2.skels)){
    return(NULL)
  }
  
  # Re-root the neuron
  for(id in neuron.ids){
    df <- subset(banc.ids,banc.ids$root_id==id)
    if(!is.na(df$nucleus_position_nm[1]) & df$nucleus_id[1]!="0"){
      soma <- nat::xyzmatrix(df$nucleus_position_nm)[1,]
      banc.l2.skels[id][[1]] <- nat::reroot(x = banc.l2.skels[id][[1]], point = c(soma))
      banc.l2.skels[id][[1]]$tags$soma <- nat::rootpoints(banc.l2.skels[id][[1]] )
    }else{ # As best we can
      leaves <- nat::endpoints(banc.l2.skels[id][[1]])
      npoints1 <- nat::xyzmatrix(banc.l2.skels[id][[1]])[leaves,]
      npoints <- npoints1
      if(nrow(npoints1)&ncol(npoints1)==3){npoints=npoints1}
      pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neuropil.surf)
      npoints2 <- data.frame(npoints[!pin,])
      if(nrow(npoints2)&ncol(npoints2)==3){npoints=npoints2}
      pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neck_connective.surf)
      npoints3 <- data.frame(npoints[!pin,])
      if(nrow(npoints3)&ncol(npoints3)==3){npoints=npoints3}
      npoints$nucleus_id <- 0
      npoints$root_id <- id
      nearest <- nabor::knn(query = nat::xyzmatrix(npoints), data = rbind(xyzmatrix(banc_neuropil.surf),xyzmatrix(banc_neck_connective.surf)), k = 1)
      soma <-nat::xyzmatrix(npoints)[which.max(nearest$nn.dists),]
      banc.l2.skels[id][[1]] <- nat::reroot(x = banc.l2.skels[id][[1]], point = c(soma))
    }
  }
  
  # Write neuron
  success <- nat::write.neurons(banc.l2.skels, file = neuron.ids, dir = banc.l2swc.save.path, format='swc', Force = TRUE)
  
  # Get summary data and save metrics
  banc.skels.stats <- summary(banc.l2.skels)
  colnames(banc.skels.stats) <- paste0("l2_",snakecase::to_snake_case(colnames(banc.skels.stats)))
  banc.skels.stats <- round(banc.skels.stats,2)
  banc.skels.stats$root_id <- rownames(banc.skels.stats)
  for(i in 1:nrow(banc.skels.stats)){
    metrics.file <- file.path(banc.metrics.save.path,paste0(banc.skels.stats[i,"root_id"],".csv"))
    if(file.exists(metrics.file)){
      id.metrics <- readr::read_csv(metrics.file, col_types = banc.col.types, show_col_types = FALSE)
      id.metrics <- id.metrics[,c("root_id",setdiff(colnames(id.metrics),colnames(banc.skels.stats)))]
      id.metrics <- dplyr::left_join(banc.skels.stats[i,], id.metrics, by = "root_id")
    }else{
      id.metrics <- banc.skels.stats[i,]
    }
    readr::write_csv(x = id.metrics, file = metrics.file)
  }
  
  # announce
  message("completed: ", length(banc.l2.skels), " L2 grabs saved as swc files")
  
  # return
  NULL
}

# Were there errors?
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}
} # end if banc.l2.todo > 0

# List new files
banc.list <- list.files(banc.l2swc.save.path, pattern = "\\.swc")
banc.l2 <- gsub("\\.swc","",basename(banc.list))
banc.l2.new <- setdiff(banc.l2, banc.l2.done)

# Announce
message("##### BANCpipeline: banc l2 skeletons updated #####")
message(sprintf("##### new L2 SWC files: %s", length(banc.l2.new)))



