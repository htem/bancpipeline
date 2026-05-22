#' banc-malecns-nblast — Pairwise BANC ↔ maleCNS (v0.9) NBLAST.
#'
#' Iterates per super_class (motor + visceral_circulatory first). Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, `malecns_09_meta.csv`, maleCNS BANC-space skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.malecns.save.path>/results/<ver>/<root_id>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

###########################
### BANC-maleCNS NBLAST ###
###########################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))
version <- banc.nblast.malecns.version

banc.meta.all <- banctable_query() %>%
  dplyr::filter(!is.na(root_id))
super_classes <- sample(unique(banc.meta.all$super_class))

# Read maleCNS meta once (not per super_class)
mcns.meta.all <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "malecns_09_meta.csv"),
                                                   col_types = hemibrainr:::sql_col_types))

# Iterate
super_classes = c("motor","visceral_circulatory",super_classes)
for(sp in super_classes){

  # Direct us to the BANC dataset
  message(sp)
  banc.nblast.malecns.swc.save.path.version <- file.path(banc.nblast.malecns.swc.save.path, version)
  dir.results <- file.path(banc.nblast.malecns.save.path, "results", version)
  dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

  # Read IDs
  # maleCNS covers the whole CNS (brain + VNC), so we NBLAST against all BANC neurons
  if(!is.na(sp)){
    banc.meta <- banc.meta.all %>%
      dplyr::filter(super_class == sp)
    mcns.meta <- mcns.meta.all %>%
      dplyr::filter(super_class == sp)
  }else{
    banc.meta <- banc.meta.all
    mcns.meta <- mcns.meta.all
  }

  # Skip if missing
  if(!nrow(mcns.meta)){
    next
  }
  mcns.ids <- unique(mcns.meta$malecns_09_id)
  banc.root.ids <- sample(unique(banc.meta$root_id))

  # NBLAST files todo
  message("##### Calculating NBLASTs we need to do #####")
  dir.create(banc.nblast.malecns.save.path, showWarnings = FALSE)
  nblast.files <- if(redo) list_files_age_sorted(dir.results, pattern = "\\.csv") else list.files(dir.results, pattern = "\\.csv")

  # Fix files with missing supervoxel_id in filename
  na_svid_files <- nblast.files[grepl("supervoxel_id_NA_", basename(nblast.files))]
  if (length(na_svid_files)) {
    if (!grepl(.Platform$file.sep, na_svid_files[1])) na_svid_files <- file.path(dir.results, na_svid_files)
    message(sprintf("  Fixing %d files with supervoxel_id_NA in filename...", length(na_svid_files)))
    na_root_ids <- gsub(".*_root_id_|\\.csv", "", basename(na_svid_files))
    svid_lookup <- banc.meta.all$supervoxel_id[match(na_root_ids, banc.meta.all$root_id)]
    can_fix <- !is.na(svid_lookup)
    if (any(can_fix)) {
      new_names <- file.path(dirname(na_svid_files[can_fix]),
                             paste0("supervoxel_id_", svid_lookup[can_fix],
                                    "_root_id_", na_root_ids[can_fix], ".csv"))
      renamed <- mapply(file.rename, na_svid_files[can_fix], new_names)
      message(sprintf("  Renamed %d/%d files", sum(renamed), sum(can_fix)))
    }
    nblast.files <- if(redo) list_files_age_sorted(dir.results, pattern = "\\.csv") else list.files(dir.results, pattern = "\\.csv")
  }

  nblast.done <- gsub(".*_root_id_|\\.csv", "", basename(nblast.files))
  if(redo){
    nblast.todo <- banc.root.ids
    undone <- setdiff(banc.root.ids, nblast.done)
    nblast.todo <- union(undone, nblast.done)
  }else{
    nblast.todo <- setdiff(banc.root.ids, nblast.done)
  }
  # Honor test mode if active
  if (exists("banc.test.ids", envir = .GlobalEnv))
    nblast.todo <- intersect(nblast.todo, banc.test.ids)
  message(length(nblast.todo), " neurons to NBLAST, done: ", length(nblast.done))
  if (length(nblast.todo) == 0) {
    message("All maleCNS NBLASTs up to date for super_class: ", sp, ". Skipping.")
    next
  }

  # Read root nodes
  banc.soma.positions <- arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather"))

  # Preload maleCNS data
  message("##### Converting maleCNS skeletons to vector clouds #####")
  mcns.swc.full.files <- file.path(banc.nblast.malecns.swc.save.path.version, paste0(mcns.ids, ".swc"))
  mcns.swc.full.files <- mcns.swc.full.files[file.exists(mcns.swc.full.files)]
  mcns.skels <- nat::read.neurons(mcns.swc.full.files)
  mcns.dps <- nat::dotprops(mcns.skels/1000,
                            k = 20,
                            topo = FALSE,
                            resample = 3,
                            UseAlpha = 1,
                            OmitFailures = TRUE,
                            .parallel = FALSE)
  save(mcns.dps, file = file.path(banc.nblast.malecns.save.path, "banc_microns_space_malecns_dps.rda"))
  rm('mcns.skels')
  # load(file = file.path(banc.nblast.malecns.save.path, "banc_microns_space_malecns_dps.rda"))
  
  # Batch for parallel processing
  message("##### Setting up parallel processing #####")
  multiplier <- 1000
  upper <- ifelse((numCores * multiplier) < length(nblast.todo), numCores * multiplier, length(nblast.todo))
  batches <- split(sample(nblast.todo), round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
  batches <- sample(batches)

  # Register cores
  cl <- setup_parallel()

  # NBLAST!
  message("##### Running BANC-maleCNS NBLAST #####")
  by.query <- foreach::foreach(batch = seq_along(batches),
                               .combine = 'c',
                               .init = list(),
                               .errorhandling = 'pass') %dopar% {

                                 # Select IDs
                                 neuron.ids <- batches[[batch]]

                                 # Conduct near-soma search
                                 success <- 0
                                 for(neuron.id in neuron.ids){

                                   # Look up supervoxel_id early for file existence check
                                   supervoxel.id <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_id)]
                                   query.soma <- subset(banc.soma.positions,
                                                        banc.soma.positions$root_id == neuron.id)
                                   use.alpha <- FALSE
                                   if(!is.null(query.soma) & nrow(query.soma)){
                                     nucleus.id <- query.soma$nucleus_id[1]
                                     if("supervoxel_id" %in% names(query.soma) && !is.na(query.soma$supervoxel_id[1])){
                                       supervoxel.id <- query.soma$supervoxel_id[1]
                                     }
                                     if(is.na(nucleus.id) | nucleus.id != "0"){
                                       use.alpha = FALSE
                                     }else{
                                       use.alpha = TRUE
                                     }
                                   }else{
                                     nucleus.id <- 0
                                   }

                                   # Check if already done — before expensive skeleton read
                                   file <- file.path(dir.results, paste0("supervoxel_id_", supervoxel.id, "_root_id_", neuron.id, ".csv"))
                                   if(file.exists(file) & !redo){
                                     next
                                   }

                                   # Read L2 skeleton directly by path (not directory scan)
                                   swc_file <- file.path(banc.l2swc.save.path, paste0(neuron.id, ".swc"))
                                   if(!file.exists(swc_file)){
                                     next
                                   }
                                   banc.query <- tryCatch(
                                     nat::read.neurons(swc_file, neuronnames = function(f) sub("\\.swc$", "", basename(f))),
                                     error = function(e) NULL
                                   )
                                   if(is.null(banc.query) | !length(banc.query)){
                                     next
                                   }

                                   banc.query.dps <- nat::dotprops(banc.query/1000,
                                                                   k = 20,
                                                                   topo = FALSE,
                                                                   resample = ifelse(summary(banc.query)$cable.length < 30000, 1, 2),
                                                                   UseAlpha = 1,
                                                                   OmitFailures = TRUE)
                                   if(is.null(banc.query.dps) | !length(banc.query.dps)){
                                     next
                                   }

                                   # Forward NBLAST: 1 query vs all targets
                                   nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = mcns.dps,
                                                                      UseAlpha = use.alpha, normalised = TRUE,
                                                                      smat = nat.nblast::smat_alpha.fcwb))

                                   # Reverse NBLAST on top 100 only
                                   top_k <- min(100, length(nb.fwd))
                                   top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                   nb.rev <- drop(nat.nblast::nblast(query = mcns.dps[top_names], target = banc.query.dps,
                                                                      UseAlpha = use.alpha, normalised = TRUE,
                                                                      smat = nat.nblast::smat_alpha.fcwb))
                                   nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                   # Make data frame (top 100 only)
                                   nb.df <- data.frame(nb = nb) %>%
                                     dplyr::arrange(dplyr::desc(nb))
                                   nb.df$bodyid <- rownames(nb.df)
                                   nb.df <- dplyr::left_join(nb.df, mcns.meta[, c("malecns_09_id", "cell_type")],
                                                             by = c("bodyid" = "malecns_09_id"))

                                   # Save NBLAST
                                   readr::write_csv(nb.df, file = file)
                                   success <- success + 1

                                 }

                                 # Return success ratio
                                 success / length(neuron.ids)
                               }

  # Stop cores
  stop_parallel(cl)

  # Were there errors?
  message("##### Displaying any errors from foreach loop #####")
  for(i in 1:length(by.query)){
    if(!is.null(by.query[[i]])){
      if(!is.numeric(by.query[[i]])){
        message(by.query[[i]])
      }
    }
  }

}

###########################
### Clean low-score files ##
###########################

# # Delete NBLAST result CSVs where no score exceeds 0
# # These are uninformative matches that waste disk space and slow compilation.
# message("##### Cleaning low-score NBLAST files #####")
# dir.results.all <- file.path(banc.nblast.malecns.save.path, "results", version)
# all_result_files <- list.files(dir.results.all, pattern = "\\.csv$", full.names = TRUE)
# message(sprintf("  Checking %d result files...", length(all_result_files)))
# 
# n_deleted <- 0L
# for (f in all_result_files) {
#   nb_df <- tryCatch(
#     readr::read_csv(f, col_types = readr::cols(.default = "c", nb = "d"),
#                     show_col_types = FALSE),
#     error = function(e) NULL
#   )
#   if (is.null(nb_df) || nrow(nb_df) == 0 || max(nb_df$nb, na.rm = TRUE) < 0) {
#     file.remove(f)
#     n_deleted <- n_deleted + 1L
#   }
# }
# message(sprintf("  Deleted %d/%d files with no score >= 0.3", n_deleted, length(all_result_files)))
