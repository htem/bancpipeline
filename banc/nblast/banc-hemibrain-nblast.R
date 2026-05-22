#' banc-hemibrain-nblast — Pairwise BANC ↔ hemibrain (v1.2.1) NBLAST (native + mirrored).
#'
#' Non-VNC BANC neurons only. Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, hemibrain BANC-space skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.hemibrain.save.path>/results_with_mirrored/<ver>/<root_id>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

#############################
### BANC-Hemibrain NBLAST ###
#############################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))
version <- banc.nblast.version

# Direct us to the BANC dataset
banc.nblast.hemibrain.swc.save.path.version <- file.path(banc.nblast.hemibrain.swc.save.path,version)
dir.results <- file.path(banc.nblast.hemibrain.save.path,"results_with_mirrored", version)
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# Read IDs
banc.meta <- banctable_query() %>%
  dplyr::filter(!is.na(root_id),
                !region %in% c("ventral_nerve_cord"))  %>%
  dplyr::filter(!grepl("glia|trachea", as.character(super_class)),
                !grepl("glia|trachea", as.character(cell_class)))
banc.root.ids <- unique(banc.meta$root_id)

# Get hemibrain meta data
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
hb.ids <- unique(hb.meta$bodyid)

# Read BANC root nodes
banc.soma.positions <- arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather"))
has_pos <- !is.na(banc.soma.positions$root_position_nm) & banc.soma.positions$root_position_nm != ""
banc.soma.positions[has_pos, c("X", "Y", "Z")] <- nat::xyzmatrix(banc.soma.positions$root_position_nm[has_pos])
banc.soma.positions <- banc_decapitate(banc.soma.positions, invert = TRUE)

# NBLAST files todo
message("##### Calculating NBLASTs we need to do #####")
dir.create(banc.nblast.hemibrain.save.path, showWarnings = FALSE)
nblast.files <- if(redo) list_files_age_sorted(dir.results, pattern = "\\.csv") else list.files(dir.results, pattern = "\\.csv")

# Fix files with missing supervoxel_id in filename
na_svid_files <- nblast.files[grepl("supervoxel_id_NA_", basename(nblast.files))]
if (length(na_svid_files)) {
  if (!grepl(.Platform$file.sep, na_svid_files[1])) na_svid_files <- file.path(dir.results, na_svid_files)
  message(sprintf("  Fixing %d files with supervoxel_id_NA in filename...", length(na_svid_files)))
  na_root_ids <- gsub(".*_root_id_|\\.csv", "", basename(na_svid_files))
  svid_lookup <- banc.meta$supervoxel_id[match(na_root_ids, banc.meta$root_id)]
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
  message("All hemibrain NBLASTs up to date. Nothing to do.")
} else {

# Preload hemibrain data
# message("##### Converting hemibrain skeletons to vector clouds #####")
# hb.swc.full.files <- unique(unlist(lapply(hb.ids, function(id) file.path(banc.nblast.hemibrain.swc.save.path.version,paste0(id,".swc")))))
# hb.swc.full.files <- hb.swc.full.files[file.exists(hb.swc.full.files)]
# hb.skels <- nat::read.neurons(hb.swc.full.files)
# hb.skels.m <- banc_mirror(hb.skels)
# names(hb.skels.m) <- paste0("m",names(hb.skels.m))
# hb.skels <- c(hb.skels,hb.skels.m)
# hb.dps <- nat::dotprops(hb.skels/1000,
#                         k = 20,
#                         topo = FALSE,
#                         resample = 3,
#                         UseAlpha = 1,
#                         OmitFailures = TRUE,
#                         .parallel = FALSE)
# save(hb.dps,file=file.path(banc.nblast.hemibrain.save.path,"banc_microns_space_hemibrain_dps.rda"))
# rm('hb.skels')

# Load pre-made neuron object
load(file.path(banc.nblast.hemibrain.save.path,"banc_microns_space_hemibrain_dps.rda"))

# Batch for parallel processing
message("##### Setting up parallel processing #####")
multiplier <- 1000
upper <- ifelse((numCores*multiplier)<length(nblast.todo),numCores*multiplier,length(nblast.todo))
batches <- split(sample(nblast.todo), round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
batches <- sample(batches)

# Register cores
cl <- setup_parallel()

# NBLAST!
message("##### Running BANC-hemibrain NBLAST #####")
by.query <- foreach::foreach(batch = seq_along(batches),
                             .combine = 'c',
                             .init = list(),
                             .errorhandling = 'pass') %dopar% {
                               
                               # Select IDs
                               neuron.ids <- batches[[batch]]
                               message("working on:", paste0(neuron.ids,collapse=", "))
                               
                               # Conduct near-soma search
                               success <- 0
                               for(neuron.id in neuron.ids){
                                 
                                 # Look up supervoxel_id early for file existence check
                                 supervoxel.id <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_id)]

                                 # Check if already done — before expensive skeleton read
                                 file <- file.path(dir.results, paste0("supervoxel_id_",supervoxel.id,"_root_id_",neuron.id,".csv"))
                                 if(file.exists(file)&!redo){
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
                                 banc.query <- tryCatch({banc_decapitate(banc.query, invert = TRUE)},
                                                        error = function(e) NULL)
                                 if(is.null(banc.query)|!length(banc.query)){
                                   message("not in hemibrain, skipping:", neuron.id)
                                   next
                                 }
                                 banc.query.pruned <- tryCatch(
                                   prune_in_volume(banc.query,
                                                   surf = banc_hemibrain.surf,
                                                   invert = TRUE),
                                   error = function(e) {
                                     warning("prune_in_volume failed for ", neuron.id, ": ", e$message)
                                     NULL
                                   }
                                 )
                                 if(is.null(banc.query.pruned)|!length(banc.query.pruned)){
                                   message("not in right hemibrain, attempting mirror ..")
                                   banc.query.pruned <- tryCatch(
                                     prune_in_volume(banc_mirror(banc.query),
                                                     surf = banc_hemibrain.surf,
                                                     invert = TRUE),
                                     error = function(e) {
                                       warning("prune_in_volume (mirror) failed for ", neuron.id, ": ", e$message)
                                       NULL
                                     }
                                   )
                                   if(is.null(banc.query.pruned)|!length(banc.query.pruned)){
                                     message("not in hemibrain, skipping:", neuron.id)
                                     next
                                   }
                                   message("in left mirrored-hemibrain")
                                   banc.query.pruned <- banc_mirror(banc.query.pruned)
                                 }
                                 message("points in hemibrain, NBLASTing:", neuron.id)
                                 banc.query.dps <- nat::dotprops(banc.query.pruned/1000,
                                                                 k = 20,
                                                                 topo = FALSE,
                                                                 resample = ifelse(summary(banc.query.pruned)$cable.length<30000,1,2),
                                                                 UseAlpha = 1,
                                                                 OmitFailures = TRUE)
                                 if(is.null(banc.query.dps)|!length(banc.query.dps)){
                                   next
                                 }
                                 
                                 # Conduct near-soma search
                                 query.soma <- subset(banc.soma.positions, 
                                                      banc.soma.positions$root_id == neuron.id)
                                 nucleus.id <- query.soma$nucleus_id[1]
                                 if(is.na(nucleus.id)|nucleus.id!="0"){
                                   use.alpha = FALSE
                                 }else{
                                   use.alpha = TRUE
                                 }
                                 if(is.null(query.soma)|!nrow(query.soma)){
                                   query.soma <- nat::xyzmatrix(banc.query)[nat::endpoints(banc.query[[1]]),]
                                   nucleus.id <- 0
                                 }
                                 
                                 # Forward NBLAST: 1 query vs all targets
                                 nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = hb.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))

                                 # Reverse NBLAST on top 100 only
                                 top_k <- min(100, length(nb.fwd))
                                 top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                 nb.rev <- drop(nat.nblast::nblast(query = hb.dps[top_names], target = banc.query.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))
                                 nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                 # Make data frame (top 100 only)
                                 nb.df <- data.frame(nb = nb) %>%
                                   dplyr::arrange(dplyr::desc(nb))
                                 nb.df$bodyid <- rownames(nb.df)
                                 nb.df <- dplyr::left_join(nb.df, hb.meta[,c("bodyid","cell_type")], by = "bodyid")
                                 
                                 # Save NBLAST
                                 readr::write_csv(nb.df, file=file)
                                 message("Saved NBLAST result: ", file)
                                 success <- success+1
                                 
                               }
                               
                               # Announce
                               message("completed: ", length(neuron.ids), " neurons NBLASTed against hemibrain neurons with nearby root")

                               # Return success ratio
                               success/length(neuron.ids)
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
} # end if nblast.todo > 0

