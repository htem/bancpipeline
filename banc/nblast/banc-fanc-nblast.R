#' banc-fanc-nblast — Pairwise BANC ↔ FANC (v1.116) NBLAST.
#'
#' VNC + ascending + descending + sensory + afferent + efferent neurons.
#' Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, FANC BANC-space skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.fanc.save.path>/results/<ver>/<root_id>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

########################
### BANC-FANC NBLAST ###
########################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))
version <- banc.nblast.version

# Direct us to the BANC dataset
banc.nblast.fanc.swc.save.path.version <- file.path(banc.nblast.fanc.swc.save.path,version)
dir.results <- file.path(banc.nblast.fanc.save.path,"results", version)
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# Read IDs
banc.meta <- banctable_query() %>%
  dplyr::filter(!is.na(root_id),
                region == "ventral_nerve_cord" | grepl("ascending|descending", super_class),
                grepl("sensory|afferent|efferent", super_class),
                !grepl("glia", super_class),
                !grepl("glia", cell_class))
banc.root.ids <- unique(banc.meta$root_id)

# Get fanc meta data
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(fc.meta$root_id)

# NBLAST files todo
message("##### Calculating NBLASTs we need to do #####")
dir.create(banc.nblast.fanc.save.path, showWarnings = FALSE)
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
  message("All FANC NBLASTs up to date. Nothing to do.")
} else {

# Read root nodes
y.cut <- 325000
banc.soma.positions <- arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather"))
has_pos <- !is.na(banc.soma.positions$root_position_nm) & banc.soma.positions$root_position_nm != ""
banc.soma.positions[has_pos, c("X", "Y", "Z")] <- nat::xyzmatrix(banc.soma.positions$root_position_nm[has_pos])
banc.soma.positions <- banc_decapitate(banc.soma.positions, invert = FALSE)

# Preload fanc data
message("##### Converting fanc skeletons to vector clouds #####")
fc.swc.full.files <- unique(unlist(lapply(fc.ids, function(id) file.path(banc.nblast.fanc.swc.save.path.version,paste0(id,".swc")))))
fc.swc.full.files <- fc.swc.full.files[file.exists(fc.swc.full.files)]
# fc.skels <- nat::read.neurons(fc.swc.full.files)
# fc.skels.m <- banc_mirror(fc.skels)
# names(fc.skels.m) <- paste0("m",names(fc.skels.m))
# fc.skels <- c(fc.skels,fc.skels.m)
# fc.dps <- nat::dotprops(fc.skels/1000,
#                         k = 20,
#                         topo = FALSE,
#                         resample = 3,
#                         UseAlpha = 1,
#                         OmitFailures = TRUE,
#                         .parallel = FALSE)
# save(fc.dps,file=file.path(banc.nblast.fanc.save.path,"banc_microns_space_fanc_dps.rda"))
# rm('fc.skels')
load(file.path(banc.nblast.fanc.save.path,"banc_microns_space_fanc_dps.rda"))

# Batch for parallel processing
message("##### Setting up parallel processing #####")
multiplier <- 10000
upper <- ifelse((numCores*multiplier)<length(nblast.todo),numCores*multiplier,length(nblast.todo))
batches <- split(sample(nblast.todo), round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
batches <- sample(batches)

# Register cores
cl <- setup_parallel()

# NBLAST!
message("##### Running BANC-MANC NBLAST #####")
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
                                 svid <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_id)]

                                 # Check if already done — before expensive skeleton read
                                 file <- file.path(dir.results, paste0("supervoxel_id_",svid,"_root_id_",neuron.id,".csv"))
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
                                 banc.query <- tryCatch({banc_decapitate(banc.query, invert = FALSE)},
                                                        error = function(e) NULL)
                                 if(is.null(banc.query)|!length(banc.query)){
                                   message("not in vnc, skipping:", neuron.id)
                                   next
                                 }
                                 banc.query.pruned <- tryCatch(
                                   prune_in_volume(banc.query,
                                                   surf = banc_vnc_neuropil.surf,
                                                   invert = TRUE),
                                   error = function(e) {
                                     warning("prune_in_volume failed for ", neuron.id, ": ", e$message)
                                     NULL
                                   }
                                 )
                                 if(is.null(banc.query.pruned)|!length(banc.query.pruned)){
                                   message("not in vnc, skipping:", neuron.id)
                                   next
                                 }
                                 message("points in vnc, NBLASTing:", neuron.id)
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
                                 if(is.null(query.soma)|!nrow(query.soma)){
                                   query.soma <- nat::xyzmatrix(banc.query)[nat::endpoints(banc.query[[1]]),]
                                }
                                 
                                 # Forward NBLAST: 1 query vs all targets
                                 nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = fc.dps,
                                                                    UseAlpha = FALSE, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))

                                 # Reverse NBLAST on top 100 only
                                 top_k <- min(100, length(nb.fwd))
                                 top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                 nb.rev <- drop(nat.nblast::nblast(query = fc.dps[top_names], target = banc.query.dps,
                                                                    UseAlpha = FALSE, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))
                                 nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                 # Make data frame (top 100 only)
                                 nb.df <- data.frame(nb = nb) %>%
                                   dplyr::arrange(dplyr::desc(nb))
                                 nb.df$query <- names(banc.query.dps)
                                 nb.df$supervoxel_id <- svid
                                 nb.df$fanc_id <- rownames(nb.df)
                                 nb.df$fanc_match <- gsub("m","",rownames(nb.df))
                                 nb.df <- dplyr::left_join(nb.df,
                                                           fc.meta %>%
                                                             dplyr::select(fanc_root_id = root_id,
                                                                           cell_id,
                                                                           fanc_supervoxel_id = supervoxel_id,
                                                                           fanc_root_position = root_position,
                                                                           cell_type),
                                                           by = c("fanc_match"="fanc_root_id"))
                                 
                                 # Save NBLAST
                                 readr::write_csv(nb.df, file=file)
                                 message("Saved NBLAST result: ", file)
                                 success <- success+1
                                 
                               }
                               
                               # Announce
                               message("completed: ", length(neuron.ids), " neurons NBLASTed against FANC neurons with nearby root")

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

