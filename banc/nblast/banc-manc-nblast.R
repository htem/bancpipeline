#' banc-manc-nblast — Pairwise BANC ↔ MANC (v1.2.1) NBLAST.
#'
#' Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, `manc_meta.csv`, MANC BANC-space skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.manc.save.path>/results/<ver>/<root_id>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

########################
### BANC-MANC NBLAST ###
########################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))
version <- banc.nblast.version

# Direct us to the BANC dataset
banc.nblast.manc.swc.save.path.version <- file.path(banc.nblast.manc.swc.save.path,version)
dir.results <- file.path(banc.nblast.manc.save.path,"results", version)
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# Read IDs
banc.meta <- banctable_query() %>%
  dplyr::filter(!is.na(root_id),
                !grepl("glia",super_class),
                !grepl("glia",cell_class))
banc.root.ids <- unique(banc.meta$root_id)

# Get manc meta data
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)

# NBLAST files todo
message("##### Calculating NBLASTs we need to do #####")
dir.create(banc.nblast.manc.save.path, showWarnings = FALSE)
nblast.files <- if(redo) list_files_age_sorted(dir.results, pattern = "\\.csv") else list.files(dir.results, pattern = "\\.csv")

# Fix files with missing supervoxel_id in filename
na_svid_files <- nblast.files[grepl("supervoxel_id_NA_", basename(nblast.files))]
if (length(na_svid_files)) {
  # list.files returns basenames only; need full paths for rename
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
  # Refresh file list after renames
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
  message("All MANC NBLASTs up to date. Nothing to do.")
} else {

# Read root nodes
y.cut <- 325000
banc.soma.positions <- arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather"))
has_pos <- !is.na(banc.soma.positions$root_position_nm) & banc.soma.positions$root_position_nm != ""
banc.soma.positions[has_pos, c("X", "Y", "Z")] <- nat::xyzmatrix(banc.soma.positions$root_position_nm[has_pos])
banc.soma.positions <- banc_decapitate(banc.soma.positions, invert = FALSE)

# # manc soma positions in banc space
# mc.positions <- mc.meta %>%
#   dplyr::mutate(root_pos = dplyr::case_when(
#     soma_location!="" ~ soma_location,
#     tosoma_location!="" ~ tosoma_location,
#     root_location!="" ~ root_location,
#     TRUE ~ ""
#   )) %>%
#   dplyr::mutate(root_pos=gsub("list\\(|\\).*","",root_pos)) %>%
#   dplyr::distinct(bodyid, root_pos) %>%
#   dplyr::filter(root_pos!="")
# mc.positions[,c("X",'Y',"Z")] <- nat::xyzmatrix(mc.positions$root_pos)
# mc.meta.no.soma.ids <- setdiff(mc.ids, mc.meta$bodyid[!mc.meta$soma])
  
# # Transforms into BANC space
# message("##### Transforming MANC roots to BANC space #####")
# mc.positions.xyz <- nat::xyzmatrix(mc.positions$root_pos)
# mc.positions.xyz <- mc.positions.xyz*rep(8/1000, 3)
# mc.positions.xyz.jrcvnc2018f <- nat.templatebrains::xform_brain(mc.positions.xyz, 
#                                                              sample = "MANC", 
#                                                              reference = "JRCVNC2018F")
# mc.positions.xyz.banc <- bancr::banc_to_JRC2018F(mc.positions.xyz.jrcvnc2018f, 
#                                                  region="vnc", 
#                                                  method="tpsreg", 
#                                                  banc.units = "nm", 
#                                                  inverse = TRUE)
# nat::xyzmatrix(mc.positions) <- mc.positions.xyz.banc

# # Preload MANC data
# message("##### Converting MANC skeletons to vector clouds #####")
# mc.swc.full.files <- unique(unlist(lapply(mc.ids, function(id) file.path(banc.nblast.manc.swc.save.path.version,paste0(id,".swc")))))
# mc.swc.full.files <- mc.swc.full.files[file.exists(mc.swc.full.files)]
# mc.skels <- nat::read.neurons(mc.swc.full.files)
# mc.skels <- manc_reroot(mc.skels,id=names(mc.skels),mc.positions=mc.positions)
# mc.skels <- prune_in_volume(mc.skels, surf = banc_vnc_neuropil.surf, invert = TRUE, OmitFailures = TRUE)
# mc.dps <- nat::dotprops(mc.skels/1000,
#                         k = 20,
#                         topo = FALSE,
#                         resample = 3,
#                         UseAlpha = 1,
#                         OmitFailures = TRUE,
#                         .parallel = FALSE)
# save(mc.dps,file=file.path(banc.nblast.manc.save.path,"banc_microns_space_manc_dps.rda"))
# rm('mc.skels')

# Load MANC neurons
load(file.path(banc.nblast.manc.save.path,"banc_microns_space_manc_dps.rda"))

# Batch for parallel processing
message("##### Setting up parallel processing #####")
multiplier <- 10
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
                                 
                                 # Look up soma info early for file existence check
                                 query.soma <- subset(banc.soma.positions,
                                                      banc.soma.positions$root_id == neuron.id)
                                 nucleus.id <- query.soma$nucleus_id[1]
                                 supervoxel.id <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_id)]
                                 if("supervoxel_id" %in% names(query.soma) && !is.na(query.soma$supervoxel_id[1])){
                                   supervoxel.id <- query.soma$supervoxel_id[1]
                                 }
                                 if(is.na(nucleus.id)|nucleus.id!="0"){
                                   use.alpha = FALSE
                                 }else{
                                   use.alpha = TRUE
                                 }

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
                                 if(is.null(query.soma)|!nrow(query.soma)){
                                   query.soma <- nat::xyzmatrix(banc.query)[nat::endpoints(banc.query[[1]]),]
                                   nucleus.id <- 0
                                 }
                                 
                                 # # Search for nearest rootpoints
                                 # search.targets <- nabor::knn(query = nat::xyzmatrix(query.soma), 
                                 #                              data = nat::xyzmatrix(mc.positions), 
                                 #                              radius = 150000, 
                                 #                              k = nrow(mc.positions)) # Must have a soma within 50 microns
                                 # if(nrow(search.targets$nn.idx)){
                                 #   chosen.targets <- unique(unlist(apply(search.targets$nn.idx, 1, function(row) row[row!=0])))
                                 # }else{
                                 #   chosen.targets <- search.targets$nn.idx[search.targets$nn.idx!=0]
                                 # }
                                 # chosen.targets <- unique(mc.positions[chosen.targets,]$bodyid)
                                 # if(!length(chosen.targets)){
                                 #   warning("no near neighbours for: ", neuron.id, " trying all with no somas")
                                 #   # This is FAFB, so find neurons with no soma, severed ANs/sensories
                                 #   chosen.targets <- mc.meta.no.soma.ids
                                 # }
                                 
                                 # # Run NBLAST
                                 # mc.dps.target <- mc.dps[names(mc.dps) %in% chosen.targets]
                                 # if(!length(mc.dps.target)){
                                 #   warning("hits not yet skeletonised: ", length(chosen.targets))
                                 #   message("hits not yet skeletonised: ", length(chosen.targets))
                                 #   next
                                 # }
                                 # nb.1 <- nat.nblast::nblast(query = banc.query.dps, target = mc.dps.target, 
                                 #                            UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 # nb.2 <- nat.nblast::nblast(query = mc.dps.target, target = banc.query.dps, 
                                 #                            UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 # nb <- (nb.1+nb.2)/2
                                 
                                 # # Is max score too low?
                                 # max.score <- max(nb, na.rm = TRUE)
                                 # if(max.score<0.65){
                                 #   # Run NBLAST
                                 #   mc.dps.target <- mc.dps
                                 #   nb.1 <- nat.nblast::nblast(query = banc.query.dps, target = mc.dps.target, 
                                 #                              UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 #   nb.2 <- nat.nblast::nblast(query = mc.dps.target, target = banc.query.dps, 
                                 #                              UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 #   nb <- (nb.1+nb.2)/2
                                 # }
                                 # Forward NBLAST: 1 query vs all targets
                                 nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = mc.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))

                                 # Reverse NBLAST on top 100 only
                                 top_k <- min(100, length(nb.fwd))
                                 top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                 nb.rev <- drop(nat.nblast::nblast(query = mc.dps[top_names], target = banc.query.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))
                                 nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                 # Make data frame (top 100 only)
                                 nb.df <- data.frame(nb = nb) %>%
                                   dplyr::arrange(dplyr::desc(nb))
                                 nb.df$bodyid <- rownames(nb.df)
                                 nb.df <- dplyr::left_join(nb.df, mc.meta[,c("bodyid","root_location","cell_type")], by = "bodyid")
                                 
                                 # Save NBLAST
                                 readr::write_csv(nb.df, file=file)
                                 message("Saved NBLAST result: ", file)
                                 success <- success+1
                                 
                               }
                               
                               # Announce
                               message("completed: ", length(neuron.ids), " neurons NBLASTed against MANC neurons with nearby root")

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

