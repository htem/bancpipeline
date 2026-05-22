#' banc-fafb-nblast — Pairwise BANC ↔ FAFB (v783) NBLAST.
#'
#' Resumable; skips done queries unless `BANC_NBLAST_REDO=TRUE`.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - FAFB BANC-space skeletons under `<banc.nblast.fafb.swc.save.path>/<ver>`
#'
#' @section Writes:
#'   - `<banc.nblast.fafb.save.path>/results/<ver>/<root_id>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.
#'
#' @section Used by:
#'   `banc-nblast-compile.R` (per-query CSVs → consolidated feather).
#'
#' @section Paper:
#'   Methods §"NBLAST cross-dataset matching".

########################
### BANC-FAFB NBLAST ###
########################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))
version <- banc.nblast.version

# Direct us to the BANC dataset
banc.nblast.fafb.swc.save.path.version <- file.path(banc.nblast.fafb.swc.save.path,version)
dir.results <- file.path(banc.nblast.fafb.save.path,"results", version) 
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# Read IDs
banc.meta <- tryCatch(banctable_query(), error = function(e) {
  warning("banctable_query failed, falling back to CSV: ", e$message)
  readr::read_csv(file.path(banc.meta.save.path, "banc_meta.csv"),
                  col_types = banc.col.types, show_col_types = FALSE)
})
banc.meta <- banc.meta %>%
  dplyr::filter(!is.na(root_id))
  #dplyr::filter(!region %in% c("vnc","ventral_nerve_cord","optic_lobe","optic_lobes"))
banc.root.ids <- unique(banc.meta$root_id)

# Fill missing supervoxel_ids from CAVE cell_info table
missing_svids <- is.na(banc.meta$supervoxel_id) | banc.meta$supervoxel_id == ""
if (any(missing_svids)) {
  message(sprintf("  %d/%d neurons missing supervoxel_id — looking up from CAVE...",
                  sum(missing_svids), nrow(banc.meta)))
  cave_svids <- tryCatch({
    dplyr::bind_rows(
      banc_cell_info(rawcoords = TRUE) %>%
        dplyr::distinct(.data$pt_root_id, .data$pt_supervoxel_id),
      banc_nuclei(rawcoords = TRUE) %>%
        dplyr::distinct(pt_root_id = .data$root_id, .data$pt_supervoxel_id)
    ) %>%
      dplyr::distinct(.data$pt_root_id, .keep_all = TRUE)
  }, error = function(e) {
    warning("CAVE lookup for supervoxel_ids failed: ", e$message)
    NULL
  })
  if (!is.null(cave_svids)) {
    idx <- match(banc.meta$root_id[missing_svids], cave_svids$pt_root_id)
    filled <- !is.na(idx)
    banc.meta$supervoxel_id[missing_svids][filled] <- cave_svids$pt_supervoxel_id[idx[filled]]
    message(sprintf("  Filled %d/%d missing supervoxel_ids from CAVE",
                    sum(filled), sum(missing_svids)))
  }
}

# Get fafb meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)
# fw.meta.no.soma <- fw.meta %>%
#   dplyr::filter(is.na(soma_x), 
#                 super_class%in%c("ascending","sensory_ascending","sensory") | flow%in%c("afferent"))
# fw.meta.no.soma.ids <- unique(fw.meta.no.soma$root_783)

# NBLAST files todo
message("##### Calculating NBLASTs we need to do #####")
dir.create(banc.nblast.fafb.save.path, showWarnings = FALSE)
nblast.files <- list_files_age_sorted(dir.results, pattern = "\\.csv")

# Fix files with missing supervoxel_id in filename
na_svid_files <- nblast.files[grepl("supervoxel_id_NA_", basename(nblast.files))]
if (length(na_svid_files)) {
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
  if (any(!can_fix)) {
    message(sprintf("  %d files still have no supervoxel_id mapping", sum(!can_fix)))
  }
  # Refresh file list after renames
  nblast.files <- list_files_age_sorted(dir.results, pattern = "\\.csv")
}

nblast.done <- gsub(".*_root_id_|\\.csv", "", basename(nblast.files))
if(redo){
  nblast.todo <- banc.root.ids
  undone <- setdiff(banc.root.ids, nblast.done)
  nblast.todo <- union(undone, nblast.done)
}else{
  nblast.todo <- sample(setdiff(banc.root.ids, nblast.done))
}
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  nblast.todo <- intersect(nblast.todo, banc.test.ids)
message(length(nblast.todo), " neurons to NBLAST, done: ", length(nblast.done))

if (length(nblast.todo) == 0) {
  message("All FAFB NBLASTs up to date. Nothing to do.")
} else {

# Read root nodes
y.cut <- 325000
banc.soma.positions <- arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather"))
has_pos <- !is.na(banc.soma.positions$root_position_nm) & banc.soma.positions$root_position_nm != ""
banc.soma.positions[has_pos, c("X", "Y", "Z")] <- nat::xyzmatrix(banc.soma.positions$root_position_nm[has_pos])
banc.soma.positions <- banc_decapitate(banc.soma.positions, invert = TRUE)

# fafb soma positions in banc space
fw.positions <- fw.meta %>%
  dplyr::mutate(X = dplyr::case_when(
    !is.na(soma_x) ~ soma_x,
    TRUE ~ pos_x
  ),
  Y = dplyr::case_when(
    !is.na(soma_y) ~ soma_y,
    TRUE ~ pos_y
  ),
  Z = dplyr::case_when(
    !is.na(soma_z) ~ soma_z,
    TRUE ~ pos_z
  )) %>%
  dplyr::distinct(root_783, nucleus_id, X, Y, Z)

# # Transforms into BANC space
# message("##### Transforming FAFB roots to BANC space #####")
# fw.positions.xyz <- flywire_raw2nm(nat::xyzmatrix(fw.positions))
# fw.positions.xyz.jrc2018f <- nat.templatebrains::xform_brain(fw.positions.xyz, 
#                                                              sample = "FAFB14", 
#                                                              reference = "JRC2018F")
# fw.positions.xyz.jrc2018f.m <- nat.templatebrains::mirror_brain(fw.positions.xyz.jrc2018f, 
#                                                                 brain = nat.flybrains::JRC2018F, 
#                                                                 transform = "flip")
# fw.positions.xyz.banc <- bancr::banc_to_JRC2018F(fw.positions.xyz.jrc2018f.m, region = "brain", method = "tpsreg", inverse = TRUE)
# nat::xyzmatrix(fw.positions) <- fw.positions.xyz.banc

# Preload FAFB data
# message("##### Converting FAFB skeletons to vector clouds #####")
# fw.swc.full.files <- unique(unlist(lapply(fw.ids, function(id) file.path(banc.nblast.fafb.swc.save.path.version,paste0(id,".swc")))))
# fw.swc.full.files <- fw.swc.full.files[file.exists(fw.swc.full.files)]
# fw.skels <- nat::read.neurons(fw.swc.full.files)
# fw.dps <- nat::dotprops(fw.skels/1000, 
#                           k = 20,
#                           topo = FALSE, 
#                           resample = 2, 
#                           UseAlpha = 1,
#                           OmitFailures = TRUE,
#                           .parallel = FALSE)
# save(fw.dps,file.path(banc.nblast.fafb.save.path,"banc_microns_space_fafb_dps.rda"))
# save(fw.skels,file.path(banc.nblast.fafb.save.path,"banc_microns_space_fafb_skels.rda"))
# rm('fw.skels')

# Remove unneeded objects from workspace
load(file.path(banc.nblast.fafb.save.path,"banc_microns_space_fafb_dps.rda"))

# Batch for parallel processing
message("##### Setting up parallel processing #####")
multiplier <- 10
upper <- ifelse((numCores*multiplier)<length(nblast.todo),numCores*multiplier,length(nblast.todo))
batches <- split(nblast.todo, round(seq(from = 1, to = upper, length.out = length(nblast.todo))))

# Register cores
# doMC::registerDoMC(1)
cl <- setup_parallel()

# NBLAST!
message("##### Running BANC-FAFB NBLAST #####")
by.query <- foreach::foreach(batch = sample(seq_along(batches)),
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
                                 file <- file.path(dir.results, paste0("supervoxel_id_",supervoxel.id, "_root_id_",neuron.id,".csv"))
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
                                   next
                                 }
                                 banc.query.dps <- nat::dotprops(banc.query/1000,
                                                                 k = 20,
                                                                 topo = FALSE,
                                                                 resample = ifelse(summary(banc.query)$cable.length<30000,1,2),
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
                                 #                              data = nat::xyzmatrix(fw.positions), 
                                 #                              radius = 150000, 
                                 #                              k = nrow(fw.positions)) # Must have a soma within 50 microns
                                 # if(nrow(search.targets$nn.idx)){
                                 #   chosen.targets <- unique(unlist(apply(search.targets$nn.idx, 1, function(row) row[row!=0])))
                                 # }else{
                                 #   chosen.targets <- search.targets$nn.idx[search.targets$nn.idx!=0]
                                 # }
                                 # chosen.targets <- unique(fw.positions[chosen.targets,]$root_783)
                                 # if(!length(chosen.targets)){
                                 #   warning("no near neighbours for: ", neuron.id, " trying all with no somas")
                                 #   # This is FAFB, so find neurons with no soma, severed ANs/sensories
                                 #   chosen.targets <- fw.meta.no.soma.ids
                                 # }
                                 # 
                                 # # Run NBLAST
                                 # fw.dps.target <- fw.dps[names(fw.dps) %in% chosen.targets]
                                 # missing.hits <- sum(!chosen.targets %in% names(fw.dps))
                                 # if(missing.hits){
                                 #   warning("hits not yet skeletonised: ", missing.hits)
                                 #   message("hits not yet skeletonised: ", missing.hits)
                                 # }
                                 # if(!length(fw.dps.target)){
                                 #   next
                                 # }
                                 # nb.1 <- nat.nblast::nblast(query = banc.query.dps, target = fw.dps.target, 
                                 #                            UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 # nb.2 <- nat.nblast::nblast(query = fw.dps.target, target = banc.query.dps, 
                                 #                            UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 # nb <- (nb.1+nb.2)/2
                                 # 
                                 # # Is max score too low?
                                 # max.score <- max(nb, na.rm = TRUE)
                                 # if(max.score<0.65){
                                 #   fw.dps.target <- fw.dps
                                 #   nb.1 <- nat.nblast::nblast(query = banc.query.dps, target = fw.dps.target, 
                                 #                              UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 #   nb.2 <- nat.nblast::nblast(query = fw.dps.target, target = banc.query.dps, 
                                 #                              UseAlpha = use.alpha, normalised = TRUE, smat = nat.nblast::smat_alpha.fcwb)
                                 #   nb <- (nb.1+nb.2)/2
                                 # }
                                 
                                 # Forward NBLAST: 1 query vs all targets
                                 nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = fw.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))

                                 # Reverse NBLAST on top 100 only (~100x faster than full reverse)
                                 top_k <- min(100, length(nb.fwd))
                                 top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                 nb.rev <- drop(nat.nblast::nblast(query = fw.dps[top_names], target = banc.query.dps,
                                                                    UseAlpha = use.alpha, normalised = TRUE,
                                                                    smat = nat.nblast::smat_alpha.fcwb))
                                 nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                 # Make data frame (top 100 only)
                                 nb.df <- data.frame(nb = nb) %>%
                                   dplyr::arrange(dplyr::desc(nb))
                                 nb.df$root_783 <- rownames(nb.df)
                                 nb.df <- dplyr::left_join(nb.df, fw.meta[,c("root_783","nucleus_id","cell_type")], by = "root_783")
                                 
                                 # Save NBLAST
                                 readr::write_csv(nb.df, file=file)
                                 message("Saved NBLAST result: ", file)
                                 success <- success+1
                                 
                               }
                               
                               # Announce
                               message("completed: ", length(neuron.ids), " neurons NBLASTed against FAFB neurons with nearby root")
                               
                               # Return success ratio
                               success/length(neuron.ids)
                             }

# Stop cores
stop_parallel(cl)

# Were there errors?
message("##### NBLAST run, displaying any errors from foreach loop #####")
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    if(!is.numeric(by.query[[i]])){
      message(by.query[[i]])
    }
  }
}
} # end if nblast.todo > 0
