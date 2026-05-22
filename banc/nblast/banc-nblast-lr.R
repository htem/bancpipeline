#' banc-nblast-lr — BANC self-NBLAST after thin-plate-spline left-right mirror.
#'
#' Surfaces left ↔ right partner candidates. Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, BANC L2 skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.mirror.save.path>/results/<root_888>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

##############################
### Left-right BANC NBLAST ###
##############################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))

# Direct us to the BANC dataset
bancr::choose_banc()
dir.results <- file.path(banc.nblast.mirror.save.path,"results")
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# Read IDs — use root_888 as stable pool identifier.
# Wrap banctable_query in tryCatch with snapshot fallback: bancr's banctable_query
# can throw "Error in !nrow(bc) : invalid argument type" when SeaTable returns
# no rows (transient API failure), which kills the script before the pipe even
# runs. Fall back to the latest dated snapshot written by banc-update-seatable.R
# (via the .banctable_read_snapshot helper defined in banc-startup.R).
banc.meta <- tryCatch(banctable_query(), error = function(e) {
  warning("banctable_query failed: ", e$message, " — falling back to latest snapshot")
  .banctable_read_snapshot()
})
if (is.null(banc.meta) || nrow(banc.meta) == 0)
  stop("banctable_query and snapshot fallback both returned no rows")
banc.meta <- banc.meta %>%
  banc_filter_neurons() %>%
  dplyr::filter(!is.na(root_888))
banc.root.ids <- unique(banc.meta$root_888)

# Build root_888 → current root_id mapping for file lookup
# L2 skeleton files are stored by current root_id, not root_888
banc.id.map <- banc.meta %>%
  dplyr::select(root_888, root_id) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)

# NBLAST files todo
dir.create(banc.nblast.mirror.save.path, showWarnings = FALSE)
nblast.files <- list.files(dir.results, pattern = "\\.csv")
nblast.files <- nblast.files[grepl("super",nblast.files)]

# Fix files with missing supervoxel_id in filename
na_svid_files <- nblast.files[grepl("supervoxel_id_NA_", basename(nblast.files))]
if (length(na_svid_files)) {
  if (!grepl(.Platform$file.sep, na_svid_files[1])) na_svid_files <- file.path(dir.results, na_svid_files)
  message(sprintf("  Fixing %d files with supervoxel_id_NA in filename...", length(na_svid_files)))
  na_root_ids <- gsub(".*_root_id_|\\.csv", "", basename(na_svid_files))
  svid_lookup <- banc.meta$supervoxel_id[match(na_root_ids, banc.meta$root_888)]
  can_fix <- !is.na(svid_lookup)
  if (any(can_fix)) {
    new_names <- file.path(dirname(na_svid_files[can_fix]),
                           paste0("supervoxel_id_", svid_lookup[can_fix],
                                  "_root_id_", na_root_ids[can_fix], ".csv"))
    renamed <- mapply(file.rename, na_svid_files[can_fix], new_names)
    message(sprintf("  Renamed %d/%d files", sum(renamed), sum(can_fix)))
  }
  nblast.files <- list.files(dir.results, pattern = "\\.csv")
  nblast.files <- nblast.files[grepl("super",nblast.files)]
}

nblast.done <- gsub(".*_root_id_|\\.csv", "", basename(nblast.files))
if(redo){
  nblast.todo <- banc.root.ids
}else{
  nblast.todo <- setdiff(banc.root.ids, nblast.done)
}
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  nblast.todo <- intersect(nblast.todo, banc.test.ids)
message(length(nblast.todo), " neurons to NBLAST")

if (length(nblast.todo) == 0) {
  message("All mirror NBLASTs up to date. Nothing to do.")
} else {

# Preload data — read L2 skeletons by current root_id, name by root_888
message("##### Reading BANC L2 skeletons #####")
banc.file.map <- banc.id.map %>%
  dplyr::mutate(swc_file = file.path(banc.l2swc.save.path, paste0(root_id, ".swc"))) %>%
  dplyr::filter(file.exists(swc_file)) %>%
  # Multiple root_888 can map to the same current root_id (post-merge), so dedup
  # on root_id to avoid nat::read.neurons "duplicate names" error. The survivor's
  # root_888 is paired with the skeleton; siblings sharing its root_id are the
  # same neuron post-merge and skipped.
  dplyr::distinct(root_id, .keep_all = TRUE)
# Some L2 SWCs fail to parse (e.g. embedded nulls in line 1). Use OmitFailures
# and re-name survivors by filename match so the names<- alignment can't fail.
banc.skels <- nat::read.neurons(banc.file.map$swc_file, OmitFailures = TRUE)
.read_rids   <- sub("\\.swc$", "", basename(names(banc.skels)))
.matched_888 <- banc.file.map$root_888[match(.read_rids, banc.file.map$root_id)]
.keep        <- !is.na(.matched_888)
banc.skels   <- banc.skels[.keep]
names(banc.skels) <- .matched_888[.keep]
.n_dropped   <- nrow(banc.file.map) - sum(.keep)
if (.n_dropped > 0L)
  message(sprintf("  Dropped %d/%d skeletons (parse failures or unmapped)",
                  .n_dropped, nrow(banc.file.map)))

# Write mirrored SWCs for neurons that don't have them yet
mirror.swc.dir <- "/n/data1/hms/neurobio/wilson/banc/banc_mirrored_swc"
for(id in names(banc.skels)){
  try({
    if(!file.exists(file.path(mirror.swc.dir, paste0(id, ".swc")))){
      nat::write.neurons(banc_mirror(banc.skels[id]), dir = mirror.swc.dir, Force = TRUE)
    }
  })
}

# Create dotprops in microns
banc.dps <- nat::dotprops(banc.skels/1000,
                          k = 20,
                          topo = TRUE,
                          resample = 2,
                          UseAlpha = 1,
                          OmitFailures = TRUE,
                          .parallel = FALSE)

# Mirror
banc.dps.m <- nat::nlapply(banc.dps, function(n){
  n <- n*1000
  new.points <- banc_mirror(nat::xyzmatrix(n), banc.units = "nm", method = "tpsreg")
  nat::xyzmatrix(n) <- new.points/1000
  n
},
OmitFailures = TRUE,
.parallel = FALSE)

# Remove from workspace
rm('banc.dps')
rm('banc.skels')

# Batch for parallel processing
multiplier <- 10
upper <- ifelse((numCores*multiplier)<length(nblast.todo),numCores*multiplier,length(nblast.todo))
batches <- split(nblast.todo, round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
batches <- sample(batches)

# Register cores
cl <- setup_parallel()

# NBLAST!
message("##### NBLASTing BANC against BANC mirror #####")
by.query <- foreach::foreach(batch = seq_along(batches),
                             .combine = 'c',
                             .multicombine = TRUE,
                             .init = list(),
                             .errorhandling = 'pass') %dopar% {

    # Select IDs (root_888)
    neuron.ids <- batches[[batch]]

    for(neuron.id in neuron.ids){

      tryCatch({

      # Map root_888 to current root_id for file lookup
      current.id <- banc.id.map$root_id[banc.id.map$root_888 == neuron.id]

      # Look up supervoxel_id early for file existence check
      sp.id <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_888)]
      file <- file.path(dir.results, paste0("supervoxel_id_",sp.id,"_root_id_",neuron.id,".csv"))
      if(file.exists(file)&!redo){
        next
      }

      # Read L2 skeleton directly by path (not directory scan)
      swc_file <- file.path(banc.l2swc.save.path, paste0(current.id, ".swc"))
      if(!file.exists(swc_file)){
        next
      }
      banc.query <- tryCatch(
        nat::read.neurons(swc_file, neuronnames = function(f) sub("\\.swc$", "", basename(f))),
        error = function(e) NULL
      )
      if(is.null(banc.query)|!length(banc.query)){
        next
      }
      # Rename to root_888 for consistent NBLAST naming
      names(banc.query) <- neuron.id
      banc.query.dps <- nat::dotprops(banc.query/1000,
                                      k = 20,
                                      topo = TRUE,
                                      resample = 2,
                                      UseAlpha = 1,
                                      OmitFailures = TRUE)
      if(is.null(banc.query.dps)|!length(banc.query.dps)){
        next
      }

      # Forward NBLAST: 1 query vs all mirrored targets
      use.alpha <- TRUE
      nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = banc.dps.m,
                                         UseAlpha = use.alpha, normalised = TRUE,
                                         smat = nat.nblast::smat_alpha.fcwb))

      # Reverse NBLAST on top 100 only
      top_k <- min(100, length(nb.fwd))
      top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
      nb.rev <- drop(nat.nblast::nblast(query = banc.dps.m[top_names], target = banc.query.dps,
                                         UseAlpha = use.alpha, normalised = TRUE,
                                         smat = nat.nblast::smat_alpha.fcwb))
      nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

      # Make data frame (top 100 only) — target names are root_888
      nb.df <- data.frame(nb = nb) %>%
        dplyr::arrange(dplyr::desc(nb))
      nb.df$root_888 <- rownames(nb.df)
      nb.df <- dplyr::left_join(nb.df, banc.meta[,c("root_888","root_id","nucleus_id","cell_type")], by = "root_888")

      # Save NBLAST
      readr::write_csv(nb.df, file=file)

      # Catch
      }, error = function(e) {
        warning(e)
      })

      # Return nothing
      NULL
    }

    # Announce
    message("completed: ", length(neuron.ids), " neurons NBLASTed against mirrored neurons")
}

# Stop cores
stop_parallel(cl)

# Were there errors?
message("##### Displaying any errors from foreach loop #####")
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}
} # end if mirror nblast.todo > 0

##########################
### native BANC NBLAST ###
##########################
dir.native.results <- file.path(banc.nblast.native.save.path,"results")
dir.create(dir.native.results, showWarnings = FALSE, recursive = TRUE)

# Remove from workspace
rm('banc.dps.m')

# NBLAST files todo
if(redo){
  nblast.todo <- banc.root.ids
}else{
  dir.create(banc.nblast.native.save.path, showWarnings = FALSE)
  nblast.files <- list.files(dir.native.results, pattern = "\\.csv")
  nblast.files <- nblast.files[grepl("super",nblast.files)]
  nblast.done <-  gsub(".*_root_id_|\\.csv","",basename(nblast.files))
  nblast.todo <- setdiff(banc.root.ids,nblast.done)
}
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  nblast.todo <- intersect(nblast.todo, banc.test.ids)
message(length(nblast.todo), " neurons to NBLAST")

if (length(nblast.todo) == 0) {
  message("All native NBLASTs up to date. Nothing to do.")
} else {

# Preload data — read L2 skeletons by current root_id, name by root_888
message("##### Reading BANC L2 skeletons #####")
banc.file.map <- banc.id.map %>%
  dplyr::mutate(swc_file = file.path(banc.l2swc.save.path, paste0(root_id, ".swc"))) %>%
  dplyr::filter(file.exists(swc_file)) %>%
  # Multiple root_888 can map to the same current root_id (post-merge), so dedup
  # on root_id to avoid nat::read.neurons "duplicate names" error. The survivor's
  # root_888 is paired with the skeleton; siblings sharing its root_id are the
  # same neuron post-merge and skipped.
  dplyr::distinct(root_id, .keep_all = TRUE)
# Same robust pattern as the mirror phase above.
banc.skels <- nat::read.neurons(banc.file.map$swc_file, OmitFailures = TRUE)
.read_rids   <- sub("\\.swc$", "", basename(names(banc.skels)))
.matched_888 <- banc.file.map$root_888[match(.read_rids, banc.file.map$root_id)]
.keep        <- !is.na(.matched_888)
banc.skels   <- banc.skels[.keep]
names(banc.skels) <- .matched_888[.keep]
.n_dropped   <- nrow(banc.file.map) - sum(.keep)
if (.n_dropped > 0L)
  message(sprintf("  Dropped %d/%d skeletons (parse failures or unmapped)",
                  .n_dropped, nrow(banc.file.map)))

# Create dotprops objects, in microns
banc.dps <- nat::dotprops(banc.skels/1000,
                          k = 20,
                          topo = TRUE,
                          resample = 2,
                          UseAlpha = 1,
                          OmitFailures = TRUE,
                          .parallel = FALSE)

# Remove from workspace
rm('banc.skels')

# Batch for parallel processing
multiplier <- 10
upper <- ifelse((numCores*multiplier)<length(nblast.todo),numCores*multiplier,length(nblast.todo))
batches <- split(nblast.todo, round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
batches <- sample(batches)

# Register cores
cl <- setup_parallel()

# NBLAST!
message("##### NBLASTing BANC against BANC native #####")
by.query <- foreach::foreach(batch = seq_along(batches),
                             .combine = 'c',
                             .multicombine = TRUE,
                             .init = list(),
                             .errorhandling = 'pass') %dopar% {

                               # Select IDs (root_888)
                               neuron.ids <- batches[[batch]]

                               for(neuron.id in neuron.ids){

                                 tryCatch({

                                   # Map root_888 to current root_id for file lookup
                                   current.id <- banc.id.map$root_id[banc.id.map$root_888 == neuron.id]

                                   # Look up supervoxel_id early for file existence check
                                   sp.id <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_888)]
                                   file <- file.path(dir.native.results, paste0("supervoxel_id_",sp.id,"_root_id_",neuron.id,".csv"))
                                   if(file.exists(file)&!redo){
                                     next
                                   }

                                   # Read L2 skeleton directly by path (not directory scan)
                                   swc_file <- file.path(banc.l2swc.save.path, paste0(current.id, ".swc"))
                                   if(!file.exists(swc_file)){
                                     next
                                   }
                                   banc.query <- tryCatch(
                                     nat::read.neurons(swc_file, neuronnames = function(f) sub("\\.swc$", "", basename(f))),
                                     error = function(e) NULL
                                   )
                                   if(is.null(banc.query)|!length(banc.query)){
                                     next
                                   }
                                   # Rename to root_888 for consistent NBLAST naming
                                   names(banc.query) <- neuron.id
                                   banc.query.dps <- nat::dotprops(banc.query/1000,
                                                                   k = 20,
                                                                   topo = TRUE,
                                                                   resample = 2,
                                                                   UseAlpha = 1,
                                                                   OmitFailures = TRUE)
                                   if(is.null(banc.query.dps)|!length(banc.query.dps)){
                                     next
                                   }

                                   # Forward NBLAST: 1 query vs all native targets
                                   use.alpha <- TRUE
                                   nb.fwd <- drop(nat.nblast::nblast(query = banc.query.dps, target = banc.dps,
                                                                      UseAlpha = use.alpha, normalised = TRUE,
                                                                      smat = nat.nblast::smat_alpha.fcwb))

                                   # Reverse NBLAST on top 100 only
                                   top_k <- min(100, length(nb.fwd))
                                   top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
                                   nb.rev <- drop(nat.nblast::nblast(query = banc.dps[top_names], target = banc.query.dps,
                                                                      UseAlpha = use.alpha, normalised = TRUE,
                                                                      smat = nat.nblast::smat_alpha.fcwb))
                                   nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

                                   # Make data frame (top 100 only) — names are root_888
                                   nb.df <- data.frame(nb = nb) %>%
                                     dplyr::arrange(dplyr::desc(nb))
                                   nb.df$root_888 <- rownames(nb.df)
                                   nb.df <- dplyr::left_join(nb.df, banc.meta[,c("root_888","root_id","nucleus_id","cell_type")], by = "root_888")

                                   # Save NBLAST
                                   readr::write_csv(nb.df, file=file)

                                   # Catch
                                 }, error = function(e) {
                                   warning(e)
                                 })

                                 # Return nothing
                                 NULL
                               }

                               # Announce
                               message("completed: ", length(neuron.ids), " neurons NBLASTed against native neurons")
                             }

# Stop cores
stop_parallel(cl)

# Were there errors?
message("##### Displaying any errors from foreach loop #####")
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}
} # end if native nblast.todo > 0
