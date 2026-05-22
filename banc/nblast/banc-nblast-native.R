#' banc-nblast-native — BANC self-NBLAST in native space (no mirroring).
#'
#' Surfaces ipsilateral sister cells, partial duplicates, and within-side
#' morphology clusters that `banc-nblast-lr.R` (mirror) cannot find.
#' Resumable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, BANC L2 skeletons
#'
#' @section Writes:
#'   - `<banc.nblast.native.save.path>/results/supervoxel_id_<sv>_root_id_<root_888>.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_native.sh`, `o2/production/o2_banc_native_array.sh`.

##############################################
### BANC self NBLAST (native, no mirroring) ###
##############################################
###
### Standalone script for the within-BANC native (non-mirrored)
### NBLAST. Each query is NBLASTed against every other BANC neuron
### in its native left/right configuration — useful for finding
### within-side ipsilateral matches (sister cells, partial duplicates,
### anatomical co-clusters) that L-R mirror NBLAST can't surface.
###
### Mirror NBLAST (banc-nblast-lr.R) covers contralateral/ipsi-mirror
### partners; this script covers the native-side neighbourhood.
###
### Output: per-neuron CSVs at
###   banc.nblast.native.save.path/results/
###   filename = supervoxel_id_<sv>_root_id_<root_888>.csv
###
### Resumable: skips any (root_888) that already has a CSV unless
### redo = TRUE.
##############################################
source("banc/banc-startup.R")
redo <- isTRUE(as.logical(Sys.getenv("BANC_NBLAST_REDO", "FALSE")))

bancr::choose_banc()
dir.results <- file.path(banc.nblast.native.save.path, "results")
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

# Universe — current proofread BANC neurons keyed by root_888 (stable across
# segmentation versions for this run).
banc.meta <- banctable_query() %>%
  banc_filter_neurons() %>%
  dplyr::filter(!is.na(root_888))
banc.root.ids <- unique(banc.meta$root_888)

# root_888 -> current root_id mapping for L2 SWC file lookup
banc.id.map <- banc.meta %>%
  dplyr::select(root_888, root_id) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)

# Resume — skip neurons that already have an output CSV
nblast.files <- list.files(dir.results, pattern = "\\.csv")
nblast.files <- nblast.files[grepl("super", nblast.files)]
nblast.done  <- gsub(".*_root_id_|\\.csv", "", basename(nblast.files))
if (redo) {
  nblast.todo <- banc.root.ids
} else {
  nblast.todo <- setdiff(banc.root.ids, nblast.done)
}
if (exists("banc.test.ids", envir = .GlobalEnv))
  nblast.todo <- intersect(nblast.todo, banc.test.ids)

# SLURM array sharding — split nblast.todo across array tasks. Each task
# processes a deterministic, disjoint stride of the todo list. The
# resume-by-existing-CSV logic above is double protection against any race.
shard_id    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",    unset = "-1"))
shard_total <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "-1"))
if (shard_id >= 0L && shard_total > 1L) {
  # Sort once so partitioning is deterministic across array tasks
  nblast.todo <- sort(nblast.todo)
  keep <- (seq_along(nblast.todo) - 1L) %% shard_total == shard_id
  nblast.todo <- nblast.todo[keep]
  message(sprintf("  SLURM array shard %d/%d: %d neurons assigned",
                  shard_id, shard_total, length(nblast.todo)))
} else {
  message(length(nblast.todo), " neurons to NBLAST (native, single job)")
}

if (length(nblast.todo) == 0L) {
  message("All native NBLASTs up to date. Nothing to do.")
} else {

# Preload all L2 skeletons in native space and dotprops them ONCE.
# Cache banc.dps to avoid the ~5h sequential build on every run/shard. The
# cache is keyed on the banc.l2swc.save.path mtime so it stays fresh when
# new L2 SWCs land (a single mtime check is cheap; rebuilding only happens
# if the cache is stale or absent).
dps_cache_dir  <- file.path(banc.save.path, "cache")
dir.create(dps_cache_dir, showWarnings = FALSE, recursive = TRUE)
dps_cache_file <- file.path(dps_cache_dir, "banc_native_dps.rds")
l2_dir_mtime   <- file.info(banc.l2swc.save.path)$mtime
cache_fresh    <- file.exists(dps_cache_file) &&
                  difftime(file.info(dps_cache_file)$mtime,
                           l2_dir_mtime, units = "secs") > 0

if (cache_fresh) {
  message(sprintf("##### Loading cached banc.dps from %s #####", dps_cache_file))
  banc.dps <- readRDS(dps_cache_file)
  message(sprintf("  loaded %d dotprops", length(banc.dps)))
} else {
  message("##### Reading BANC L2 skeletons #####")
  banc.file.map.full <- banc.id.map %>%
    dplyr::mutate(swc_file = file.path(banc.l2swc.save.path, paste0(root_id, ".swc"))) %>%
    dplyr::filter(file.exists(swc_file))
  # Multiple root_888 can map to the same current root_id (post-v888 merges).
  # That gives duplicate swc_file entries → nat::read.neurons aborts with
  # "Neurons cannot have duplicate names" (OmitFailures only catches per-file
  # parse errors, not the dispatcher-level dedup check). Dedupe by swc_file
  # for the read; nblast.todo is unchanged so each root_888 still gets its own
  # output CSV (the foreach below reads the SWC fresh per query — duplicates
  # produce identical NBLAST results under different root_888-keyed filenames,
  # which is the desired behavior for downstream compile).
  banc.file.map <- banc.file.map.full %>%
    dplyr::distinct(swc_file, .keep_all = TRUE)
  n_shared_root <- nrow(banc.file.map.full) - nrow(banc.file.map)
  if (n_shared_root > 0L)
    message(sprintf("  %d root_888 share a current root_id with another; will reuse the same skeleton",
                    n_shared_root))

  # Some L2 SWCs fail to parse (e.g. embedded nulls in line 1). OmitFailures
  # returns only the successful reads; we then map the survivors back to
  # root_888 by filename so the names alignment can never mismatch.
  banc.skels <- nat::read.neurons(banc.file.map$swc_file, OmitFailures = TRUE)
  read_rids   <- sub("\\.swc$", "", basename(names(banc.skels)))
  matched_888 <- banc.file.map$root_888[match(read_rids, banc.file.map$root_id)]
  keep        <- !is.na(matched_888)
  banc.skels  <- banc.skels[keep]
  names(banc.skels) <- matched_888[keep]
  n_dropped <- length(read_rids) - sum(keep) +
               (nrow(banc.file.map) - length(read_rids))
  if (n_dropped > 0L)
    message(sprintf("  Dropped %d/%d skeletons (parse failures or unmapped)",
                    n_dropped, nrow(banc.file.map)))

  banc.dps <- nat::dotprops(banc.skels / 1000,
                            k          = 20,
                            topo       = TRUE,
                            resample   = 2,
                            UseAlpha   = 1,
                            OmitFailures = TRUE,
                            .parallel  = FALSE)
  rm(banc.skels); gc(verbose = FALSE)

  # Atomically write the cache so multiple shards racing don't tear the file.
  tmp_cache <- paste0(dps_cache_file, ".tmp.", Sys.getpid())
  saveRDS(banc.dps, tmp_cache, compress = FALSE)
  file.rename(tmp_cache, dps_cache_file)
  message(sprintf("  cached banc.dps -> %s (%d dotprops)",
                  dps_cache_file, length(banc.dps)))
}

# Batch + parallel
multiplier <- 10
upper   <- min(numCores * multiplier, length(nblast.todo))
batches <- split(nblast.todo,
                 round(seq(from = 1, to = upper, length.out = length(nblast.todo))))
batches <- sample(batches)

cl <- setup_parallel()
on.exit(stop_parallel(cl), add = TRUE)

message("##### NBLASTing BANC against BANC native #####")
by.query <- foreach::foreach(batch = seq_along(batches),
                             .combine     = "c",
                             .multicombine = TRUE,
                             .init        = list(),
                             .errorhandling = "pass") %dopar% {
  neuron.ids <- batches[[batch]]
  for (neuron.id in neuron.ids) {
    tryCatch({
      current.id <- banc.id.map$root_id[banc.id.map$root_888 == neuron.id]
      sp.id      <- banc.meta$supervoxel_id[match(neuron.id, banc.meta$root_888)]
      out_file   <- file.path(dir.results,
                              paste0("supervoxel_id_", sp.id,
                                     "_root_id_",      neuron.id, ".csv"))
      if (file.exists(out_file) && !redo) next

      swc_file <- file.path(banc.l2swc.save.path, paste0(current.id, ".swc"))
      if (!file.exists(swc_file)) next

      banc.query <- tryCatch(
        nat::read.neurons(swc_file,
                          neuronnames = function(f) sub("\\.swc$", "", basename(f))),
        error = function(e) NULL
      )
      if (is.null(banc.query) || !length(banc.query)) next
      names(banc.query) <- neuron.id

      banc.query.dps <- nat::dotprops(banc.query / 1000,
                                      k          = 20,
                                      topo       = TRUE,
                                      resample   = 2,
                                      UseAlpha   = 1,
                                      OmitFailures = TRUE)
      if (is.null(banc.query.dps) || !length(banc.query.dps)) next

      use.alpha <- TRUE
      nb.fwd <- drop(nat.nblast::nblast(
        query      = banc.query.dps,
        target     = banc.dps,
        UseAlpha   = use.alpha,
        normalised = TRUE,
        smat       = nat.nblast::smat_alpha.fcwb))

      top_k     <- min(100L, length(nb.fwd))
      top_names <- names(sort(nb.fwd, decreasing = TRUE))[1:top_k]
      nb.rev <- drop(nat.nblast::nblast(
        query      = banc.dps[top_names],
        target     = banc.query.dps,
        UseAlpha   = use.alpha,
        normalised = TRUE,
        smat       = nat.nblast::smat_alpha.fcwb))
      nb <- (nb.fwd[top_names] + nb.rev[top_names]) / 2

      nb.df <- data.frame(nb = nb) %>%
        dplyr::arrange(dplyr::desc(nb))
      nb.df$root_888 <- rownames(nb.df)
      nb.df <- dplyr::left_join(
        nb.df,
        banc.meta[, c("root_888", "root_id", "nucleus_id", "cell_type")],
        by = "root_888")

      readr::write_csv(nb.df, file = out_file)
    },
    error = function(e) warning(conditionMessage(e)))
    NULL
  }
  message("completed batch: ", length(neuron.ids), " neurons (native)")
}

stop_parallel(cl)

message("##### Errors from foreach loop #####")
for (i in seq_along(by.query)) {
  if (!is.null(by.query[[i]])) message(by.query[[i]])
}

}  # end if length(nblast.todo) > 0
