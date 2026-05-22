#' banc-synapses-v3-optimised — Bounding-box-accelerated v3 synapse classification.
#'
#' Same pipeline as `banc-synapses-v3.R` but `classify_chunk` adds an XYZ
#' bbox pre-filter before each `alphashape3d::inashape3d()` call, cutting
#' point-in-surface tests from ~5M to ~10⁴–10⁵ per neuropil.
#'
#' @section Reads:
#'   - GCS v3 raw synapse parquet shards, SeaTable `banc_meta`, CAVE
#'
#' @section Writes:
#'   - per-batch parquet caches, capture-rate CSVs, v2↔v3 comparison CSV
#'
#' @section Notes:
#'   - Local v3 parquet had ~30% spurious "0" root_ids; the trusted v3
#'     source is now the CAVE-ingested GCS export (see
#'     `banc-calculate-connectivity.R`). This script remains for spatial
#'     neuropil annotation only.
#'   - Test mode: `TEST_MODE = TRUE`.

###########################################################
### banc-synapses-v3-optimised.R
###
### Process the v3 synapse predictions (260326_assignment)
### and compare capture rates against v2.
###
### OPTIMISATION (2026-04-10):
###   classify_chunk now applies a bounding-box pre-filter before each
###   alphashape3d::inashape3d() call. For each alpha shape (volume or
###   neuropil), the XYZ bounding box of the shape's vertices is computed
###   once upfront and stored alongside the shape. At classification time,
###   only points that fall within a neuropil's bounding box are tested
###   against the full alpha shape. Because most neuropils are small
###   relative to the whole CNS, this dramatically reduces the number of
###   expensive point-in-surface tests (from ~5M per neuropil to typically
###   a few thousand to a few hundred thousand).
###
### Pipeline stages (all resumable):
###   0. Setup + test mode flag
###   1. Download raw parquet shards from GCS
###   2. Region + neuropil assignment via alpha shape point-in-surface
###      tests, in parallel chunks. Cached per-batch parquet output.
###   3. Supervoxel -> root_id via banc.meta cache + CAVE fallback.
###   4. Classify pre/post_status as neuron/fragment and compute
###      4 capture rate CSVs (gross, inout, region, neuropil).
###   5. v2 <-> v3 comparison CSV.
###   6. Push outputs to GCS (optional, guarded by PUSH_TO_GCS).
###
### Test mode: set TEST_MODE=TRUE at top to download only a few
### shards, sample 1000 synapses, and run end-to-end in seconds.
###########################################################
source("banc/banc-startup.R")

##########################################
### CONFIG                            ###
##########################################

TEST_MODE      <- FALSE     # TRUE = sample 1k synapses from a few shards
TEST_N         <- 1000L
TEST_SHARDS    <- 20L       # how many raw shards to pull in test mode (each is small)

V3_NAME        <- "synapses_v3"
V3_RELEASE     <- "260326_assignment"
V3_GCS_SRC     <- sprintf("gs://zetta_lee_fly_cns_001_synapse/%s/seg/metadata_scored",
                           V3_RELEASE)

# Target version: normally banc.version, but can be overridden via
# BANC_V3_TARGET_VERSION env var to run a one-off finalization for a prior
# version (e.g. produce banc_850_synapses_v3.parquet from the v888-default
# codebase, so it can serve as the basis for a subsequent v888 reuse run).
# Affects ONLY the paths/version literals in this script — banc.version
# (sourced from banc-startup.R) is left untouched so other code paths still
# see the project's current version.
target.version <- {
  v <- Sys.getenv("BANC_V3_TARGET_VERSION", unset = "")
  if (!nzchar(v)) banc.version else v
}

V3_GCS_DST     <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/v%s/%s",
                           target.version, V3_NAME)

# v3 centroid units: 16x16x45 nm per voxel (from info file "resolution")
V3_VOXEL_NM    <- c(16, 16, 45)

BATCH_ROWS     <- 5e6L      # rows per region/neuropil classification batch
CAVE_CHUNKSIZE <- 1e5L      # svids per CAVE call (fafbseg default)
# Size threshold removed (2026-04-19): was MIN_SIZE = 5L. The synapse-id ->
# neuropil correspondence should be COMPLETE (no size filter) so downstream
# analyses can choose any size threshold. Filtering happens in the consumer
# (see SIZE_THRESHOLD in banc/metrics/banc-calculate-completion.R).
MIN_SIZE       <- 0L        # 0 = no filter; preserved for backward-compat

PUSH_TO_GCS    <- TRUE

# Reuse mode: when targeting a new BANC version (e.g. v888) and a completed
# v3 run already exists for a prior version (e.g. v850), skip Stages 1+2
# (GCS download + alpha-shape region/neuropil classification) and just re-run
# Stage 3 (svid -> root_id) against CAVE for the new version. Stable across
# versions: syn_id, X/Y/Z coords, region, neuropil, side.
# Controlled by BANC_V3_REUSE_BASIS env var (e.g. "850"); unset = no reuse.
V3_REUSE_BASIS <- {
  v <- Sys.getenv("BANC_V3_REUSE_BASIS", unset = "")
  if (!nzchar(v)) NA_character_ else v
}

##########################################
### PATHS                              ###
##########################################

# In test mode use a stable local scratch dir (banc.save.path is on O2 and
# doesn't exist locally). In production use the cluster path.
if (TEST_MODE) {
  banc.synapses.v3.save.path <- normalizePath("~/banc_synapses_v3_test", mustWork = FALSE)
} else {
  banc.synapses.v3.save.path <- file.path(banc.save.path, V3_NAME)
}
raw_dir       <- file.path(banc.synapses.v3.save.path, "raw")
# Stage-2 outputs are versioned by target.version so a reuse-based re-run for a
# new BANC version doesn't clobber the previous version's root_id-annotated
# batch files.
processed_dir <- file.path(banc.synapses.v3.save.path, "processed",
                           paste0("v", target.version))
cache_dir     <- file.path(banc.synapses.v3.save.path, "cache")
rates_dir     <- file.path(banc.synapses.v3.save.path, "capture_rates")
for (d in c(banc.synapses.v3.save.path, raw_dir, processed_dir, cache_dir, rates_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

final_parquet      <- file.path(banc.synapses.v3.save.path,
                                 sprintf("banc_%s_synapses_v3.parquet", target.version))
neuropil_lookup    <- file.path(banc.synapses.v3.save.path,
                                 sprintf("banc_%s_synapses_v3_neuropil_lookup.parquet",
                                         target.version))
# Per-target-version svid->root_id cache. Was a single shared file
# `svid_rootid_cache.parquet` until 2026-04-19, when we discovered that a v850
# run filled the cache with v850 root_ids and a subsequent v888 reuse run read
# the same cache and never queried CAVE — silently producing
# `banc_888_synapses_v3.parquet` with v850 root_ids mislabeled as v888.
# Version-keying makes each (target.version, svid) lookup independent.
svid_cache_file    <- file.path(cache_dir,
                                 sprintf("svid_rootid_cache_v%s.parquet",
                                         target.version))
# Legacy unsuffixed cache file — useful for one-time migration / inspection
# but never read going forward (different versions need fresh CAVE queries).
svid_cache_file_legacy <- file.path(cache_dir, "svid_rootid_cache.parquet")

# Legacy (pre-versioned-processed-dir) location. Some earlier runs wrote batch
# files directly under `processed/` without a version subdir. If the current
# version's processed/ is empty and V3_REUSE_BASIS isn't set, but legacy files
# exist, we'll pick them up automatically.
legacy_processed_dir <- file.path(banc.synapses.v3.save.path, "processed")

# Reuse-basis file paths — resolved when V3_REUSE_BASIS is set.
if (!is.na(V3_REUSE_BASIS)) {
  basis_final_parquet <- file.path(banc.synapses.v3.save.path,
                                    sprintf("banc_%s_synapses_v3.parquet",
                                            V3_REUSE_BASIS))
  basis_neuropil_lookup <- file.path(banc.synapses.v3.save.path,
                                      sprintf("banc_%s_synapses_v3_neuropil_lookup.parquet",
                                              V3_REUSE_BASIS))
}

##########################################
### 0. MODE BANNER                     ###
##########################################

message(sprintf("\n### banc-synapses-v3 | TEST_MODE=%s | target_version=%s ###",
                TEST_MODE, target.version))
if (target.version != banc.version) {
  message(sprintf("  TARGET OVERRIDE: target.version=%s, banc.version=%s",
                  target.version, banc.version))
}
message(sprintf("  Working directory: %s", banc.synapses.v3.save.path))
if (!is.na(V3_REUSE_BASIS)) {
  message(sprintf("  REUSE MODE: basing off v%s final parquet (Stages 1+2 skipped)",
                  V3_REUSE_BASIS))
}

# If processed_dir for the target version is empty but a legacy flat
# processed/ has batch files (pre-versioning refactor), adopt those in place.
# Skip migration in REUSE mode -- there Stage 2 synthesises a fresh batch from
# the basis parquet, and stealing the legacy batches would conflate Stage-2
# outputs across versions.
if (length(list.files(processed_dir, pattern = "^batch_.*\\.parquet$")) == 0L &&
    length(list.files(legacy_processed_dir,
                      pattern = "^batch_.*\\.parquet$")) > 0L &&
    is.na(V3_REUSE_BASIS)) {
  legacy_batches <- list.files(legacy_processed_dir,
                                pattern = "^batch_.*\\.parquet$",
                                full.names = TRUE)
  message(sprintf("  migrating %d legacy batch files from %s -> %s",
                  length(legacy_batches), legacy_processed_dir, processed_dir))
  for (lf in legacy_batches) {
    file.rename(lf, file.path(processed_dir, basename(lf)))
  }
}

##########################################
### 1. DOWNLOAD RAW SHARDS FROM GCS   ###
##########################################

if (!is.na(V3_REUSE_BASIS)) {
  message("\n### Stage 1: skipped (reuse mode) ###")
} else {

message("\n### Stage 1: download raw shards ###")

# List all shards on GCS
t0 <- Sys.time()
gcs_shards <- system2("gsutil", c("ls", paste0(V3_GCS_SRC, "/data/")),
                       stdout = TRUE, stderr = FALSE)
gcs_shards <- gcs_shards[grepl("\\.parquet$", gcs_shards)]
message(sprintf("  %d shards on GCS", length(gcs_shards)))

if (TEST_MODE) {
  # Random sample shards across the whole volume so the test is representative
  set.seed(42)
  gcs_shards <- sample(gcs_shards, size = min(TEST_SHARDS, length(gcs_shards)))
  message(sprintf("  TEST_MODE: using %d randomly sampled shards", length(gcs_shards)))
}

# Determine which shards are missing locally
local_existing <- list.files(raw_dir, pattern = "\\.parquet$", full.names = FALSE)
gcs_basenames <- basename(gcs_shards)
missing_idx <- !(gcs_basenames %in% local_existing)
n_missing <- sum(missing_idx)
message(sprintf("  %d/%d shards already local; need to download %d",
                length(gcs_shards) - n_missing, length(gcs_shards), n_missing))

if (n_missing > 0) {
  # Write URI list to a temp file and use `gsutil -m cp -I` for parallel download
  uri_file <- tempfile(fileext = ".txt")
  writeLines(gcs_shards[missing_idx], uri_file)
  message("  gsutil -m cp downloading...")
  cmd_status <- system2("gsutil",
                         c("-q", "-m", "cp", "-I", raw_dir),
                         stdin = uri_file, stdout = "", stderr = "")
  if (cmd_status != 0) warning("gsutil cp returned non-zero status: ", cmd_status)
  unlink(uri_file)
}
message(sprintf("  Stage 1 done in %s",
                format(round(difftime(Sys.time(), t0, units = "mins"), 1))))

}  # end !V3_REUSE_BASIS gate for Stage 1

##########################################
### 2. LOAD + REGION/NEUROPIL          ###
##########################################

if (!is.na(V3_REUSE_BASIS)) {
  message("\n### Stage 2: skipped (reuse mode — synthesising batches from basis) ###")

  # Ensure basis outputs exist.
  if (!file.exists(basis_final_parquet)) {
    stop(sprintf(paste0("V3_REUSE_BASIS=%s requested but basis parquet not ",
                        "found at %s"), V3_REUSE_BASIS, basis_final_parquet))
  }

  # If the target version's processed_dir already has batch files (e.g. a
  # prior partial reuse run), skip synthesis.
  existing_reuse_batches <- list.files(processed_dir,
                                        pattern = "^batch_.*\\.parquet$",
                                        full.names = TRUE)
  if (length(existing_reuse_batches) > 0L) {
    message(sprintf("  target processed_dir already has %d batch files; reusing",
                    length(existing_reuse_batches)))
  } else {
    message(sprintf("  loading basis parquet: %s", basis_final_parquet))
    # Read as Arrow Table + integer64 to preserve uint64 sv_id precision in
    # newly-produced basis files. (Pre-fix basis files have already-corrupted
    # sv_ids stored as double; reading via integer64 doesn't restore them but
    # also doesn't make them worse.)
    basis_tbl <- arrow::read_parquet(basis_final_parquet, as_data_frame = FALSE)
    basis_df  <- as.data.frame(basis_tbl, integer64 = "integer64")
    rm(basis_tbl)
    # Drop stale root_id columns — Stage 3 re-populates these for target.version.
    drop_cols <- intersect(c("pre_root_id", "post_root_id", "pre_status",
                              "post_status"), names(basis_df))
    if (length(drop_cols) > 0L) {
      basis_df <- basis_df %>% dplyr::select(-dplyr::all_of(drop_cols))
    }
    message(sprintf("  basis rows: %s (dropped cols: %s)",
                    format(nrow(basis_df), big.mark = ","),
                    paste(drop_cols, collapse = ", ")))
    # Write as a single batch file that Stage 3 can consume.
    out_f <- file.path(processed_dir, "batch_000001.parquet")
    arrow::write_parquet(basis_df, out_f)
    message(sprintf("  wrote %s (%s rows)", out_f,
                    format(nrow(basis_df), big.mark = ",")))
    rm(basis_df)
    gc()
  }
} else {

message("\n### Stage 2: region + neuropil assignment ###")

# Open all downloaded shards as one virtual dataset
ds <- arrow::open_dataset(raw_dir, format = "parquet")
n_total_rows <- ds %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
message(sprintf("  Total rows in v3 raw dataset: %s",
                format(n_total_rows, big.mark = ",")))

# Columns we need downstream
keep_cols <- c("syn_id", "size",
               "centroid_x", "centroid_y", "centroid_z",
               "presyn_x", "presyn_y", "presyn_z",
               "postsyn_x", "postsyn_y", "postsyn_z",
               "presyn_sv_id", "postsyn_sv_id",
               "mean_score", "median_score")

# Test mode: collect a sample and process as a single batch
if (TEST_MODE) {
  message(sprintf("  TEST_MODE: sampling %d rows", TEST_N))
  set.seed(42)
  # Cast uint64 -> int64 in arrow before collect so sv_ids stay as
  # bit64::integer64 in R (see comment near rb_to_df below).
  uint_to_int64_cols <- intersect(
    c("syn_id", "size", "presyn_sv_id", "postsyn_sv_id"),
    keep_cols
  )
  v3_tbl <- ds %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    arrow::as_arrow_table()
  if (length(uint_to_int64_cols) > 0L) {
    sch <- v3_tbl$schema
    for (cn in uint_to_int64_cols) sch[[cn]] <- arrow::int64()
    v3_tbl <- v3_tbl$cast(sch)
  }
  v3 <- as.data.frame(v3_tbl, integer64 = "integer64")
  v3 <- v3 %>% dplyr::slice_sample(n = min(TEST_N, nrow(v3)))
  message(sprintf("  Sampled %d rows", nrow(v3)))
}

##########################################
### Alpha-shape surfaces (match v2)   ###
##########################################

alpha_value <- 50000  # nm — matches banc-calculate-neuropil-inclusion.R:45

message("  Building alpha shapes (one-time)...")
t_alpha <- Sys.time()

volume_surfs <- list(
  neck          = banc_neck_connective.surf,
  brain         = banc_brain_neuropil.surf,
  optic_lobes   = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf, "optic"))),
  sez           = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,
                                              "GNG|CAN|FLA|AMMC|SAD|PRW"))),
  central_brain = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf, "midbrain"))),
  vnc           = banc_vnc_neuropil.surf
)

volume_ashapes <- lapply(names(volume_surfs), function(nm) {
  message(sprintf("    %s...", nm))
  alphashape3d::ashape3d(nat::xyzmatrix(volume_surfs[[nm]]), alpha = alpha_value)
})
names(volume_ashapes) <- names(volume_surfs)

brain_nps <- sort(banc_brain_neuropils.surf$RegionList)
vnc_nps   <- sort(banc_vnc_neuropils.surf$RegionList)
np_region_map <- c(
  setNames(ifelse(grepl("^LO|^ME|^AME|^LOP", brain_nps), "optic_lobes",
           ifelse(grepl("^CAN|^GNG|^FLA|^AMMC|^SAD|^PRW", brain_nps),
                  "suboesophageal_zone", "central_brain")),
           brain_nps),
  setNames(rep("vnc", length(vnc_nps)), vnc_nps)
)

all_nps <- c(brain_nps, vnc_nps)
all_np_surfs <- c(
  lapply(brain_nps, function(np) subset(banc_brain_neuropils.surf, np)),
  lapply(vnc_nps,   function(np) subset(banc_vnc_neuropils.surf, np))
)
names(all_np_surfs) <- all_nps

np_ashapes <- vector("list", length(all_nps))
names(np_ashapes) <- all_nps
for (i in seq_along(all_nps)) {
  np_ashapes[[i]] <- alphashape3d::ashape3d(nat::xyzmatrix(all_np_surfs[[i]]),
                                              alpha = alpha_value)
}
message(sprintf("  Alpha shapes ready in %s (%d volumes, %d neuropils)",
                format(round(difftime(Sys.time(), t_alpha, units = "mins"), 1)),
                length(volume_ashapes), length(np_ashapes)))

##########################################
### Bounding-box pre-filter setup     ###
##########################################
###
### OPTIMISATION: compute the axis-aligned bounding box (AABB) of each
### alpha shape's vertices once. In classify_chunk, we first test which
### points fall inside the AABB (a cheap vectorised comparison) and only
### pass those to inashape3d(). For small neuropils this can reduce the
### number of point-in-surface tests from ~5M to a few thousand.

compute_ashape_bbox <- function(as3d) {
  # The alpha shape stores its input points in as3d$x (an Nx3 matrix)
  pts <- as3d$x
  list(
    xmin = min(pts[, 1]), xmax = max(pts[, 1]),
    ymin = min(pts[, 2]), ymax = max(pts[, 2]),
    zmin = min(pts[, 3]), zmax = max(pts[, 3])
  )
}

message("  Pre-computing bounding boxes for alpha shapes...")

volume_bboxes <- lapply(volume_ashapes, compute_ashape_bbox)
np_bboxes     <- lapply(np_ashapes, compute_ashape_bbox)

message(sprintf("  Bounding boxes ready (%d volumes, %d neuropils)",
                length(volume_bboxes), length(np_bboxes)))

##########################################
### Chunk classifier (optimised)      ###
##########################################
###
### OPTIMISATION vs original: before calling inashape3d() for each
### neuropil/volume, we check which points lie within the alpha shape's
### axis-aligned bounding box (AABB). Only those candidate points are
### passed to inashape3d(), which is the expensive O(N * faces) operation.
### Since most neuropils occupy a small fraction of the full CNS volume,
### this typically eliminates 90-99% of points before the expensive test.
### The result is identical to the original -- the AABB is a superset of
### the alpha shape, so no true positives are lost.

classify_chunk <- function(chunk_df) {
  chunk_df$neuropil <- ""
  chunk_df$region   <- ""
  chunk_df$side     <- ""

  # v3 centroid is in 16x16x45 nm voxel units -- convert to nm
  points <- cbind(
    X = chunk_df$centroid_x * V3_VOXEL_NM[1],
    Y = chunk_df$centroid_y * V3_VOXEL_NM[2],
    Z = chunk_df$centroid_z * V3_VOXEL_NM[3]
  )

  n_pts <- nrow(points)
  px <- points[, "X"]
  py <- points[, "Y"]
  pz <- points[, "Z"]

  # Side
  lrdiffs <- bancr:::banc_lr_position(points, units = "nm")
  chunk_df$side <- ifelse(lrdiffs > 0, "right", "left")

  # Helper: test points against an alpha shape with bounding-box pre-filter.
  # Returns a logical vector of length n_pts (TRUE = inside the shape).
  test_with_bbox <- function(as3d, bbox) {
    # Step 1: cheap AABB filter (vectorised, ~instant)
    in_box <- (px >= bbox$xmin & px <= bbox$xmax &
               py >= bbox$ymin & py <= bbox$ymax &
               pz >= bbox$zmin & pz <= bbox$zmax)
    candidates <- which(in_box)

    if (length(candidates) == 0L) return(rep(FALSE, n_pts))

    # Step 2: expensive alpha-shape test only on candidate points
    inside_candidates <- alphashape3d::inashape3d(
      points     = points[candidates, , drop = FALSE],
      as3d       = as3d,
      indexAlpha  = "ALL"
    )

    result <- rep(FALSE, n_pts)
    result[candidates] <- inside_candidates
    result
  }

  # Regional volumes
  for (vol_name in names(volume_ashapes)) {
    inside <- test_with_bbox(volume_ashapes[[vol_name]], volume_bboxes[[vol_name]])
    if (any(inside)) chunk_df$region[which(inside)] <- vol_name
  }

  # Individual neuropils (one synapse can be in multiple -> comma-joined)
  for (np_name in names(np_ashapes)) {
    inside <- test_with_bbox(np_ashapes[[np_name]], np_bboxes[[np_name]])
    if (any(inside)) {
      idx <- which(inside)
      chunk_df$neuropil[idx] <- vapply(chunk_df$neuropil[idx], function(x) {
        paste(unique(c(strsplit(x, ",", fixed = TRUE)[[1]], np_name)), collapse = ",")
      }, character(1))
      unassigned <- idx[chunk_df$region[idx] == ""]
      if (length(unassigned)) chunk_df$region[unassigned] <- np_region_map[np_name]
    }
  }

  chunk_df$neuropil[chunk_df$neuropil == ""] <- "outside"
  chunk_df$region[chunk_df$region == ""]   <- "outside"
  chunk_df$neuropil <- gsub("^,", "", chunk_df$neuropil)
  chunk_df$region   <- gsub("^,", "", chunk_df$region)

  # Cast X/Y/Z (nm) into a tidy location
  chunk_df$X <- points[, "X"]
  chunk_df$Y <- points[, "Y"]
  chunk_df$Z <- points[, "Z"]
  chunk_df
}

##########################################
### Parallel cluster (Stage 2 only)   ###
##########################################
###
### foreach + doSNOW: workers each classify a disjoint slice of one
### macro-batch and the master concatenates the results into a single
### batch_NNNNNN.parquet — workers do not write files. This is
### concurrency-safe *within* a single SLURM job. It does NOT make the
### script safe to run as multiple parallel SLURM jobs: the batch
### numbering race in stage 2 and the svid-cache race in stage 3 still
### apply across processes.

cl <- NULL
if (!TEST_MODE && numCores > 1L) {
  message(sprintf("  Setting up parallel backend (%d workers) for stage 2...",
                  numCores))
  cl <- parallel::makeCluster(numCores)
  doSNOW::registerDoSNOW(cl)
  parallel::clusterExport(
    cl,
    varlist = c("volume_ashapes", "np_ashapes", "np_region_map",
                "volume_bboxes", "np_bboxes",
                "V3_VOXEL_NM", "classify_chunk"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(alphashape3d)
      library(bancr)
      library(nat)
    })
  })
}

classify_chunk_par <- function(df) {
  if (is.null(cl) || nrow(df) < 2L * numCores) {
    return(classify_chunk(df))
  }
  sub_idx <- split(seq_len(nrow(df)),
                   cut(seq_len(nrow(df)), numCores, labels = FALSE))
  parts <- foreach::foreach(i = seq_along(sub_idx),
                            .packages = c("alphashape3d", "bancr", "nat"),
                            .errorhandling = "stop") %dopar% {
    classify_chunk(df[sub_idx[[i]], , drop = FALSE])
  }
  dplyr::bind_rows(parts)
}

##########################################
### Batched processing w/ resume      ###
##########################################

# Inventory already-processed batches so we can skip them on restart
existing_batches <- list.files(processed_dir, pattern = "^batch_.*\\.parquet$",
                                full.names = TRUE)
batch_index <- length(existing_batches) + 1L
message(sprintf("  %d batch files already processed — next batch index = %d",
                length(existing_batches), batch_index))

process_and_save_batch <- function(df, idx) {
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df %>% dplyr::filter(size >= MIN_SIZE)
  if (nrow(df) == 0) return(invisible(NULL))

  df <- classify_chunk_par(df)

  out_file <- file.path(processed_dir, sprintf("batch_%06d.parquet", idx))
  arrow::write_parquet(df, out_file)
  message(sprintf("    saved %s (%s rows)",
                  basename(out_file), format(nrow(df), big.mark = ",")))
  invisible(out_file)
}

if (TEST_MODE) {
  message("  TEST_MODE: single batch")
  process_and_save_batch(v3, batch_index)
} else {
  # Iterate over the dataset in chunks. Scanner$create() works across arrow
  # versions; the older NewScan()/UseAsync() builder was removed in arrow >= 17.
  scanner <- arrow::Scanner$create(
    ds,
    projection = keep_cols,
    batch_size = 1e5L
  )
  reader <- scanner$ToRecordBatchReader()
  accumulator <- list()
  acc_rows <- 0L
  batch_count <- 0L
  rows_seen <- 0L

  # Seek past already-processed rows (simpler: if N batches exist, skip first N*BATCH_ROWS rows)
  rows_to_skip <- (batch_index - 1L) * BATCH_ROWS
  if (rows_to_skip > 0L) {
    message(sprintf("  Skipping ~%s rows already processed",
                    format(rows_to_skip, big.mark = ",")))
  }

  # Preserve uint64 precision when crossing into R. Raw shards store syn_id,
  # size, presyn_sv_id and postsyn_sv_id as uint64; sv_ids in BANC are ~7e16,
  # well above 2^53 (~9e15), so a plain as.data.frame(rb) silently rounds
  # them to the nearest representable double. Cast uint64 -> int64 in arrow
  # first, then convert to bit64::integer64 in R.
  uint_to_int64_cols <- intersect(
    c("syn_id", "size", "presyn_sv_id", "postsyn_sv_id"),
    keep_cols
  )
  rb_to_df <- function(rb) {
    tbl <- arrow::as_arrow_table(rb)
    if (length(uint_to_int64_cols) > 0L) {
      sch <- tbl$schema
      for (cn in uint_to_int64_cols) sch[[cn]] <- arrow::int64()
      tbl <- tbl$cast(sch)
    }
    as.data.frame(tbl, integer64 = "integer64")
  }

  repeat {
    rb <- reader$read_next_batch()
    if (is.null(rb)) break
    df_piece <- rb_to_df(rb)

    # Skip already-processed rows
    if (rows_seen + nrow(df_piece) <= rows_to_skip) {
      rows_seen <- rows_seen + nrow(df_piece)
      next
    }
    if (rows_seen < rows_to_skip) {
      drop_n <- rows_to_skip - rows_seen
      df_piece <- df_piece[(drop_n + 1L):nrow(df_piece), , drop = FALSE]
      rows_seen <- rows_to_skip
    }
    rows_seen <- rows_seen + nrow(df_piece)

    accumulator[[length(accumulator) + 1L]] <- df_piece
    acc_rows <- acc_rows + nrow(df_piece)

    if (acc_rows >= BATCH_ROWS) {
      df_batch <- dplyr::bind_rows(accumulator)
      process_and_save_batch(df_batch, batch_index)
      batch_index <- batch_index + 1L
      batch_count <- batch_count + 1L
      accumulator <- list()
      acc_rows <- 0L
      message(sprintf("  progress: %s / %s rows scanned",
                      format(rows_seen, big.mark = ","),
                      format(n_total_rows, big.mark = ",")))
    }
  }
  # Flush the final partial batch
  if (acc_rows > 0L) {
    df_batch <- dplyr::bind_rows(accumulator)
    process_and_save_batch(df_batch, batch_index)
    batch_index <- batch_index + 1L
    batch_count <- batch_count + 1L
  }
  message(sprintf("  Stage 2 done: %d new batches processed", batch_count))
}

}  # end !V3_REUSE_BASIS gate for Stage 2

# Tear down the stage-2 worker pool before stage 3 starts.
if (exists("cl") && !is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- NULL
}

##########################################
### 3. SUPERVOXEL -> ROOT_ID via CAVE  ###
##########################################

message("\n### Stage 3: supervoxel -> root_id ###")

# Load banc.meta for local cache
banc.meta <- tryCatch(banctable_query(), error = function(e) {
  warning("banctable_query failed, falling back to CSV: ", e$message)
  readr::read_csv(file.path(banc.meta.save.path, "banc_meta.csv"),
                  col_types = banc.col.types, show_col_types = FALSE)
}) %>% dplyr::filter(!is.na(root_id))

# Load persistent svid -> root_id cache if it exists
svid_cache <- if (file.exists(svid_cache_file)) {
  arrow::read_parquet(svid_cache_file) %>%
    dplyr::mutate(supervoxel_id = as.character(supervoxel_id),
                  root_id = as.character(root_id))
} else {
  data.frame(supervoxel_id = character(0), root_id = character(0),
             stringsAsFactors = FALSE)
}
message(sprintf("  loaded svid cache: %s entries",
                format(nrow(svid_cache), big.mark = ",")))

# Collect all unique svids across processed batches
batch_files <- sort(list.files(processed_dir, pattern = "^batch_.*\\.parquet$",
                                full.names = TRUE))
message(sprintf("  scanning %d processed batch files for unique svids...",
                length(batch_files)))
read_sv_cols <- function(bf, cols) {
  # arrow::read_parquet() drops int64 -> R double, which silently rounds
  # sv_ids (~7e16, > 2^53). Read as Arrow Table and pull cols as integer64.
  tbl <- arrow::read_parquet(bf, col_select = cols, as_data_frame = FALSE)
  as.data.frame(tbl, integer64 = "integer64")
}
all_svids <- character(0)
for (bf in batch_files) {
  df <- read_sv_cols(bf, c("presyn_sv_id", "postsyn_sv_id"))
  all_svids <- unique(c(all_svids,
                        as.character(df$presyn_sv_id),
                        as.character(df$postsyn_sv_id)))
}
all_svids <- setdiff(all_svids, c("0", "", NA_character_))
message(sprintf("  %s unique svids", format(length(all_svids), big.mark = ",")))

# Step 1: fill from banc.meta (local, instant)
meta_map <- setNames(as.character(banc.meta$root_id),
                     as.character(banc.meta$supervoxel_id))
svids_from_meta <- all_svids[all_svids %in% names(meta_map)]
cache_from_meta <- data.frame(
  supervoxel_id = svids_from_meta,
  root_id       = unname(meta_map[svids_from_meta]),
  stringsAsFactors = FALSE
)
message(sprintf("  found %s svids in banc.meta",
                format(nrow(cache_from_meta), big.mark = ",")))

# Step 2: fill from existing cache
already_cached <- svid_cache$supervoxel_id
need_cave <- setdiff(setdiff(all_svids, svids_from_meta), already_cached)
message(sprintf("  %s svids already in persistent cache",
                format(length(intersect(all_svids, already_cached)), big.mark = ",")))
message(sprintf("  %s svids to query via CAVE",
                format(length(need_cave), big.mark = ",")))

# Step 3: CAVE fallback in chunks, appending to persistent cache
if (length(need_cave) > 0L) {
  n_cave <- length(need_cave)
  n_chunks <- ceiling(n_cave / CAVE_CHUNKSIZE)
  message(sprintf("  querying CAVE in %d chunks of %s",
                  n_chunks, format(CAVE_CHUNKSIZE, big.mark = ",")))
  for (ci in seq_len(n_chunks)) {
    from <- (ci - 1L) * CAVE_CHUNKSIZE + 1L
    to   <- min(ci * CAVE_CHUNKSIZE, n_cave)
    chunk_svids <- need_cave[from:to]
    res <- tryCatch(
      banc_rootid(chunk_svids, version = target.version),
      error = function(e) {
        warning(sprintf("CAVE chunk %d/%d failed: %s", ci, n_chunks, e$message))
        rep(NA_character_, length(chunk_svids))
      }
    )
    svid_cache <- dplyr::bind_rows(svid_cache, data.frame(
      supervoxel_id = chunk_svids,
      root_id       = as.character(res),
      stringsAsFactors = FALSE
    ))
    # Persist cache after each chunk so we can resume if killed
    arrow::write_parquet(svid_cache, svid_cache_file)
    message(sprintf("    CAVE chunk %d/%d done (%s svids)",
                    ci, n_chunks, format(length(chunk_svids), big.mark = ",")))
  }
}

# Build the full lookup table combining banc.meta + persistent cache
full_lookup <- dplyr::bind_rows(cache_from_meta, svid_cache) %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE)
svid_to_root <- setNames(full_lookup$root_id, full_lookup$supervoxel_id)
message(sprintf("  full svid -> root lookup: %s entries",
                format(length(svid_to_root), big.mark = ",")))

##########################################
### Re-write each batch with root_ids  ###
##########################################

message("\n  joining root_ids into processed batches...")
for (bf in batch_files) {
  # Read as Arrow Table + integer64 to preserve uint64 sv_id precision
  # (see read_sv_cols comment above).
  tbl <- arrow::read_parquet(bf, as_data_frame = FALSE)
  df  <- as.data.frame(tbl, integer64 = "integer64")
  if (all(c("pre_root_id", "post_root_id") %in% names(df))) next  # already done
  df$pre_root_id  <- unname(svid_to_root[as.character(df$presyn_sv_id)])
  df$post_root_id <- unname(svid_to_root[as.character(df$postsyn_sv_id)])
  arrow::write_parquet(df, bf)
}
message("  done.")

##########################################
### 4. MERGE + CAPTURE RATES          ###
##########################################

message("\n### Stage 4: merge + capture rates ###")

# Two prior failures here:
#   1. arrow::open_dataset(processed_dir) — SLURM-only "Invalid: Unrecognized
#      filesystem type in URI: file:///_"
#   2. dplyr::bind_rows(lapply(read_parquet)) — Floating Point Exception during
#      bind_rows on the int64 syn_id column across 52 × ~5M-row tibbles
# Stay in arrow C++ space until write: read each batch as a lazy Arrow Table
# (as_data_frame=FALSE), concat in arrow, write final_parquet IMMEDIATELY so
# downstream consumers (v888 reuse, completion) have the basis even if the
# capture-rate sections trip.
batch_files_stage4 <- sort(list.files(
  processed_dir, pattern = "^batch_.*\\.parquet$", full.names = TRUE
))
message(sprintf("  reading %d batch files as Arrow tables...",
                length(batch_files_stage4)))
batch_tables <- lapply(batch_files_stage4, function(f) {
  arrow::read_parquet(f, as_data_frame = FALSE)
})
combined_arrow <- do.call(arrow::concat_tables, batch_tables)
rm(batch_tables); gc(verbose = FALSE)
message(sprintf("  concat: %s rows in arrow Table",
                format(combined_arrow$num_rows, big.mark = ",")))

# Write the consolidated parquet first — this is the artifact the v888 reuse
# run needs as its basis. Capture-rate failures downstream don't block J2.
# (neuropil_lookup write was previously here but the file has no actual
# consumers — Stage 2 reuse only reads basis_final_parquet — so it's dropped.)
arrow::write_parquet(combined_arrow, final_parquet)
message(sprintf("  saved %s", basename(final_parquet)))

# Now collect to R for capture-rate stats. Dedup by syn_id (defensive — covers
# Stage-2 resume overlap and the cleaned batch_000003 boundary).
combined <- dplyr::collect(combined_arrow)
rm(combined_arrow); gc(verbose = FALSE)
message(sprintf("  collected %s rows to R",
                format(nrow(combined), big.mark = ",")))
n_pre_dedup <- nrow(combined)
combined <- combined %>% dplyr::distinct(syn_id, .keep_all = TRUE)
n_removed <- n_pre_dedup - nrow(combined)
if (n_removed > 0L) {
  message(sprintf("  de-duplicated: removed %s rows (%.3f%% of merged)",
                  format(n_removed, big.mark = ","),
                  100 * n_removed / n_pre_dedup))
}

# Classify pre/post as neuron vs fragment
neuron.ids <- unique(as.character(banc.meta$root_id))
combined <- combined %>%
  dplyr::mutate(
    pre_status  = ifelse(pre_root_id  %in% neuron.ids, "neuron", "fragment"),
    post_status = ifelse(post_root_id %in% neuron.ids, "neuron", "fragment")
  )

write_rate <- function(df, name) {
  f <- file.path(rates_dir,
                 sprintf("banc_%s_v3_%s_capture_rates.csv", target.version, name))
  readr::write_csv(df, f)
  message(sprintf("  saved %s (%d rows)", basename(f), nrow(df)))
  f
}

# 1. Gross rates
summary_gross <- combined %>%
  dplyr::count(pre_status, post_status) %>%
  dplyr::mutate(prop = round(n / sum(n), 4))
write_rate(summary_gross, "gross")

# 2. Inside/outside
summary_inout <- combined %>%
  dplyr::mutate(in_mesh = ifelse(region == "outside", "outside", "inside")) %>%
  dplyr::group_by(in_mesh) %>%
  dplyr::count(pre_status, post_status, in_mesh) %>%
  dplyr::mutate(prop = round(n / sum(n), 4)) %>%
  dplyr::ungroup()
write_rate(summary_inout, "inout")

# 3. By region x side
summary_region <- combined %>%
  dplyr::group_by(region, side) %>%
  dplyr::count(pre_status, post_status, side, region) %>%
  dplyr::mutate(prop = round(n / sum(n), 4)) %>%
  dplyr::ungroup()
write_rate(summary_region, "region")

# 4. By neuropil
summary_neuropil <- combined %>%
  dplyr::group_by(region, side, neuropil) %>%
  dplyr::count(pre_status, post_status, side, region, neuropil) %>%
  dplyr::mutate(prop = round(n / sum(n), 4)) %>%
  dplyr::ungroup()
write_rate(summary_neuropil, "neuropil")

# Overall summary to console
n_total       <- nrow(combined)
n_pre_neuron  <- sum(combined$pre_status  == "neuron")
n_post_neuron <- sum(combined$post_status == "neuron")
n_both_neuron <- sum(combined$pre_status == "neuron" & combined$post_status == "neuron")
n_outside     <- sum(combined$region == "outside")
message(sprintf("\n  === v3 overall ===  total=%s", format(n_total, big.mark = ",")))
message(sprintf("  pre neuron:  %s (%.1f%%)",
                format(n_pre_neuron,  big.mark = ","), 100 * n_pre_neuron  / n_total))
message(sprintf("  post neuron: %s (%.1f%%)",
                format(n_post_neuron, big.mark = ","), 100 * n_post_neuron / n_total))
message(sprintf("  both:        %s (%.1f%%)",
                format(n_both_neuron, big.mark = ","), 100 * n_both_neuron / n_total))
message(sprintf("  outside:     %s (%.1f%%)",
                format(n_outside,     big.mark = ","), 100 * n_outside     / n_total))

##########################################
### 5. v2 vs v3 COMPARISON            ###
##########################################

message("\n### Stage 5: v2 vs v3 comparison ###")

v2_region_file <- file.path(banc.connectivity.save.path,
                             sprintf("banc_%s_region_capture_rates.csv", target.version))
v2_neuropil_file <- file.path(banc.connectivity.save.path,
                               sprintf("banc_%s_neuropil_capture_rates.csv", target.version))

compare_rates <- function(v2_file, v3_df, by_cols, out_name) {
  if (!file.exists(v2_file)) {
    message(sprintf("  v2 file not found, skipping: %s", v2_file))
    return(invisible(NULL))
  }
  v2 <- readr::read_csv(v2_file, show_col_types = FALSE)
  cmp <- dplyr::full_join(
    v2    %>% dplyr::rename(n_v2 = n, prop_v2 = prop),
    v3_df %>% dplyr::rename(n_v3 = n, prop_v3 = prop),
    by = by_cols
  ) %>%
    dplyr::mutate(
      delta_n    = n_v3 - n_v2,
      delta_prop = prop_v3 - prop_v2
    )
  out <- file.path(rates_dir,
                   sprintf("banc_%s_v2_vs_v3_%s.csv", target.version, out_name))
  readr::write_csv(cmp, out)
  message(sprintf("  saved %s (%d rows)", basename(out), nrow(cmp)))
}

compare_rates(v2_region_file, summary_region,
              c("region", "side", "pre_status", "post_status"),
              "region_comparison")
compare_rates(v2_neuropil_file, summary_neuropil,
              c("region", "side", "neuropil", "pre_status", "post_status"),
              "neuropil_comparison")

##########################################
### 6. PUSH TO GCS                     ###
##########################################

if (PUSH_TO_GCS && !TEST_MODE) {
  message("\n### Stage 6: push to GCS ###")
  to_push <- c(
    final_parquet,
    neuropil_lookup,
    list.files(rates_dir, pattern = "\\.csv$", full.names = TRUE)
  )
  for (f in to_push) {
    dest <- paste0(V3_GCS_DST, "/", basename(f))
    message(sprintf("  gsutil cp %s -> %s", basename(f), dest))
    system2("gsutil", c("-q", "cp", f, dest), stdout = "", stderr = "")
  }
  message("  GCS push complete.")
} else if (TEST_MODE) {
  message("\n### Stage 6: skipping GCS push in TEST_MODE ###")
}

message("\n### banc-synapses-v3 done ###")
