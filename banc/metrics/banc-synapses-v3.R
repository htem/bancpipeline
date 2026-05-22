#' banc-synapses-v3 — Process v3 synapse predictions and compare capture vs v2.
#'
#' Resumable: download v3 shards → region/neuropil via alpha-shape →
#' svid→root_id (cache + CAVE) → capture rates → v2↔v3 comparison.
#'
#' @section Reads:
#'   - GCS v3 raw synapse parquet shards
#'   - SeaTable `banc_meta`, CAVE (svid fallbacks)
#'
#' @section Writes:
#'   - per-batch parquet caches, capture-rate CSVs, v2↔v3 comparison CSV
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_synapses_v3.sh`.
#'
#' @section Notes:
#'   - See `banc-synapses-v3-optimised.R` for the bbox-prefilter variant.
#'   - Test mode: `TEST_MODE = TRUE`.

###########################################################
### banc-synapses-v3.R
###
### Process the v3 synapse predictions (260326_assignment)
### and compare capture rates against v2.
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
V3_GCS_DST     <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/v%s/%s",
                           banc.version, V3_NAME)

# v3 centroid units: 16x16x45 nm per voxel (from info file "resolution")
V3_VOXEL_NM    <- c(16, 16, 45)

BATCH_ROWS     <- 5e6L      # rows per region/neuropil classification batch
CAVE_CHUNKSIZE <- 1e5L      # svids per CAVE call (fafbseg default)
MIN_SIZE       <- 5L        # size threshold, matches v2

PUSH_TO_GCS    <- TRUE

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
processed_dir <- file.path(banc.synapses.v3.save.path, "processed")
cache_dir     <- file.path(banc.synapses.v3.save.path, "cache")
rates_dir     <- file.path(banc.synapses.v3.save.path, "capture_rates")
for (d in c(banc.synapses.v3.save.path, raw_dir, processed_dir, cache_dir, rates_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

final_parquet      <- file.path(banc.synapses.v3.save.path,
                                 sprintf("banc_%s_synapses_v3.parquet", banc.version))
neuropil_lookup    <- file.path(banc.synapses.v3.save.path,
                                 sprintf("banc_%s_synapses_v3_neuropil_lookup.parquet",
                                         banc.version))
svid_cache_file    <- file.path(cache_dir, "svid_rootid_cache.parquet")

##########################################
### 0. MODE BANNER                     ###
##########################################

message(sprintf("\n### banc-synapses-v3 | TEST_MODE=%s | version=%s ###",
                TEST_MODE, banc.version))
message(sprintf("  Working directory: %s", banc.synapses.v3.save.path))

##########################################
### 1. DOWNLOAD RAW SHARDS FROM GCS   ###
##########################################

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

##########################################
### 2. LOAD + REGION/NEUROPIL          ###
##########################################

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
  v3 <- ds %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    dplyr::collect()
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
### Chunk classifier (matches v2)     ###
##########################################

classify_chunk <- function(chunk_df) {
  chunk_df$neuropil <- ""
  chunk_df$region   <- ""
  chunk_df$side     <- ""

  # v3 centroid is in 16x16x45 nm voxel units — convert to nm
  points <- cbind(
    X = chunk_df$centroid_x * V3_VOXEL_NM[1],
    Y = chunk_df$centroid_y * V3_VOXEL_NM[2],
    Z = chunk_df$centroid_z * V3_VOXEL_NM[3]
  )

  # Side
  lrdiffs <- bancr:::banc_lr_position(points, units = "nm")
  chunk_df$side <- ifelse(lrdiffs > 0, "right", "left")

  # Regional volumes
  for (vol_name in names(volume_ashapes)) {
    inside <- alphashape3d::inashape3d(points = points,
                                        as3d = volume_ashapes[[vol_name]],
                                        indexAlpha = "ALL")
    if (any(inside)) chunk_df$region[which(inside)] <- vol_name
  }

  # Individual neuropils (one synapse can be in multiple -> comma-joined)
  for (np_name in names(np_ashapes)) {
    inside <- alphashape3d::inashape3d(points = points,
                                        as3d = np_ashapes[[np_name]],
                                        indexAlpha = "ALL")
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

  repeat {
    rb <- reader$read_next_batch()
    if (is.null(rb)) break
    df_piece <- as.data.frame(rb)

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

# Tear down the stage-2 worker pool before stage 3 starts.
if (!is.null(cl)) {
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
all_svids <- character(0)
for (bf in batch_files) {
  df <- arrow::read_parquet(bf, col_select = c("presyn_sv_id", "postsyn_sv_id"))
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
      banc_rootid(chunk_svids, version = banc.version),
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
  df <- arrow::read_parquet(bf)
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

combined <- arrow::open_dataset(processed_dir, format = "parquet") %>%
  dplyr::collect()
message(sprintf("  merged %s rows", format(nrow(combined), big.mark = ",")))

# Write the unified parquet + lookup
arrow::write_parquet(combined, final_parquet)
arrow::write_parquet(combined %>%
                       dplyr::select(syn_id, X, Y, Z, region, neuropil, side,
                                     pre_root_id, post_root_id),
                     neuropil_lookup)
message(sprintf("  saved %s", basename(final_parquet)))
message(sprintf("  saved %s", basename(neuropil_lookup)))

# Classify pre/post as neuron vs fragment
neuron.ids <- unique(as.character(banc.meta$root_id))
combined <- combined %>%
  dplyr::mutate(
    pre_status  = ifelse(pre_root_id  %in% neuron.ids, "neuron", "fragment"),
    post_status = ifelse(post_root_id %in% neuron.ids, "neuron", "fragment")
  )

write_rate <- function(df, name) {
  f <- file.path(rates_dir,
                 sprintf("banc_%s_v3_%s_capture_rates.csv", banc.version, name))
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
                             sprintf("banc_%s_region_capture_rates.csv", banc.version))
v2_neuropil_file <- file.path(banc.connectivity.save.path,
                               sprintf("banc_%s_neuropil_capture_rates.csv", banc.version))

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
                   sprintf("banc_%s_v2_vs_v3_%s.csv", banc.version, out_name))
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
