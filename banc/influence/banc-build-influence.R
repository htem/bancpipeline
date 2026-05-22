#' banc-build-influence — Compute all-to-all BANC influence scores (sharded).
#'
#' Uses the `influencer` R package (Python PETSc/SLEPc backend) for efficient
#' sparse-matrix solves. Outputs partitioned parquet files; resumable —
#' skips already-written chunks. Supports SLURM-array sharding via
#' `BANC_INFLUENCE_SHARD_{IDX,TOTAL}` env vars.
#'
#' @section Reads:
#'   - `<banc.versioned.save.path>/banc_<ver>_edgelist_simple_<src>.feather`
#'   - `<banc.versioned.save.path>/banc_<ver>_meta.feather`
#'   - env vars `BANC_INFLUENCE_SHARD_IDX`, `BANC_INFLUENCE_SHARD_TOTAL`
#'
#' @section Writes:
#'   - `<banc.versioned.save.path>/influence/all_to_all/chunk_NNNN.parquet`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_influence_array.sh`,
#'   `o2/production/o2_banc_influence.sh`
#'
#' @section Schema:
#'   BANC-project/manuscript/print/dataverse/documentation/influence_all_to_all.md
#'
#' @section Paper:
#'   Methods §"Influence".
#'
#' @section Notes:
#'   - GCS sync is handled separately by `banc-sync-influence.R`.

###############################################
### Calculate all-to-all BANC influence scores
### using the influencer package (Python backend).
###
### Uses Python InfluenceCalculator with PETSc/SLEPc
### for efficient sparse matrix solves.
###
### Output: partitioned parquet files in
###   {versioned save.path}/influence/all_to_all/
### Columns: upstream_id, downstream_id, raw_influence
###
### Resumable: skips already-written chunks.
###############################################
source("banc/banc-startup.R")

local({

message("### banc: calculating all-to-all influence scores ###")

library(influencer)

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R

# Versioned output directory (matches banc-data.R save.path pattern)
save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                        paste0("banc_", banc.version))
dir.create(save.path, recursive = TRUE, showWarnings = FALSE)

# Output directory for partitioned parquet files
output_dir <- file.path(save.path, "influence", "all_to_all")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# To delete corrupt/stale output and rebuild from scratch, uncomment:
# unlink(list.files(output_dir, pattern = "\\.parquet$", full.names = TRUE))
# unlink(file.path(save.path, sprintf("banc_%s_proofread_influence.sqlite", banc.version)))

# Overwrite existing chunks? Set TRUE to recompute all (e.g., after bugfix).
# When FALSE (default), already-written chunks are skipped (resumable).
overwrite_chunks <- FALSE

# Batch size: number of seed neurons per parquet chunk
# Each chunk = batch_size * n_neurons rows (~500 * 155K = ~77M rows)
batch_size <- 500

# ---- Sharding (parallel SLURM jobs) ----
#
# To run multiple jobs concurrently, set:
#   BANC_INFLUENCE_SHARD_TOTAL = N   (total number of shards)
#   BANC_INFLUENCE_SHARD_IDX   = i   (this job's shard index, 0..N-1)
#
# Each shard owns batches where (b - 1) %% N == i, so shards write disjoint
# chunk_NNNN.parquet files — no write races. Shard 0 builds the SQLite if it
# doesn't exist; other shards wait for it. The GCS sync at the end is only
# performed for un-sharded runs (TOTAL == 1) — for sharded runs, run
# banc-sync-influence.R after all shards finish.
shard_idx   <- as.integer(Sys.getenv("BANC_INFLUENCE_SHARD_IDX",   "0"))
shard_total <- as.integer(Sys.getenv("BANC_INFLUENCE_SHARD_TOTAL", "1"))
if (is.na(shard_idx) || is.na(shard_total) || shard_total < 1L ||
    shard_idx < 0L || shard_idx >= shard_total) {
  stop(sprintf("Invalid shard env: BANC_INFLUENCE_SHARD_IDX=%s, BANC_INFLUENCE_SHARD_TOTAL=%s",
               Sys.getenv("BANC_INFLUENCE_SHARD_IDX", "<unset>"),
               Sys.getenv("BANC_INFLUENCE_SHARD_TOTAL", "<unset>")))
}
if (shard_total > 1L) {
  message(sprintf("  SHARDED RUN: shard %d of %d (will process batches where (b-1) %%%% %d == %d)",
                  shard_idx, shard_total, shard_total, shard_idx))
}

# Synapse count threshold — only include pre->post edges with count >= this.
# (1 = un-thresholded, include all connections.)
count_thresh <- 5

# Minimum raw influence score to store (NULL = store all pairs)
# Set to exp(-24) to only store meaningful scores and reduce storage
# NULL keeps all ~24 billion pairs (~576 GB uncompressed, ~100-150 GB parquet)
min_influence <- NULL

# Adjusted influence constant
const <- 24

#################################
### LOAD PROOFREAD NEURON IDs ###
#################################

# Only include neurons that are proofread or roughly proofread
message("Loading proofread neuron IDs...")

# Backbone proofread from CAVE
proofread_ids <- tryCatch({
  pr <- banc_backbone_proofread()
  as.character(pr$pt_root_id)
}, error = function(e) {
  message("  Could not fetch backbone_proofread from CAVE: ", e$message)
  character(0)
})

# Roughly proofread from SeaTable
roughly_proofread_ids <- tryCatch({
  rp <- banctable_query_cached()
  as.character(rp$root_id)
}, error = function(e) {
  message("  Could not fetch roughly_proofread from SeaTable: ", e$message)
  character(0)
})

proofread_all <- unique(c(proofread_ids, roughly_proofread_ids))
message(sprintf("  Proofread neurons: %s (backbone: %s, roughly: %s)",
                format(length(proofread_all), big.mark = ","),
                format(length(proofread_ids), big.mark = ","),
                format(length(roughly_proofread_ids), big.mark = ",")))

#######################
### BUILD SQLITE    ###
#######################

# Persistent SQLite database for the Python InfluenceCalculator
# Created from the feather edgelist, filtered to proofread neurons only.
# Delete the existing SQLite to rebuild with updated proofread set.
sqlite_path <- file.path(save.path,
                         sprintf("banc_%s_proofread_influence.sqlite", banc.version))

if (!file.exists(sqlite_path) && shard_idx > 0L) {
  # Non-zero shards must NOT race to build the SQLite. Wait for shard 0
  # (or a previous un-sharded run) to create it. Poll every 30s for up
  # to 4 hours, then give up.
  message(sprintf("  Shard %d waiting for SQLite to appear at %s",
                  shard_idx, sqlite_path))
  wait_secs <- 0L
  max_wait  <- 4L * 60L * 60L
  while (!file.exists(sqlite_path) && wait_secs < max_wait) {
    Sys.sleep(30)
    wait_secs <- wait_secs + 30L
    if (wait_secs %% 300L == 0L) {
      message(sprintf("    still waiting (%d min elapsed)...", wait_secs %/% 60L))
    }
  }
  if (!file.exists(sqlite_path)) {
    stop(sprintf("Shard %d gave up waiting for SQLite after %d min — start shard 0 first.",
                 shard_idx, max_wait %/% 60L))
  }
  # Small delay so the file is fully fsynced before we open it
  Sys.sleep(5)
}

if (!file.exists(sqlite_path)) {
  # Default synapse source = v3. Override with BANC_SYN_SOURCE env var.
  syn_source <- Sys.getenv("BANC_SYN_SOURCE", unset = "")
  if (!nzchar(syn_source)) {
    syn_source <- if (exists("banc.synapse.source.default"))
      banc.synapse.source.default else "v3"
  }
  stopifnot(syn_source %in% c("v2", "v3"))
  # Read versioned edgelist produced by banc-calculate-connectivity.R
  edgelist_path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                             sprintf("banc_%s", banc.version),
                             sprintf("banc_%s_edgelist_simple_%s.feather",
                                     banc.version, syn_source))
  stopifnot(file.exists(edgelist_path))
  message("Reading edgelist from: ", edgelist_path)
  edgelist <- arrow::read_feather(edgelist_path)

  # Ensure IDs are character (preserves 64-bit integer precision)
  edgelist$pre <- as.character(edgelist$pre)
  edgelist$post <- as.character(edgelist$post)

  message(sprintf("Full edgelist: %s connections between %s neurons",
                  format(nrow(edgelist), big.mark = ","),
                  format(length(unique(c(edgelist$pre, edgelist$post))), big.mark = ",")))

  # Filter to proofread neurons only (both pre and post must be proofread)
  edgelist <- edgelist[edgelist$pre %in% proofread_all & edgelist$post %in% proofread_all, ]
  message(sprintf("Proofread edgelist: %s connections between %s neurons",
                  format(nrow(edgelist), big.mark = ","),
                  format(length(unique(c(edgelist$pre, edgelist$post))), big.mark = ",")))

  # Build minimal metadata from filtered neuron IDs
  all_ids <- unique(c(edgelist$pre, edgelist$post))
  meta_ic <- data.frame(root_id = all_ids, stringsAsFactors = FALSE)
  message(sprintf("Metadata: %s neurons", format(nrow(meta_ic), big.mark = ",")))

  # Create SQLite database for the Python calculator
  message("Creating SQLite database: ", sqlite_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  tryCatch({
    # Add post_count if needed
    if (!"post_count" %in% names(edgelist)) {
      edgelist$post_count <- round(edgelist$count / edgelist$norm, digits = 6)
    }
    DBI::dbWriteTable(con, "meta", meta_ic, overwrite = TRUE)
    DBI::dbWriteTable(con, "edgelist_simple", edgelist, overwrite = TRUE)
  }, finally = DBI::dbDisconnect(con))
  message("SQLite database created: ", format(structure(file.size(sqlite_path), class = "object_size"), units = "auto"))

  rm(edgelist, meta_ic, all_ids)
  gc()
} else {
  message("Using existing SQLite database: ", sqlite_path)
}

####################################
### INITIALIZE CALCULATOR        ###
####################################

message("Initializing Python influence calculator...")
message("  count_thresh = ", count_thresh,
        if (count_thresh <= 1) " (un-thresholded)" else " (edges kept where count >= count_thresh)")
message("  signed = FALSE")

t_init_start <- Sys.time()
ic <- influence_calculator_py(
  filename = sqlite_path,
  count_thresh = count_thresh,
  signed = FALSE
)

# Get neuron IDs from the Python calculator's index mapping
# NOTE: id_to_index is a bidict.bidict (not a plain dict).
# reticulate cannot auto-convert bidict — names(py_to_r()) returns Python
# method names ("clear","copy",...) instead of neuron IDs.
# Fix: use Python's list() builtin via reticulate to extract keys directly,
# avoiding r.variable references which fail inside local() scope.
builtins <- reticulate::import_builtins()
neuron_ids <- as.character(reticulate::py_to_r(
  builtins$list(ic$id_to_index$keys())
))
n_neurons <- as.integer(ic$n_neurons)

# Validate extracted neuron IDs
bidict_methods <- c("clear", "copy", "equals_order_sensitive", "forceput",
                    "forceupdate", "get", "inv", "inverse", "items", "keys",
                    "on_dup", "pop", "popitem", "put", "putall",
                    "setdefault", "update", "values")
bad_ids <- neuron_ids %in% bidict_methods | is.na(neuron_ids) | neuron_ids == ""
if (any(bad_ids)) {
  stop(sprintf(
    "neuron_ids extraction failed: %d/%d IDs are bidict method names or NA (got: %s). ",
    sum(bad_ids), length(neuron_ids),
    paste(head(neuron_ids[bad_ids], 5), collapse = ", ")),
    "Check that reticulate can convert id_to_index.keys() correctly.")
}
if (length(neuron_ids) != n_neurons) {
  stop(sprintf("neuron_ids length (%d) != n_neurons (%d). Key extraction mismatch.",
               length(neuron_ids), n_neurons))
}
if (!all(grepl("^[0-9]+$", neuron_ids))) {
  non_numeric <- head(neuron_ids[!grepl("^[0-9]+$", neuron_ids)], 5)
  stop(sprintf("neuron_ids contains non-numeric values: %s. Expected BANC root IDs.",
               paste(non_numeric, collapse = ", ")))
}

message(sprintf("Calculator ready: %s neurons in connectivity matrix [init: %s]",
                format(n_neurons, big.mark = ","),
                format(round(difftime(Sys.time(), t_init_start, units = "mins"), 1))))

gc()

# Pre-flight: show sample neuron IDs to confirm extraction worked
message(sprintf("Sample neuron IDs: %s ...",
                paste(head(neuron_ids, 5), collapse = ", ")))

####################################
### CALCULATE ALL-TO-ALL         ###
####################################

n_batches <- ceiling(n_neurons / batch_size)
message(sprintf("Processing %s neurons in %d batches of %d",
                format(n_neurons, big.mark = ","), n_batches, batch_size))
message(sprintf("Estimated time: ~%.1f hours (after initial factorization)",
                n_neurons * 0.4 / 3600))

t_start <- Sys.time()
batches_completed <- 0

for (b in seq_len(n_batches)) {
  # Sharding: only this shard's stride
  if (((b - 1L) %% shard_total) != shard_idx) next

  chunk_file <- file.path(output_dir, sprintf("chunk_%04d.parquet", b))

  # Skip already-processed chunks (resumable), unless overwrite requested
  if (file.exists(chunk_file) && !overwrite_chunks) {
    batches_completed <- batches_completed + 1
    if (b %% 50 == 0 || b == 1) {
      message(sprintf("  Skipping batch %d/%d (already exists)", b, n_batches))
    }
    next
  }

  start_idx <- (b - 1) * batch_size + 1
  end_idx <- min(b * batch_size, n_neurons)
  batch_neuron_ids <- neuron_ids[start_idx:end_idx]

  t_batch_start <- Sys.time()

  batch_results <- vector("list", length(batch_neuron_ids))

  for (i in seq_along(batch_neuron_ids)) {
    seed_id <- batch_neuron_ids[i]

    # Skip invalid seed IDs (NA or non-numeric — guards against bidict bug)
    if (is.na(seed_id) || !grepl("^[0-9]+$", seed_id)) {
      warning(sprintf("Skipping invalid seed_id at index %d: '%s'", i, seed_id))
      next
    }

    inf <- tryCatch(
      calculate_influence_py(ic, seed_ids = seed_id, const = const),
      error = function(e) {
        warning(sprintf("Error for seed %s: %s", seed_id, e$message))
        NULL
      }
    )

    if (is.null(inf) || nrow(inf) == 0) next

    # Extract raw influence scores
    raw_scores <- inf[["Influence_score_(unsigned)"]]
    downstream_ids <- inf$id

    # Optional threshold filter to reduce storage
    if (!is.null(min_influence)) {
      keep <- raw_scores >= min_influence
      raw_scores <- raw_scores[keep]
      downstream_ids <- downstream_ids[keep]
    }

    batch_results[[i]] <- data.frame(
      upstream_id = seed_id,
      downstream_id = downstream_ids,
      raw_influence = raw_scores,
      stringsAsFactors = FALSE
    )
  }

  # Combine and write
  batch_results <- batch_results[!sapply(batch_results, is.null)]
  if (length(batch_results) > 0) {
    batch_df <- dplyr::bind_rows(batch_results)
    # Validate: upstream_id must be a valid numeric root_id, not a bidict method
    valid_rows <- grepl("^[0-9]+$", batch_df$upstream_id) & !is.na(batch_df$upstream_id)
    if (sum(valid_rows) == 0) {
      warning(sprintf("Batch %d: no valid upstream_ids, skipping write", b))
    } else {
      if (sum(!valid_rows) > 0) {
        warning(sprintf("Batch %d: dropping %d rows with invalid upstream_id",
                        b, sum(!valid_rows)))
        batch_df <- batch_df[valid_rows, ]
      }
      arrow::write_parquet(batch_df, chunk_file)
    }
  } else {
    warning(sprintf("Batch %d: no results (all seeds returned NULL), skipping write", b))
  }

  batches_completed <- batches_completed + 1
  elapsed_hrs <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))
  batch_secs <- as.numeric(difftime(Sys.time(), t_batch_start, units = "secs"))
  eta_hrs <- (n_batches - batches_completed) * batch_secs / 3600

  message(sprintf("  Batch %d/%d (%.0f%%): neurons %d-%d | %.1fs/batch | elapsed: %.1fh | ETA: %.1fh",
                  b, n_batches, 100 * batches_completed / n_batches,
                  start_idx, end_idx, batch_secs, elapsed_hrs, eta_hrs))

  rm(batch_results)
  if (exists("batch_df")) rm(batch_df)
  gc()
}

t_end <- Sys.time()
message(sprintf("### banc: all-to-all influence calculation complete ###"))
message(sprintf("Total time: %s", format(round(difftime(t_end, t_start, units = "hours"), 1))))
message(sprintf("Results in: %s", output_dir))

# Report output size
chunk_files <- list.files(output_dir, pattern = "\\.parquet$", full.names = TRUE)
total_size <- sum(file.info(chunk_files)$size, na.rm = TRUE)
message(sprintf("Total output: %d chunk files, %s",
                length(chunk_files),
                format(structure(total_size, class = "object_size"), units = "auto")))

# Post-completion validation: spot-check a random chunk
if (length(chunk_files) > 0) {
  check_file <- sample(chunk_files, 1)
  check_df <- arrow::read_parquet(check_file)
  n_valid_upstream <- sum(grepl("^[0-9]+$", check_df$upstream_id), na.rm = TRUE)
  n_nonzero <- sum(check_df$raw_influence > 0, na.rm = TRUE)
  message(sprintf("Spot-check %s: %d rows, %d valid upstream_ids, %d non-zero influence",
                  basename(check_file), nrow(check_df), n_valid_upstream, n_nonzero))
  if (n_valid_upstream == 0) {
    stop("VALIDATION FAILED: spot-checked chunk has no valid upstream_ids. ",
         "The neuron ID extraction may still be broken. Aborting upload.")
  }
  if (n_nonzero == 0) {
    warning("Spot-checked chunk has all-zero influence scores — this may indicate a problem.")
  }
  rm(check_df)
}

# Sync influence directory to versioned GCS path (matches banc-data.R pattern).
# Skip when sharded — run banc-sync-influence.R after all shards complete.
if (shard_total > 1L) {
  message(sprintf("Sharded run (shard %d/%d): skipping GCS sync. ",
                  shard_idx, shard_total),
          "Run banc-sync-influence.R after all shards finish.")
} else {
  message("Syncing influence files to Google Cloud Storage...")
  gs_dest <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/influence/all_to_all/",
                     banc.version)
  upload_cmd <- sprintf("gsutil -m rsync -r %s %s", output_dir, gs_dest)
  upload_status <- system(upload_cmd)
  if (upload_status == 0) {
    message(sprintf("Influence parquet files (%d chunks) synced to: %s",
                    length(chunk_files), gs_dest))
  } else {
    warning("Failed to sync influence parquet files to GCS (exit code: ", upload_status, ")")
  }
}

})