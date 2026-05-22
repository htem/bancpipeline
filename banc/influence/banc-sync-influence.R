#' banc-sync-influence — Validate + rsync influence parquet shards to GCS.
#'
#' Run AFTER all `banc-build-influence.R` array shards have completed.
#' Repeats the chunk completeness check + spot-check + gsutil rsync that
#' the un-sharded build script does at the end.
#'
#' @section Reads:
#'   - `<banc.versioned.save.path>/influence/all_to_all/chunk_*.parquet`
#'
#' @section Writes:
#'   - GCS `gs://lee-lab_..._/compiled_data/banc_<ver>/influence/all_to_all/`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_influence.sh`,
#'   `o2/production/o2_banc_influence_array.sh`,
#'   `o2/production/o2_banc_aggregate_influence.sh`,
#'   `o2/production/o2_banc_aggregate_influence_priority.sh`

###############################################
### Validate + sync BANC all-to-all influence
### parquet chunks to GCS.
###
### Run AFTER all banc-build-influence.R shards
### have completed. Performs the same chunk
### completeness check, spot-check, and gsutil
### rsync that the un-sharded build script does
### at the end.
###############################################
source("banc/banc-startup.R")

local({

message("### banc: validate + sync influence chunks ###")

# Match paths used by banc-build-influence.R
save.path  <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                        paste0("banc_", banc.version))
output_dir <- file.path(save.path, "influence", "all_to_all")
if (!dir.exists(output_dir)) {
  stop("Influence output directory does not exist: ", output_dir)
}

chunk_files <- list.files(output_dir, pattern = "^chunk_[0-9]+\\.parquet$",
                          full.names = TRUE)
if (length(chunk_files) == 0L) {
  stop("No chunk files found in ", output_dir)
}

total_size <- sum(file.info(chunk_files)$size, na.rm = TRUE)
message(sprintf("Found %d chunk files, %s total",
                length(chunk_files),
                format(structure(total_size, class = "object_size"), units = "auto")))

# Completeness check: chunk indices should be a contiguous run 1..N
indices <- as.integer(sub("^chunk_0*([0-9]+)\\.parquet$", "\\1",
                          basename(chunk_files)))
indices <- sort(indices)
expected <- seq_len(max(indices))
missing  <- setdiff(expected, indices)
if (length(missing) > 0L) {
  warning(sprintf("Missing %d chunk(s): %s%s",
                  length(missing),
                  paste(head(missing, 10), collapse = ", "),
                  if (length(missing) > 10L) ", ..." else ""))
} else {
  message(sprintf("All %d chunk indices present (1..%d)",
                  length(indices), max(indices)))
}

# Spot-check a random chunk for valid upstream IDs and non-zero scores
check_file <- sample(chunk_files, 1)
check_df   <- arrow::read_parquet(check_file)
n_valid_upstream <- sum(grepl("^[0-9]+$", check_df$upstream_id), na.rm = TRUE)
n_nonzero        <- sum(check_df$raw_influence > 0, na.rm = TRUE)
message(sprintf("Spot-check %s: %d rows, %d valid upstream_ids, %d non-zero influence",
                basename(check_file), nrow(check_df), n_valid_upstream, n_nonzero))
if (n_valid_upstream == 0L) {
  stop("VALIDATION FAILED: spot-checked chunk has no valid upstream_ids. Aborting upload.")
}
if (n_nonzero == 0L) {
  warning("Spot-checked chunk has all-zero influence scores — this may indicate a problem.")
}
rm(check_df)

# Sync to versioned GCS path
message("Syncing influence files to Google Cloud Storage...")
gs_dest <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/influence/all_to_all/",
                   banc.version)
upload_cmd <- sprintf("gsutil -m rsync -r %s %s", output_dir, gs_dest)
upload_status <- system(upload_cmd)
if (upload_status == 0) {
  message(sprintf("Influence parquet files (%d chunks) synced to: %s",
                  length(chunk_files), gs_dest))
} else {
  stop("Failed to sync influence parquet files to GCS (exit code: ", upload_status, ")")
}

message("### banc: influence sync complete ###")

})
