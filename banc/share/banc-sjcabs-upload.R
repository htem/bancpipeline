#' banc-sjcabs-upload — Sync exported BANC/maleCNS/MANC connectomes to GCS.
#'
#' @section Reads:
#'   - `/n/data1/hms/neurobio/wilson/connectomes/{banc,malecns,manc}/`
#'
#' @section Writes:
#'   - GCS `compiled_data/`
#'
#' @section CLI:
#'   [datasets...]  subset of `banc malecns manc` (default all)

###########################################################
### Upload processed connectome data to GCS
###
### Syncs exported data for BANC, maleCNS, and MANC
### to the shared GCS bucket for SJCABS and general access.
###
### GCS bucket: gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data
###
### Dependencies (must run first):
###   BANC:    banc/share/banc-sjcabs.R
###   maleCNS: malecns/malecns-sjcabs.R
###   MANC:    manc/manc-sjcabs.R
###
### These scripts export processed files to:
###   /n/data1/hms/neurobio/wilson/connectomes/{banc,malecns,manc}/
###
### Usage: Rscript banc/sharing/sjcabs.R
###        Rscript banc/sharing/sjcabs.R banc
###        Rscript banc/sharing/sjcabs.R malecns manc
###########################################################

gcs_bucket <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
local_root <- "/n/data1/hms/neurobio/wilson/connectomes"

datasets <- list(
  banc    = file.path(local_root, "banc"),
  malecns = file.path(local_root, "malecns"),
  manc    = file.path(local_root, "manc")
)

# Parse command-line args: if specified, only sync those datasets
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  valid <- intersect(args, names(datasets))
  if (length(valid) == 0) {
    stop("Unknown dataset(s): ", paste(args, collapse = ", "),
         "\n  Valid options: ", paste(names(datasets), collapse = ", "))
  }
  datasets <- datasets[valid]
}

message("=== Uploading processed data to GCS ===")
message(sprintf("Bucket: %s", gcs_bucket))
message(sprintf("Datasets: %s", paste(names(datasets), collapse = ", ")))

for (name in names(datasets)) {
  local_path <- datasets[[name]]
  gcs_path <- file.path(gcs_bucket, name)

  if (!dir.exists(local_path)) {
    message(sprintf("\n[%s] Local path not found: %s — skipping", name, local_path))
    next
  }

  # Show what will be synced
  n_files <- length(list.files(local_path, recursive = TRUE))
  total_size <- sum(file.info(list.files(local_path, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
  message(sprintf("\n[%s] %d files (%.1f GB) -> %s",
                  name, n_files, total_size / 1024^3, gcs_path))

  # gsutil -m rsync: multi-threaded recursive sync
  # -r: recursive
  # -d: delete extra files at destination (mirror exact copy)
  cmd <- sprintf("gsutil -m rsync -r %s %s", shQuote(local_path), gcs_path)
  message(sprintf("  Running: %s", cmd))

  t0 <- Sys.time()
  exit_code <- system(cmd)
  elapsed <- difftime(Sys.time(), t0, units = "mins")

  if (exit_code == 0) {
    message(sprintf("  [%s] Sync complete [%.1f min]", name, as.numeric(elapsed)))
  } else {
    message(sprintf("  [%s] gsutil exited with code %d", name, exit_code))
  }
}

message("\n=== Upload complete ===")
