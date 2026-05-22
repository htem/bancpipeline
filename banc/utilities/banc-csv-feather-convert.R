#' banc-csv-feather-convert — One-time migration: CSV → per-metric feather files.
#'
#' Reads the legacy monolithic `banc_root_positions.csv` and splits it
#' into the per-metric feather files the modern pipeline expects.
#' Run-once migration; safe to invoke idempotently.
#'
#' @section Reads:
#'   - `<banc.save.path>/banc_root_positions.csv`
#'   - NBLAST summary CSVs under `<banc.nblast.save.path>/`
#'
#' @section Writes:
#'   - `<banc.save.path>/banc_root_positions.feather`
#'   - `<banc.save.path>/banc_regions.feather`
#'   - `<banc.save.path>/banc_detections.feather`
#'   - `<banc.save.path>/banc_l2_metrics.feather`
#'   - `<banc.save.path>/banc_volumes.feather`
#'   - NBLAST .feather siblings of the input CSVs

###########################################################
### One-time migration: CSV → per-metric feather files
###
### Reads the monolithic banc_root_positions.csv and splits
### it into per-metric feather files:
###   - banc_root_positions.feather
###   - banc_regions.feather
###   - banc_detections.feather
###   - banc_l2_metrics.feather
###   - banc_volumes.feather
###
### Also converts NBLAST summary CSVs → feather if the
### feather versions don't already exist.
###
### Usage: Rscript banc/utilities/banc-csv-feather-convert.R
###########################################################
source("banc/banc-startup.R")

message("========================================")
message("=== CSV → Feather Migration          ===")
message("========================================")

###########################
### Metrics CSV → split ###
###########################

csv_file <- file.path(banc.save.path, "banc_root_positions.csv")
if (!file.exists(csv_file)) {
  stop("banc_root_positions.csv not found at: ", csv_file)
}

message("Reading banc_root_positions.csv...")
banc.metrics <- readr::read_csv(csv_file, col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("  %d rows read", nrow(banc.metrics)))

# Helper: select only columns that exist, write feather
write_metric_feather <- function(df, cols, filename) {
  cols_present <- intersect(cols, colnames(df))
  if (length(cols_present) <= 1) {
    message(sprintf("  Skipping %s — no metric columns found", filename))
    return(invisible())
  }
  out <- df %>% dplyr::select(dplyr::all_of(cols_present))
  outpath <- file.path(banc.save.path, filename)
  arrow::write_feather(out, outpath)
  message(sprintf("  Wrote %s (%d rows, cols: %s)",
                  filename, nrow(out), paste(cols_present, collapse = ", ")))
}

# 1. Root positions
write_metric_feather(banc.metrics,
  c("root_id", "root_position", "root_position_nm", "nucleus_id", "side"),
  "banc_root_positions.feather")

# 2. Regions
write_metric_feather(banc.metrics,
  c("root_id", "region"),
  "banc_regions.feather")

# 3. Synapses
write_metric_feather(banc.metrics,
  c("root_id", "input_connections", "output_connections",
    "input_side_index", "output_side_index",
    "mitochondria", "mitochondria_volume"),
  "banc_detections.feather")

# 4. L2 metrics
write_metric_feather(banc.metrics,
  c("root_id", "l2_nodes", "l2_cable_length_um", "pd_width", "segregation_index"),
  "banc_l2_metrics.feather")

# 5. Volumes
write_metric_feather(banc.metrics,
  c("root_id", "volume_nm3"),
  "banc_volumes.feather")

###########################
### NBLAST CSVs         ###
###########################

message("\nConverting NBLAST summary CSVs → feather...")

nblast_csvs <- list(
  "banc_fafb_783_nblast",
  "banc_fanc_1116_nblast",
  "banc_manc_v1.2.1_nblast",
  "banc_hemibrain_v1.2.1_nblast",
  "banc_malecns_v0.9_nblast",
  "banc_mirror_nblast"
)

for (name in nblast_csvs) {
  feather_path <- file.path(banc.meta.save.path, paste0(name, ".feather"))
  csv_path <- file.path(banc.meta.save.path, paste0(name, ".csv"))

  if (file.exists(feather_path)) {
    message(sprintf("  %s.feather already exists — skipping", name))
    next
  }
  if (!file.exists(csv_path)) {
    message(sprintf("  %s.csv not found — skipping", name))
    next
  }

  tryCatch({
    df <- readr::read_csv(csv_path, col_types = banc.col.types, show_col_types = FALSE)
    arrow::write_feather(df, feather_path)
    message(sprintf("  Converted %s.csv → .feather (%d rows)", name, nrow(df)))
  }, error = function(e) {
    message(sprintf("  Error converting %s: %s", name, e$message))
  })
}

message("\n=== Migration complete ===")
message("You can now run the pipeline with feather-based scripts.")
message("The original CSV files have NOT been deleted — remove manually when ready.")
