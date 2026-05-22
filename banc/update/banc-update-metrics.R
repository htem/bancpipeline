#' banc-update-metrics — Join per-metric feathers and push metric columns to SeaTable.
#'
#' Reads per-metric feather files (`banc_root_positions`, `banc_regions`,
#' `banc_detections`, `banc_l2_metrics`, `banc_volumes`), joins by
#' `root_id`, writes combined `banc_metrics.feather`, then pushes metric
#' columns to SeaTable.
#'
#' @section Reads:
#'   - `<banc.save.path>/banc_{root_positions,regions,detections,l2_metrics,volumes}.feather`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `<banc.save.path>/banc_metrics.feather`
#'   - SeaTable `banc_meta`: cols `region`, `side`, `l2_nodes`,
#'     `l2_cable_length_um`, `input_connections`, `output_connections`,
#'     `input_side_index`, `output_side_index`, `mitochondria`,
#'     `mitochondria_volume`, `pd_width`, `segregation_index`, `volume_nm3`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_metrics.sh`,
#'   `o2/production/o2_banc_update.sh`

###########################################################
### Join per-metric feather files, write combined
### banc_metrics.feather, and push to BANC seatable
###
### Reads per-metric feather files (populated by calculate-*
### scripts), joins by root_id, writes combined feather,
### then pushes metric columns to seatable.
###
### Per-metric files:
###   banc_root_positions.feather, banc_regions.feather,
###   banc_detections.feather, banc_l2_metrics.feather,
###   banc_volumes.feather
###
### Columns pushed: region, side, l2_nodes, l2_cable_length_um,
###   input_connections, output_connections,
###   input_side_index, output_side_index,
###   mitochondria, mitochondria_volume,
###   pd_width, segregation_index, volume_nm3
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: joining metrics and pushing to seatable ###")
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read feather files  ###
###########################

# Helper: read feather if it exists, else return empty df
read_metric <- function(filename) {
  path <- file.path(banc.save.path, filename)
  if (file.exists(path)) {
    arrow::read_feather(path) %>%
      dplyr::mutate(root_id = as.character(root_id)) %>%
      dplyr::filter(!is.na(root_id)) %>%
      dplyr::distinct(root_id, .keep_all = TRUE)
  } else {
    message(sprintf("  %s not found — skipping", filename))
    data.frame(root_id = character(0))
  }
}

positions  <- read_metric("banc_root_positions.feather")
regions    <- read_metric("banc_regions.feather")
synapses   <- read_metric("banc_detections.feather")
l2_metrics <- read_metric("banc_l2_metrics.feather")
volumes    <- read_metric("banc_volumes.feather")

if (nrow(positions) == 0 && nrow(regions) == 0 && nrow(synapses) == 0 &&
    nrow(l2_metrics) == 0 && nrow(volumes) == 0) {
  message("No metric feather files found. Run calculate-* scripts first.")
  return(invisible())
}

# Join all by root_id
banc.metrics <- positions %>%
  dplyr::full_join(regions, by = "root_id") %>%
  dplyr::full_join(synapses, by = "root_id") %>%
  dplyr::full_join(l2_metrics, by = "root_id") %>%
  dplyr::full_join(volumes, by = "root_id") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

message(sprintf("Combined %d neurons from per-metric feather files", nrow(banc.metrics)))

combined_file <- file.path(banc.save.path, "banc_metrics.feather")

# Select metric columns for seatable push
banc.metrics2 <- banc.metrics %>%
  dplyr::select(dplyr::any_of(c("root_id", "region", "side",
    "mitochondria", "mitochondria_volume",
    "input_connections", "output_connections",
    "input_side_index", "output_side_index",
    "l2_nodes", "l2_cable_length_um",
    "pd_width", "segregation_index", "volume_nm3")))

###########################
### Join with seatable  ###
###########################

bc <- banctable_query("SELECT _id, status, root_id, super_class, region, root_region from banc_meta") %>%
  dplyr::mutate(root_id = as.character(root_id))

# Snapshot original SeaTable region for the anti-join at push time. Only
# rows whose final pushed value differs from this go to banctable_update_rows.
bc.metrics.orig <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::select(`_id`, orig_region = region) %>%
  dplyr::mutate(orig_region = dplyr::coalesce(as.character(orig_region), ""))

# IMPORTANT: invert the join — start from SeaTable rows so newly-added
# neurons (not yet in any per-metric feather) still reach the case_when
# region inference. Pre-fix the script left_joined feather → SeaTable, which
# silently dropped new SeaTable rows: they could never get a region pushed,
# even when `super_class` / `root_region` would have been enough to infer one.
metrics.update <- bc %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::select(`_id`, root_region, root_id, bc_region = region, super_class) %>%
  dplyr::left_join(banc.metrics2, by = "root_id") %>%
  # Region refinement: prefer existing seatable region, then calculated, then infer from super_class.
  # "neck_connective" is no longer a valid region value — ascending/descending
  # neurons are identified via super_class; cervical_connective lives as a
  # separate "tract" entry. Stale "neck_connective" inputs collapse to NA so
  # downstream consumers fall back to super_class.
  dplyr::mutate(region = dplyr::case_when(
    !is.na(bc_region) & !bc_region %in% c("brain", "outside", "rind", "neck_connective", "", " ", "0") ~ bc_region,
    !is.na(region) & !region %in% c("brain", "rind", "neck_connective", "", "0") ~ region,
    super_class %in% c("visual_centrifugal") ~ "central_brain",
    super_class %in% c("visual_projection") ~ "optic_lobe",
    super_class %in% c("ventral_nerve_cord_intrinsic") ~ "ventral_nerve_cord",
    super_class %in% c("central_brain_intrinsic") ~ "central_brain",
    region %in% c("midbrain", "central_brain") ~ "central_brain",
    region %in% c("vnc", "ventral_nerve_cord") ~ "ventral_nerve_cord",
    region %in% c("optic") ~ "optic_lobe",
    grepl('CAN|FLA|GNG|AMMC|SAD|PRW', root_region) | region == "sez" ~ "central_brain",
    grepl('_midbrain_', root_region) ~ "central_brain",
    grepl('optic', root_region) ~ "optic_lobe",
    grepl("midbrain", root_region) ~ "central_brain",
    grepl("vnc", root_region) ~ "ventral_nerve_cord",
    TRUE ~ NA_character_
  )) %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(`_id`)) %>%
  as.data.frame()

# Write combined feather (after region enrichment so it includes refined regions)
feather_out <- metrics.update %>%
  dplyr::select(-dplyr::any_of(c("bc_region", "super_class", "root_region", "_id")))
arrow::write_feather(feather_out, combined_file)
message(sprintf("Wrote combined banc_metrics.feather (%d rows)", nrow(feather_out)))

###########################
### Clean and push      ###
###########################

# Correct data types
if ("segregation_index" %in% colnames(metrics.update))
  metrics.update$segregation_index <- as.numeric(metrics.update$segregation_index)
if ("region" %in% colnames(metrics.update))
  metrics.update$region[is.na(metrics.update$region)] <- ""
if ("side" %in% colnames(metrics.update))
  metrics.update$side[is.na(metrics.update$side)] <- ""
if ("mitochondria" %in% colnames(metrics.update))
  metrics.update$mitochondria <- as.numeric(metrics.update$mitochondria)
if ("mitochondria_volume" %in% colnames(metrics.update))
  metrics.update$mitochondria_volume <- as.numeric(metrics.update$mitochondria_volume)

# Remove helper columns
metrics.update$bc_region <- NULL
metrics.update$super_class <- NULL
metrics.update$root_region <- NULL

# Anti-join on region only (the column the inverted join can newly affect):
# rows whose pushed region matches the SeaTable's current region don't need
# to round-trip. Other metric columns may still differ, but pre-fix the same
# update was happening for every neuron in feather every run, so we're not
# regressing on those — only avoiding the new rows the invert added that
# carry region == orig region (typically "" → "" no-ops).
.before_anti <- nrow(metrics.update)
metrics.update <- metrics.update %>%
  dplyr::left_join(bc.metrics.orig, by = "_id") %>%
  dplyr::filter(is.na(orig_region) | orig_region != ifelse(is.na(region), "", region) |
                  # keep all rows that have OTHER metric data to push
                  !is.na(l2_cable_length_um) | !is.na(input_connections) |
                  !is.na(volume_nm3)) %>%
  dplyr::select(-orig_region)

message(sprintf("Pushing metrics for %d neurons to seatable (anti-join skipped %d region-no-op rows)",
                nrow(metrics.update), .before_anti - nrow(metrics.update)))

if (nrow(metrics.update)) {
  banctable_update_rows(base = 'banc_meta',
                        table = 'banc_meta',
                        df = metrics.update,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

message(sprintf("### banc: metrics push complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
