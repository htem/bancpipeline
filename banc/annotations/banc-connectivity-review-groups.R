#' banc-connectivity-review-groups — Assemble cross-matched neuron groups for review.
#'
#' Three predefined groups (`antennal_lobe`, `central_complex`, `brain_motor`):
#' joins BANC↔FAFB matches and Codex alignment scores into a manual-review CSV.
#'
#' @section Reads:
#'   - GCS `compiled_data/banc_<ver>/*.feather`, SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `data/codex/connectivity_review_<group>.csv`

##############################################################
### GROUPS OF CROSS MATCHED NEURONS FOR CONNECIVITY REVIEW ###
##############################################################
source("banc/banc-startup.R")
redo <- FALSE

# DEFINE GROUPS
groups <- c("antennal_lobe","central_complex", "brain_motor")

# GCS data paths
gcs_base <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
cache_dir <- file.path(tempdir(), "gcs_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: download from GCS if not cached locally
gcs_cache <- function(gcs_path) {
  local_file <- file.path(cache_dir, basename(gcs_path))
  if (!file.exists(local_file)) {
    message("  Downloading ", basename(gcs_path), " from GCS...")
    system2("gsutil", c("cp", gcs_path, local_file),
            stdout = FALSE, stderr = FALSE)
    if (!file.exists(local_file)) stop("gsutil cp failed for ", gcs_path)
  } else {
    message("  Using cached ", basename(local_file))
  }
  local_file
}

# DATA
## meta
message("=== Loading metadata ===")
banc.meta <- banctable_query("SELECT root_id, root_626, supervoxel_id, region, side, root_position, position, hemilineage, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, hemibrain_cell_type, manc_cell_type, malecns_cell_type, fafb_match, hemibrain_match, manc_match, malecns_match FROM banc_meta")
banc.meta <- banc.meta %>%
  dplyr::mutate(root_626 = as.character(root_626),
                root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id))

# Fill missing side from root_position (voxels), falling back to position (voxels)
missing_side <- is.na(banc.meta$side) | banc.meta$side == ""
if (any(missing_side)) {
  # Try root_position first
  has_root_pos <- missing_side & !is.na(banc.meta$root_position) & banc.meta$root_position != ""
  if (any(has_root_pos)) {
    rp_nm <- bancr::banc_raw2nm(nat::xyzmatrix(banc.meta$root_position[has_root_pos]))
    lr <- bancr:::banc_lr_position(rp_nm, units = "nm")
    banc.meta$side[has_root_pos] <- ifelse(lr > 0, "right", "left")
  }
  # Fall back to position for any still missing
  still_missing <- is.na(banc.meta$side) | banc.meta$side == ""
  has_pos <- still_missing & !is.na(banc.meta$position) & banc.meta$position != ""
  if (any(has_pos)) {
    pos_nm <- bancr::banc_raw2nm(nat::xyzmatrix(banc.meta$position[has_pos]))
    lr <- bancr:::banc_lr_position(pos_nm, units = "nm")
    banc.meta$side[has_pos] <- ifelse(lr > 0, "right", "left")
  }
  n_filled <- sum(missing_side) - sum(is.na(banc.meta$side) | banc.meta$side == "")
  message(sprintf("  Filled side for %d / %d neurons from root_position/position", n_filled, sum(missing_side)))
}

# Get root_746 (banc_746_id) by joining on supervoxel_id
message("  Mapping root_626 -> root_746 via supervoxel_id...")
banc_746_meta <- arrow::read_feather(gcs_cache(
  file.path(gcs_base, "banc_746", "banc_746_meta.feather")))
banc_746_ids <- banc_746_meta %>%
  dplyr::select(banc_746_id, supervoxel_id) %>%
  dplyr::filter(!is.na(supervoxel_id)) %>%
  dplyr::mutate(supervoxel_id = as.character(supervoxel_id),
                banc_746_id = as.character(banc_746_id)) %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE)
banc.meta <- banc.meta %>%
  dplyr::left_join(banc_746_ids, by = "supervoxel_id") %>%
  dplyr::rename(root_746 = banc_746_id)
message(sprintf("  BANC: %d neurons, %d with root_746",
                nrow(banc.meta), sum(!is.na(banc.meta$root_746))))
rm(banc_746_meta, banc_746_ids); gc()

franken.meta <- franken_meta("SELECT fafb_id, manc_id, region, side, hemilineage, super_class, cell_class, cell_sub_class, cell_type FROM franken_meta",
                                base = "cns_meta")
fafb.meta <- dplyr::filter(franken.meta, !is.na(fafb_id)) %>%
  dplyr::mutate(fafb_id = as.character(fafb_id))
manc.meta <- dplyr::filter(franken.meta, !is.na(manc_id)) %>%
  dplyr::mutate(manc_id = as.character(manc_id))
malecns.meta <- franken_meta("SELECT malecns_09_id, region, side, hemilineage, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, hemibrain_cell_type, manc_cell_type FROM malecns",
                                base = "cns_meta") %>%
  dplyr::mutate(malecns_09_id = as.character(malecns_09_id))
# message(sprintf("  FAFB: %d, MANC: %d, maleCNS: %d neurons",
#                 nrow(fafb.meta), nrow(manc.meta), nrow(malecns.meta)))

## edgelists
message("\n=== Loading edgelists ===")
banc.el <- arrow::read_feather(gcs_cache(
  file.path(gcs_base, "banc_746", "banc_746_simple_edgelist.feather")))
banc.el <- banc.el %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
message(sprintf("  BANC edgelist: %d edges", nrow(banc.el)))

fafb.el <- arrow::read_feather(gcs_cache(
  file.path(gcs_base, "fafb_783", "fafb_783_simple_edgelist.feather")))
fafb.el <- fafb.el %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
message(sprintf("  FAFB edgelist: %d edges", nrow(fafb.el)))

manc.el <- arrow::read_feather(gcs_cache(
  file.path(gcs_base, "manc_121", "manc_121_simple_edgelist.feather")))
manc.el <- manc.el %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
message(sprintf("  MANC edgelist: %d edges", nrow(manc.el)))

malecns.el <- arrow::read_feather(gcs_cache(
  file.path(gcs_base, "malecns_09", "malecns_09_simple_edgelist.feather")))
malecns.el <- malecns.el %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
message(sprintf("  maleCNS edgelist: %d edges", nrow(malecns.el)))

# Save path
save.path <- "data/connectivity_review/"

# Helper: clean a cell type vector (remove NA, "", and "unknown"-containing)
clean_cts <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & x != ""]
  x[!grepl("unknown", x, ignore.case = TRUE)]
}

# GRAB
for(group in groups){

  message(sprintf("\n=== Processing group: %s ===", group))

  # Select meta data
  if(group == "central_complex"){
    fafb.meta.chosen <- fafb.meta %>%
      dplyr::filter(grepl("central_complex|CX",cell_class)|grepl("central_complex|CX",cell_sub_class))
    manc.meta.chosen <- manc.meta %>%
      dplyr::filter(grepl("central_complex|CX",cell_class)|grepl("central_complex|CX",cell_sub_class))
    malecns.meta.chosen <- malecns.meta %>%
      dplyr::filter(grepl("central_complex|CX",cell_class)|grepl("central_complex|CX",cell_sub_class))
    all.cts <- clean_cts(c(fafb.meta.chosen$cell_type,
                            manc.meta.chosen$cell_type,
                            malecns.meta.chosen$cell_type))
    banc.meta.chosen <- banc.meta %>%
      dplyr::filter(grepl("central_complex|CX",cell_class)|grepl("central_complex|CX",cell_sub_class)|cell_type%in%all.cts)
  }else if (group == "antennal_lobe"){
    fafb.meta.chosen <- fafb.meta %>%
      dplyr::filter(grepl("antennal_lobe|^AL",cell_class)|grepl("antennal_lobe|^AL",cell_sub_class))
    manc.meta.chosen <- manc.meta %>%
      dplyr::filter(grepl("antennal_lobe|^AL",cell_class)|grepl("antennal_lobe|^AL",cell_sub_class))
    malecns.meta.chosen <- malecns.meta %>%
      dplyr::filter(grepl("antennal_lobe|^AL",cell_class)|grepl("antennal_lobe|^AL",cell_sub_class))
    all.cts <- clean_cts(c(fafb.meta.chosen$cell_type,
                            manc.meta.chosen$cell_type,
                            malecns.meta.chosen$cell_type))
    banc.meta.chosen <- banc.meta %>%
      dplyr::filter(grepl("antennal_lobe|^AL",cell_class)|grepl("antennal_lobe|^AL",cell_sub_class)|cell_type%in%all.cts)
  }else if(group == "brain_motor"){
    fafb.meta.chosen <- fafb.meta %>%
      dplyr::filter(grepl("motor",super_class)&grepl("brain",region))
    manc.meta.chosen <- manc.meta %>%
      dplyr::filter(grepl("motor",super_class)&grepl("brain",region))
    malecns.meta.chosen <- malecns.meta %>%
      dplyr::filter(grepl("motor",super_class)&grepl("brain",region))
    all.cts <- clean_cts(c(fafb.meta.chosen$cell_type,
                            manc.meta.chosen$cell_type,
                            malecns.meta.chosen$cell_type))
    banc.meta.chosen <- banc.meta %>%
      dplyr::filter(grepl("motor",super_class)&grepl("brain",region)|cell_type%in%all.cts)
  }

  # Pool sizes
  message(sprintf("  Pool sizes:"))
  message(sprintf("    BANC:    %d neurons (%d with root_746)",
                  nrow(banc.meta.chosen), sum(!is.na(banc.meta.chosen$root_746))))
  message(sprintf("    FAFB:    %d neurons", nrow(fafb.meta.chosen)))
  message(sprintf("    MANC:    %d neurons", nrow(manc.meta.chosen)))
  message(sprintf("    maleCNS: %d neurons", nrow(malecns.meta.chosen)))
  message(sprintf("    TOTAL:   %d neurons across all datasets",
                  nrow(banc.meta.chosen) + nrow(fafb.meta.chosen) +
                    nrow(manc.meta.chosen) + nrow(malecns.meta.chosen)))

  # Get edgelists: pre AND post both in selected neuron pool
  banc_ids <- na.omit(unique(banc.meta.chosen$root_746))
  banc.el.chosen <- banc.el %>%
    dplyr::filter(pre %in% banc_ids, post %in% banc_ids) %>%
    dplyr::select(pre, post, count, norm)

  fafb_ids <- na.omit(unique(fafb.meta.chosen$fafb_id))
  fafb.el.chosen <- fafb.el %>%
    dplyr::filter(pre %in% fafb_ids, post %in% fafb_ids) %>%
    dplyr::select(pre, post, count, norm)

  manc_ids <- na.omit(unique(manc.meta.chosen$manc_id))
  manc.el.chosen <- manc.el %>%
    dplyr::filter(pre %in% manc_ids, post %in% manc_ids) %>%
    dplyr::select(pre, post, count, norm)

  malecns_ids <- na.omit(unique(malecns.meta.chosen$malecns_09_id))
  malecns.el.chosen <- malecns.el %>%
    dplyr::filter(pre %in% malecns_ids, post %in% malecns_ids) %>%
    dplyr::select(pre, post, count, norm)

  message(sprintf("  Edgelist sizes (pre & post in pool):"))
  message(sprintf("    BANC:    %d edges", nrow(banc.el.chosen)))
  message(sprintf("    FAFB:    %d edges", nrow(fafb.el.chosen)))
  message(sprintf("    MANC:    %d edges", nrow(manc.el.chosen)))
  message(sprintf("    maleCNS: %d edges", nrow(malecns.el.chosen)))

  # Expand edgelists via maleCNS neighbour cell types.
  # 1. Find maleCNS neurons connected to chosen neurons with count >= 5
  # 2. Collect their cell types (+ cross-dataset type columns) into one vector
  # 3. Use that vector to expand all datasets by matching on cell_type

  # maleCNS neighbours of chosen neurons with count >= 5
  malecns_all_ids <- na.omit(unique(malecns.meta$malecns_09_id))
  malecns_neighbour_edges <- malecns.el %>%
    dplyr::filter(count >= 10,
                  (pre %in% malecns_ids | post %in% malecns_ids),
                  pre %in% malecns_all_ids, 
                  post %in% malecns_all_ids)
  malecns_neighbour_ids <- setdiff(
    unique(c(malecns_neighbour_edges$pre, malecns_neighbour_edges$post)),
    malecns_ids
  )

  # Look up cell types for those neighbour neurons in maleCNS meta
  malecns_neighbour_meta <- malecns.meta %>%
    dplyr::filter(malecns_09_id %in% malecns_neighbour_ids)
  expansion_cts <- clean_cts(c(
    malecns_neighbour_meta$cell_type,
    malecns_neighbour_meta$fafb_cell_type,
    malecns_neighbour_meta$hemibrain_cell_type,
    malecns_neighbour_meta$manc_cell_type
  ))
  message(sprintf("  Expansion: %d maleCNS neighbour neurons -> %d unique cell types",
                  length(malecns_neighbour_ids), length(expansion_cts)))

  # Expand: edges where BOTH pre and post are in chosen + expansion neurons
  malecns_expanded_ids <- na.omit(unique(c(
    malecns_ids,
    malecns.meta$malecns_09_id[malecns.meta$cell_type %in% expansion_cts]
  )))
  malecns.el.chosen <- malecns.el %>%
    dplyr::filter(pre %in% malecns_expanded_ids & post %in% malecns_ids| 
                  post %in% malecns_expanded_ids & pre %in% malecns_ids) %>%
    dplyr::select(pre, post, count, norm)

  banc_expanded_ids <- na.omit(unique(c(
    banc_ids,
    banc.meta$root_746[banc.meta$cell_type %in% expansion_cts]
  )))
  banc.el.chosen <- banc.el %>%
    dplyr::filter(pre %in% banc_expanded_ids & post %in% banc_ids|
                  post %in% banc_expanded_ids & pre %in% banc_ids) %>%
    dplyr::select(pre, post, count, norm)

  fafb_expanded_ids <- na.omit(unique(c(
    fafb_ids,
    fafb.meta$fafb_id[fafb.meta$cell_type %in% expansion_cts]
  )))
  fafb.el.chosen <- fafb.el %>%
    dplyr::filter(pre %in% fafb_expanded_ids & post %in% fafb_ids|
                  post %in% fafb_expanded_ids & pre %in% fafb_ids) %>%
    dplyr::select(pre, post, count, norm)

  manc_expanded_ids <- na.omit(unique(c(
    manc_ids,
    manc.meta$manc_id[manc.meta$cell_type %in% expansion_cts]
  )))
  manc.el.chosen <- manc.el %>%
    dplyr::filter(pre %in% manc_expanded_ids & post %in% manc_ids| 
                  post %in% manc_expanded_ids & pre %in% manc_ids) %>%
    dplyr::select(pre, post, count, norm)

  # Count unique neurons in expanded edgelists
  n_banc_exp <- length(unique(c(banc.el.chosen$pre, banc.el.chosen$post)))
  n_fafb_exp <- length(unique(c(fafb.el.chosen$pre, fafb.el.chosen$post)))
  n_manc_exp <- length(unique(c(manc.el.chosen$pre, manc.el.chosen$post)))
  n_malecns_exp <- length(unique(c(malecns.el.chosen$pre, malecns.el.chosen$post)))

  message(sprintf("  Expanded edgelist sizes (+ neighbours with count >= 2):"))
  message(sprintf("    BANC:    %d edges, %d neurons", nrow(banc.el.chosen), n_banc_exp))
  message(sprintf("    FAFB:    %d edges, %d neurons", nrow(fafb.el.chosen), n_fafb_exp))
  message(sprintf("    MANC:    %d edges, %d neurons", nrow(manc.el.chosen), n_manc_exp))
  message(sprintf("    maleCNS: %d edges, %d neurons", nrow(malecns.el.chosen), n_malecns_exp))
  message(sprintf("    TOTAL:   %d neurons across all datasets",
                  n_banc_exp + n_fafb_exp + n_manc_exp + n_malecns_exp))

  # Save
  save.path.group <- file.path(save.path, group)
  dir.create(save.path.group, recursive = TRUE, showWarnings = FALSE)

  # Meta
  if (nrow(banc.meta.chosen) > 0)
    readr::write_csv(banc.meta.chosen, file.path(save.path.group, sprintf("banc_%s_meta.csv", group)))
  if (nrow(fafb.meta.chosen) > 0)
    readr::write_csv(fafb.meta.chosen, file.path(save.path.group, sprintf("fafb_%s_meta.csv", group)))
  #if (nrow(manc.meta.chosen) > 0)
  #  readr::write_csv(manc.meta.chosen, file.path(save.path.group, sprintf("manc_%s_meta.csv", group)))
  if (nrow(malecns.meta.chosen) > 0)
    readr::write_csv(malecns.meta.chosen, file.path(save.path.group, sprintf("malecns_%s_meta.csv", group)))

  # Edgelists
  if (nrow(banc.el.chosen) > 0)
    arrow::write_feather(banc.el.chosen, file.path(save.path.group, sprintf("banc_%s_edgelist.feather", group)))
  if (nrow(fafb.el.chosen) > 0)
    arrow::write_feather(fafb.el.chosen, file.path(save.path.group, sprintf("fafb_%s_edgelist.feather", group)))
  #if (nrow(manc.el.chosen) > 0)
  #  arrow::write_feather(manc.el.chosen, file.path(save.path.group, sprintf("manc_%s_edgelist.feather", group)))
  if (nrow(malecns.el.chosen) > 0)
    arrow::write_feather(malecns.el.chosen, file.path(save.path.group, sprintf("malecns_%s_edgelist.feather", group)))

  message(sprintf("  Saved to %s", save.path.group))
}

message("\n=== Done ===")
