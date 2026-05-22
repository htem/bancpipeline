#' banc-data — Assemble the versioned BANC release: meta, metrics, synapses, edgelists.
#'
#' @section Reads:
#'   - `banc_<ver>_synapses_<src>.parquet`, `banc_<ver>_edgelist_simple_<src>.feather`
#'   - `banc_metrics.feather`, NT-prediction feathers, per-neuron split CSVs
#'   - SeaTable `banc_meta`
#'
#' @section Writes (under versioned save path):
#'   - `banc_<ver>_meta.feather`, `banc_<ver>_metrics.feather`
#'   - `banc_<ver>_synapses_{v2,v3}_enriched.parquet` (v2: NT+compartment; v3: CAVE root_ids, no NT)
#'   - `banc_<ver>_edgelist_split.feather`, L2 skeletons, region cut-outs, neuropil meshes
#'
#' @section CLI:
#'   --source {v2,v3}
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   BANC-project/R/startup/banc-meta.R, banc-meta-live.R (meta feather → dated parquet
#'   snapshot consumed by every figure); banc-edgelist.R (loads simple edgelist via
#'   `construct_path`).
#'
#' @section Schema:
#'   banc_888_meta.md, banc_888_metrics.md, banc_888_synapses_v2_enriched.md,
#'   banc_888_edgelist_split_v2.md.
#'
#' @section Paper:
#'   Methods §"Annotation taxonomy".

###########################################################
### Compile versioned BANC data for distribution
###
### Consolidates banc-build.R, banc-data.R, and the data
### sections of banc-sjcabs.R into a single script.
###
### Dependencies (must run before this script):
###   banc-calculate-connectivity.R → synapses parquet,
###     simple edgelist, basic meta
###   banc-calculate-neuropil-inclusion.R → neuropil parquet
###   banc-calculate-ntpred.R → NT prediction feathers
###   banc-calculate-split.R → per-neuron synapse CSVs
###   banc-update-metrics.R → combined metrics feather
###
### Outputs (in versioned save.path):
###   banc_{version}_meta.feather
###   banc_{version}_metrics.feather
###   banc_{version}_synapses_v2_enriched.parquet  (--source v2; NT predictions + compartment labels)
###   banc_{version}_synapses_v3_enriched.parquet  (--source v3; spatial + CAVE root_ids; no NT)
###   banc_{version}_edgelist_simple_{v2|v3}.feather
###   banc_{version}_edgelist_split.feather
###   banc_banc_space_l2_swc/   (L2 skeletons)
###   {cut_out}/                 (7 brain regions)
###   obj/                       (neuropil meshes)
###
### Source for the synapse parquet AND simple edgelist is chosen via --source
### v2|v3 CLI arg or BANC_SYN_SOURCE env var; defaults to
### banc.synapse.source.default.
###
### Usage: Rscript banc/share/banc-data.R [--source v2|v3]
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: compiling versioned data ###")
t_start <- Sys.time()
bancr::choose_banc()

#######################
### CONFIGURATION   ###
#######################

# BANC export version — must match banc-calculate-connectivity.R
# banc.version set in banc-startup.R

# Versioned output directory
save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                        paste0("banc_", banc.version))
dir.create(save.path, recursive = TRUE, showWarnings = FALSE)

# Synapse size threshold for the enriched parquet. Matches
# banc-calculate-connectivity.R (=5) so the enriched table and the
# simple/split edgelists share the same per-synapse inclusion criterion.
# (Previously 2 — produced an enriched parquet that was a superset of the
# edgelists; consumers had to guess which to trust.)
banc.size.threshold <- 5

# Synapse source for Sections 3 + 4 — chosen via --source v2|v3 CLI arg, then
# BANC_SYN_SOURCE env var, then banc.synapse.source.default. Lifted up here
# so Section 3 (synapse parquet) can branch on it; Section 4 (edgelist) reuses
# the same value.
.syn_source <- {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == "--source")
  if (length(i) == 1 && length(args) >= i + 1) {
    tolower(args[i + 1])
  } else {
    env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
    if (!is.na(env) && nzchar(env)) tolower(env)
    else if (exists("banc.synapse.source.default")) tolower(banc.synapse.source.default)
    else "v3"
  }
}
if (!.syn_source %in% c("v2", "v3")) {
  stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", .syn_source))
}
message(sprintf("synapse source: %s", .syn_source))

# Per-synapse NT prediction parquet — paired with .syn_source
nt.pred.path <- switch(.syn_source,
  v2 = "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v2/banc_nt_prediction_w_sizethresh_5_11102025.parquet",
  v3 = "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v3/banc_nt_prediction_v3_w_sizethresh_10_05042026.parquet"
)

# Helper: extract leading letters for seed columns
extract_three_letters <- function(text) {
  sapply(text, function(t) {
    three_letters <- stringr::str_extract(t, "^[A-Za-z]{3}")
    if (!is.na(three_letters)) return(three_letters)
    two_letters <- stringr::str_extract(t, "^[A-Za-z]{2}")
    if (!is.na(two_letters)) return(two_letters)
    stringr::str_extract(t, "^[A-Za-z]{1}")
  })
}

version_id_col <- paste0("banc_", banc.version, "_id")

##############################
### SECTION 1: META        ###
##############################

message("\n=== Section 1: Meta ===")

# # # Query SeaTable
banc.meta.orig <- banctable_query() 
  #  %>%
  # dplyr::filter(!grepl("glia|trachea|not_a_neuron|merge|orphan|tadpole", super_class),
  #               !grepl("GLIA|TRACHEA|NOT_A_NEURON|DEBRIS|MERGE|TADPOLE", status),
  #               proofread=="TRUE"|roughly_proofread=="TRUE")
banc.meta <- banc.meta.orig

# Convert supervoxel IDs → version-specific root IDs
message(sprintf("Converting %d supervoxel IDs to v%s root IDs...",
                nrow(banc.meta), banc.version))
banc.meta[[version_id_col]] <- banc_rootid(banc.meta$supervoxel_id,
                                            version = banc.version)
banc.meta <- banc.meta %>%
  dplyr::mutate(!!version_id_col := ifelse(
    is.na(.data[[version_id_col]]) | .data[[version_id_col]] == "0",
    root_id, .data[[version_id_col]])) %>%
  dplyr::filter(!is.na(.data[[version_id_col]])|!is.na(.data[[version_id_col]])=="0")

# # Get CAVE functions table for cell_function enrichment
# cns.functions <- banctable_query(sql = "SELECT * FROM functions") %>%
#   dplyr::select(-starts_with("_")) %>%
#   dplyr::filter(!is.na(cell_type))
# cns.functions.singular <- cns.functions %>%
#   dplyr::mutate(cell_function_cave = dplyr::case_when(
#     !is.na(response) & response != "" ~ response,
#     !is.na(behaviour) & behaviour != "" ~ behaviour,
#     !is.na(valence) & valence != "" ~ valence,
#     !is.na(modality) & modality != "" ~ modality,
#     TRUE ~ NA
#   )) %>%
#   dplyr::distinct(cell_type, .keep_all = TRUE) %>%
#   dplyr::distinct(cell_type, cell_function_cave)

# # Join CAVE functions (prefer existing SeaTable cell_function)
# banc.meta <- banc.meta %>%
#   dplyr::left_join(cns.functions.singular, by = "cell_type") %>%
#   dplyr::mutate(cell_function = dplyr::case_when(
#     !is.na(cell_function) & cell_function != "" ~ cell_function,
#     TRUE ~ cell_function_cave
#   )) %>%
#   dplyr::select(-cell_function_cave) %>%
#   dplyr::mutate(cell_function = gsub("\\/|\\,", "_", cell_function))

# # Join per-metric feather (pipeline-computed metrics)
# metrics_file <- file.path(banc.save.path, "banc_metrics.feather")
# if (file.exists(metrics_file)) {
#   metrics <- arrow::read_feather(metrics_file) %>%
#     dplyr::distinct(root_id, .keep_all = TRUE) %>%
#     dplyr::select(dplyr::any_of(c("root_id",
#       "l2_nodes", "l2_cable_length_um",
#       "input_connections", "output_connections",
#       "input_side_index", "output_side_index",
#       "mitochondria", "mitochondria_volume",
#       "pd_width", "segregation_index", "volume_nm3")))
#   # Remove metric columns already in banc.meta before joining
#   overlap_cols <- setdiff(intersect(colnames(metrics), colnames(banc.meta)), "root_id")
#   if (length(overlap_cols)) banc.meta <- banc.meta %>% dplyr::select(-dplyr::all_of(overlap_cols))
#   banc.meta <- banc.meta %>%
#     dplyr::left_join(metrics, by = "root_id")
#   message(sprintf("  Joined pipeline metrics (%d metric columns)", ncol(metrics) - 1))
# } else {
#   message("  banc_metrics.feather not found — skipping metric join")
# }

# # Join per-neuron NT predictions
# ntpred_file <- file.path(banc.save.path, "banc_ntpred.feather")
# if (file.exists(ntpred_file)) {
#   ntpred <- arrow::read_feather(ntpred_file) %>%
#     dplyr::distinct(root_id, .keep_all = TRUE) %>%
#     dplyr::select(root_id, neurotransmitter_predicted, neurotransmitter_score)
#   # Only join if columns not already present from SeaTable
#   nt_overlap <- intersect(c("neurotransmitter_predicted", "neurotransmitter_score"),
#                           colnames(banc.meta))
#   if (length(nt_overlap)) banc.meta <- banc.meta %>% dplyr::select(-dplyr::all_of(nt_overlap))
#   banc.meta <- banc.meta %>%
#     dplyr::left_join(ntpred, by = "root_id")
#   message(sprintf("  Joined NT predictions (%d neurons with predictions)",
#                   sum(!is.na(ntpred$neurotransmitter_predicted))))
# }

# # Enrich from franken.meta (MANC/FAFB cross-dataset annotations via SeaTable)
# message("  Enriching from franken.meta...")
# franken.meta <- tryCatch({
#   bancr::franken_meta()
# }, error = function(e) {
#   message("  WARNING: franken.meta not available: ", e$message)
#   NULL
# })
# if (!is.null(franken.meta) && nrow(franken.meta) > 0) {
#   enrich_cols <- c("neurotransmitter_verified", "neuropeptide_verified", "neuromere")
# 
#   # Join via manc_match → manc_id
#   fm_manc <- franken.meta %>%
#     dplyr::filter(!is.na(manc_id)) %>%
#     dplyr::select(manc_id, cell_type, dplyr::any_of(enrich_cols)) %>%
#     dplyr::rename(fm_manc_cell_type = cell_type) %>%
#     dplyr::distinct(manc_id, .keep_all = TRUE)
# 
#   # Join via fafb_match → fafb_id
#   fm_fafb <- franken.meta %>%
#     dplyr::filter(!is.na(fafb_id)) %>%
#     dplyr::select(fafb_id, cell_type, dplyr::any_of(enrich_cols)) %>%
#     dplyr::rename(fm_fafb_cell_type = cell_type) %>%
#     dplyr::distinct(fafb_id, .keep_all = TRUE)
# 
#   # Remove existing columns to avoid .x/.y duplication
#   for (.fc in enrich_cols) {
#     if (.fc %in% colnames(banc.meta)) banc.meta[[.fc]] <- NULL
#   }
# 
#   banc.meta <- banc.meta %>%
#     dplyr::left_join(fm_manc, by = c("manc_match" = "manc_id"), suffix = c("", ".manc")) %>%
#     dplyr::left_join(fm_fafb, by = c("fafb_match" = "fafb_id"), suffix = c("", ".fafb"))
# 
#   # Coalesce enriched columns from both join paths
#   for (.fc in enrich_cols) {
#     manc_col <- paste0(.fc, ".manc")
#     fafb_col <- paste0(.fc, ".fafb")
#     banc.meta[[.fc]] <- dplyr::coalesce(
#       if (.fc %in% colnames(banc.meta)) banc.meta[[.fc]] else NA_character_,
#       if (manc_col %in% colnames(banc.meta)) banc.meta[[manc_col]] else NA_character_,
#       if (fafb_col %in% colnames(banc.meta)) banc.meta[[fafb_col]] else NA_character_
#     )
#     banc.meta[[manc_col]] <- NULL
#     banc.meta[[fafb_col]] <- NULL
#   }
# 
#   # Build composite_cell_type: franken cell_type is manc_cell_type or fafb_cell_type
#   banc.meta <- banc.meta %>%
#     dplyr::mutate(composite_cell_type = dplyr::case_when(
#       is.na(cell_type) ~ root_id,
#       grepl("ascending", super_class) & !is.na(fm_manc_cell_type) ~ fm_manc_cell_type,
#       grepl("descending", super_class) & !is.na(fm_fafb_cell_type) ~ fm_fafb_cell_type,
#       TRUE ~ cell_type
#     )) %>%
#     dplyr::select(-dplyr::any_of(c("fm_manc_cell_type", "fm_fafb_cell_type")))
# 
#   message(sprintf("  Joined franken.meta: %d via manc_match, %d via fafb_match",
#                   sum(!is.na(banc.meta$manc_match)), sum(!is.na(banc.meta$fafb_match))))
# } else {
#   # Fallback: composite_cell_type without franken.meta
#   banc.meta <- banc.meta %>%
#     dplyr::mutate(composite_cell_type = dplyr::case_when(
#       is.na(cell_type) ~ root_id,
#       TRUE ~ cell_type
#     ))
# }

# Compute derived columns and seed columns
# banc.meta <- banc.meta %>%
#   dplyr::mutate(
#     cell_function = dplyr::case_when(
#       is.na(cell_function) | cell_function == "" ~ "unknown",
#       TRUE ~ cell_function),
#     body_part_sensory = dplyr::case_when(
#       grepl("head|^frontal|^frontoorbital|^orbital|^interocellar|^vibrissa|^interommatidial|^occipital_dorsal|^occipital_ventral|^postorbital_dorsal|^postorbital_ventral|^vertical|^postocellar|^supracervical", body_part_sensory) ~ "head",
#       is.na(body_part_sensory) | body_part_sensory == "" ~ "unknown",
#       TRUE ~ body_part_sensory),
#     cell_function_detailed = dplyr::case_when(
#       !is.na(cell_function_detailed) & !is.na(cell_function) ~ paste0(cell_function, "_", cell_function_detailed),
#       !is.na(cell_function_detailed) & !is.na(cell_function) ~ cell_function,
#       is.na(cell_function) ~ "unknown",
#       TRUE ~ cell_function_detailed),
#     cell_function_nerve = dplyr::case_when(
#       !is.na(cell_function_detailed) ~ cell_function_detailed,
#       !is.na(nerve) ~ gsub("_r$|_l$|_left$|_right$|_R$|_L$|^right_|^left_", "", nerve),
#       TRUE ~ cell_function),
#     .is_sensory = grepl("sensory", super_class) | (!is.na(body_part_sensory) & body_part_sensory != "" & body_part_sensory != "unknown"),
#     seed_00 = dplyr::case_when(
#       is.na(side) | !side %in% c("left", "right") ~ NA,
#       .is_sensory ~ cell_type,
#       grepl("central_complex_output|mushroom_body_output_neuron", cell_class) ~ cell_type,
#       grepl("visual_projection", super_class) ~ cell_type,
#       grepl("ascending|descending", super_class) ~ cell_type,
#       TRUE ~ NA),
#     seed_01 = dplyr::case_when(
#       grepl("efferent", flow) ~ NA,
#       is.na(cell_function) | cell_function == "unknown" ~ NA,
#       TRUE ~ paste0(super_class, "_", cell_function)),
#     seed_02 = dplyr::case_when(
#       is.na(side) | !side %in% c("left", "right") | is.na(cell_sub_class) ~ NA,
#       .is_sensory ~ cell_sub_class,
#       TRUE ~ NA),
#     seed_03 = dplyr::case_when(
#       !(.is_sensory | grepl("visual_projection", super_class)) | is.na(cell_function_detailed) ~ NA,
#       grepl("visual_projection", super_class) ~ paste0("visual_projection_", cell_function_detailed),
#       .is_sensory ~ paste0(body_part_sensory, "_", cell_function_detailed),
#       TRUE ~ NA),
#     seed_04 = dplyr::case_when(
#       is.na(side) | !side %in% c("left", "right") | is.na(cell_sub_class) ~ NA,
#       .is_sensory ~ paste0(cell_sub_class, "_", side),
#       TRUE ~ NA),
#     seed_05 = dplyr::case_when(
#       grepl("central_complex_output", cell_class) & !is.na(cell_sub_class) & cell_sub_class != "" ~ cell_sub_class,
#       grepl("ascending|descending", super_class) & !is.na(cluster) & cluster != "" ~ cluster,
#       grepl("mushroom_body_output_neuron", cell_class) & !is.na(cell_class) & cell_class != "" ~ cell_class,
#       grepl("visual_projection", super_class) & !is.na(cell_type) & cell_type != "" ~ extract_three_letters(cell_type),
#       TRUE ~ NA),
#     seed_06 = dplyr::case_when(
#       is.na(side) | !side %in% c("left", "right", "midline", "center") ~ NA,
#       grepl("central_complex_output", cell_class) & !is.na(cell_sub_class) & cell_sub_class != "" ~ paste0(cell_sub_class, "_", side),
#       grepl("ascending|descending", super_class) & !is.na(cluster) & cluster != "" ~ paste0(cluster, "_", side),
#       grepl("mushroom_body_output_neuron", cell_class) & !is.na(cell_class) & cell_class != "" ~ paste0(cell_class, "_", side),
#       grepl("visual_projection", super_class) & !is.na(cell_type) & cell_type != "" ~ paste0(extract_three_letters(cell_type), "_", side),
#       is.na(nerve) | nerve == "" ~ NA,
#       TRUE ~ NA),
#     seed_07 = dplyr::case_when(
#       grepl("central_complex_output|mushroom_body_output_neuron", cell_class) ~ cell_type,
#       grepl("^EPG|^EL", cell_type) ~ cell_type,
#       grepl("sensory_ascending", super_class) & grepl("SA", cell_type) ~ cell_type,
#       grepl("ascending|descending", cell_class) & !is.na(cell_type) & cell_type != "" ~ cell_type,
#       grepl("visual_projection", super_class) ~ cell_type,
#       TRUE ~ NA),
#     seed_08 = dplyr::case_when(
#       !.is_sensory ~ NA,
#       grepl("proprio|tactile|contract|vib", cell_function) ~ paste0(cell_function, "_", body_part_sensory, "_", cell_function_nerve, "_", side),
#       TRUE ~ NA),
#     seed_09 = dplyr::case_when(
#       grepl("efferent", flow) ~ NA,
#       is.na(peripheral_target_type) ~ NA,
#       TRUE ~ peripheral_target_type),
#     seed_10 = dplyr::case_when(
#       !.is_sensory | is.na(peripheral_target_type) ~ NA,
#       TRUE ~ paste0(body_part_sensory, "_", peripheral_target_type)),
#     seed_11 = dplyr::case_when(
#       !is.na(cluster) ~ cluster,
#       !is.na(cns_network) ~ cns_network,
#       TRUE ~ NA),
#     seed_12 = dplyr::case_when(
#       grepl("ascending|descending", super_class) ~ paste0(cell_type, "_", root_id),
#       TRUE ~ NA),
#     seed_13 = root_id,
#     seed_14 = cns_network
#   ) %>%
#   # Sanitize seed columns (replace semicolons)
#   dplyr::mutate(dplyr::across(starts_with("seed_"), ~ gsub(";", "_", .x)))

# Select and deduplicate
banc.meta <- banc.meta %>%
  dplyr::distinct(.data[[version_id_col]], .keep_all = TRUE) %>%
  dplyr::select(dplyr::any_of(c(
    version_id_col, "root_id", "supervoxel_id", "position", "root_626", "root_850", "root_888", "root_890", "nucleus_id",
    "proofread", "roughly_proofread", "status",
    "side", "root_position", "root_position_nm", "root_region", "region",
    "hemilineage", "nerve", "tract", "neuromere", "flow",
    "super_class", "cell_class", "cell_sub_class", "cell_type",
    "fafb_cell_type", "fafb_alignment_cell_type", "manc_cell_type", "malecns_cell_type", "hemibrain_cell_type", "fanc_cell_type",
    "fafb_match","manc_match", "malecns_match", "hemibrain_match", "fanc_match",
    "fafb_nblast_match", "fafn_alignment_match", "manc_nblast_match", "hemibrain_nblast_match", "fanc_nblast_match","malecns_nblast_match",
    "sexually_dimorphic", "cluster", "manual_cluster", "super_cluster", "cns_network",
    "body_part_sensory", "body_part_effector", "peripheral_target_type",
    "cell_function", "cell_function_detailed",
    "composite_cell_type",
    "neurotransmitter_predicted", "neurotransmitter_score",
    "neurotransmitter_verified", "neuropeptide_verified",
    "l2_nodes", "l2_cable_length_um", "volume_nm3",
    "input_connections", "output_connections",
    "input_side_index", "output_side_index",
    "mitochondria", "mitochondria_volume",
    "pd_width", "segregation_index",
    paste0("seed_", sprintf("%02d", 0:14))
  ))) %>%
  as.data.frame()

banc.meta[banc.meta == ""] <- NA
root.ids <- banc.meta[[version_id_col]]

# Save meta
meta_file <- file.path(save.path,
                        sprintf("banc_%s_meta.feather", banc.version))
arrow::write_feather(banc.meta, meta_file)
message(sprintf("Meta saved: %s (%d neurons, %d columns)",
                meta_file, nrow(banc.meta), ncol(banc.meta)))

system(sprintf("gsutil cp %s gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/",
               meta_file, banc.version))

##############################
### SECTION 2: METRICS     ###
##############################

message("\n=== Section 2: Metrics feather ===")

metric_cols <- c("region", "side",
                 "l2_nodes", "l2_cable_length_um",
                 "input_connections", "output_connections",
                 "input_side_index", "output_side_index",
                 "mitochondria", "mitochondria_volume",
                 "pd_width", "segregation_index", "volume_nm3")

# Prefer the combined per-metric feather written by banc-update-metrics.R.
# It's keyed on current root_id and is denser than the SeaTable round-trip
# (SeaTable lags whenever a metric was just (re)calculated for new neurons).
combined_metrics_file <- file.path(banc.save.path, "banc_metrics.feather")
if (file.exists(combined_metrics_file)) {
  message(sprintf("  Reading combined metrics feather: %s", combined_metrics_file))
  banc.metrics.combined <- arrow::read_feather(combined_metrics_file) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  avail_metric_cols <- intersect(metric_cols, colnames(banc.metrics.combined))

  banc.metrics <- banc.meta %>%
    dplyr::select(dplyr::all_of(c(version_id_col, "root_id"))) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::left_join(
      banc.metrics.combined %>% dplyr::select(dplyr::all_of(c("root_id", avail_metric_cols))),
      by = "root_id") %>%
    dplyr::select(-root_id) %>%
    dplyr::distinct(.data[[version_id_col]], .keep_all = TRUE)
} else {
  message("  Combined metrics feather not found — falling back to SeaTable columns")
  avail_metric_cols <- intersect(metric_cols, colnames(banc.meta))
  banc.metrics <- banc.meta %>%
    dplyr::select(dplyr::all_of(c(version_id_col, avail_metric_cols)))
}

if (length(avail_metric_cols)) {
  metrics_out <- file.path(save.path,
                            sprintf("banc_%s_metrics.feather", banc.version))
  arrow::write_feather(banc.metrics, metrics_out)
  n_dense <- sum(rowSums(!is.na(banc.metrics[, avail_metric_cols, drop = FALSE])) > 0)
  message(sprintf("Metrics saved: %s (%d neurons, %d metric columns, %d/%d with at least one non-NA metric)",
                  metrics_out, nrow(banc.metrics), length(avail_metric_cols),
                  n_dense, nrow(banc.metrics)))
} else {
  message("  No metric columns found — skipping metrics feather")
}

##############################
### SECTION 3: SYNAPSES    ###
##############################

message(sprintf("\n=== Section 3: Synapse table (--source %s) ===", .syn_source))

if (.syn_source == "v3") {
  # v3 path: trust CAVE-curated root_ids from the GCS export, take neuropil /
  # region / side from the local spatial parquet (whose own root_id columns are
  # corrupt and ignored). No NT predictions and no compartment labels.
  v3_export_local <- file.path("/n/data1/hms/neurobio/wilson/banc/synapses_v3",
                               "cache", sprintf("v%s_inputs", banc.version),
    "synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet")
  v3_spatial <- file.path("/n/data1/hms/neurobio/wilson/banc/synapses_v3",
                          sprintf("banc_%s_synapses_v3.parquet", banc.version))
  if (!file.exists(v3_export_local)) {
    message(sprintf("v3 GCS export not cached at %s — run banc-calculate-connectivity.R --source v3 first",
                    v3_export_local))
    return(invisible())
  }
  if (!file.exists(v3_spatial)) {
    message(sprintf("v3 spatial parquet not found at %s — run banc-synapses-v3-optimised.R first",
                    v3_spatial))
    return(invisible())
  }

  message("Reading v3 GCS export (trustworthy root_ids)...")
  v3_export <- arrow::read_parquet(v3_export_local,
    col_select = c("id", "size", "pre_root_id", "post_root_id"))
  v3_export$id <- as.character(as.numeric(v3_export$id))
  v3_export$pre_root_id <- as.character(v3_export$pre_root_id)
  v3_export$post_root_id <- as.character(v3_export$post_root_id)

  message("Reading v3 spatial parquet (neuropil / region / side + scores + endpoint coords)...")
  # Widened from 7 to 15 cols: include mean_score + median_score (v3's quality
  # proxies) and presyn_/postsyn_xyz (cleft endpoints, V3 voxel = 16x16x45 nm).
  # The 15-col one-shot read has historically tripped the SLURM-only arrow C++
  # URI bug ("Invalid: Unrecognized filesystem type in URI: file:///_") on this
  # 15 GB parquet — fall back to two narrow reads joined on syn_id.
  v3_sp_cols_meta    <- c("syn_id", "neuropil", "region", "side",
                           "X", "Y", "Z", "mean_score", "median_score")
  v3_sp_cols_spatial <- c("syn_id",
                           "presyn_x", "presyn_y", "presyn_z",
                           "postsyn_x", "postsyn_y", "postsyn_z")
  v3_sp <- tryCatch(
    arrow::read_parquet(v3_spatial,
      col_select = dplyr::all_of(c(v3_sp_cols_meta, v3_sp_cols_spatial[-1]))),
    error = function(e) {
      message(sprintf("  one-pass 15-col read failed (%s); 2-pass narrow read",
                      conditionMessage(e)))
      a <- arrow::read_parquet(v3_spatial, col_select = dplyr::all_of(v3_sp_cols_meta))
      b <- arrow::read_parquet(v3_spatial, col_select = dplyr::all_of(v3_sp_cols_spatial))
      dplyr::inner_join(a, b, by = "syn_id")
    }
  )
  v3_sp <- dplyr::rename(v3_sp, id = "syn_id")
  v3_sp$id <- as.character(v3_sp$id)
  v3_sp <- v3_sp[!duplicated(v3_sp$id), , drop = FALSE]

  message(sprintf("  GCS export: %s rows | spatial: %s rows",
                  format(nrow(v3_export), big.mark = ","),
                  format(nrow(v3_sp), big.mark = ",")))

  banc.syns <- dplyr::inner_join(v3_export, v3_sp, by = "id")
  n_join <- nrow(banc.syns)
  if (n_join < 0.99 * nrow(v3_export)) {
    stop(sprintf("v3 spatial join coverage too low: %.4f%% (%s of %s)",
                 100 * n_join / nrow(v3_export),
                 format(n_join, big.mark = ","),
                 format(nrow(v3_export), big.mark = ",")))
  }
  message(sprintf("  Joined: %s rows (%.4f%% of GCS export)",
                  format(n_join, big.mark = ","),
                  100 * n_join / nrow(v3_export)))
  rm(v3_export, v3_sp); gc()

  # Filter: size threshold, no autapses, at least one end in our neuron set.
  # NOTE: synapses labelled "outside" (or with NA / empty neuropil) are
  # deliberately retained — neuropil is not part of the inclusion criterion.
  banc.syns <- banc.syns %>%
    dplyr::filter(size >= banc.size.threshold,
                  pre_root_id != post_root_id,
                  pre_root_id %in% root.ids | post_root_id %in% root.ids) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  message(sprintf("  Filtered synapses: %s rows", format(nrow(banc.syns), big.mark = ",")))

  message("  v3 compartment labels: derived via spatial NN on v2 split CSVs in Section 5")

} else {
  # v2 path: read neuropil-labelled synapses from banc-calculate-neuropil-inclusion.R
  neuropil_parquet <- file.path(banc.connectivity.save.path,
                                 sprintf("banc_%s_synapses_v2_neuropils.parquet", banc.version))
  raw_parquet <- file.path(banc.connectivity.save.path,
                            sprintf("banc_%s_synapses.parquet", banc.version))

  if (file.exists(neuropil_parquet)) {
    message("Reading neuropil-labelled synapse parquet...")
    banc.syns <- arrow::read_parquet(neuropil_parquet)
  } else if (file.exists(raw_parquet)) {
    message("Neuropil parquet not found — using raw synapse parquet")
    banc.syns <- arrow::read_parquet(raw_parquet)
  } else {
    message("No synapse parquet found — run banc-calculate-connectivity.R first")
    return(invisible())
  }

  message(sprintf("  Raw synapses: %s rows", format(nrow(banc.syns), big.mark = ",")))

  # Filter: size threshold, no autapses, at least one end in our neuron set.
  # NOTE: synapses labelled "outside" (or with NA / empty neuropil) are
  # deliberately retained — neuropil is not part of the inclusion criterion.
  banc.syns <- banc.syns %>%
    dplyr::filter(size >= banc.size.threshold,
                  pre_root_id != post_root_id,
                  pre_root_id %in% root.ids | post_root_id %in% root.ids) %>%
    dplyr::distinct(id, .keep_all = TRUE)

  message(sprintf("  Filtered synapses: %s rows", format(nrow(banc.syns), big.mark = ",")))
}

# Join per-synapse NT predictions (works for both v2 and v3 sources;
# nt.pred.path is selected to match .syn_source above)
if (file.exists(nt.pred.path)) {
  message("  Joining per-synapse NT predictions...")
  banc.nt <- arrow::read_parquet(nt.pred.path)
  banc.syns <- banc.syns %>%
    dplyr::mutate(id = as.character(id)) %>%
    dplyr::left_join(banc.nt %>%
                       dplyr::mutate(id = as.character(id)) %>%
                       dplyr::rename(syn_top_nt = predicted_nt,
                                     syn_top_p = probability),
                     by = "id")
  n_with_nt <- sum(!is.na(banc.syns$syn_top_nt))
  message(sprintf("  NT predictions joined: %s/%s synapses (%.1f%%)",
                  format(n_with_nt, big.mark = ","),
                  format(nrow(banc.syns), big.mark = ","),
                  100 * n_with_nt / nrow(banc.syns)))
} else {
  message(sprintf("  NT prediction parquet not found at %s — skipping", nt.pred.path))
}

# v3 compartment labels via per-neuron spatial NN against v2 split CSVs.
# Persisted into enriched parquet so consumers (incl. Section 5 split edgelist)
# can reuse them instead of re-running the NN.
if (.syn_source == "v3") {
  message("  v3 spatial-NN labels: loading v2 split CSVs...")
  .syn_folder <- file.path(banc.l2split.save.path, "synapses")
  if (!dir.exists(.syn_folder) || length(list.files(.syn_folder)) == 0) {
    message("    No v2 split CSVs found — pre_label/post_label will be 'unknown'.")
    banc.syns$pre_label  <- "unknown"
    banc.syns$post_label <- "unknown"
  } else {
    .bc_syn <- read_synapse_csvs(.syn_folder)
    colnames(.bc_syn) <- snakecase::to_snake_case(colnames(.bc_syn))
    .bc_syn <- .bc_syn[, !duplicated(colnames(.bc_syn))]
    if ("pre_id"  %in% colnames(.bc_syn)) .bc_syn$pre_id  <- as.character(.bc_syn$pre_id)
    if ("post_id" %in% colnames(.bc_syn)) .bc_syn$post_id <- as.character(.bc_syn$post_id)
    .bc_syn <- .bc_syn %>%
      dplyr::filter(size >= banc.size.threshold) %>%
      dplyr::select(dplyr::any_of(c("x","y","z","prepost","pre_id","post_id","label"))) %>%
      dplyr::mutate(label = hemibrainr:::standard_compartments(label))
    v2_pre  <- .bc_syn %>% dplyr::filter(prepost == 0) %>%
      dplyr::select(neuron_id = pre_id,  x, y, z, label) %>%
      dplyr::filter(!is.na(neuron_id), !is.na(label), !is.na(x), !is.na(y), !is.na(z))
    v2_post <- .bc_syn %>% dplyr::filter(prepost == 1) %>%
      dplyr::select(neuron_id = post_id, x, y, z, label) %>%
      dplyr::filter(!is.na(neuron_id), !is.na(label), !is.na(x), !is.na(y), !is.na(z))
    rm(.bc_syn); gc()
    message(sprintf("  v2 labeled points: %s pre, %s post",
                    format(nrow(v2_pre),  big.mark = ","),
                    format(nrow(v2_post), big.mark = ",")))

    # Per-neuron NN. Operate on plain matrices/vectors extracted ONCE.
    # Crucially: split BOTH v3 and v2 by index (cheap vector split), then
    # slice extracted matrices inside the loop — base R's
    # split.data.frame on 100M+ rows is ~1000× slower than index split.
    V3_VOXEL_NM <- c(16, 16, 45)

    assign_label_by_nn_v3 <- function(v3_xyz_mat, v3_neuron_vec,
                                       v2_xyz_mat, v2_label_vec, v2_neuron_vec,
                                       side_name) {
      message(sprintf("  Assigning %s_label via per-neuron NN...", side_name))
      n <- nrow(v3_xyz_mat)
      out <- rep(NA_character_, n)

      # Build index lists. v2 has ~165k neurons; v3 may have millions of
      # one-side-only "garbage" root_ids that won't have v2 entries —
      # iterate ONLY over the intersection to avoid the O(n) [[ lookup
      # cost across the full v3 universe.
      message(sprintf("    Splitting v2 index by neuron (n=%s rows, ~%s neurons)...",
                      format(length(v2_neuron_vec), big.mark = ","),
                      format(length(unique(v2_neuron_vec)), big.mark = ",")))
      v2_grp <- split(seq_along(v2_neuron_vec), v2_neuron_vec)
      v2_neurons <- names(v2_grp)

      message(sprintf("    Splitting v3 index by neuron (n=%s rows)...",
                      format(n, big.mark = ",")))
      v3_grp <- split(seq_len(n), v3_neuron_vec)
      v3_neurons <- names(v3_grp)

      common <- intersect(v3_neurons, v2_neurons)
      message(sprintf("    v3 neurons: %s | v2 neurons: %s | intersect: %s — iterating intersect only",
                      format(length(v3_neurons), big.mark = ","),
                      format(length(v2_neurons), big.mark = ","),
                      format(length(common), big.mark = ",")))

      n_total <- length(common); .n_hit <- 0L
      t0 <- Sys.time()
      for (i in seq_along(common)) {
        rid <- common[i]
        v2_idx <- v2_grp[[rid]]
        v3_idx <- v3_grp[[rid]]
        qry <- v3_xyz_mat[v3_idx, , drop = FALSE]
        dat <- v2_xyz_mat[v2_idx, , drop = FALSE]
        nn  <- tryCatch(nabor::knn(dat, qry, k = 1L), error = function(e) NULL)
        if (is.null(nn)) next
        out[v3_idx] <- v2_label_vec[v2_idx[nn$nn.idx[, 1]]]
        .n_hit <- .n_hit + length(v3_idx)
        if (i %% 5000L == 0L) {
          dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          message(sprintf("    %s: %d/%d neurons (%.0fs, %s labeled)",
                          side_name, i, n_total, dt,
                          format(.n_hit, big.mark = ",")))
        }
      }
      message(sprintf("  %s_label assigned for %s / %s v3 rows",
                      side_name, format(.n_hit, big.mark = ","),
                      format(n, big.mark = ",")))
      out
    }

    # Extract everything to plain matrices/vectors ONCE.
    pre_xyz <- cbind(banc.syns$presyn_x  * V3_VOXEL_NM[1],
                     banc.syns$presyn_y  * V3_VOXEL_NM[2],
                     banc.syns$presyn_z  * V3_VOXEL_NM[3])
    v2_pre_xyz <- as.matrix(v2_pre[, c("x","y","z"), drop = FALSE])
    v2_pre_lab <- v2_pre$label
    v2_pre_nid <- v2_pre$neuron_id
    rm(v2_pre); gc()
    banc.syns$pre_label <- assign_label_by_nn_v3(
      pre_xyz, banc.syns$pre_root_id,
      v2_pre_xyz, v2_pre_lab, v2_pre_nid, "pre")
    rm(pre_xyz, v2_pre_xyz, v2_pre_lab, v2_pre_nid); gc()

    post_xyz <- cbind(banc.syns$postsyn_x * V3_VOXEL_NM[1],
                      banc.syns$postsyn_y * V3_VOXEL_NM[2],
                      banc.syns$postsyn_z * V3_VOXEL_NM[3])
    v2_post_xyz <- as.matrix(v2_post[, c("x","y","z"), drop = FALSE])
    v2_post_lab <- v2_post$label
    v2_post_nid <- v2_post$neuron_id
    rm(v2_post); gc()
    banc.syns$post_label <- assign_label_by_nn_v3(
      post_xyz, banc.syns$post_root_id,
      v2_post_xyz, v2_post_lab, v2_post_nid, "post")
    rm(post_xyz, v2_post_xyz, v2_post_lab, v2_post_nid); gc()
  }

  # Afferent/efferent override (mirrors v2 path + Section 5)
  afferent_ids <- banc.meta[[version_id_col]][
    !is.na(banc.meta$flow) & banc.meta$flow == "afferent"]
  efferent_ids <- banc.meta[[version_id_col]][
    !is.na(banc.meta$flow) & banc.meta$flow == "efferent"]
  banc.syns <- banc.syns %>%
    dplyr::mutate(
      pre_label  = dplyr::case_when(
        pre_root_id  %in% afferent_ids ~ "axon",
        pre_root_id  %in% efferent_ids ~ "dendrite",
        is.na(pre_label) ~ "unknown",
        TRUE ~ pre_label),
      post_label = dplyr::case_when(
        post_root_id %in% afferent_ids ~ "axon",
        post_root_id %in% efferent_ids ~ "dendrite",
        is.na(post_label) ~ "unknown",
        TRUE ~ post_label))
  n_pre_lab  <- sum(banc.syns$pre_label  != "unknown", na.rm = TRUE)
  n_post_lab <- sum(banc.syns$post_label != "unknown", na.rm = TRUE)
  message(sprintf("  v3 labels persisted: pre %s/%s (%.1f%%), post %s/%s (%.1f%%)",
                  format(n_pre_lab,  big.mark = ","),
                  format(nrow(banc.syns), big.mark = ","),
                  100 * n_pre_lab  / nrow(banc.syns),
                  format(n_post_lab, big.mark = ","),
                  format(nrow(banc.syns), big.mark = ","),
                  100 * n_post_lab / nrow(banc.syns)))
}

# Compartment labels (v2 path only — split CSVs are v2-keyed)
if (.syn_source == "v2") {
# Join compartment labels from flow centrality splits (banc-calculate-split.R)
# Label mapping (hemibrainr:::standard_compartments):
#   0=unknown, 1=soma, 2=axon, 3=dendrite, 4=primary.dendrite, 7=primary.neurite
# Prefer detailed splits; fill missing with L2 splits.
message("  Joining compartment labels from split CSVs...")
detailed_split_syn_dir <- file.path(banc.split.save.path, "synapses")
l2_split_syn_dir <- file.path(banc.l2split.save.path, "synapses")

read_split_labels <- function(dir) {
  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(csvs)) return(data.table::data.table(connector_id = character(0),
                                                     Label = integer(0)))
  message(sprintf("    Reading %d split CSVs from %s...", length(csvs), basename(dir)))
  pb <- utils::txtProgressBar(max = length(csvs), style = 3)
  result <- data.table::rbindlist(lapply(seq_along(csvs), function(i) {
    utils::setTxtProgressBar(pb, i)
    tryCatch(
      data.table::fread(csvs[i], select = c("connector_id", "Label"),
                        colClasses = list(character = "connector_id", integer = "Label")),
      error = function(e) NULL
    )
  }))
  close(pb)
  result
}

label_lookup <- data.table::data.table(connector_id = character(0), Label = integer(0))

# Detailed splits first (higher priority)
if (dir.exists(detailed_split_syn_dir)) {
  detailed_labels <- read_split_labels(detailed_split_syn_dir)
  detailed_labels <- detailed_labels[!is.na(connector_id) & !is.na(Label)]
  detailed_labels <- detailed_labels[!duplicated(connector_id)]
  label_lookup <- detailed_labels
  message(sprintf("    Detailed splits: %s labels", format(nrow(label_lookup), big.mark = ",")))
}

# L2 splits to fill gaps
if (dir.exists(l2_split_syn_dir)) {
  l2_labels <- read_split_labels(l2_split_syn_dir)
  l2_labels <- l2_labels[!is.na(connector_id) & !is.na(Label)]
  l2_labels <- l2_labels[!duplicated(connector_id)]
  # Only keep L2 labels for connector_ids not already covered by detailed
  l2_labels <- l2_labels[!connector_id %in% label_lookup$connector_id]
  label_lookup <- data.table::rbindlist(list(label_lookup, l2_labels))
  message(sprintf("    L2 splits (gap-fill): %s labels", format(nrow(l2_labels), big.mark = ",")))
}

if (nrow(label_lookup) > 0) {
  data.table::setnames(label_lookup, "Label", "label")
  banc.syns <- banc.syns %>%
    dplyr::left_join(as.data.frame(label_lookup),
                     by = c("id" = "connector_id"))
  n_with_label <- sum(!is.na(banc.syns$label))
  message(sprintf("  Compartment labels joined: %s/%s synapses (%.1f%%)",
                  format(n_with_label, big.mark = ","),
                  format(nrow(banc.syns), big.mark = ","),
                  100 * n_with_label / nrow(banc.syns)))
} else {
  banc.syns$label <- NA_integer_
  message("  No split CSVs found — label column set to NA")
}

rm(label_lookup)
if (exists("detailed_labels")) rm(detailed_labels)
if (exists("l2_labels")) rm(l2_labels)

}  # end if (.syn_source == "v2") for compartment labels

# Save enriched synapse table (source-suffixed: _v2 or _v3)
synapses_file <- file.path(save.path,
                            sprintf("banc_%s_synapses_%s_enriched.parquet",
                                    banc.version, .syn_source))
write_connectome_data(banc.syns, synapses_file, format = "parquet")
message(sprintf("Synapse table saved: %s (%s rows)",
                synapses_file, format(nrow(banc.syns), big.mark = ",")))

##############################
### SECTION 4: SIMPLE ELIST ##
##############################

message("\n=== Section 4: Simple edgelist ===")

message(sprintf("  edgelist source: %s", .syn_source))

# Read from pipeline output if it exists. For v3 we deliberately bypass the
# legacy source feather and always build from banc.syns — the standalone v3
# connectivity script isn't re-run as part of the data push, so its source
# feather can sit stale (it shipped as 802 bytes / 0 rows on 2026-05-17).
edgelist_src <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_edgelist_simple_%s.feather",
                                   banc.version, .syn_source))
edgelist_file <- file.path(save.path,
                            sprintf("banc_%s_edgelist_simple_%s.feather",
                                    banc.version, .syn_source))
if (.syn_source != "v3" && file.exists(edgelist_src)) {
  banc.elist.simp <- arrow::read_feather(edgelist_src)
  arrow::write_feather(banc.elist.simp, edgelist_file)
  message(sprintf("Simple edgelist copied: %s (%s connections)",
                  edgelist_file, format(nrow(banc.elist.simp), big.mark = ",")))
} else {
  # Build from enriched synapses
  message("Simple edgelist not found — building from synapse table...")
  banc.elist.simp <- banc.syns %>%
    dplyr::group_by(pre_root_id, post_root_id) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(pre = pre_root_id, post = post_root_id) %>%
    dplyr::mutate(pre = as.character(pre), post = as.character(post)) %>%
    dplyr::group_by(post) %>%
    dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pre) %>%
    dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(norm = signif(count / post_count, 4)) %>%
    dplyr::select(pre, post, count, norm, post_count, pre_count)
  arrow::write_feather(banc.elist.simp, edgelist_file)
  message(sprintf("Simple edgelist saved: %s (%s connections)",
                  edgelist_file, format(nrow(banc.elist.simp), big.mark = ",")))
}

##############################
### SECTION 5: SPLIT ELIST ###
##############################

message("\n=== Section 5: Split edgelist ===")

# Split CSVs (banc.l2split.save.path/synapses/) are keyed on v2 connector_ids
# and carry v2 synapse positions with compartment labels from flow centrality.
#   --source v2  → connector_id join (exact)
#   --source v3  → per-neuron spatial NN onto v2 split points (see v3 branch
#                  below); inherits Label from the closest same-neuron v2
#                  pre/post point. Falls back to "unknown" if the neuron has
#                  no v2 split CSV at all.
synapses.folder <- file.path(banc.l2split.save.path, "synapses")
split_elist_file <- file.path(save.path,
                               sprintf("banc_%s_edgelist_split_v2.feather",
                                       banc.version))

if (!dir.exists(synapses.folder) || length(list.files(synapses.folder)) == 0) {
  message("Split pipeline synapse CSVs not found — skipping split edgelist")
} else if (.syn_source == "v2") {
  message("Reading split pipeline synapse CSVs...")
  bc.synapses <- read_synapse_csvs(synapses.folder)
  colnames(bc.synapses) <- snakecase::to_snake_case(colnames(bc.synapses))
  bc.synapses <- bc.synapses[, !duplicated(colnames(bc.synapses))]
  # Ensure connector_id is character (fread may infer integer64 for large IDs)
  if ("connector_id" %in% colnames(bc.synapses))
    bc.synapses$connector_id <- as.character(bc.synapses$connector_id)
  if ("pre_id" %in% colnames(bc.synapses))
    bc.synapses$pre_id <- as.character(bc.synapses$pre_id)
  if ("post_id" %in% colnames(bc.synapses))
    bc.synapses$post_id <- as.character(bc.synapses$post_id)

  # Clean and filter
  message("Building split edgelist...")
  bc.synapses <- bc.synapses %>%
    dplyr::filter(size >= banc.size.threshold,
                  pre_id != post_id) %>%
    dplyr::select(dplyr::any_of(c("connector_id", "x", "y", "z", "size",
                                   "prepost", "pre_id", "post_id",
                                   "inside", "status", "label"))) %>%
    dplyr::mutate(label = hemibrainr:::standard_compartments(label))

  # Get per-neuron NT predictions for join
  ntpred_file <- file.path(banc.save.path, "banc_ntpred.feather")
  bc.nt <- data.frame(root_id = character(0),
                      conf_nt = character(0),
                      conf_nt_p = numeric(0))
  if (file.exists(ntpred_file)) {
    bc.nt <- arrow::read_feather(ntpred_file) %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::select(root_id,
                    conf_nt = neurotransmitter_predicted,
                    conf_nt_p = neurotransmitter_score)
    message(sprintf("  Loaded NT predictions for %d neurons", nrow(bc.nt)))
  } else {
    message(sprintf("  NT prediction file not found at %s — proceeding without NT join", ntpred_file))
  }

  # Separate pre- and post-synapses
  bc.presynapses <- bc.synapses %>% dplyr::filter(prepost == 0) %>% dplyr::ungroup()
  bc.postsynapses <- bc.synapses %>% dplyr::filter(prepost == 1) %>% dplyr::ungroup()
  rm(bc.synapses); gc()

  # Build pre_label lookup from presynapses (vectorized join, not per-row match)
  pre_label_lookup <- bc.presynapses %>%
    dplyr::distinct(connector_id, .keep_all = TRUE) %>%
    dplyr::select(connector_id, pre_label = label)

  # Build edgelist from postsynapses (both endpoints must be in banc meta)
  bc.ids <- root.ids

  # Dataset-boundary correction: sensory (flow == "afferent") neurons only
  # have axons inside the dataset; effector (flow == "efferent") neurons only
  # have dendrites inside. Flow-centrality splits on truncated morphology are
  # unreliable for these, so override both endpoints' compartment labels
  # before aggregation so the edgelist reflects biological reality.
  afferent_ids <- banc.meta[[version_id_col]][
    !is.na(banc.meta$flow) & banc.meta$flow == "afferent"]
  efferent_ids <- banc.meta[[version_id_col]][
    !is.na(banc.meta$flow) & banc.meta$flow == "efferent"]
  message(sprintf("  Label override: %s afferent -> axon, %s efferent -> dendrite",
                  format(length(afferent_ids), big.mark = ","),
                  format(length(efferent_ids), big.mark = ",")))

  bc.elist.split <- bc.postsynapses %>%
    dplyr::filter(post_id %in% bc.ids, pre_id %in% bc.ids) %>%
    dplyr::rename(post = post_id, pre = pre_id, post_label = label) %>%
    dplyr::left_join(pre_label_lookup, by = "connector_id") %>%
    dplyr::mutate(pre_label = ifelse(is.na(pre_label), "unknown", pre_label)) %>%
    dplyr::mutate(
      pre_label = dplyr::case_when(
        pre %in% afferent_ids ~ "axon",
        pre %in% efferent_ids ~ "dendrite",
        TRUE ~ pre_label),
      post_label = dplyr::case_when(
        post %in% afferent_ids ~ "axon",
        post %in% efferent_ids ~ "dendrite",
        TRUE ~ post_label)
    ) %>%
    # Count synapses per compartment-level connection
    dplyr::group_by(post, pre, post_label, pre_label) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(count > 0) %>%
    # Norm = fraction of total input to post neuron (across all compartments)
    dplyr::group_by(post) %>%
    dplyr::mutate(post_count = sum(count),
                  norm = signif(count / post_count, 4)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pre) %>%
    dplyr::mutate(pre_count = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(connection = paste(pre_label, post_label, sep = "-")) %>%
    # Join NT predictions
    dplyr::left_join(bc.nt, by = c("post" = "root_id")) %>%
    dplyr::mutate(conf_nt = ifelse(is.na(conf_nt), "unknown", conf_nt)) %>%
    dplyr::rename(post_conf_nt = conf_nt, post_conf_nt_p = conf_nt_p) %>%
    dplyr::left_join(bc.nt, by = c("pre" = "root_id")) %>%
    dplyr::mutate(conf_nt = ifelse(is.na(conf_nt), "unknown", conf_nt)) %>%
    dplyr::rename(pre_conf_nt = conf_nt, pre_conf_nt_p = conf_nt_p) %>%
    as.data.frame(stringsAsFactors = FALSE)

  colnames(bc.elist.split) <- snakecase::to_snake_case(colnames(bc.elist.split))

  arrow::write_feather(bc.elist.split, split_elist_file)
  message(sprintf("Split edgelist saved: %s (%s connections, %s pre, %s post)",
                  split_elist_file,
                  format(nrow(bc.elist.split), big.mark = ","),
                  format(length(unique(bc.elist.split$pre)), big.mark = ","),
                  format(length(unique(bc.elist.split$post)), big.mark = ",")))
} else if (.syn_source == "v3") {
  ############################################################################
  # v3 source: labels are pre-computed and persisted in the enriched parquet
  # (see Section 3 v3 spatial-NN block). We read them back and aggregate to
  # edge level — no re-NN here.
  ############################################################################
  split_elist_file_v3 <- file.path(save.path,
                                    sprintf("banc_%s_edgelist_split_v3.feather",
                                            banc.version))

  message("Reading v3 enriched (with persisted pre_label/post_label)...")
  v3_path <- file.path(save.path,
                       sprintf("banc_%s_synapses_v3_enriched.parquet", banc.version))
  v3 <- arrow::read_parquet(
    v3_path,
    col_select = c("id", "size", "pre_root_id", "post_root_id",
                   "pre_label", "post_label"))
  v3$pre_root_id  <- as.character(v3$pre_root_id)
  v3$post_root_id <- as.character(v3$post_root_id)
  v3 <- v3 %>%
    dplyr::filter(size >= banc.size.threshold,
                  pre_root_id != post_root_id,
                  pre_root_id %in% root.ids, post_root_id %in% root.ids)
  message(sprintf("  v3 synapses (filtered): %s rows",
                  format(nrow(v3), big.mark = ",")))

  # Belt-and-braces: ensure any NA labels fall through to "unknown"
  v3$pre_label  <- ifelse(is.na(v3$pre_label),  "unknown", v3$pre_label)
  v3$post_label <- ifelse(is.na(v3$post_label), "unknown", v3$post_label)

  # NT predictions (same join as v2 path)
  ntpred_file <- file.path(banc.save.path, "banc_ntpred.feather")
  bc.nt <- data.frame(root_id = character(0),
                      conf_nt = character(0),
                      conf_nt_p = numeric(0))
  if (file.exists(ntpred_file)) {
    bc.nt <- arrow::read_feather(ntpred_file) %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::select(root_id,
                    conf_nt = neurotransmitter_predicted,
                    conf_nt_p = neurotransmitter_score)
  }

  v3_elist <- v3 %>%
    dplyr::rename(pre = pre_root_id, post = post_root_id) %>%
    dplyr::group_by(pre, post, pre_label, post_label) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(count > 0) %>%
    dplyr::group_by(post) %>%
    dplyr::mutate(post_count = sum(count),
                  norm = signif(count / post_count, 4)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pre) %>%
    dplyr::mutate(pre_count = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(connection = paste(pre_label, post_label, sep = "-")) %>%
    dplyr::left_join(bc.nt, by = c("post" = "root_id")) %>%
    dplyr::mutate(conf_nt = ifelse(is.na(conf_nt), "unknown", conf_nt)) %>%
    dplyr::rename(post_conf_nt = conf_nt, post_conf_nt_p = conf_nt_p) %>%
    dplyr::left_join(bc.nt, by = c("pre" = "root_id")) %>%
    dplyr::mutate(conf_nt = ifelse(is.na(conf_nt), "unknown", conf_nt)) %>%
    dplyr::rename(pre_conf_nt = conf_nt, pre_conf_nt_p = conf_nt_p) %>%
    as.data.frame(stringsAsFactors = FALSE)

  colnames(v3_elist) <- snakecase::to_snake_case(colnames(v3_elist))
  arrow::write_feather(v3_elist, split_elist_file_v3)
  message(sprintf("Split edgelist (v3, spatial-NN) saved: %s (%s connections, %s pre, %s post)",
                  split_elist_file_v3,
                  format(nrow(v3_elist), big.mark = ","),
                  format(length(unique(v3_elist$pre)), big.mark = ","),
                  format(length(unique(v3_elist$post)), big.mark = ",")))
}

##############################
### SECTION 6: SKELETONS   ###
##############################

# message("\n=== Section 6: L2 Skeletons ===")
# 
# skel_folder <- file.path(save.path, "banc_banc_space_l2_swc")
# dir.create(skel_folder, showWarnings = FALSE)
# skel_files <- list.files(banc.l2swc.save.path, full.names = TRUE, pattern = "\\.swc$")
# if (length(skel_files)) {
#   file.copy(from = skel_files, to = skel_folder, overwrite = FALSE)
#   # Remove skeletons not in this version's neuron set
#   existing <- list.files(skel_folder, full.names = TRUE, pattern = "\\.swc$")
#   toremove <- existing[!gsub("\\.swc$", "", basename(existing)) %in% root.ids]
#   if (length(toremove)) file.remove(toremove)
#   message(sprintf("  Copied %d L2 skeletons (removed %d not in v%s)",
#                   length(list.files(skel_folder, pattern = "\\.swc$")),
#                   length(toremove), banc.version))
# } else {
#   message("  No SWC files found in ", banc.l2swc.save.path)
# }

##############################
### SECTION 7: CUT-OUTS    ###
##############################

# message("\n=== Section 7: Cut-outs ===")
# 
# cut.outs <- c("mushroom_body", "antennal_lobe", "central_complex",
#               "optic", "suboesophageal_zone", "front_leg", "abdominal_neuromere")
# 
# for (cut.out in cut.outs) {
#   cut.out.good <- snakecase::to_snake_case(cut.out)
#   save.path.cut.out <- file.path(save.path, cut.out.good)
#   dir.create(save.path.cut.out, showWarnings = FALSE)
# 
#   if (cut.out == "optic") {
#     cut.out <- "^LO|^LOP|^AME|^ME"
#   }
#   if (cut.out == "antennal_lobe") {
#     cut.out <- "antennal_lobe|olfactory_receptor|thermosensory_receptor|hygrosensory_receptor|CSD"
#   }
#   if (cut.out == "suboesophageal_zone") {
#     cut.out <- "^FLA|^SEZ|^GNG|^SAD|^AMMC|^PRW"
#   }
#   if (cut.out == "front_leg") {
#     cut.out <- "^LegNp\\(T1\\)|T1|^ProNM-T1|^LNp_T1"
#   }
#   if (cut.out == "abdominal_neuromere") {
#     cut.out <- "^ANm|^ABDNM"
#   }
# 
#   # Select neurons for this cut-out
#   if (cut.out == "mushroom_body") {
#     cut.out <- "mushroom_body|kenyon_cell|APL|DPM|LHMB1|OA-VPM3"
#     kc.ids <- banc.meta %>%
#       dplyr::filter(cell_class == "kenyon_cell", side == "right") %>%
#       dplyr::pull(.data[[version_id_col]])
#     kc.elist.simp <- banc.elist.simp %>%
#       dplyr::filter(pre %in% kc.ids | post %in% kc.ids) %>%
#       dplyr::mutate(pre = ifelse(pre %in% kc.ids, "KC", pre),
#                     post = ifelse(post %in% kc.ids, "KC", post)) %>%
#       dplyr::group_by(pre, post) %>%
#       dplyr::mutate(count = sum(count)) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count >= 100)
#     chosen.ids <- setdiff(unique(c(kc.elist.simp$pre, kc.elist.simp$post)), kc.ids)
#     chosen.ids <- setdiff(chosen.ids, "KC")
#     banc.meta.cutout <- banc.meta %>%
#       dplyr::filter(side == "right" & (grepl(cut.out, super_class) |
#                                          grepl(cut.out, cell_class) |
#                                          grepl(cut.out, cell_sub_class) |
#                                          grepl(cut.out, cell_type)) |
#                       .data[[version_id_col]] %in% chosen.ids)
#   } else if (cut.out.good %in% c("optic", "front_leg")) {
#     if (nrow(banc.syns)) {
#       banc.chosen.pre <- banc.syns %>%
#         dplyr::filter(grepl(cut.out, neuropil), side == "right") %>%
#         dplyr::group_by(pre_root_id) %>%
#         dplyr::mutate(count = dplyr::n()) %>%
#         dplyr::ungroup() %>%
#         dplyr::filter(count >= 100) %>%
#         dplyr::pull(pre_root_id)
#       banc.chosen.post <- banc.syns %>%
#         dplyr::filter(grepl(cut.out, neuropil), side == "right") %>%
#         dplyr::group_by(post_root_id) %>%
#         dplyr::mutate(count = dplyr::n()) %>%
#         dplyr::ungroup() %>%
#         dplyr::filter(count >= 100) %>%
#         dplyr::pull(post_root_id)
#       chosen.ids <- unique(c(banc.chosen.pre, banc.chosen.post))
#     } else {
#       chosen.ids <- character(0)
#     }
#     banc.meta.cutout <- banc.meta %>%
#       dplyr::filter(.data[[version_id_col]] %in% chosen.ids)
#   } else if (cut.out.good %in% c("suboesophageal_zone", "abdominal_neuromere")) {
#     if (nrow(banc.syns)) {
#       banc.chosen.pre <- banc.syns %>%
#         dplyr::filter(grepl(cut.out, neuropil)) %>%
#         dplyr::group_by(pre_root_id) %>%
#         dplyr::mutate(count = dplyr::n()) %>%
#         dplyr::ungroup() %>%
#         dplyr::filter(count >= 100) %>%
#         dplyr::pull(pre_root_id)
#       banc.chosen.post <- banc.syns %>%
#         dplyr::filter(grepl(cut.out, neuropil)) %>%
#         dplyr::group_by(post_root_id) %>%
#         dplyr::mutate(count = dplyr::n()) %>%
#         dplyr::ungroup() %>%
#         dplyr::filter(count >= 100) %>%
#         dplyr::pull(post_root_id)
#       chosen.ids <- unique(c(banc.chosen.pre, banc.chosen.post))
#     } else {
#       chosen.ids <- character(0)
#     }
#     banc.meta.cutout <- banc.meta %>%
#       dplyr::filter(.data[[version_id_col]] %in% chosen.ids)
#   } else {
#     banc.meta.cutout <- banc.meta %>%
#       dplyr::filter(grepl(cut.out, super_class) |
#                       grepl(cut.out, cell_class) |
#                       grepl(cut.out, cell_sub_class) |
#                       grepl(cut.out, cell_type))
#   }
# 
#   root.ids.cut.out <- na.omit(unique(banc.meta.cutout[[version_id_col]]))
# 
#   # Subset edgelist and synapses
#   banc.elist.simp.cut.out <- banc.elist.simp %>%
#     dplyr::filter(pre %in% root.ids.cut.out & post %in% root.ids.cut.out)
#   banc.syns.cut.out <- banc.syns %>%
#     dplyr::filter(pre_root_id %in% root.ids.cut.out | post_root_id %in% root.ids.cut.out)
# 
#   # Save
#   arrow::write_feather(banc.meta.cutout,
#     file.path(save.path.cut.out, sprintf("banc_%s_%s_meta.feather", banc.version, cut.out.good)))
#   arrow::write_feather(banc.elist.simp.cut.out,
#     file.path(save.path.cut.out, sprintf("banc_%s_%s_simple_edgelist.feather", banc.version, cut.out.good)))
#   arrow::write_feather(banc.syns.cut.out,
#     file.path(save.path.cut.out, sprintf("banc_%s_%s_synapses.feather", banc.version, cut.out.good)))
#   message(sprintf("  %s: %d neurons, %d connections, %d synapses",
#                   cut.out.good, nrow(banc.meta.cutout),
#                   nrow(banc.elist.simp.cut.out), nrow(banc.syns.cut.out)))
# }

##############################
### SECTION 8: OBJ MESHES  ###
##############################

# message("\n=== Section 8: OBJ meshes ===")
# 
# save.path.obj <- file.path(save.path, "obj")
# dir.create(save.path.obj, showWarnings = FALSE)
# save.path.obj.np <- file.path(save.path.obj, "neuropils")
# dir.create(save.path.obj.np, showWarnings = FALSE)
# 
# Rvcg::vcgObjWrite(as.mesh3d(bancr::banc.surf),
#                    filename = file.path(save.path.obj, "banc_volume_nm.obj"))
# Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_vnc_neuropil.surf),
#                    filename = file.path(save.path.obj, "banc_vnc_neuropil_nm.obj"))
# Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_brain_neuropil.surf),
#                    filename = file.path(save.path.obj, "banc_brain_neuropil_nm.obj"))
# regs <- bancr::banc_vnc_neuropils.surf$RegionList
# for (reg in regs) {
#   Rvcg::vcgObjWrite(as.mesh3d(subset(banc_vnc_neuropils.surf, reg)),
#                      filename = file.path(save.path.obj.np, sprintf("banc_neuropil_%s_nm.obj", reg)))
# }
# regs <- bancr::banc_brain_neuropils.surf$RegionList
# for (reg in regs) {
#   Rvcg::vcgObjWrite(as.mesh3d(subset(banc_brain_neuropils.surf, reg)),
#                      filename = file.path(save.path.obj.np, sprintf("banc_neuropil_%s_nm.obj", reg)))
# }
# message("  OBJ meshes exported")

##############################
### SECTION 9: GCS SYNC    ###
##############################

message("\n=== Section 9: GCS bucket sync ===")
system(sprintf("gsutil -m rsync -r %s gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s",
               save.path, banc.version))
message("  Sync complete")

##############################
### SECTION 10: PNG ZIPs   ###
##############################
# Bundle the manually-reviewed cross-match PNGs (per dataset) into flat ZIPs
# and push each to gs://lee-lab.../neuron_annotations/, replacing whatever
# was there before. Only the verified buckets (1_perfect..5_possible) are
# included — rejected buckets (6_no_match, *_wrong, etc.) are skipped.
# Mirror is intra-BANC and not a "cross-match", so excluded.

message("\n=== Section 10: Verified PNG match ZIPs -> GCS ===")

png_zip_targets <- list(
  fafb       = banc.fafb.correct.match.path,
  manc       = banc.manc.correct.match.path,
  hemibrain  = banc.hemibrain.correct.match.path,
  fanc       = banc.fanc.correct.match.path,
  malecns    = banc.malecns.correct.match.path
)
png_zip_levels <- c("1_perfect", "2_confident", "3_good", "4_likely", "5_possible")
png_zip_gcs_base <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_annotations"

for (dataset in names(png_zip_targets)) {
  base_dir <- png_zip_targets[[dataset]]
  if (!dir.exists(base_dir)) {
    message(sprintf("  %s: correct/ folder not found at %s — skipping", dataset, base_dir))
    next
  }
  png_files <- unlist(lapply(png_zip_levels, function(lvl) {
    folder <- file.path(base_dir, lvl)
    if (!dir.exists(folder)) return(character(0))
    list.files(folder, pattern = "\\.png$", full.names = TRUE)
  }))
  png_files <- unique(png_files)
  if (!length(png_files)) {
    message(sprintf("  %s: no PNGs in verified buckets — skipping ZIP", dataset))
    next
  }

  zip_name <- sprintf("banc_to_%s_verified_png_matches.zip", dataset)
  zip_local <- file.path(tempdir(), zip_name)
  if (file.exists(zip_local)) unlink(zip_local)

  # Use `zip -j -@` (junk paths -> flat structure; read filenames from stdin)
  # to avoid argv length limits with tens of thousands of files.
  files_list <- tempfile(fileext = ".txt")
  writeLines(png_files, files_list)
  on.exit(unlink(files_list), add = TRUE)
  st <- system(sprintf("zip -q -j -@ %s < %s",
                        shQuote(zip_local), shQuote(files_list)))
  unlink(files_list)
  if (st != 0L || !file.exists(zip_local)) {
    message(sprintf("  %s: zip failed (status %d) — skipping push", dataset, st))
    next
  }
  zip_size_mb <- round(file.info(zip_local)$size / 1024^2, 1)

  gcs_path <- sprintf("%s/%s", png_zip_gcs_base, zip_name)
  message(sprintf("  %s: %d PNGs -> %s (%.1f MB) -> %s",
                  dataset, length(png_files), basename(zip_local),
                  zip_size_mb, gcs_path))
  st <- system2("gsutil", c("cp", zip_local, gcs_path),
                 stdout = FALSE, stderr = FALSE)
  if (st != 0L) {
    message(sprintf("    gsutil cp failed (status %d)", st))
  }
  unlink(zip_local)
}

##############################
### SUMMARY                ###
##############################

message(sprintf("\n### banc: data compilation complete for v%s [total: %s] ###",
                banc.version,
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

# Report output files
output_files <- list.files(save.path, recursive = TRUE, full.names = TRUE)
message(sprintf("Output directory: %s (%d files)", save.path, length(output_files)))
top_files <- list.files(save.path, full.names = TRUE, recursive = FALSE)
top_files <- top_files[!file.info(top_files)$isdir]
for (f in top_files) {
  size_mb <- round(file.info(f)$size / 1024^2, 1)
  message(sprintf("  %s (%.1f MB)", basename(f), size_mb))
}

})

# Hard exit: bypasses R's normal teardown which otherwise triggers
# `free(): invalid pointer` from reticulate/arrow destructor ordering.
# Job 38651583 wrote every output successfully then SIGABRT'd on cleanup,
# which `set -e` propagated to SLURM as FAILED and cancelled the chain.
quit(save = "no", status = 0, runLast = FALSE)
