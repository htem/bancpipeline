#' prep (optic-lobe preset) — Build BANC + target meta / edgelist / NBLAST / seeds for optic-lobe alignment.
#'
#' Loaded by `alignment/banc-alignment-prep.R --region optic-lobe`. Prepares
#' BANC and target optic-lobe neuron sets, edgelists, NBLAST scores, and
#' seed labels for the iterative connectivity alignment algorithm. Excludes
#' R1-6 and Lai (not reconstructed in BANC). Includes DNc01/DNc02 because
#' they dendrite in the optic lobe but aren't tagged optic.
#'
#' @section Reads:
#'   - `alignment-data-sources.R`-resolved BANC edgelist + target edgelist + NBLAST feather
#'   - SeaTable `banc_meta` (live curation)
#'   - target metadata via `load_target_meta_for_alignment()`
#'
#' @section Writes (under `data/optic_lobe/`):
#'   - `{banc,target}_optic_{side}_meta.csv`
#'   - `{banc,target}_optic_{side}_edgelist.feather`
#'   - `banc_target_optic_{side}_nblast.csv`
#'   - `banc_optic_{side}_seeds.csv`
#'
#' @section CLI (via dispatcher):
#'   --region optic-lobe [side]    side = right | left | both (default right)
#'   --source {local|gcs}          required
#'   --local-root PATH             optional
#'   --syn-source {v2|v3}          default banc.synapse.source.default
#'
#' @section Invoked by:
#'   `o2/alignment/o2_banc_alignment.sh`,
#'   `o2/alignment/o2_banc_optic_prep.sh`,
#'   `o2/oneshots/o2_banc_optic_sweep.sh`
#'
#' @section Paper:
#'   Methods §"Optic lobe cell-type alignment".
#'
#' @section Notes:
#'   - Output target columns are named `target_*` (target-dataset-agnostic).
#'     The companion `update-seatable.R` maps them onto the existing
#'     `fafb_alignment_*` SeaTable schema columns — when a non-FAFB target
#'     is wired in, those SeaTable columns also need parallel naming.

###########################################################
### Optic Lobe Cross-Dataset Matching: Data Preparation
###
### Prepares BANC and FAFB optic lobe neuron sets,
### edgelists, NBLAST scores, and seed labels for the
### iterative connectivity alignment algorithm.
###
### Outputs to data/optic_lobe/:
###   - {banc,fafb}_optic_{side}_meta.csv
###   - {banc,fafb}_optic_{side}_edgelist.feather
###   - banc_fafb_optic_{side}_nblast.csv
###   - banc_optic_{side}_seeds.csv
###########################################################
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

local({

######################
### CONFIGURATION  ###
######################

# CLI: Rscript banc-alignment-prep.R --region optic-lobe [side] --source {local|gcs} [--local-root PATH]
#                                     [--syn-source {v2|v3}]
parsed <- parse_alignment_data_args(commandArgs(trailingOnly = TRUE))
side_to_process <- if (length(parsed$positional) > 0) parsed$positional[1] else "right"
data_source     <- parsed$source
data_local_root <- parsed$local_root
syn_source      <- parsed$syn_source
min_edge_count  <- 1        # no threshold — use all connections
seed_holdout    <- 0.2      # fraction of seeds held out for validation
nblast_version  <- "783"    # FAFB NBLAST version

# Output directory
output_dir <- "data/optic_lobe"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Neuron pool filters
optic_regions      <- c("optic_lobe", "optic_lobes")
optic_super_classes <- c("visual_centrifugal", "visual_projection", "optic_lobe_intrinsic")
exclude_super_classes <- c("glia", "not_a_neuron", "trachea")
# Exclude R1-6 (not reconstructed in BANC) and Lai
lamina_types <- c("R1-6", "Lai")
# Extra cell types to include (dendrite in optic lobes but not classified as optic)
extra_optic_types <- c("DNc01", "DNc02")

######################
### DATA LOADING   ###
######################

# Resolve data paths from --source {local|gcs} [+ --local-root]
data_paths <- resolve_alignment_paths(
  source        = data_source,
  local_root    = data_local_root,
  banc_version  = banc.version,
  nblast_version = nblast_version,
  syn_source    = syn_source
)

# Load edgelists
message("=== Loading edgelists ===")
banc.el <- arrow::read_feather(data_paths$banc_edgelist)
target.el <- arrow::read_feather(data_paths$fafb_edgelist)
banc.el <- banc.el %>% dplyr::mutate(pre = as.character(pre), post = as.character(post))
target.el <- target.el %>% dplyr::mutate(pre = as.character(pre), post = as.character(post))
message(sprintf("  BANC edgelist: %s edges", format(nrow(banc.el), big.mark = ",")))
message(sprintf("  FAFB edgelist: %s edges", format(nrow(target.el), big.mark = ",")))

##############################
### STEP 0.1: NEURON SETS  ###
##############################

message("\n=== Step 0.1: Building optic lobe neuron sets ===")

# BANC metadata — always pull fresh from SeaTable (not cached)
banc.meta <- banctable_query()
banc.meta <- banc.meta %>%
  dplyr::mutate(root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id))

# Fill missing side from position (same method as banc-ids.R)
banc_missing_side <- is.na(banc.meta$side) | banc.meta$side == ""
banc_has_pos <- !is.na(banc.meta$root_position_nm) & banc.meta$root_position_nm != ""
fill_side <- banc_missing_side & banc_has_pos
if (any(fill_side)) {
  message(sprintf("  Computing side for %d neurons from root_position_nm...",
                  sum(fill_side)))
  side_xyz <- nat::xyzmatrix(banc.meta$root_position_nm[fill_side])
  lrdiffs <- bancr:::banc_lr_position(side_xyz, units = "nm")
  banc.meta$side[fill_side] <- ifelse(lrdiffs > 0, "right", "left")
  message(sprintf("  Filled: %d right, %d left",
                  sum(banc.meta$side[fill_side] == "right"),
                  sum(banc.meta$side[fill_side] == "left")))
}

# Target metadata. Default target = "fafb" gives the paper-run behaviour.
# To align against a different reference connectome, change `target_name` and
# extend `load_target_meta()` in alignment-data-sources.R.
target_name <- "fafb"
target.meta <- load_target_meta(target_name)

# BANC filter: region OR super_class OR root_region neuropil
# root_region is only used for neurons WITHOUT a super_class — avoids pulling in
# descending/central_brain neurons whose soma happens to be near an optic neuropil
# Side filter: "both" includes all sides
side_filter_banc <- if (side_to_process == "both") {
  rep(TRUE, nrow(banc.meta))
} else {
  tolower(banc.meta$side) == tolower(side_to_process)
}

banc.optic <- banc.meta %>%
  dplyr::filter(
    side_filter_banc,
    (tolower(region) %in% optic_regions |
       super_class %in% optic_super_classes |
       (grepl("_(ME|LO|LOP|AME)_", root_region) &
          (is.na(super_class) | super_class == "")) |
       cell_type %in% extra_optic_types),
    !super_class %in% exclude_super_classes | cell_type %in% extra_optic_types,
    !cell_type %in% lamina_types,
    proofread == TRUE | roughly_proofread == TRUE
  ) %>%
  dplyr::filter(!is.na(root_id)) %>%
  # GT cell_type: prefer fafb_cell_type, fall back to cell_type
  # Drop "auto:" prefixed, "unknown", "fragment", "glia" labels
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(fafb_cell_type) & fafb_cell_type != "" &
        !grepl("^auto:", fafb_cell_type) &
        !grepl("unknown|fragment|glia", tolower(fafb_cell_type)) ~ fafb_cell_type,
      !is.na(cell_type) & cell_type != "" &
        !grepl("^auto:", cell_type) &
        !grepl("unknown|fragment|glia", tolower(cell_type)) ~ cell_type,
      TRUE ~ NA_character_
    )
  )

# Target filter: region OR super_class only (no root_region column on target)
side_filter_target <- if (side_to_process == "both") {
  rep(TRUE, nrow(target.meta))
} else {
  tolower(target.meta$target_side) == tolower(side_to_process)
}

target.pool <- target.meta %>%
  dplyr::filter(
    side_filter_target,
    (tolower(target_region) %in% optic_regions |
       target_super_class %in% optic_super_classes |
       target_cell_type %in% extra_optic_types),
    !target_super_class %in% exclude_super_classes | target_cell_type %in% extra_optic_types,
    !target_cell_type %in% lamina_types
  ) %>%
  dplyr::filter(!is.na(target_id)) %>%
  # Override: R7/R8 neurons get histamine as NT (FAFB-specific photoreceptor handling)
  dplyr::mutate(target_top_nt = ifelse(target_cell_type %in% c("R7", "R8"), "histamine", target_top_nt))

message(sprintf("  query (BANC) optic %s: %s neurons (%s with cell_type)",
                side_to_process,
                format(nrow(banc.optic), big.mark = ","),
                sum(!is.na(banc.optic$cell_type) & banc.optic$cell_type != "")))
message(sprintf("  target (%s) optic %s: %s neurons (%s with target_cell_type)",
                toupper(target_name),
                side_to_process,
                format(nrow(target.pool), big.mark = ","),
                sum(!is.na(target.pool$target_cell_type) & target.pool$target_cell_type != "")))

# Target types summary
target_types <- target.pool %>%
  dplyr::filter(!is.na(target_cell_type), target_cell_type != "") %>%
  dplyr::count(target_cell_type, sort = TRUE)
message(sprintf("  target (%s) unique types: %d", toupper(target_name), nrow(target_types)))

# Save metadata. Filenames keep their region+side+dataset infix; for the paper
# run target_name == "fafb" so output naming matches the published intermediate.
readr::write_csv(banc.optic,
  alignment_path("prep-banc-meta", query = "banc", target = target_name,
                 region = "optic-lobe", side = side_to_process,
                 vq = banc.version, syn = syn_source,
                 ext = "csv", dir = output_dir))
readr::write_csv(target.pool,
  alignment_path("prep-target-meta", query = "banc", target = target_name,
                 region = "optic-lobe", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = output_dir))

################################
### STEP 0.2: EDGELISTS      ###
################################

message("\n=== Step 0.2: Extracting optic lobe edgelists ===")

# Filter edgelists: at least one endpoint in the optic pool, count >= threshold
banc_optic_ids <- unique(banc.optic$root_id)
banc.el.optic <- banc.el %>%
  dplyr::filter(count >= min_edge_count,
                (pre %in% banc_optic_ids | post %in% banc_optic_ids)) %>%
  dplyr::select(pre, post, count, norm)

# Exclude lamina types from target edgelist (not reconstructed in query)
target_lamina_ids <- target.meta %>%
  dplyr::filter(target_cell_type %in% lamina_types) %>%
  dplyr::pull(target_id) %>%
  unique()
message(sprintf("  Excluding %d target lamina neurons from edgelist", length(target_lamina_ids)))

target_optic_ids <- unique(target.pool$target_id)
target.el.optic <- target.el %>%
  dplyr::filter(count >= min_edge_count,
                !pre %in% target_lamina_ids, !post %in% target_lamina_ids,
                (pre %in% target_optic_ids | post %in% target_optic_ids)) %>%
  dplyr::select(pre, post, count, norm)

message(sprintf("  BANC optic edgelist: %s edges (>=%d synapses, at least one endpoint in pool)",
                format(nrow(banc.el.optic), big.mark = ","), min_edge_count))
message(sprintf("  target (%s) optic edgelist: %s edges (>=%d synapses, at least one endpoint in pool)",
                toupper(target_name),
                format(nrow(target.el.optic), big.mark = ","), min_edge_count))

arrow::write_feather(banc.el.optic,
  alignment_path("prep-banc-edges", query = "banc", target = target_name,
                 region = "optic-lobe", side = side_to_process,
                 vq = banc.version, syn = syn_source,
                 ext = "feather", dir = output_dir))
arrow::write_feather(target.el.optic,
  alignment_path("prep-target-edges", query = "banc", target = target_name,
                 region = "optic-lobe", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "feather", dir = output_dir))

################################
### STEP 0.3: NBLAST SCORES  ###
################################

message("\n=== Step 0.3: Extracting NBLAST scores ===")

nblast_file <- data_paths$nblast

if (!is.na(nblast_file) && file.exists(nblast_file)) {
  nblast.all <- arrow::read_feather(nblast_file)
  message(sprintf("  Loaded %s NBLAST rows from %s",
                  format(nrow(nblast.all), big.mark = ","), basename(nblast_file)))

  # Filter to BANC neurons in optic pool (keep all FAFB matches)
  nblast.optic <- nblast.all %>%
    dplyr::filter(as.character(pt_root_id) %in% banc_optic_ids) %>%
    dplyr::select(banc_id = pt_root_id, target_id = match_id, nblast_score = score) %>%
    dplyr::mutate(banc_id = as.character(banc_id),
                  target_id = as.character(target_id))

  message(sprintf("  NBLAST scores for optic pool: %s rows covering %s/%s BANC neurons",
                  format(nrow(nblast.optic), big.mark = ","),
                  format(length(unique(nblast.optic$banc_id)), big.mark = ","),
                  format(length(banc_optic_ids), big.mark = ",")))

  readr::write_csv(nblast.optic,
    alignment_path("prep-nblast", query = "banc", target = target_name,
                   region = "optic-lobe", side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = output_dir))
  rm(nblast.all); gc()
} else {
  message("  No NBLAST file resolved — proceeding without NBLAST (Tier 2 empty)")
  nblast.optic <- data.frame(banc_id = character(), target_id = character(),
                              nblast_score = numeric())
  readr::write_csv(nblast.optic,
    alignment_path("prep-nblast", query = "banc", target = target_name,
                   region = "optic-lobe", side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = output_dir))
}

################################
### STEP 0.4: SEED LABELS    ###
################################

message("\n=== Step 0.4: Identifying seed labels ===")

# Tier 1: BANC neurons with cell_type in SeaTable
# optic_lobe_intrinsic neurons with cell_type are ALL held out as validation
# (the algorithm should discover their types from connectivity alone)
# Non-intrinsic typed neurons (visual_projection, visual_centrifugal) can serve as anchors
typed_banc <- banc.optic %>%
  dplyr::filter(!is.na(cell_type), cell_type != "") %>%
  dplyr::select(root_888 = root_id, cell_type, super_class)

intrinsic_typed <- typed_banc %>%
  dplyr::filter(super_class == "optic_lobe_intrinsic") %>%
  dplyr::mutate(tier = 1L, is_holdout = TRUE) %>%
  dplyr::select(-super_class)

non_intrinsic_typed <- typed_banc %>%
  dplyr::filter(super_class != "optic_lobe_intrinsic") %>%
  dplyr::mutate(tier = 1L, is_holdout = FALSE) %>%
  dplyr::select(-super_class)

tier1 <- dplyr::bind_rows(non_intrinsic_typed, intrinsic_typed)

message(sprintf("  Tier 1 (SeaTable cell_type): %d neurons (%d anchor [non-intrinsic], %d holdout [intrinsic])",
                nrow(tier1), sum(!tier1$is_holdout), sum(tier1$is_holdout)))

# Tier 2: BANC neurons with NBLAST matches but no cell_type
# Assign type of best-scoring FAFB match
tier1_ids <- tier1$root_888
if (nrow(nblast.optic) > 0) {
  # Get best NBLAST match per BANC neuron
  best_nblast <- nblast.optic %>%
    dplyr::filter(!banc_id %in% tier1_ids) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::slice_max(nblast_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # Look up target cell_type for the best match
  best_nblast <- best_nblast %>%
    dplyr::left_join(target.pool %>%
      dplyr::select(target_id, target_cell_type) %>%
      dplyr::distinct(target_id, .keep_all = TRUE),
      by = "target_id") %>%
    dplyr::filter(!is.na(target_cell_type), target_cell_type != "")

  tier2 <- best_nblast %>%
    dplyr::select(root_888 = banc_id, cell_type = target_cell_type) %>%
    dplyr::mutate(tier = 2L, is_holdout = FALSE)
} else {
  tier2 <- data.frame(root_888 = character(), cell_type = character(),
                       tier = integer(), is_holdout = logical())
}
message(sprintf("  Tier 2 (NBLAST-initialized): %d neurons", nrow(tier2)))

# Tier 3: everything else (no initial type)
tier3_ids <- setdiff(banc_optic_ids, c(tier1$root_888, tier2$root_888))
tier3 <- data.frame(root_888 = tier3_ids, cell_type = NA_character_,
                     tier = 3L, is_holdout = FALSE)
message(sprintf("  Tier 3 (unassigned): %d neurons", nrow(tier3)))

seeds <- dplyr::bind_rows(tier1, tier2, tier3)
readr::write_csv(seeds,
  alignment_path("prep-seeds", query = "banc", target = target_name,
                 region = "optic-lobe", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = output_dir))

################################
### STEP 0.5: TYPE CAPACITY  ###
################################

message("\n=== Step 0.5: Computing per-type capacity from FAFB L-R variability ===")

# Get FAFB optic neurons for BOTH sides (not just the side being processed)
# to compute natural left-right count variability per type
target.pool.both <- target.meta %>%
  dplyr::filter(
    (tolower(region) %in% optic_regions |
       super_class %in% optic_super_classes),
    !super_class %in% exclude_super_classes,
    !cell_type %in% lamina_types,
    !grepl("photoreceptor|lamina_intrinsic", cell_class, ignore.case = TRUE),
    !is.na(cell_type), cell_type != ""
  )

type_counts_by_side <- target.pool.both %>%
  dplyr::count(cell_type, side, name = "count") %>%
  tidyr::pivot_wider(names_from = side, values_from = count,
                     values_fill = 0, names_prefix = "target_")

# Capacity depends on whether we're running one side or both
if (side_to_process == "both") {
  # Bilateral: capacity = left + right + slack (both hemispheres compete)
  type_capacity <- type_counts_by_side %>%
    dplyr::mutate(
      target_total = target_left + target_right,
      target_max = pmax(target_left, target_right, na.rm = TRUE),
      lr_diff = abs(target_left - target_right),
      give = pmax(lr_diff, ceiling(0.1 * target_total), 2),
      capacity = target_total + give
    ) %>%
    dplyr::select(cell_type, target_left, target_right, target_max, lr_diff, give, capacity)
  message("  Bilateral capacity: left + right + slack")
} else {
  # Single side: capacity = max(left, right) + slack
  type_capacity <- type_counts_by_side %>%
    dplyr::mutate(
      target_max = pmax(target_left, target_right, na.rm = TRUE),
      lr_diff = abs(target_left - target_right),
      give = pmax(lr_diff, ceiling(0.1 * target_max), 1),
      capacity = target_max + give
    ) %>%
    dplyr::select(cell_type, target_left, target_right, target_max, lr_diff, give, capacity)
}

capacity_path <- alignment_path("prep-capacity", query = "banc", target = target_name,
                                 region = "optic-lobe", side = side_to_process,
                                 vq = banc.version, vt = nblast_version,
                                 ext = "csv", dir = output_dir)
readr::write_csv(type_capacity, capacity_path)

message(sprintf("  Types with capacity: %d", nrow(type_capacity)))
if (side_to_process == "both") {
  target_total <- type_capacity$target_left + type_capacity$target_right
  message(sprintf("  Mean capacity/(L+R) ratio: %.2f",
                  mean(type_capacity$capacity / pmax(target_total, 1))))
} else {
  message(sprintf("  Mean capacity/max ratio: %.2f",
                  mean(type_capacity$capacity / pmax(type_capacity$target_max, 1))))
}
message(sprintf("  Saved: %s", capacity_path))

################################
### SUMMARY                  ###
################################

message("\n=== Preparation complete ===")
message(sprintf("Side: %s", side_to_process))
message(sprintf("BANC neurons: %s", format(nrow(banc.optic), big.mark = ",")))
message(sprintf("FAFB neurons: %s", format(nrow(target.pool), big.mark = ",")))
message(sprintf("target unique types: %d", nrow(target_types)))
message(sprintf("BANC edgelist: %s edges", format(nrow(banc.el.optic), big.mark = ",")))
message(sprintf("FAFB edgelist: %s edges", format(nrow(target.el.optic), big.mark = ",")))
message(sprintf("Seeds: %d Tier1 + %d Tier2 + %d Tier3 = %d total",
                nrow(tier1), nrow(tier2), nrow(tier3), nrow(seeds)))
message(sprintf("Output: %s/", output_dir))

})
