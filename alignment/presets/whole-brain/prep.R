#' prep (whole-brain preset) — Build BANC + target meta / edgelist / NBLAST / seeds for whole-brain alignment.
#'
#' Loaded by `alignment/banc-alignment-prep.R --region whole-brain`. Prepares
#' BANC and target (default FAFB) brain neuron sets (optic_lobe +
#' central_brain + neck_connective), edgelists, NBLAST scores, and seed
#' labels for the iterative connectivity alignment algorithm. Includes
#' ascending / descending classes (unlike the optic-lobe preset).
#'
#' @section Reads:
#'   - `alignment-data-sources.R`-resolved BANC edgelist + target edgelist + NBLAST feather
#'   - SeaTable `banc_meta` (live curation)
#'   - target metadata via `load_target_meta_for_alignment()`
#'   - env var `BANC_WB_OUTPUT_DIR` — override output directory
#'
#' @section Writes (under `data/whole_brain_alignment/` by default):
#'   - `{banc,target}_brain_{side}_meta.csv`
#'   - `{banc,target}_brain_{side}_edgelist.feather`
#'   - `banc_target_brain_{side}_nblast.csv`
#'   - `banc_brain_{side}_seeds.csv`
#'   - `target_type_capacity.csv`
#'
#' @section CLI (via dispatcher):
#'   --region whole-brain [side]   side = right | left | both (default both)
#'   --source {local|gcs}          required
#'   --local-root PATH             optional
#'   --syn-source {v2|v3}          default banc.synapse.source.default
#'
#' @section Invoked by:
#'   `o2/alignment/o2_banc_alignment.sh`,
#'   `o2/alignment/o2_banc_wb_prep_v888v2.sh`,
#'   `o2/alignment/o2_banc_wb_prep_v888v2_priority.sh`,
#'   `o2/alignment/o2_banc_wb_cosine_v888v2.sh`,
#'   `o2/alignment/o2_banc_wb_cosine_v888v3.sh`,
#'   `o2/alignment/o2_banc_wb_production_v888v2.sh`,
#'   `o2/oneshots/o2_banc_wb_probe.sh`
#'
#' @section Notes:
#'   - Output target columns are named `target_*` (target-dataset-agnostic).
#'     The companion `update-seatable.R` maps them onto the existing
#'     `fafb_alignment_*` SeaTable schema columns — when a non-FAFB target
#'     is wired in, those SeaTable columns also need parallel naming.

###########################################################
### Whole-Brain Cross-Dataset Matching: Data Preparation
###
### Prepares BANC and FAFB brain neuron sets (optic_lobe +
### central_brain), edgelists, NBLAST scores, and seed labels
### for the iterative connectivity alignment algorithm.
###
### Excludes ascending/descending neurons (neck_connective
### region or ascending/descending super_class).
###
### Usage: Rscript alignment/banc-alignment-prep.R --region whole-brain [right|left]
###
### Outputs to data/whole_brain_alignment/:
###   - {banc,fafb}_brain_{side}_meta.csv
###   - {banc,fafb}_brain_{side}_edgelist.feather
###   - banc_fafb_brain_{side}_nblast.csv
###   - banc_brain_{side}_seeds.csv
###   - fafb_type_capacity.csv
###########################################################
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

local({

######################
### CONFIGURATION  ###
######################

# CLI: Rscript banc-alignment-prep.R --region whole-brain [side] --source {local|gcs} [--local-root PATH]
#                                     [--syn-source {v2|v3}]
parsed <- parse_alignment_data_args(commandArgs(trailingOnly = TRUE))
side_to_process <- if (length(parsed$positional) > 0) parsed$positional[1] else "both"
data_source     <- parsed$source
data_local_root <- parsed$local_root
syn_source      <- parsed$syn_source
min_edge_count  <- 1
seed_holdout    <- 0.2
nblast_version  <- "783"

output_dir <- Sys.getenv("BANC_WB_OUTPUT_DIR",
                         unset = "data/whole_brain_alignment")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Neuron pool filters — whole brain (optic + central + ascending/descending)
include_regions      <- c("optic_lobe", "optic_lobes", "central_brain", "neck_connective")
exclude_super_classes <- c("glia", "not_a_neuron", "trachea")
exclude_regions      <- character(0)  # no region exclusions
# Types not reconstructed in BANC
exclude_types        <- c("R1-6", "Lai")
# Super_classes whose typed neurons serve as anchors (not holdout).
# Well-characterised projection/motor/sensory types. BANC uses both bare
# ("ascending", "descending", "sensory") and compound forms
# ("sensory_ascending", "sensory_descending", "ascending_visceral_circulatory")
# in the super_class column — list both so %in% catches all variants.
# A bug omitting the compound forms put ~530 well-typed neurons (mostly
# sensory_ascending) in holdout instead of anchor; restored 2026-05-02.
anchor_super_classes <- c("visual_centrifugal", "visual_projection",
                          "sensory", "sensory_ascending", "sensory_descending",
                          "motor", "endocrine", "efferent",
                          "ascending", "ascending_visceral_circulatory",
                          "descending")

######################
### DATA LOADING   ###
######################

data_paths <- resolve_alignment_paths(
  source         = data_source,
  local_root     = data_local_root,
  banc_version   = banc.version,
  nblast_version = nblast_version,
  syn_source     = syn_source
)

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

message("\n=== Step 0.1: Building whole-brain neuron sets ===")

# BANC metadata — prefer the curated feather written by banc/share/banc-data.R
# (it has match-quality fixes baked in via banc_filter_neurons + curated columns)
# and fall back to live SeaTable if the feather is missing. Override with
# BANC_META_SOURCE=seatable to force the live query.
.meta_src <- tolower(Sys.getenv("BANC_META_SOURCE", unset = "feather"))
.meta_feather <- file.path(banc.versioned.save.path,
                            sprintf("banc_%s_meta.feather", banc.version))
if (.meta_src == "feather" && file.exists(.meta_feather)) {
  message(sprintf("  Loading BANC meta from feather: %s (mtime %s)",
                  basename(.meta_feather),
                  format(file.info(.meta_feather)$mtime, "%Y-%m-%d %H:%M")))
  banc.meta <- arrow::read_feather(.meta_feather)
  # banc-data.R writes the current root_id as banc_<version>_id (not root_id).
  # Rename so downstream code can refer to root_id uniformly.
  .id_col <- sprintf("banc_%s_id", banc.version)
  if (.id_col %in% colnames(banc.meta) && !"root_id" %in% colnames(banc.meta))
    banc.meta <- banc.meta %>% dplyr::rename(root_id = !!.id_col)
} else {
  if (.meta_src == "feather")
    message("  BANC meta feather not found; falling back to banctable_query()")
  banc.meta <- banctable_query()
}
banc.meta <- banc.meta %>%
  dplyr::mutate(root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id))

# Fill missing side from position
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
target.meta <- load_target_meta(target_name) %>%
  dplyr::mutate(
    # Treat "unknown" as missing
    target_cell_type = ifelse(tolower(target_cell_type) == "unknown",
                              NA_character_, target_cell_type),
    # Neurons without cell_type use their super_class as label
    target_cell_type = ifelse(is.na(target_cell_type) | target_cell_type == "",
                              ifelse(!is.na(target_super_class) & target_super_class != "",
                                     target_super_class, NA_character_),
                              target_cell_type)
  )

# Target filter first: we need its type vocabulary to validate BANC GT labels.
# Mirrors BANC's region OR asc/desc rescue so AN_* types living in FAFB
# region == "ventral_nerve_cord" stay in vocab.
side_filter_target <- if (side_to_process == "both") {
  rep(TRUE, nrow(target.meta))
} else {
  tolower(target.meta$target_side) == tolower(side_to_process)
}

target.pool <- target.meta %>%
  dplyr::filter(
    side_filter_target,
    tolower(target_region) %in% include_regions |
      grepl("ascending|descending", target_super_class),
    !target_super_class %in% exclude_super_classes,
    !target_cell_type %in% exclude_types
  ) %>%
  dplyr::filter(!is.na(target_id)) %>%
  # R7/R8 NT override (FAFB-specific photoreceptor handling)
  dplyr::mutate(target_top_nt = ifelse(target_cell_type %in% c("R7", "R8"),
                                       "histamine", target_top_nt))

target_type_set <- unique(target.pool$target_cell_type[
  !is.na(target.pool$target_cell_type) & target.pool$target_cell_type != ""])

# BANC filter: region-based with exclusions
# When side_to_process="both", include all sides (whole brain bilateral)
side_filter_banc <- if (side_to_process == "both") {
  rep(TRUE, nrow(banc.meta))
} else {
  tolower(banc.meta$side) == tolower(side_to_process)
}

banc.pool <- banc.meta %>%
  dplyr::filter(
    side_filter_banc,
    # ascending|descending super_class is the backwards-compatible substitute
    # for region == "neck_connective" (the BANC region column is being
    # restructured; cervical_connective will live as a tract entry).
    tolower(region) %in% include_regions |
      grepl("ascending|descending", super_class),
    !super_class %in% exclude_super_classes,
    !cell_type %in% exclude_types,
    proofread == TRUE | roughly_proofread == TRUE
  ) %>%
  dplyr::filter(!is.na(root_id)) %>%
  # GT cell_type: prefer SeaTable cell_type when it exists in the target
  # vocab, else fall back to fafb_cell_type (BANC's curated cross-match
  # column) when *that* is in vocab; else NA. Drops "auto:" prefix and
  # unknown/fragment/glia placeholders. `fafb_cell_type` here is a BANC
  # column (the curated FAFB match), not a target-side column.
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) & cell_type != "" &
        !grepl("^auto:", cell_type) &
        !grepl("unknown|fragment|glia", tolower(cell_type)) &
        cell_type %in% target_type_set ~ cell_type,
      !is.na(fafb_cell_type) & fafb_cell_type != "" &
        !grepl("^auto:", fafb_cell_type) &
        !grepl("unknown|fragment|glia", tolower(fafb_cell_type)) &
        fafb_cell_type %in% target_type_set ~ fafb_cell_type,
      TRUE ~ NA_character_
    )
  )

message(sprintf("  BANC brain %s: %s neurons (%s with cell_type)",
                side_to_process,
                format(nrow(banc.pool), big.mark = ","),
                sum(!is.na(banc.pool$cell_type) & banc.pool$cell_type != "")))
message(sprintf("  FAFB brain %s: %s neurons (%s with cell_type)",
                side_to_process,
                format(nrow(target.pool), big.mark = ","),
                sum(!is.na(target.pool$cell_type) & target.pool$cell_type != "")))

# Region breakdown
banc_region_counts <- banc.pool %>%
  dplyr::count(region, name = "n") %>%
  dplyr::arrange(dplyr::desc(n))
message("  BANC by region:")
for (i in seq_len(nrow(banc_region_counts))) {
  message(sprintf("    %s: %s", banc_region_counts$region[i],
                  format(banc_region_counts$n[i], big.mark = ",")))
}

target_types <- target.pool %>%
  dplyr::filter(!is.na(cell_type), cell_type != "") %>%
  dplyr::count(cell_type, sort = TRUE)
message(sprintf("  FAFB unique types: %d", nrow(target_types)))

# Save metadata
readr::write_csv(banc.pool,
  alignment_path("prep-banc-meta", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
                 vq = banc.version, syn = syn_source,
                 ext = "csv", dir = output_dir))
readr::write_csv(target.pool,
  alignment_path("prep-target-meta", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = output_dir))

################################
### STEP 0.2: EDGELISTS      ###
################################

message("\n=== Step 0.2: Extracting brain edgelists ===")

banc_pool_ids <- unique(banc.pool$root_id)
banc.el.pool <- banc.el %>%
  dplyr::filter(count >= min_edge_count,
                (pre %in% banc_pool_ids | post %in% banc_pool_ids)) %>%
  dplyr::select(pre, post, count, norm)

# Exclude lamina neurons from FAFB edgelist
target_exclude_ids <- target.meta %>%
  dplyr::filter(cell_type %in% exclude_types) %>%
  dplyr::pull(target_id) %>%
  unique()
message(sprintf("  Excluding %d FAFB lamina neurons from edgelist", length(target_exclude_ids)))

target_pool_ids <- unique(target.pool$target_id)
target.el.pool <- target.el %>%
  dplyr::filter(count >= min_edge_count,
                !pre %in% target_exclude_ids, !post %in% target_exclude_ids,
                (pre %in% target_pool_ids | post %in% target_pool_ids)) %>%
  dplyr::select(pre, post, count, norm)

message(sprintf("  BANC brain edgelist: %s edges", format(nrow(banc.el.pool), big.mark = ",")))
message(sprintf("  FAFB brain edgelist: %s edges", format(nrow(target.el.pool), big.mark = ",")))

arrow::write_feather(banc.el.pool,
  alignment_path("prep-banc-edges", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
                 vq = banc.version, syn = syn_source,
                 ext = "feather", dir = output_dir))
arrow::write_feather(target.el.pool,
  alignment_path("prep-target-edges", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
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

  nblast.pool <- nblast.all %>%
    dplyr::filter(as.character(pt_root_id) %in% banc_pool_ids) %>%
    dplyr::select(banc_id = pt_root_id, target_id = match_id, nblast_score = score) %>%
    dplyr::mutate(banc_id = as.character(banc_id),
                  target_id = as.character(target_id))

  message(sprintf("  NBLAST scores for brain pool: %s rows covering %s/%s BANC neurons",
                  format(nrow(nblast.pool), big.mark = ","),
                  format(length(unique(nblast.pool$banc_id)), big.mark = ","),
                  format(length(banc_pool_ids), big.mark = ",")))

  readr::write_csv(nblast.pool,
    alignment_path("prep-nblast", query = "banc", target = target_name,
                   region = "whole-brain", side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = output_dir))
  rm(nblast.all); gc()
} else {
  message("  No NBLAST file resolved — proceeding without NBLAST (Tier 2 empty)")
  nblast.pool <- data.frame(banc_id = character(), target_id = character(),
                              nblast_score = numeric())
  readr::write_csv(nblast.pool,
    alignment_path("prep-nblast", query = "banc", target = target_name,
                   region = "whole-brain", side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = output_dir))
}

################################
### STEP 0.4: SEED LABELS    ###
################################

message("\n=== Step 0.4: Identifying seed labels ===")

typed_banc <- banc.pool %>%
  dplyr::filter(!is.na(cell_type), cell_type != "") %>%
  dplyr::select(root_888 = root_id, cell_type, super_class, region)

# Anchors: typed neurons with anchor super_classes (well-characterised projection types)
# Holdout: everything else (intrinsic/central types the algorithm must discover)
anchor_typed <- typed_banc %>%
  dplyr::filter(super_class %in% anchor_super_classes) %>%
  dplyr::mutate(tier = 1L, is_holdout = FALSE) %>%
  dplyr::select(-super_class, -region)

holdout_typed <- typed_banc %>%
  dplyr::filter(!super_class %in% anchor_super_classes) %>%
  dplyr::mutate(tier = 1L, is_holdout = TRUE) %>%
  dplyr::select(-super_class, -region)

tier1 <- dplyr::bind_rows(anchor_typed, holdout_typed)

message(sprintf("  Tier 1 (SeaTable cell_type): %d neurons (%d anchor, %d holdout)",
                nrow(tier1), sum(!tier1$is_holdout), sum(tier1$is_holdout)))

# Tier 2: NBLAST-initialized (no cell_type but has NBLAST match)
tier1_ids <- tier1$root_888
if (nrow(nblast.pool) > 0) {
  best_nblast <- nblast.pool %>%
    dplyr::filter(!banc_id %in% tier1_ids) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::slice_max(nblast_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

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

tier3_ids <- setdiff(banc_pool_ids, c(tier1$root_888, tier2$root_888))
tier3 <- data.frame(root_888 = tier3_ids, cell_type = NA_character_,
                     tier = 3L, is_holdout = FALSE)
message(sprintf("  Tier 3 (unassigned): %d neurons", nrow(tier3)))

seeds <- dplyr::bind_rows(tier1, tier2, tier3)
readr::write_csv(seeds,
  alignment_path("prep-seeds", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = output_dir))

################################
### STEP 0.5: TYPE CAPACITY  ###
################################

message("\n=== Step 0.5: Computing per-type capacity from FAFB L-R variability ===")

target.pool.both <- target.meta %>%
  dplyr::filter(
    tolower(region) %in% include_regions,
    !super_class %in% exclude_super_classes,
    !cell_type %in% exclude_types,
    !is.na(cell_type), cell_type != ""
  )

type_counts_by_side <- target.pool.both %>%
  dplyr::count(cell_type, side, name = "count") %>%
  tidyr::pivot_wider(names_from = side, values_from = count,
                     values_fill = 0, names_prefix = "target_")

if (side_to_process == "both") {
  # Bilateral: pooled capacity = left + right + slack
  # Per-side hard caps (cap_left/cap_right) split give in half so the total
  # cap_left + cap_right >= capacity (off by at most 1 due to ceiling).
  type_capacity <- type_counts_by_side %>%
    dplyr::mutate(
      target_total = target_left + target_right,
      target_max = pmax(target_left, target_right, na.rm = TRUE),
      lr_diff = abs(target_left - target_right),
      give = pmax(lr_diff, ceiling(0.1 * target_total), 2),
      capacity = target_total + give,
      give_half = ceiling(give / 2),
      cap_left = target_left + give_half,
      cap_right = target_right + give_half
    ) %>%
    dplyr::select(cell_type, target_left, target_right, target_max,
                  lr_diff, give, capacity, cap_left, cap_right)
  message("  Bilateral capacity: left + right + slack (per-side caps included)")
} else {
  # Per-side run: capacity is single-side; cap_left/cap_right unused but emitted
  # for schema consistency (set equal to capacity).
  type_capacity <- type_counts_by_side %>%
    dplyr::mutate(
      target_max = pmax(target_left, target_right, na.rm = TRUE),
      lr_diff = abs(target_left - target_right),
      give = pmax(lr_diff, ceiling(0.1 * target_max), 1),
      capacity = target_max + give,
      cap_left = capacity,
      cap_right = capacity
    ) %>%
    dplyr::select(cell_type, target_left, target_right, target_max,
                  lr_diff, give, capacity, cap_left, cap_right)
}

readr::write_csv(type_capacity,
  alignment_path("prep-capacity", query = "banc", target = target_name,
                 region = "whole-brain", side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = output_dir))

message(sprintf("  Types with capacity: %d", nrow(type_capacity)))
message(sprintf("  Mean capacity/max ratio: %.2f",
                mean(type_capacity$capacity / pmax(type_capacity$target_max, 1))))

################################
### SUMMARY                  ###
################################

message("\n=== Preparation complete ===")
message(sprintf("Side: %s", side_to_process))
message(sprintf("BANC neurons: %s", format(nrow(banc.pool), big.mark = ",")))
message(sprintf("FAFB neurons: %s", format(nrow(target.pool), big.mark = ",")))
message(sprintf("FAFB types: %d", nrow(target_types)))
message(sprintf("BANC edgelist: %s edges", format(nrow(banc.el.pool), big.mark = ",")))
message(sprintf("FAFB edgelist: %s edges", format(nrow(target.el.pool), big.mark = ",")))
message(sprintf("Seeds: %d Tier1 + %d Tier2 + %d Tier3 = %d total",
                nrow(tier1), nrow(tier2), nrow(tier3), nrow(seeds)))
message(sprintf("Output: %s/", output_dir))

})
