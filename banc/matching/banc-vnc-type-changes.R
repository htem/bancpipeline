#' banc-vnc-type-changes — Detect VNC cell_type changes vs the new VNC typology.
#'
#' Two parts: (A) VNC intrinsic neurons — joins
#' `banc_vnc_cell_type_and_dimorphism.csv` to SeaTable on `root_626` to
#' detect cell_type / hemilineage changes and dimorphism assignments;
#' (B) Effector neurons — joins `banc_to_manc_effectors.csv` and resolves
#' new types via MANC matches in `franken_meta()`.
#'
#' @section Reads:
#'   - `data/codex/banc_vnc_cell_type_and_dimorphism.csv`
#'   - `data/codex/banc_to_manc_effectors.csv`
#'   - SeaTable `banc_meta`
#'   - `franken_meta()`
#'
#' @section Writes:
#'   - `data/codex/vnc_type_changes.csv` — proposed cell_type / dimorphism updates
#'
#' @section Notes:
#'   - banc IDs read as character — 18-digit root IDs lose precision as numeric.

###########################################################
### Compare new vs old VNC types
###
### Part A: VNC intrinsic neurons
### Reads banc_vnc_dimorphism CSV, joins to seatable
### by banc ID -> root_626 to get old cell_type.
### Computes type-changed status, hemilineage changes.
###
### Part B: Effector neurons (motor + visceral_circulatory)
### Reads banc_to_manc_effectors CSV, joins manc match ID
### to franken_meta to find "new" type, assigns dimorphism
### from median dimorph score per type.
###
### Plots are in panels-vnc-morphology-connectivity-types.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: VNC type change analysis ###")

###########################
### Read data           ###
###########################

# Read banc ID as character to avoid 64-bit integer precision loss.
# R's double type can only represent ~15.9 significant digits; BANC root IDs
# have 18 digits, so reading as numeric silently corrupts the last 2-3 digits.
# cell_type in the CSV is the "new" type from the updated VNC typology.
vnc <- readr::read_csv("data/codex/banc_vnc_cell_type_and_dimorphism.csv",
                        col_types = readr::cols(ID = "c"),
                        show_col_types = FALSE) %>%
  dplyr::rename(banc_id = ID, new_type = Type, dimorphism = Dimorphism) %>%
  dplyr::mutate(dimorphism = gsub("^known_", "", dimorphism),
                dimorphism = dplyr::case_when(
                  dimorphism == "Homologous" ~ "isomorphic",
                  dimorphism == "Sex Specific" ~ "female-specific",
                  dimorphism == "Dimorphic" ~ "dimorphic",
                  dimorphism == "female_specific" ~ "female-specific",
                  dimorphism == "male_specific" ~ "male-specific",
                  TRUE ~ dimorphism
                ))

bc <- banctable_query("SELECT _id, root_id, root_626, supervoxel_id, super_class, flow, cell_type, side, region, manc_cell_type, malecns_cell_type, manc_png_match, manc_nblast_match, manc_match, malecns_match, status, sexually_dimorphic, body_part_effector, hemilineage, banc_match, banc_match_supervoxel_id FROM banc_meta") %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::mutate(root_626 = as.character(root_626),
                root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id),
                manc_png_match = as.character(manc_png_match)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)

# Recover old cell_type via manc_png_match -> franken_meta manc_id -> cell_type.
# cell_type in banc_meta will be overwritten with new types, so we use the
# MANC PNG match body ID to look up the original type from franken_meta.
fm_types <- franken_meta("SELECT manc_id, cell_type FROM franken_meta") %>%
  dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
  dplyr::mutate(manc_id = as.character(manc_id)) %>%
  dplyr::distinct(manc_id, .keep_all = TRUE)

bc <- bc %>%
  dplyr::left_join(fm_types %>% dplyr::rename(seatable_type = cell_type),
                   by = c("manc_png_match" = "manc_id"))

# Join to get old type and identifiers
df <- vnc %>%
  dplyr::left_join(bc, by = c("banc_id" = "root_626")) %>%
  dplyr::filter(super_class == "ventral_nerve_cord_intrinsic")

# Clean dimorphism: any remaining empty -> "isomorphic"
df <- df %>%
  dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism) | dimorphism == "",
                                             "isomorphic", dimorphism))

# Propagate dimorphism within cell types: if any member of a type has a
# non-isomorphic label, all members get the strongest label.
# Priority: female-specific > male-specific > dimorphic > isomorphic
.dim_rank <- c("female-specific" = 4, "male-specific" = 3,
               "dimorphic" = 2, "isomorphic" = 1)
type_level_dim <- df %>%
  dplyr::filter(!is.na(new_type), new_type != "") %>%
  dplyr::mutate(.rank = .dim_rank[dimorphism]) %>%
  dplyr::group_by(new_type) %>%
  dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(new_type, dimorphism_type = dimorphism)
n_propagated <- sum(df$new_type %in% type_level_dim$new_type &
  df$dimorphism != type_level_dim$dimorphism_type[match(df$new_type, type_level_dim$new_type)],
  na.rm = TRUE)
df <- df %>%
  dplyr::left_join(type_level_dim, by = "new_type") %>%
  dplyr::mutate(dimorphism = dplyr::if_else(!is.na(dimorphism_type),
                                             dimorphism_type, dimorphism)) %>%
  dplyr::select(-dimorphism_type)
if (n_propagated > 0)
  message(sprintf("  Propagated dimorphism labels to %d neurons via cell type", n_propagated))

# Determine type changed: compare CSV cell_type (new_type) against seatable cell_type.
# A neuron's type "changed" if the new typology assigns a different name from seatable.
df <- df %>%
  dplyr::mutate(type_changed = dplyr::if_else(
    !is.na(seatable_type) & !is.na(new_type) & seatable_type != new_type,
    "type changed", "type did not change"
  ))

# Effective type: the CSV cell_type is the resolved type from the updated typology.
df <- df %>%
  dplyr::mutate(effective_type = new_type)

message(sprintf("  VNC intrinsic: %d neurons total, %d matched to seatable",
                nrow(df), sum(!is.na(df$seatable_type))))
message(sprintf("  Type changes: %d changed, %d unchanged",
                sum(df$type_changed == "type changed"),
                sum(df$type_changed == "type did not change")))

#########################################
### Effector neurons (motor/endocrine) ##
#########################################

# Read effector CSV: banc_to_manc_effectors.csv
# Columns: banc ID (root_626), type (BANC seatable cell_type),
#   dimorph score, twin ID, twin score, manc match ID, manc match score
eff <- readr::read_csv("data/codex/id_id_matches/banc_to_manc_effectors.csv",
                        col_types = readr::cols(`banc ID` = "c",
                                                `twin ID` = "c",
                                                `manc match ID` = "c"),
                        show_col_types = FALSE) %>%
  dplyr::rename(banc_id = `banc ID`,
                csv_type = type,
                dimorph_score = `dimorph score`,
                twin_id = `twin ID`,
                twin_score = `twin score`,
                manc_match_id = `manc match ID`,
                manc_match_score = `manc match score`)

# Join manc match ID to franken_meta to find the "new" type
eff <- eff %>%
  dplyr::left_join(fm_types, by = c("manc_match_id" = "manc_id")) %>%
  dplyr::rename(new_type = cell_type)

# Assign dimorphism from median dimorph score per type:
# >4 = sex-specific, >2 = dimorphic, <=2 = isomorphic
eff_dimorph <- eff %>%
  dplyr::filter(!is.na(csv_type), csv_type != "") %>%
  dplyr::group_by(csv_type) %>%
  dplyr::summarise(median_dimorph = stats::median(dimorph_score, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(dimorphism = dplyr::case_when(
    median_dimorph > 4 ~ "sex-specific",
    median_dimorph > 2 ~ "dimorphic",
    TRUE ~ "isomorphic"
  ))

eff <- eff %>%
  dplyr::left_join(eff_dimorph %>% dplyr::select(csv_type, dimorphism),
                   by = "csv_type") %>%
  dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism), "isomorphic", dimorphism))

# Join to seatable for root_id, supervoxel_id, and seatable_type (from manc_png_match)
eff <- eff %>%
  dplyr::left_join(bc %>% dplyr::select(root_626, root_id, supervoxel_id,
                                         super_class, manc_nblast_match,
                                         seatable_type),
                   by = c("banc_id" = "root_626"))

# Determine type changed: seatable_type (from manc_png_match) vs new_type (from manc_match_id).
# If seatable_type is NA (no manc_png_match), the neuron is likely a type change.
eff <- eff %>%
  dplyr::mutate(
    effective_type = dplyr::coalesce(new_type, csv_type),
    type_changed = dplyr::case_when(
      is.na(seatable_type) & !is.na(new_type) ~ "type changed",
      !is.na(seatable_type) & !is.na(new_type) & seatable_type != new_type ~ "type changed",
      TRUE ~ "type did not change"
    )
  )

# Propagate dimorphism within effective_type (same logic as VNC intrinsic)
eff_type_level_dim <- eff %>%
  dplyr::filter(!is.na(effective_type), effective_type != "") %>%
  dplyr::mutate(.rank = .dim_rank[dimorphism]) %>%
  dplyr::group_by(effective_type) %>%
  dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(effective_type, dimorphism_type = dimorphism)
n_eff_propagated <- sum(eff$effective_type %in% eff_type_level_dim$effective_type &
  eff$dimorphism != eff_type_level_dim$dimorphism_type[match(eff$effective_type, eff_type_level_dim$effective_type)],
  na.rm = TRUE)
eff <- eff %>%
  dplyr::left_join(eff_type_level_dim, by = "effective_type") %>%
  dplyr::mutate(dimorphism = dplyr::if_else(!is.na(dimorphism_type),
                                             dimorphism_type, dimorphism)) %>%
  dplyr::select(-dimorphism_type)

message(sprintf("  Effectors: %d neurons, %d with MANC match type, %d type changes",
                nrow(eff),
                sum(!is.na(eff$new_type)),
                sum(eff$type_changed == "type changed")))
message(sprintf("  Effector dimorphism: %s",
                paste(sprintf("%s=%d", names(table(eff$dimorphism)),
                              table(eff$dimorphism)), collapse = ", ")))
if (n_eff_propagated > 0)
  message(sprintf("  Propagated dimorphism labels to %d effector neurons via effective_type", n_eff_propagated))

###########################
### Hemilineage lookup  ###
###########################

# VNC intrinsic type names encode hemilineage directly:
#   IN19B057 -> hemilineage 19B
#   IN06A085 -> hemilineage 06A
#   IN20A.22A085 -> hemilineage 20A.22A
# For non-VNC types (brain types like SMP496), use franken_meta lookup.
extract_hemilineage <- function(type_name) {
  # Use gsub to strip the leading "IN" and trailing digits, keeping only
  # the hemilineage code. Returns NA for non-matching types.
  result <- rep(NA_character_, length(type_name))
  valid <- !is.na(type_name)
  is_vnc <- valid & grepl("^IN\\d+[A-Z]", type_name, perl = TRUE)
  # Strip trailing digit-only suffix to get "IN" + hemilineage
  # e.g. IN19B057 -> IN19B, IN20A.22A085 -> IN20A.22A
  hemi <- sub("\\d+$", "", type_name[is_vnc])
  result[is_vnc] <- sub("^IN", "", hemi)
  result
}

fm <- franken_meta()

# Franken_meta lookup for brain types (ascending, descending, etc.)
hemi_lookup <- fm %>%
  dplyr::filter(!is.na(cell_type), cell_type != "",
                !is.na(hemilineage), hemilineage != "") %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::select(cell_type, hemilineage)

# For each type: first try VNC name extraction, then fall back to franken_meta.
# Old type is the seatable cell_type; new type is the CSV cell_type.
df <- df %>%
  dplyr::mutate(
    old_hemilineage = extract_hemilineage(seatable_type),
    new_hemilineage = extract_hemilineage(new_type)
  )

# Fall back to franken_meta for types where VNC extraction failed (brain types)
old_fm <- hemi_lookup %>% dplyr::rename(old_hemilineage_fm = hemilineage)
new_fm <- hemi_lookup %>% dplyr::rename(new_hemilineage_fm = hemilineage)

df <- df %>%
  dplyr::left_join(old_fm, by = c("seatable_type" = "cell_type")) %>%
  dplyr::left_join(new_fm, by = c("new_type" = "cell_type")) %>%
  dplyr::mutate(
    old_hemilineage = dplyr::coalesce(old_hemilineage, old_hemilineage_fm),
    new_hemilineage = dplyr::coalesce(new_hemilineage, new_hemilineage_fm)
  ) %>%
  dplyr::select(-old_hemilineage_fm, -new_hemilineage_fm)

df <- df %>%
  dplyr::mutate(hemi_changed = dplyr::case_when(
    is.na(old_hemilineage) | is.na(new_hemilineage) ~ "no hemilineage",
    old_hemilineage == new_hemilineage ~ "hemilineage did not change",
    TRUE ~ "hemilineage changed"
  ))

message(sprintf("  Hemilineage changes: %d changed, %d unchanged, %d unknown",
                sum(df$hemi_changed == "hemilineage changed"),
                sum(df$hemi_changed == "hemilineage did not change"),
                sum(df$hemi_changed == "no hemilineage")))

################################################
### NBLAST scores for type-changed neurons   ###
################################################

# BANC-MANC NBLAST scores from GCS feather file.
# The feather file has columns: root_626, match_id, match_cell_type, score.
# For each type-changed neuron, look up NBLAST scores for three categories:
# nblast (top hit), seatable (current assignment), new (effective CSV type).
message("  Loading BANC-MANC NBLAST scores...")

gcs_nblast_path <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_manc_v1.2.1_nblast.feather"
nblast_cache_dir <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache_dir, showWarnings = FALSE, recursive = TRUE)
nblast_local <- file.path(nblast_cache_dir, basename(gcs_nblast_path))

if (!file.exists(nblast_local)) {
  message("  Downloading NBLAST feather from GCS...")
  system2("gsutil", c("cp", gcs_nblast_path, nblast_local),
          stdout = TRUE, stderr = TRUE)
}

type_changed_df <- df %>%
  dplyr::filter(type_changed == "type changed")

nblast_available <- file.exists(nblast_local)
if (nblast_available) {
  nb_all <- arrow::read_feather(nblast_local)
  message(sprintf("  NBLAST feather: %d rows", nrow(nb_all)))

  # MANC IDs that share the same super_class as the VNC intrinsic query neurons.
  # Used to restrict VNC intrinsic NBLAST top hits to same super_class.
  vnc_match_ids <- fm %>%
    dplyr::filter(super_class == "ventral_nerve_cord_intrinsic",
                  !is.na(manc_id)) %>%
    dplyr::pull(manc_id) %>% as.character() %>% unique()

  # Build banc_id -> root_id map: nb_all$root_626 matches bc$root_626 (= banc_id),
  # while nb_all$pt_root_id matches bc$root_id (current root_id).
  # Include both VNC intrinsic and effector neurons.
  id_map <- dplyr::bind_rows(
    df %>% dplyr::select(banc_id, root_id),
    eff %>% dplyr::select(banc_id, root_id)
  ) %>%
    dplyr::filter(!is.na(root_id), root_id != "", root_id != banc_id) %>%
    dplyr::distinct(banc_id, root_id)

  # Map pt_root_id (= current root_id) back to banc_id (= root_626)
  rootid_to_bancid <- stats::setNames(id_map$banc_id, id_map$root_id)

  # For each type-changed neuron, find NBLAST scores for three categories:
  #   "nblast"   — top MANC hit (any type), score and cell type
  #   "seatable" — best score for the seatable cell_type assignment
  #   "new"      — best score for the effective type from the CSV
  type_changed_ids <- unique(type_changed_df$banc_id)
  type_changed_rootids <- id_map$root_id[id_map$banc_id %in% type_changed_ids]

  nb_subset <- nb_all %>%
    dplyr::filter(root_626 %in% type_changed_ids |
                    pt_root_id %in% type_changed_rootids) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
    dplyr::mutate(root_626 = dplyr::case_when(
      root_626 %in% type_changed_ids ~ root_626,
      pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
      TRUE ~ root_626
    )) %>%
    dplyr::select(-pt_root_id)

  # Top NBLAST hit (any type) — restricted to same super_class as query
  top_scores <- nb_subset %>%
    dplyr::mutate(match_id = as.character(match_id)) %>%
    dplyr::filter(match_id %in% vnc_match_ids) %>%
    dplyr::group_by(root_626) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(banc_id = root_626,
                     nblast_score = score,
                     nblast_type = match_cell_type)

  # Best score for seatable type
  st_scores <- type_changed_df %>%
    dplyr::transmute(banc_id, lookup_type = seatable_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_subset, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")

  # Best score for effective type (new, filled with old if blank)
  new_scores <- type_changed_df %>%
    dplyr::transmute(banc_id, lookup_type = effective_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_subset, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

  type_changed_df <- type_changed_df %>%
    dplyr::left_join(top_scores, by = "banc_id") %>%
    dplyr::left_join(st_scores, by = "banc_id") %>%
    dplyr::left_join(new_scores, by = "banc_id")

  n_with_scores <- sum(!is.na(type_changed_df$nblast_score))
  message(sprintf("  Found NBLAST scores for %d/%d type-changed neurons",
                  n_with_scores, nrow(type_changed_df)))
} else {
  type_changed_df$nblast_score <- NA_real_
  type_changed_df$nblast_type <- NA_character_
  type_changed_df$seatable_nblast <- NA_real_
  type_changed_df$new_nblast <- NA_real_
  message("  NBLAST feather not available; skipping scores")
}

########################################
### MANC body ID lookup for bancsee  ###
########################################

# Look up MANC body IDs by cell_type for neuroglancer links.
# The NBLAST feather itself maps match_cell_type -> match_id (MANC body ID),
# so we use it directly as the authoritative source.
message("  Building MANC body ID lookup...")

if (nblast_available) {
  # Extract unique (match_id, match_cell_type) pairs from NBLAST results
  manc_lookup <- nb_all %>%
    dplyr::filter(!is.na(match_cell_type), match_cell_type != "",
                  !is.na(match_id), match_id != "") %>%
    dplyr::distinct(match_id, match_cell_type) %>%
    dplyr::rename(manc_id = match_id, cell_type = match_cell_type)
  message(sprintf("  MANC lookup: %d unique (manc_id, cell_type) pairs from NBLAST",
                  nrow(manc_lookup)))
} else {
  # Fallback: franken_meta
  manc_lookup <- fm %>%
    dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
    dplyr::distinct(manc_id, cell_type)
}

# Fallback: use manc_nblast_match (MANC body ID) from seatable to fill nblast_type
# where the feather had no hit. Look up the cell_type for the body ID.
if ("manc_nblast_match" %in% names(type_changed_df)) {
  nblast_match_types <- manc_lookup %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::rename(nblast_type_fb = cell_type)
  type_changed_df <- type_changed_df %>%
    dplyr::left_join(nblast_match_types,
                     by = c("manc_nblast_match" = "manc_id")) %>%
    dplyr::mutate(nblast_type = dplyr::coalesce(nblast_type, nblast_type_fb)) %>%
    dplyr::select(-nblast_type_fb)
}

###########################################
### Neuroglancer URLs for type changes  ###
###########################################

# Build neuroglancer links for each type-changed neuron.
# Each link shows:
#   1. The BANC neuron (root_626 / static ID) in the BANC segmentation layer
#   2. MANC neurons of the seatable (old) cell type in one MANC layer
#   3. MANC neurons of the new (CSV) cell type in another MANC layer
# This lets the reviewer compare old vs new MANC morphological matches.
message(sprintf("  Building neuroglancer links for %d type-changed neurons...",
                nrow(type_changed_df)))

# Decode the base URL once — this URL has two MANC imported-mesh layers
ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6384965498437632"
ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                    return = "text", cache = TRUE)
ngl_base <- fafbseg::ngl_decode_scene(
  fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

# Find layer indices by name
ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))
banc_layer_idx <- match("v626 neurons", ngl_ls$name)
if (is.na(banc_layer_idx)) banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)

# MANC layers: seatable type, new type, and NBLAST type
manc_layer_idxs <- grep("manc", ngl_ls$name, ignore.case = TRUE)
manc_old_idx <- manc_layer_idxs[1]
manc_new_idx <- if (length(manc_layer_idxs) >= 2) manc_layer_idxs[2] else manc_layer_idxs[1]
manc_nblast_idx <- grep("nblast", ngl_ls$name, ignore.case = TRUE)
manc_nblast_idx <- if (length(manc_nblast_idx) > 0) manc_nblast_idx[1] else NA_integer_
message(sprintf("  BANC layer: '%s' [%d]", ngl_ls$name[banc_layer_idx], banc_layer_idx))
message(sprintf("  MANC layers: '%s' [%d] (seatable type), '%s' [%d] (new type)",
                ngl_ls$name[manc_old_idx], manc_old_idx,
                ngl_ls$name[manc_new_idx], manc_new_idx))
if (!is.na(manc_nblast_idx)) {
  message(sprintf("  NBLAST layer: '%s' [%d]", ngl_ls$name[manc_nblast_idx], manc_nblast_idx))
}

first_error <- NULL
type_changed_df$neuroglancer_url <- vapply(seq_len(nrow(type_changed_df)), function(i) {
  row <- type_changed_df[i, ]

  # BANC neuron: use banc_id (= root_626, static ID)
  banc_rid <- row$banc_id

  # MANC body IDs for seatable (old) type
  st_manc <- if (!is.na(row$seatable_type) && row$seatable_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$seatable_type) %>% dplyr::pull(manc_id)
  } else character(0)

  # MANC body IDs for new (CSV) type
  new_manc <- if (!is.na(row$effective_type) && row$effective_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$effective_type) %>% dplyr::pull(manc_id)
  } else character(0)

  # MANC body IDs for NBLAST top-hit type
  nblast_manc <- if (!is.na(row$nblast_type) && row$nblast_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$nblast_type) %>% dplyr::pull(manc_id)
  } else character(0)

  tryCatch({
    sc <- ngl_base
    # Set BANC segment
    sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
    sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
    # Set MANC seatable-type segments
    sc[["layers"]][[manc_old_idx]][["segments"]] <- as.character(st_manc)
    sc[["layers"]][[manc_old_idx]][["hiddenSegments"]] <- NULL
    # Set MANC new-type segments
    sc[["layers"]][[manc_new_idx]][["segments"]] <- as.character(new_manc)
    sc[["layers"]][[manc_new_idx]][["hiddenSegments"]] <- NULL
    # Set MANC NBLAST-type segments (if layer exists)
    if (!is.na(manc_nblast_idx) && length(nblast_manc) > 0) {
      sc[["layers"]][[manc_nblast_idx]][["segments"]] <- as.character(nblast_manc)
      sc[["layers"]][[manc_nblast_idx]][["hiddenSegments"]] <- NULL
    }
    as.character(sc)
  }, error = function(e) {
    if (is.null(first_error)) first_error <<- conditionMessage(e)
    NA_character_
  })
}, character(1))
if (!is.null(first_error)) message(sprintf("  bancsee first error: %s", first_error))

n_with_links <- sum(!is.na(type_changed_df$neuroglancer_url))
message(sprintf("  Generated %d/%d neuroglancer links", n_with_links, nrow(type_changed_df)))

#############################################
### Save CSVs: all changes + hemilineage  ###
#############################################

# CSV 1: All type changes (3440 neurons)
all_changes_csv <- type_changed_df %>%
  dplyr::transmute(
    banc_id,
    root_id,
    seatable_type,
    effective_type,
    nblast_type,
    old_hemilineage,
    new_hemilineage,
    hemi_changed,
    dimorphism,
    nblast_score,
    seatable_nblast,
    new_nblast,
    neuroglancer_url
  )

all_changes_csv <- all_changes_csv %>%
  dplyr::arrange(dplyr::desc(seatable_nblast - new_nblast))

csv_path_all <- "data/codex/vnc_type_changes.csv"
readr::write_csv(all_changes_csv, csv_path_all)
message(sprintf("  Saved all type changes CSV: %s (%d rows)", csv_path_all, nrow(all_changes_csv)))

# CSV 2: Hemilineage changes only (subset), sorted by new_nblast descending
hemi_csv <- all_changes_csv %>%
  dplyr::filter(hemi_changed == "hemilineage changed") %>%
  dplyr::select(-hemi_changed) %>%
  dplyr::arrange(dplyr::desc(new_nblast))

csv_path_hemi <- "data/codex/vnc_hemilineage_changes.csv"
readr::write_csv(hemi_csv, csv_path_hemi)
message(sprintf("  Saved hemilineage changes CSV: %s (%d rows)", csv_path_hemi, nrow(hemi_csv)))

# CSV 3: Effector type changes assessment
# For each effector neuron, compare seatable_type (from manc_png_match) vs
# MANC-derived new_type. Includes NBLAST scores, neuroglancer links,
# dimorphism, and manc_match_id/score for reference.
eff_type_changed <- eff %>%
  dplyr::filter(type_changed == "type changed")

if (nblast_available) {
  # Look up NBLAST scores for effector type-changed neurons
  eff_tc_ids <- unique(eff_type_changed$banc_id)
  eff_tc_rootids <- id_map$root_id[id_map$banc_id %in% eff_tc_ids]

  nb_eff <- nb_all %>%
    dplyr::filter(root_626 %in% eff_tc_ids |
                    pt_root_id %in% eff_tc_rootids) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
    dplyr::mutate(root_626 = dplyr::case_when(
      root_626 %in% eff_tc_ids ~ root_626,
      pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
      TRUE ~ root_626
    )) %>%
    dplyr::select(-pt_root_id)

  # Top NBLAST hit
  eff_top <- nb_eff %>%
    dplyr::group_by(root_626) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(banc_id = root_626,
                     nblast_score = score,
                     nblast_type = match_cell_type)

  # Best score for seatable type
  eff_st <- eff_type_changed %>%
    dplyr::transmute(banc_id, lookup_type = seatable_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_eff, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")

  # Best score for new type (from MANC match)
  eff_new <- eff_type_changed %>%
    dplyr::transmute(banc_id, lookup_type = new_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_eff, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

  eff_type_changed <- eff_type_changed %>%
    dplyr::left_join(eff_top, by = "banc_id") %>%
    dplyr::left_join(eff_st, by = "banc_id") %>%
    dplyr::left_join(eff_new, by = "banc_id")
} else {
  eff_type_changed$nblast_score <- NA_real_
  eff_type_changed$nblast_type <- NA_character_
  eff_type_changed$seatable_nblast <- NA_real_
  eff_type_changed$new_nblast <- NA_real_
}

# Fallback: use manc_nblast_match from seatable to fill nblast_type
if ("manc_nblast_match" %in% names(eff_type_changed)) {
  eff_nblast_fb <- manc_lookup %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::rename(nblast_type_fb = cell_type)
  eff_type_changed <- eff_type_changed %>%
    dplyr::left_join(eff_nblast_fb,
                     by = c("manc_nblast_match" = "manc_id")) %>%
    dplyr::mutate(nblast_type = dplyr::coalesce(nblast_type, nblast_type_fb)) %>%
    dplyr::select(-nblast_type_fb)
}

# Build neuroglancer links for effector type-changed neurons
message(sprintf("  Building neuroglancer links for %d effector type-changed neurons...",
                nrow(eff_type_changed)))

eff_first_error <- NULL
eff_type_changed$neuroglancer_url <- vapply(seq_len(nrow(eff_type_changed)), function(i) {
  row <- eff_type_changed[i, ]
  banc_rid <- row$banc_id

  # MANC body IDs for seatable (old) type
  st_manc <- if (!is.na(row$seatable_type) && row$seatable_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$seatable_type) %>% dplyr::pull(manc_id)
  } else character(0)

  # MANC body IDs for effective (new) type
  new_manc <- if (!is.na(row$effective_type) && row$effective_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$effective_type) %>% dplyr::pull(manc_id)
  } else character(0)

  # MANC body IDs for NBLAST top-hit type
  nblast_manc <- if (!is.na(row$nblast_type) && row$nblast_type != "") {
    manc_lookup %>% dplyr::filter(cell_type == row$nblast_type) %>% dplyr::pull(manc_id)
  } else character(0)

  tryCatch({
    sc <- ngl_base
    sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
    sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
    sc[["layers"]][[manc_old_idx]][["segments"]] <- as.character(st_manc)
    sc[["layers"]][[manc_old_idx]][["hiddenSegments"]] <- NULL
    sc[["layers"]][[manc_new_idx]][["segments"]] <- as.character(new_manc)
    sc[["layers"]][[manc_new_idx]][["hiddenSegments"]] <- NULL
    if (!is.na(manc_nblast_idx) && length(nblast_manc) > 0) {
      sc[["layers"]][[manc_nblast_idx]][["segments"]] <- as.character(nblast_manc)
      sc[["layers"]][[manc_nblast_idx]][["hiddenSegments"]] <- NULL
    }
    as.character(sc)
  }, error = function(e) {
    if (is.null(eff_first_error)) eff_first_error <<- conditionMessage(e)
    NA_character_
  })
}, character(1))
if (!is.null(eff_first_error)) message(sprintf("  Effector bancsee first error: %s", eff_first_error))
message(sprintf("  Generated %d/%d effector neuroglancer links",
                sum(!is.na(eff_type_changed$neuroglancer_url)), nrow(eff_type_changed)))

eff_changes_csv <- eff_type_changed %>%
  dplyr::transmute(
    banc_id,
    root_id,
    seatable_type,
    effective_type,
    nblast_type,
    manc_match_id,
    manc_match_score,
    dimorphism,
    dimorph_score,
    nblast_score,
    seatable_nblast,
    new_nblast,
    neuroglancer_url
  ) %>%
  dplyr::arrange(dplyr::desc(new_nblast))

csv_path_eff <- "data/codex/vnc_effector_type_changes.csv"
readr::write_csv(eff_changes_csv, csv_path_eff)
message(sprintf("  Saved effector type changes CSV: %s (%d rows)",
                csv_path_eff, nrow(eff_changes_csv)))

#########################################################
### CSV 4: maleCNS effector type changes assessment   ###
#########################################################

# For each BANC effector, compare current seatable malecns_cell_type
# vs the maleCNS cell type inferred from banc_to_mcns_effectors.csv.
# The CSV maps banc ID -> mcns match ID; we look up the match's
# cell_type in maleCNS metadata to get the "new" malecns type.
# Neuroglancer links use a separate scene with maleCNS layers.
message("  === maleCNS effector type changes ===")

mcns_eff <- readr::read_csv("data/codex/id_id_matches/banc_to_mcns_effectors.csv",
                             col_types = readr::cols(`banc ID` = "c",
                                                     `twin ID` = "c",
                                                     `mcns match ID` = "c"),
                             show_col_types = FALSE) %>%
  dplyr::rename(banc_id = `banc ID`,
                csv_type = type,
                dimorph_score = `dimorph score`,
                twin_id = `twin ID`,
                twin_score = `twin score`,
                mcns_match_id = `mcns match ID`,
                mcns_match_score = `mcns match score`)

# Load maleCNS metadata to map mcns_match_id -> cell_type
mcns_meta_gcs <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09/malecns_09_meta.feather"
mcns_meta_local <- file.path(nblast_cache_dir, basename(mcns_meta_gcs))
if (!file.exists(mcns_meta_local)) {
  message("  Downloading maleCNS metadata from GCS...")
  system2("gsutil", c("cp", mcns_meta_gcs, mcns_meta_local),
          stdout = TRUE, stderr = TRUE)
}
mcns_meta <- arrow::read_feather(mcns_meta_local,
                                  col_select = c("malecns_09_id", "cell_type"))

# Join mcns_match_id to maleCNS metadata to get the "new" malecns cell type
mcns_eff <- mcns_eff %>%
  dplyr::left_join(mcns_meta %>%
                     dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
                     dplyr::distinct(malecns_09_id, .keep_all = TRUE),
                   by = c("mcns_match_id" = "malecns_09_id")) %>%
  dplyr::rename(new_malecns_type = cell_type)

# Assign dimorphism from median dimorph score per type
mcns_eff_dimorph <- mcns_eff %>%
  dplyr::filter(!is.na(csv_type), csv_type != "") %>%
  dplyr::group_by(csv_type) %>%
  dplyr::summarise(median_dimorph = stats::median(dimorph_score, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(dimorphism = dplyr::case_when(
    median_dimorph > 4 ~ "sex-specific",
    median_dimorph > 2 ~ "dimorphic",
    TRUE ~ "isomorphic"
  ))

mcns_eff <- mcns_eff %>%
  dplyr::left_join(mcns_eff_dimorph %>% dplyr::select(csv_type, dimorphism),
                   by = "csv_type") %>%
  dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism), "isomorphic", dimorphism))

# Join to seatable for current malecns_cell_type
mcns_eff <- mcns_eff %>%
  dplyr::left_join(bc %>% dplyr::select(root_626, root_id, supervoxel_id,
                                         super_class,
                                         current_malecns_type = malecns_cell_type),
                   by = c("banc_id" = "root_626"))

# Detect type changes: current malecns_cell_type vs CSV-inferred
mcns_eff <- mcns_eff %>%
  dplyr::mutate(
    type_changed = dplyr::case_when(
      is.na(current_malecns_type) & !is.na(new_malecns_type) ~ "type changed",
      !is.na(current_malecns_type) & !is.na(new_malecns_type) &
        current_malecns_type != new_malecns_type ~ "type changed",
      TRUE ~ "type did not change"
    )
  )

message(sprintf("  maleCNS effectors: %d neurons, %d with maleCNS match type, %d type changes",
                nrow(mcns_eff),
                sum(!is.na(mcns_eff$new_malecns_type)),
                sum(mcns_eff$type_changed == "type changed")))

mcns_type_changed <- mcns_eff %>%
  dplyr::filter(type_changed == "type changed")

# NBLAST scores from banc_malecns_v0.9 feather
mcns_nblast_gcs <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_malecns_v0.9_nblast.feather"
mcns_nblast_local <- file.path(nblast_cache_dir, basename(mcns_nblast_gcs))
if (!file.exists(mcns_nblast_local)) {
  message("  Downloading maleCNS NBLAST feather from GCS...")
  system2("gsutil", c("cp", mcns_nblast_gcs, mcns_nblast_local),
          stdout = TRUE, stderr = TRUE)
}

mcns_nblast_available <- file.exists(mcns_nblast_local)
if (mcns_nblast_available && nrow(mcns_type_changed) > 0) {
  nb_mcns <- arrow::read_feather(mcns_nblast_local)
  message(sprintf("  maleCNS NBLAST feather: %d rows", nrow(nb_mcns)))

  mcns_tc_ids <- unique(mcns_type_changed$banc_id)
  mcns_tc_rootids <- if (exists("id_map")) {
    id_map$root_id[id_map$banc_id %in% mcns_tc_ids]
  } else character(0)

  nb_mcns_sub <- nb_mcns %>%
    dplyr::filter(root_626 %in% mcns_tc_ids |
                    pt_root_id %in% mcns_tc_rootids) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score)

  # Remap pt_root_id hits back to banc_id
  if (exists("rootid_to_bancid")) {
    nb_mcns_sub <- nb_mcns_sub %>%
      dplyr::mutate(root_626 = dplyr::case_when(
        root_626 %in% mcns_tc_ids ~ root_626,
        pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
        TRUE ~ root_626
      ))
  }
  nb_mcns_sub <- nb_mcns_sub %>% dplyr::select(-pt_root_id)

  # Top NBLAST hit (any type)
  mcns_top <- nb_mcns_sub %>%
    dplyr::group_by(root_626) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(banc_id = root_626,
                     nblast_score = score,
                     nblast_type = match_cell_type)

  # Best score for current malecns_cell_type
  mcns_st <- mcns_type_changed %>%
    dplyr::transmute(banc_id, lookup_type = current_malecns_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_mcns_sub, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")

  # Best score for new maleCNS type (from CSV match)
  mcns_new <- mcns_type_changed %>%
    dplyr::transmute(banc_id, lookup_type = new_malecns_type) %>%
    dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
    dplyr::inner_join(nb_mcns_sub, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == lookup_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

  mcns_type_changed <- mcns_type_changed %>%
    dplyr::left_join(mcns_top, by = "banc_id") %>%
    dplyr::left_join(mcns_st, by = "banc_id") %>%
    dplyr::left_join(mcns_new, by = "banc_id")

  rm(nb_mcns, nb_mcns_sub); gc()
} else {
  mcns_type_changed$nblast_score <- NA_real_
  mcns_type_changed$nblast_type <- NA_character_
  mcns_type_changed$seatable_nblast <- NA_real_
  mcns_type_changed$new_nblast <- NA_real_
}

# Build neuroglancer links for maleCNS type-changed effectors.
# Uses a separate scene template with maleCNS segmentation layers.
message(sprintf("  Building neuroglancer links for %d maleCNS effector type-changed neurons...",
                nrow(mcns_type_changed)))

# maleCNS body ID lookup by cell_type (for populating neuroglancer layers)
mcns_lookup <- mcns_meta %>%
  dplyr::filter(!is.na(cell_type), cell_type != "") %>%
  dplyr::mutate(mcns_id = as.character(malecns_09_id)) %>%
  dplyr::select(mcns_id, cell_type)

# Decode maleCNS neuroglancer scene template
mcns_ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6343573275410432"
mcns_ngl_url2 <- sub("#!middleauth+", "?", mcns_ngl_url, fixed = TRUE)
mcns_ngl_parts <- unlist(strsplit(mcns_ngl_url2, "?", fixed = TRUE))
mcns_ngl_json <- fafbseg::flywire_fetch(mcns_ngl_parts[2], token = bancr:::banc_token(),
                                          return = "text", cache = TRUE)
mcns_ngl_base <- fafbseg::ngl_decode_scene(
  fafbseg::ngl_encode_url(mcns_ngl_json, baseurl = mcns_ngl_parts[1]))

# Find layer indices in maleCNS scene (3 maleCNS layers: old, new, nblast)
mcns_ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(mcns_ngl_base))
mcns_banc_layer_idx <- match("v626 neurons", mcns_ngl_ls$name)
if (is.na(mcns_banc_layer_idx)) mcns_banc_layer_idx <- match("segmentation proofreading", mcns_ngl_ls$name)

mcns_layer_idxs <- grep("malecns|male.?cns", mcns_ngl_ls$name, ignore.case = TRUE)
mcns_old_idx <- if (length(mcns_layer_idxs) >= 1) mcns_layer_idxs[1] else NA_integer_
mcns_new_idx <- if (length(mcns_layer_idxs) >= 2) mcns_layer_idxs[2] else mcns_old_idx
mcns_nblast_layer_idx <- if (length(mcns_layer_idxs) >= 3) mcns_layer_idxs[3] else NA_integer_

message(sprintf("  maleCNS scene layers: %s", paste(mcns_ngl_ls$name, collapse = ", ")))
if (!is.na(mcns_banc_layer_idx))
  message(sprintf("  BANC layer: '%s' [%d]", mcns_ngl_ls$name[mcns_banc_layer_idx], mcns_banc_layer_idx))
if (!is.na(mcns_old_idx))
  message(sprintf("  maleCNS layers: '%s' [%d] (current), '%s' [%d] (new), '%s' [%d] (nblast)",
                  mcns_ngl_ls$name[mcns_old_idx], mcns_old_idx,
                  mcns_ngl_ls$name[mcns_new_idx], mcns_new_idx,
                  if (!is.na(mcns_nblast_layer_idx)) mcns_ngl_ls$name[mcns_nblast_layer_idx] else "none",
                  ifelse(is.na(mcns_nblast_layer_idx), -1L, mcns_nblast_layer_idx)))

mcns_first_error <- NULL
mcns_type_changed$neuroglancer_url <- vapply(seq_len(nrow(mcns_type_changed)), function(i) {
  row <- mcns_type_changed[i, ]
  banc_rid <- row$banc_id

  # maleCNS body IDs for current malecns_cell_type
  st_mcns <- if (!is.na(row$current_malecns_type) && row$current_malecns_type != "") {
    mcns_lookup %>% dplyr::filter(cell_type == row$current_malecns_type) %>% dplyr::pull(mcns_id)
  } else character(0)

  # maleCNS body IDs for new maleCNS type
  new_mcns <- if (!is.na(row$new_malecns_type) && row$new_malecns_type != "") {
    mcns_lookup %>% dplyr::filter(cell_type == row$new_malecns_type) %>% dplyr::pull(mcns_id)
  } else character(0)

  # maleCNS body IDs for NBLAST top-hit type
  nblast_mcns <- if (!is.na(row$nblast_type) && row$nblast_type != "") {
    mcns_lookup %>% dplyr::filter(cell_type == row$nblast_type) %>% dplyr::pull(mcns_id)
  } else character(0)

  tryCatch({
    sc <- mcns_ngl_base
    if (!is.na(mcns_banc_layer_idx)) {
      sc[["layers"]][[mcns_banc_layer_idx]][["segments"]] <- as.character(banc_rid)
      sc[["layers"]][[mcns_banc_layer_idx]][["hiddenSegments"]] <- NULL
    }
    if (!is.na(mcns_old_idx)) {
      sc[["layers"]][[mcns_old_idx]][["segments"]] <- as.character(st_mcns)
      sc[["layers"]][[mcns_old_idx]][["hiddenSegments"]] <- NULL
    }
    if (!is.na(mcns_new_idx)) {
      sc[["layers"]][[mcns_new_idx]][["segments"]] <- as.character(new_mcns)
      sc[["layers"]][[mcns_new_idx]][["hiddenSegments"]] <- NULL
    }
    if (!is.na(mcns_nblast_layer_idx) && length(nblast_mcns) > 0) {
      sc[["layers"]][[mcns_nblast_layer_idx]][["segments"]] <- as.character(nblast_mcns)
      sc[["layers"]][[mcns_nblast_layer_idx]][["hiddenSegments"]] <- NULL
    }
    as.character(sc)
  }, error = function(e) {
    if (is.null(mcns_first_error)) mcns_first_error <<- conditionMessage(e)
    NA_character_
  })
}, character(1))
if (!is.null(mcns_first_error)) message(sprintf("  maleCNS bancsee first error: %s", mcns_first_error))
message(sprintf("  Generated %d/%d maleCNS effector neuroglancer links",
                sum(!is.na(mcns_type_changed$neuroglancer_url)), nrow(mcns_type_changed)))

mcns_changes_csv <- mcns_type_changed %>%
  dplyr::transmute(
    banc_id,
    root_id,
    current_malecns_type,
    new_malecns_type,
    nblast_type,
    mcns_match_id,
    mcns_match_score,
    dimorphism,
    dimorph_score,
    nblast_score,
    seatable_nblast,
    new_nblast,
    neuroglancer_url
  ) %>%
  dplyr::arrange(dplyr::desc(new_nblast))

csv_path_mcns <- "data/codex/mcns_effector_type_changes.csv"
readr::write_csv(mcns_changes_csv, csv_path_mcns)
message(sprintf("  Saved maleCNS effector type changes CSV: %s (%d rows)",
                csv_path_mcns, nrow(mcns_changes_csv)))

rm(mcns_eff, mcns_type_changed, mcns_changes_csv, mcns_meta); gc()

###############################################################
### Update manc_cell_type, manc_match, malecns_cell_type    ###
### for reviewed VNC type changes (intrinsic + effector).    ###
### Only neurons with accept_new == "T" in reviewed CSV.    ###
### COMMENTED OUT — uncomment to run the seatable update.    ###
###############################################################

reviewed_csv_path_update <- "data/codex/vnc_type_changes_reveiwed.csv"
if (file.exists(reviewed_csv_path_update) && nblast_available) {
  message("=== SeaTable update: reviewed VNC type changes ===")

  # Read reviewed CSV and filter to accepted type changes
  accepted <- readr::read_csv(reviewed_csv_path_update,
                               col_types = readr::cols(banc_id = "c", root_id = "c",
                                                        .default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "T") %>%
    dplyr::select(banc_id, effective_type)

  # Join side, region, super_class, _id, current cell_type/manc_cell_type/malecns_cell_type
  accepted <- accepted %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`, super_class, flow, side, region,
                                           current_cell_type = cell_type,
                                           current_manc_cell_type = manc_cell_type,
                                           current_malecns_cell_type = malecns_cell_type,
                                           seatable_type),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Filter: VNC super_classes only (intrinsic, motor, visceral_circulatory)
  accepted <- accepted %>%
    dplyr::filter(grepl("ventral_nerve_cord_intrinsic|motor|visceral_circulatory",
                        super_class, ignore.case = TRUE))

  # Filter: exclude central_brain neurons
  accepted <- accepted %>%
    dplyr::filter(is.na(region) | !grepl("central_brain", region))

  # # TEST FILTER: effectors only — comment out to run all VNC neurons
  # accepted <- accepted %>%
  #   dplyr::filter(flow == "efferent", region == "ventral_nerve_cord")

  message(sprintf("  %d accepted neurons after super_class + region filters", nrow(accepted)))

  # --- manc_match: Tier 1 — best per-neuron NBLAST hit for effective_type ---
  accepted_ids <- unique(accepted$banc_id)
  accepted_rootids <- id_map$root_id[id_map$banc_id %in% accepted_ids]

  nb_accepted <- nb_all %>%
    dplyr::filter(root_626 %in% accepted_ids |
                    pt_root_id %in% accepted_rootids) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
    dplyr::mutate(
      match_id = as.character(match_id),
      root_626 = dplyr::case_when(
        root_626 %in% accepted_ids ~ root_626,
        pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
        TRUE ~ root_626
      )) %>%
    dplyr::select(-pt_root_id)

  nblast_match <- accepted %>%
    dplyr::select(banc_id, effective_type) %>%
    dplyr::inner_join(nb_accepted, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == effective_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(banc_id, manc_match_nblast = match_id)

  message(sprintf("  manc_match: %d/%d from NBLAST",
                  nrow(nblast_match), nrow(accepted)))

  # --- manc_match: Tier 2 — franken_meta fallback (side-aware) ---
  fm_pool <- fm %>%
    dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
    dplyr::mutate(manc_id = as.character(manc_id)) %>%
    dplyr::select(manc_id, cell_type, fm_side = side)

  no_nblast <- accepted %>%
    dplyr::filter(!banc_id %in% nblast_match$banc_id)

  # Try side-matched first
  fm_side <- no_nblast %>%
    dplyr::inner_join(fm_pool,
                      by = c("effective_type" = "cell_type", "side" = "fm_side")) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(banc_id, manc_match_fm = manc_id)

  # Any side for remaining
  still_no <- no_nblast %>%
    dplyr::filter(!banc_id %in% fm_side$banc_id)
  fm_any <- still_no %>%
    dplyr::inner_join(fm_pool %>% dplyr::distinct(cell_type, .keep_all = TRUE),
                      by = c("effective_type" = "cell_type")) %>%
    dplyr::transmute(banc_id, manc_match_fm = manc_id)

  fm_fallback <- dplyr::bind_rows(fm_side, fm_any)
  message(sprintf("  manc_match: %d/%d from franken_meta fallback (%d side-matched)",
                  nrow(fm_fallback), nrow(no_nblast), nrow(fm_side)))

  # Combine manc_match sources
  accepted <- accepted %>%
    dplyr::left_join(nblast_match, by = "banc_id") %>%
    dplyr::left_join(fm_fallback, by = "banc_id") %>%
    dplyr::mutate(manc_match = dplyr::coalesce(manc_match_nblast, manc_match_fm))

  n_no_match <- sum(is.na(accepted$manc_match))
  if (n_no_match > 0) {
    message(sprintf("  WARNING: %d neurons have no manc_match (NBLAST or franken_meta)",
                    n_no_match))
  }

  # --- malecns_cell_type from maleCNS metadata ---
  malecns_meta_gcs <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09/malecns_09_meta.feather"
  malecns_meta_local <- file.path(nblast_cache_dir, basename(malecns_meta_gcs))
  if (!file.exists(malecns_meta_local)) {
    message("  Downloading maleCNS metadata from GCS...")
    system2("gsutil", c("cp", malecns_meta_gcs, malecns_meta_local),
            stdout = TRUE, stderr = TRUE)
  }
  malecns_meta <- arrow::read_feather(malecns_meta_local,
                                       col_select = c("malecns_09_id", "cell_type",
                                                       "manc_cell_type"))

  # Type-level mapping: MANC cell_type name -> maleCNS cell_type name
  malecns_type_lookup <- malecns_meta %>%
    dplyr::filter(!is.na(manc_cell_type), manc_cell_type != "",
                  !is.na(cell_type), cell_type != "") %>%
    dplyr::distinct(manc_cell_type, .keep_all = TRUE) %>%
    dplyr::select(manc_cell_type, malecns_ct = cell_type)

  # All known maleCNS cell_type names (for direct-match fallback)
  all_malecns_types <- unique(malecns_meta$cell_type[!is.na(malecns_meta$cell_type)])

  accepted <- accepted %>%
    dplyr::left_join(malecns_type_lookup,
                     by = c("effective_type" = "manc_cell_type")) %>%
    dplyr::mutate(malecns_cell_type = dplyr::case_when(
      # If effective_type is itself a known maleCNS cell_type, use directly
      effective_type %in% all_malecns_types ~ effective_type,
      # Otherwise use the manc_cell_type -> maleCNS cell_type mapping
      !is.na(malecns_ct) ~ malecns_ct,
      TRUE ~ NA_character_
    )) %>%
    dplyr::select(-malecns_ct)

  n_malecns <- sum(!is.na(accepted$malecns_cell_type))
  message(sprintf("  malecns_cell_type resolved for %d/%d neurons",
                  n_malecns, nrow(accepted)))

  # --- Build push data frame ---
  # Conditional column updates:
  #   manc_cell_type: always overwrite with effective_type
  #   manc_match: always overwrite with new match
  #   malecns_cell_type: only update if currently blank/NA in seatable
  #   cell_type: overwrite with effective_type only if cell_type == seatable_type (old)
  push_df <- accepted %>%
    dplyr::transmute(
      `_id`,
      manc_match = as.character(manc_match),
      manc_cell_type = effective_type,
      # Only fill malecns_cell_type if it is currently blank/NA
      malecns_cell_type = dplyr::if_else(
        is.na(current_malecns_cell_type) | current_malecns_cell_type == "",
        malecns_cell_type,
        current_malecns_cell_type
      ),
      # Only update cell_type if it matches the old manc_cell_type (seatable_type)
      cell_type = dplyr::case_when(
        !is.na(seatable_type) & !is.na(current_cell_type) &
          current_cell_type == seatable_type ~ effective_type,
        TRUE ~ current_cell_type
      )
    ) %>%
    as.data.frame()

  # Replace NAs with empty string for SeaTable text columns
  push_df$manc_match[is.na(push_df$manc_match)] <- ""
  push_df$malecns_cell_type[is.na(push_df$malecns_cell_type)] <- ""
  push_df$cell_type[is.na(push_df$cell_type)] <- ""

  n_ct_update <- sum(push_df$cell_type != "" &
    push_df$cell_type != accepted$current_cell_type, na.rm = TRUE)
  n_malecns_update <- sum(
    is.na(accepted$current_malecns_cell_type) | accepted$current_malecns_cell_type == "",
    na.rm = TRUE)

  message(sprintf("  Pushing update for %d neurons:", nrow(push_df)))
  message(sprintf("    manc_match: %d non-empty", sum(push_df$manc_match != "")))
  message(sprintf("    manc_cell_type: %d (all overwritten)", nrow(push_df)))
  message(sprintf("    cell_type: %d updated (where old == seatable_type)", n_ct_update))
  message(sprintf("    malecns_cell_type: %d newly filled (was blank)",
                  n_malecns_update))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_df,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(accepted, nb_accepted, nblast_match, no_nblast, fm_pool,
     fm_side, fm_any, still_no, fm_fallback, push_df,
     malecns_meta, malecns_type_lookup); gc()
} else {
  if (!file.exists(reviewed_csv_path_update))
    message("  Skipping SeaTable update: reviewed CSV not found")
  if (!nblast_available)
    message("  Skipping SeaTable update: NBLAST data not available")
}

###############################################################
### Update manc_match for accept_new == "F" (old type better)
### Set manc_match = manc_png_match (keep existing assignment)
### COMMENTED OUT — uncomment to run the seatable update.
###############################################################

if (file.exists(reviewed_csv_path_update)) {
  message("=== SeaTable update: accept_new == 'F' (manc_match <- manc_png_match) ===")

  rejected <- readr::read_csv(reviewed_csv_path_update,
                               col_types = readr::cols(banc_id = "c", root_id = "c",
                                                        .default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "F") %>%
    dplyr::select(banc_id)

  # Join to get _id and current manc_png_match from seatable
  rejected <- rejected %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`, manc_png_match, flow, region),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Build push: set manc_match = manc_png_match
  push_reject <- rejected %>%
    dplyr::transmute(
      `_id`,
      manc_png_match,
      manc_match = dplyr::if_else(is.na(manc_png_match) | manc_png_match == "",
                                  "", as.character(manc_png_match))
    ) %>%
    dplyr::filter(manc_match != "") %>%
    as.data.frame()

  message(sprintf("  %d neurons: setting manc_match = manc_png_match", nrow(push_reject)))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_reject,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(rejected, push_reject); gc()
}

###############################################################
### Update status for accept_new == "A" (all wrong)
### Append MANC_PNG_MATCH_WRONG to the status column.
### COMMENTED OUT — uncomment to run the seatable update.
###############################################################

if (file.exists(reviewed_csv_path_update)) {
  message("=== SeaTable update: accept_new == 'A' (append MANC_PNG_MATCH_WRONG) ===")

  all_wrong <- readr::read_csv(reviewed_csv_path_update,
                                col_types = readr::cols(banc_id = "c", root_id = "c",
                                                         .default = "c"),
                                show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "A") %>%
    dplyr::select(banc_id)

  # Join to get _id and current status from seatable
  all_wrong <- all_wrong %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`, status, flow, region),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Append MANC_PNG_MATCH_WRONG to status (skip if already present)
  push_wrong <- all_wrong %>%
    dplyr::filter(!grepl("MANC_PNG_MATCH_WRONG", status, fixed = TRUE) | is.na(status)) %>%
    dplyr::mutate(status = append_status(status, "MANC_PNG_MATCH_WRONG")) %>%
    dplyr::select(`_id`, status) %>%
    as.data.frame()

  message(sprintf("  %d neurons: appending MANC_PNG_MATCH_WRONG to status", nrow(push_wrong)))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_wrong,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(all_wrong, push_wrong); gc()
}

###############################################################
### Newly typed VNC intrinsic: neurons with new_type from
### CSV but no prior seatable_type (no manc_png_match).
### These were not in the reviewed CSV because type_changed
### requires !is.na(seatable_type).
### Sets cell_type, manc_cell_type, manc_match, malecns_cell_type.
### COMMENTED OUT — uncomment to run the seatable update.
###############################################################

if (nblast_available) {
  message("=== SeaTable update: newly typed VNC intrinsic (no prior manc_png_match) ===")

  newly_typed <- df %>%
    dplyr::filter(!is.na(new_type), new_type != "",
                  is.na(seatable_type) | seatable_type == "") %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  message(sprintf("  %d VNC intrinsic neurons with new type but no prior manc_png_match",
                  nrow(newly_typed)))

  if (nrow(newly_typed) > 0) {
    # --- manc_match: Tier 1 — best NBLAST hit for new_type ---
    nt_ids <- unique(newly_typed$banc_id)
    nt_rootids <- id_map$root_id[id_map$banc_id %in% nt_ids]

    nb_nt <- nb_all %>%
      dplyr::filter(root_626 %in% nt_ids |
                      pt_root_id %in% nt_rootids) %>%
      dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
      dplyr::mutate(
        match_id = as.character(match_id),
        root_626 = dplyr::case_when(
          root_626 %in% nt_ids ~ root_626,
          pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
          TRUE ~ root_626
        )) %>%
      dplyr::select(-pt_root_id)

    nt_nblast_match <- newly_typed %>%
      dplyr::select(banc_id, new_type) %>%
      dplyr::inner_join(nb_nt, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == new_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(banc_id, manc_match_nblast = match_id)

    message(sprintf("  manc_match: %d/%d from NBLAST",
                    nrow(nt_nblast_match), nrow(newly_typed)))

    # --- manc_match: Tier 2 — franken_meta fallback (side-aware) ---
    nt_fm_pool <- fm %>%
      dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
      dplyr::mutate(manc_id = as.character(manc_id)) %>%
      dplyr::select(manc_id, cell_type, fm_side = side)

    no_nblast_nt <- newly_typed %>%
      dplyr::filter(!banc_id %in% nt_nblast_match$banc_id)

    fm_side_nt <- no_nblast_nt %>%
      dplyr::inner_join(nt_fm_pool,
                        by = c("new_type" = "cell_type", "side" = "fm_side")) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(banc_id, manc_match_fm = manc_id)

    still_no_nt <- no_nblast_nt %>%
      dplyr::filter(!banc_id %in% fm_side_nt$banc_id)
    fm_any_nt <- still_no_nt %>%
      dplyr::inner_join(nt_fm_pool %>% dplyr::distinct(cell_type, .keep_all = TRUE),
                        by = c("new_type" = "cell_type")) %>%
      dplyr::transmute(banc_id, manc_match_fm = manc_id)

    nt_fm_fallback <- dplyr::bind_rows(fm_side_nt, fm_any_nt)
    message(sprintf("  manc_match: %d/%d from franken_meta fallback (%d side-matched)",
                    nrow(nt_fm_fallback), nrow(no_nblast_nt), nrow(fm_side_nt)))

    # Combine manc_match sources
    newly_typed <- newly_typed %>%
      dplyr::left_join(nt_nblast_match, by = "banc_id") %>%
      dplyr::left_join(nt_fm_fallback, by = "banc_id") %>%
      dplyr::mutate(new_manc_match = dplyr::coalesce(manc_match_nblast, manc_match_fm))

    n_no_match <- sum(is.na(newly_typed$new_manc_match))
    if (n_no_match > 0) {
      message(sprintf("  WARNING: %d neurons have no manc_match", n_no_match))
    }

    # --- malecns_cell_type from maleCNS metadata ---
    nt_malecns_local <- file.path(nblast_cache_dir, "malecns_09_meta.feather")
    if (!file.exists(nt_malecns_local)) {
      message("  Downloading maleCNS metadata from GCS...")
      system2("gsutil", c("cp",
              "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09/malecns_09_meta.feather",
              nt_malecns_local), stdout = TRUE, stderr = TRUE)
    }
    nt_malecns_meta <- arrow::read_feather(nt_malecns_local,
                                            col_select = c("malecns_09_id", "cell_type",
                                                            "manc_cell_type"))

    nt_malecns_lookup <- nt_malecns_meta %>%
      dplyr::filter(!is.na(manc_cell_type), manc_cell_type != "",
                    !is.na(cell_type), cell_type != "") %>%
      dplyr::distinct(manc_cell_type, .keep_all = TRUE) %>%
      dplyr::select(manc_cell_type, malecns_ct = cell_type)

    nt_all_malecns_types <- unique(nt_malecns_meta$cell_type[!is.na(nt_malecns_meta$cell_type)])

    newly_typed <- newly_typed %>%
      dplyr::left_join(nt_malecns_lookup,
                       by = c("new_type" = "manc_cell_type")) %>%
      dplyr::mutate(new_malecns_cell_type = dplyr::case_when(
        new_type %in% nt_all_malecns_types ~ new_type,
        !is.na(malecns_ct) ~ malecns_ct,
        TRUE ~ NA_character_
      )) %>%
      dplyr::select(-malecns_ct)

    n_malecns <- sum(!is.na(newly_typed$new_malecns_cell_type))
    message(sprintf("  malecns_cell_type resolved for %d/%d neurons",
                    n_malecns, nrow(newly_typed)))

    # --- Build push data frame ---
    # Conditional updates:
    #   manc_cell_type: always set to new_type
    #   manc_match: always set to NBLAST/franken_meta match
    #   cell_type: set to new_type only if currently blank
    #   malecns_cell_type: set only if currently blank
    push_new <- newly_typed %>%
      dplyr::transmute(
        `_id`,
        manc_match = as.character(new_manc_match),
        manc_cell_type = new_type,
        cell_type = dplyr::if_else(
          is.na(cell_type) | cell_type == "",
          new_type, cell_type),
        malecns_cell_type = dplyr::if_else(
          is.na(malecns_cell_type) | malecns_cell_type == "",
          new_malecns_cell_type,
          malecns_cell_type)
      ) %>%
      as.data.frame()

    push_new$manc_match[is.na(push_new$manc_match)] <- ""
    push_new$malecns_cell_type[is.na(push_new$malecns_cell_type)] <- ""
    push_new$cell_type[is.na(push_new$cell_type)] <- ""

    n_ct_fill <- sum(
      (is.na(newly_typed$cell_type) | newly_typed$cell_type == "") &
        push_new$cell_type != "")
    n_malecns_fill <- sum(
      (is.na(newly_typed$malecns_cell_type) | newly_typed$malecns_cell_type == "") &
        push_new$malecns_cell_type != "")

    message(sprintf("  Pushing update for %d newly typed neurons:", nrow(push_new)))
    message(sprintf("    manc_match: %d non-empty", sum(push_new$manc_match != "")))
    message(sprintf("    manc_cell_type: %d (all set to new_type)", nrow(push_new)))
    message(sprintf("    cell_type: %d filled (was blank)", n_ct_fill))
    message(sprintf("    malecns_cell_type: %d filled (was blank)", n_malecns_fill))

    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = push_new,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  SeaTable update complete")

    rm(newly_typed, nb_nt, nt_nblast_match, no_nblast_nt, fm_side_nt,
       fm_any_nt, still_no_nt, nt_fm_fallback, nt_fm_pool, push_new,
       nt_malecns_meta, nt_malecns_lookup); gc()
  }
}

###############################################################
### maleCNS effector review: accept_new == "T"
### Update malecns_match, malecns_cell_type, and cell_type.
###############################################################

mcns_reviewed_csv_path <- "data/codex/mcns_effector_type_changes_reviewed.csv"
if (file.exists(mcns_reviewed_csv_path)) {
  message("=== SeaTable update: maleCNS effector reviewed type changes (T) ===")

  mcns_accepted <- readr::read_csv(mcns_reviewed_csv_path,
                                    col_types = readr::cols(banc_id = "c", root_id = "c",
                                                             mcns_match_id = "c",
                                                             .default = "c"),
                                    show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "T") %>%
    dplyr::select(banc_id, new_malecns_type, mcns_match_id)

  # Join to get _id, current cell_type and malecns_cell_type from seatable
  mcns_accepted <- mcns_accepted %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`,
                                           current_cell_type = cell_type,
                                           current_malecns_type = malecns_cell_type),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  message(sprintf("  %d accepted maleCNS effector neurons", nrow(mcns_accepted)))

  # Build push data frame:
  #   malecns_match: always set to mcns_match_id
  #   malecns_cell_type: always set to new_malecns_type
  #   cell_type: update only if blank/NA or == current_malecns_type (old assignment)
  push_mcns_t <- mcns_accepted %>%
    dplyr::transmute(
      `_id`,
      malecns_match = as.character(mcns_match_id),
      malecns_cell_type = new_malecns_type,
      cell_type = dplyr::case_when(
        is.na(current_cell_type) | current_cell_type == "" ~ new_malecns_type,
        !is.na(current_malecns_type) & current_cell_type == current_malecns_type ~ new_malecns_type,
        TRUE ~ current_cell_type
      )
    ) %>%
    as.data.frame()

  # Replace NAs with empty string for SeaTable text columns
  push_mcns_t$malecns_match[is.na(push_mcns_t$malecns_match)] <- ""
  push_mcns_t$malecns_cell_type[is.na(push_mcns_t$malecns_cell_type)] <- ""
  push_mcns_t$cell_type[is.na(push_mcns_t$cell_type)] <- ""

  n_ct_update <- sum(push_mcns_t$cell_type != mcns_accepted$current_cell_type, na.rm = TRUE)
  message(sprintf("  Pushing update for %d neurons:", nrow(push_mcns_t)))
  message(sprintf("    malecns_match: %d non-empty", sum(push_mcns_t$malecns_match != "")))
  message(sprintf("    malecns_cell_type: %d (all overwritten)", nrow(push_mcns_t)))
  message(sprintf("    cell_type: %d updated (was blank or == old malecns type)", n_ct_update))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_mcns_t,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(mcns_accepted, push_mcns_t); gc()
}

###############################################################
### maleCNS effector review: accept_new == "A"
### Append MALECNS_PNG_MATCH_WRONG to the status column.
### COMMENTED OUT — uncomment to run the seatable update.
###############################################################

if (file.exists(mcns_reviewed_csv_path)) {
  message("=== SeaTable update: maleCNS effector accept_new == 'A' (append MALECNS_PNG_MATCH_WRONG) ===")

  mcns_all_wrong <- readr::read_csv(mcns_reviewed_csv_path,
                                     col_types = readr::cols(banc_id = "c", root_id = "c",
                                                              .default = "c"),
                                     show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "A") %>%
    dplyr::select(banc_id)

  # Join to get _id and current status from seatable
  mcns_all_wrong <- mcns_all_wrong %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`, status),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Append MALECNS_PNG_MATCH_WRONG to status (skip if already present)
  push_mcns_wrong <- mcns_all_wrong %>%
    dplyr::filter(!grepl("MALECNS_PNG_MATCH_WRONG", status, fixed = TRUE) | is.na(status)) %>%
    dplyr::mutate(status = append_status(status, "MALECNS_PNG_MATCH_WRONG")) %>%
    dplyr::select(`_id`, status) %>%
    as.data.frame()

  message(sprintf("  %d neurons: appending MALECNS_PNG_MATCH_WRONG to status", nrow(push_mcns_wrong)))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_mcns_wrong,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(mcns_all_wrong, push_mcns_wrong); gc()
}

###############################################################
### Efferent VNC hemilineage update
### Infer hemilineage from manc_cell_type / malecns_cell_type
### using VNC name extraction + franken_meta lookup.
### Priority: manc_cell_type first, then malecns_cell_type.
###############################################################

message("=== Efferent VNC hemilineage update ===")

# Efferent VNC neurons with missing or TBD hemilineage
eff_vnc <- bc %>%
  dplyr::filter(flow == "efferent",
                grepl("ventral_nerve_cord", region)) %>%
  dplyr::filter(is.na(hemilineage) | hemilineage %in% c("", "TBD"))

message(sprintf("  %d efferent VNC neurons with missing/TBD hemilineage", nrow(eff_vnc)))

if (nrow(eff_vnc) > 0) {
  # franken_meta lookup table: cell_type -> hemilineage
  if (!exists("fm")) fm <- franken_meta()
  hemi_lookup_fm <- fm %>%
    dplyr::filter(!is.na(cell_type), cell_type != "",
                  !is.na(hemilineage), !hemilineage %in% c("", "TBD")) %>%
    dplyr::distinct(cell_type, .keep_all = TRUE) %>%
    dplyr::select(cell_type, hemilineage)

  # Try MANC cell type first
  eff_vnc <- eff_vnc %>%
    dplyr::mutate(
      hemi_manc_extract = extract_hemilineage(manc_cell_type)
    ) %>%
    dplyr::left_join(hemi_lookup_fm %>% dplyr::rename(hemi_manc_fm = hemilineage),
                     by = c("manc_cell_type" = "cell_type")) %>%
    dplyr::mutate(
      hemi_manc = dplyr::coalesce(hemi_manc_extract, hemi_manc_fm),
      hemi_manc = dplyr::if_else(!is.na(hemi_manc) & hemi_manc != "TBD",
                                  hemi_manc, NA_character_)
    )

  # Fall back to maleCNS cell type
  eff_vnc <- eff_vnc %>%
    dplyr::mutate(
      hemi_mcns_extract = extract_hemilineage(malecns_cell_type)
    ) %>%
    dplyr::left_join(hemi_lookup_fm %>% dplyr::rename(hemi_mcns_fm = hemilineage),
                     by = c("malecns_cell_type" = "cell_type")) %>%
    dplyr::mutate(
      hemi_mcns = dplyr::coalesce(hemi_mcns_extract, hemi_mcns_fm),
      hemi_mcns = dplyr::if_else(!is.na(hemi_mcns) & hemi_mcns != "TBD",
                                  hemi_mcns, NA_character_)
    )

  # Final: MANC first, then maleCNS
  eff_vnc <- eff_vnc %>%
    dplyr::mutate(new_hemilineage = dplyr::coalesce(hemi_manc, hemi_mcns))

  push_hemi <- eff_vnc %>%
    dplyr::filter(!is.na(new_hemilineage), new_hemilineage != "") %>%
    dplyr::select(`_id`, root_id, hemilineage = new_hemilineage) %>%
    as.data.frame()

  message(sprintf("  %d neurons can be assigned a hemilineage (%d from MANC, %d from maleCNS)",
                  nrow(push_hemi),
                  sum(!is.na(eff_vnc$hemi_manc)),
                  sum(is.na(eff_vnc$hemi_manc) & !is.na(eff_vnc$hemi_mcns))))

  if (nrow(push_hemi) > 0) {
    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = push_hemi,
    #                       append_allowed = FALSE,
    #                       chunksize = 1000)
    # message("  SeaTable hemilineage update complete")

    # Refresh bc with updated hemilineage values for downstream plots
    bc <- bc %>%
      dplyr::rows_update(push_hemi %>% dplyr::select(`_id`, hemilineage),
                         by = "_id")
  }

  rm(eff_vnc, push_hemi); gc()
}


###########################################################
### Reviewed type changes                                ###
###########################################################

# Read the manually reviewed type changes CSV.
# accept_new: T = new type better, F = old type better, A = all wrong
reviewed_csv_path <- "data/codex/vnc_type_changes_reveiwed.csv"
if (file.exists(reviewed_csv_path)) {
  message("  Reading reviewed type changes...")
  reviewed <- readr::read_csv(reviewed_csv_path,
                               col_types = readr::cols(banc_id = "c", root_id = "c",
                                                        .default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(
      accept_new = trimws(accept_new),
      seatable_nblast = as.numeric(seatable_nblast),
      new_nblast = as.numeric(new_nblast),
      nblast_score = as.numeric(nblast_score)
    ) %>%
    dplyr::filter(accept_new %in% c("T", "F", "A"))

  # Join super_class and sexually_dimorphic from seatable data
  reviewed <- reviewed %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, super_class, sexually_dimorphic),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::mutate(super_class = dplyr::if_else(is.na(super_class), "unknown", super_class))

  # Effector super_classes
  effector_classes <- c("motor", "visceral_circulatory", "endocrine")

  # For effectors that are female-specific in BANC, override to "All wrong"
  # (BANC is female — a female-specific effector type means the MANC match is wrong)
  n_eff_override <- sum(reviewed$super_class %in% effector_classes &
                          reviewed$sexually_dimorphic == "female-specific" &
                          reviewed$accept_new != "A", na.rm = TRUE)
  reviewed <- reviewed %>%
    dplyr::mutate(accept_new = dplyr::if_else(
      super_class %in% effector_classes & sexually_dimorphic == "female-specific",
      "A", accept_new))

  reviewed$accept_new <- factor(reviewed$accept_new,
    levels = c("T", "F", "A"),
    labels = c("Connectivity match better", "Morphology match better", "All wrong"))

  message(sprintf("  %d reviewed neurons: %s",
                  nrow(reviewed),
                  paste(table(reviewed$accept_new), names(table(reviewed$accept_new)),
                        collapse = ", ")))
  if (n_eff_override > 0) {
    message(sprintf("  %d effector neurons overridden to 'All wrong' (female-specific in BANC)",
                    n_eff_override))
  }
  # Effector summary
  eff_reviewed <- reviewed %>% dplyr::filter(super_class %in% effector_classes)
  if (nrow(eff_reviewed) > 0) {
    message(sprintf("  Effectors: %d total — %s",
                    nrow(eff_reviewed),
                    paste(table(eff_reviewed$accept_new),
                          names(table(eff_reviewed$accept_new)), collapse = ", ")))
  }

} else {
  message("  Reviewed CSV not found at ", reviewed_csv_path)
}

###########################################################
### Left-right mirror matches (banc_match)
###########################################################

message("### Left-right mirror matches ###")

# Read left-right match CSV: Ego ID -> Twin ID (both are root_626 values)
lr_matches <- readr::read_csv("data/codex/banc_vnc_intrinsic_left_right_matches.csv",
                               col_types = readr::cols(`Ego ID` = "c", `Twin ID` = "c",
                                                        Conf = "d"),
                               show_col_types = FALSE) %>%
  dplyr::rename(ego_id = `Ego ID`, twin_id = `Twin ID`, conf = Conf)

message(sprintf("  LR matches CSV: %d pairs", nrow(lr_matches)))

# Build lookup: root_626 -> supervoxel_id from bc (already loaded)
r626_to_svid <- bc %>%
  dplyr::filter(!is.na(supervoxel_id), supervoxel_id != "") %>%
  dplyr::select(root_626, supervoxel_id) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)

# Map: ego_id (root_626) -> _id for the row to update
# twin_id (root_626) -> banc_match (root_id for that row) + supervoxel_id
ego_lookup <- bc %>%
  dplyr::select(`_id`, root_626, root_id) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)

twin_lookup <- bc %>%
  dplyr::select(root_626, root_id, supervoxel_id) %>%
  dplyr::distinct(root_626, .keep_all = TRUE) %>%
  dplyr::rename(twin_root_id = root_id, twin_svid = supervoxel_id)

lr_update <- lr_matches %>%
  dplyr::inner_join(ego_lookup, by = c("ego_id" = "root_626")) %>%
  dplyr::inner_join(twin_lookup, by = c("twin_id" = "root_626")) %>%
  dplyr::select(`_id`, root_id, banc_match = twin_root_id,
                banc_match_supervoxel_id = twin_svid) %>%
  dplyr::filter(!is.na(banc_match), banc_match != "",
                !is.na(banc_match_supervoxel_id), banc_match_supervoxel_id != "") %>%
  as.data.frame()

message(sprintf("  %d/%d pairs mapped to seatable rows with valid banc_match + supervoxel_id",
                nrow(lr_update), nrow(lr_matches)))

# SeaTable upload (commented out — review first)
# banctable_update_rows(base = 'banc_meta', table = "banc_meta",
#                       df = lr_update, append_allowed = FALSE, chunksize = 1000)
# message("  Pushed banc_match + banc_match_supervoxel_id to SeaTable")


###########################################################
### Summary statistics                                   ###
###########################################################

message("\n=== Summary: VNC intrinsic neurons ===")

bc_vnc <- bc %>% dplyr::filter(super_class == "ventral_nerve_cord_intrinsic")
n_vnc_st <- nrow(bc_vnc)

# Q: VNC intrinsic with morphology match (manc_png_match)
n_Q <- bc_vnc %>%
  dplyr::filter(!is.na(manc_png_match), manc_png_match != "") %>% nrow()

# P: of those with morphology match, how many differ from top NBLAST hit
p_df <- bc_vnc %>%
  dplyr::filter(!is.na(seatable_type), seatable_type != "",
                !is.na(manc_nblast_match), manc_nblast_match != "") %>%
  dplyr::mutate(manc_nblast_match = as.character(manc_nblast_match)) %>%
  dplyr::left_join(fm_types %>% dplyr::rename(nblast_match_type = cell_type),
                   by = c("manc_nblast_match" = "manc_id"))
n_P <- sum(!is.na(p_df$nblast_match_type) &
           p_df$seatable_type != p_df$nblast_match_type, na.rm = TRUE)

# J: VNC intrinsic in connectivity-matching CSV
n_J <- nrow(df)

# K: type agreed (had morphology match AND connectivity confirmed same type)
n_K <- df %>%
  dplyr::filter(type_changed == "type did not change",
                !is.na(seatable_type), seatable_type != "") %>% nrow()

# L: type change suggestions
n_L <- sum(df$type_changed == "type changed")

# N: new matches (no prior morphology match)
n_N <- df %>%
  dplyr::filter(!is.na(new_type), new_type != "",
                is.na(seatable_type) | seatable_type == "") %>% nrow()

# X, Y, R: from manual review of type changes (VNC intrinsic only)
vnc_changed_ids <- df %>%
  dplyr::filter(type_changed == "type changed") %>%
  dplyr::pull(banc_id)

if (exists("reviewed") && is.data.frame(reviewed)) {
  vnc_rev <- reviewed %>%
    dplyr::filter(banc_id %in% vnc_changed_ids)
  n_X <- sum(vnc_rev$accept_new == "Connectivity match better")
  n_Y <- sum(vnc_rev$accept_new == "Morphology match better")
  n_R <- sum(vnc_rev$accept_new == "All wrong")
} else {
  n_X <- n_Y <- n_R <- NA_integer_
}

# H: agreement rate (among morphology-matched neurons: connectivity agreed or was better)
H <- round(100 * (n_K + n_X) / (n_K + n_L), 1)

# G: improvement from connectivity matching (additional correct matches)
G <- round(100 * (n_X + n_N) / n_J, 1)

message(sprintf(paste0(
  "\nIn our review of ventral nerve cord intrinsic neurons, we started by ",
  "matching neurons by morphology using NBLAST. Potential matches were then ",
  "reviewed by a human annotator. %d/%d (%.0f%%) ventral nerve cord intrinsic ",
  "neurons were matched to a MANC cell type. Of these, %d differed from ",
  "their top NBLAST hit (manc_png_match versus manc_nblast_match).\n\n",
  "Our connectivity matching suggested a match for %d ventral nerve cord ",
  "intrinsic neurons. Of these, %d already had a reviewed morphology match, ",
  "%d were type change suggestions and %d were new matches.\n\n",
  "We re-reviewed the %d type change suggestions, confirming %d as better ",
  "than the original morphology assignment, %d as too different to be a ",
  "cell type match, and %d as a poor match in both cases.\n\n",
  "This means that the total agreement rate of our connectivity matching ",
  "is %.1f%% [100*(K+X)/(K+L) = 100*(%d+%d)/(%d+%d)], and that our cell ",
  "type matching was improved by %.1f%% [100*(X+N)/J = 100*(%d+%d)/%d]. ",
  "We estimate that the human review took a total of ~50 person-hours."),
  n_Q, n_vnc_st, 100 * n_Q / n_vnc_st,
  n_P,
  n_J, n_K, n_L, n_N,
  n_L, n_X, n_Y, n_R,
  H, n_K, n_X, n_K, n_L,
  G, n_X, n_N, n_J))

message("\n=== Summary: Effector neurons ===")

# Effector neurons (motor + visceral_circulatory)
n_eff_total <- nrow(eff)
n_eff_png <- eff %>%
  dplyr::filter(!is.na(seatable_type), seatable_type != "") %>% nrow()
n_eff_manc_match <- eff %>%
  dplyr::filter(!is.na(new_type), new_type != "") %>% nrow()
n_eff_changed <- sum(eff$type_changed == "type changed")
n_eff_unchanged <- eff %>%
  dplyr::filter(type_changed == "type did not change",
                !is.na(seatable_type), seatable_type != "") %>% nrow()
n_eff_new <- eff %>%
  dplyr::filter(!is.na(new_type), new_type != "",
                is.na(seatable_type) | seatable_type == "") %>% nrow()

# Effector reviews: check separate effector review CSV first, then main reviewed CSV
eff_changed_ids <- eff %>%
  dplyr::filter(type_changed == "type changed") %>%
  dplyr::pull(banc_id)

eff_rev <- NULL
# Try separate effector review CSV (matches maleCNS pattern)
eff_reviewed_csv <- "data/codex/vnc_effector_type_changes_reveiwed.csv"
if (!file.exists(eff_reviewed_csv)) {
  # Fall back: user may have added accept_new column directly to the generated CSV
  eff_reviewed_csv_alt <- "data/codex/vnc_effector_type_changes.csv"
  if (file.exists(eff_reviewed_csv_alt)) {
    eff_rev_test <- readr::read_csv(eff_reviewed_csv_alt,
                                     col_types = readr::cols(banc_id = "c", .default = "c"),
                                     show_col_types = FALSE)
    if ("accept_new" %in% names(eff_rev_test)) {
      eff_reviewed_csv <- eff_reviewed_csv_alt
      message("  Reading effector reviews from: ", eff_reviewed_csv)
    }
    rm(eff_rev_test)
  }
}

if (file.exists(eff_reviewed_csv)) {
  eff_rev_all <- readr::read_csv(eff_reviewed_csv,
                                  col_types = readr::cols(banc_id = "c", root_id = "c",
                                                           .default = "c"),
                                  show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new %in% c("T", "F", "A"),
                  banc_id %in% eff_changed_ids)
  message(sprintf("  Effector reviews from CSV: %d/%d type-changed effectors reviewed",
                  nrow(eff_rev_all), length(eff_changed_ids)))
  if (nrow(eff_rev_all) > 0) {
    eff_rev_all$accept_new <- factor(eff_rev_all$accept_new,
      levels = c("T", "F", "A"),
      labels = c("Connectivity match better", "Morphology match better", "All wrong"))
    eff_rev <- eff_rev_all
  }
  rm(eff_rev_all)
} else if (exists("reviewed") && is.data.frame(reviewed)) {
  # Fall back to main reviewed CSV (in case effector reviews were combined there)
  eff_rev <- reviewed %>%
    dplyr::filter(banc_id %in% eff_changed_ids)
  message(sprintf("  Effector reviews from main reviewed CSV: %d/%d type-changed effectors found",
                  nrow(eff_rev), length(eff_changed_ids)))
}

if (!is.null(eff_rev) && nrow(eff_rev) > 0) {
  n_eff_X <- sum(eff_rev$accept_new == "Connectivity match better")
  n_eff_Y <- sum(eff_rev$accept_new == "Morphology match better")
  n_eff_R <- sum(eff_rev$accept_new == "All wrong")
  eff_H <- round(100 * (n_eff_unchanged + n_eff_X) / (n_eff_unchanged + n_eff_changed), 1)
} else {
  n_eff_X <- n_eff_Y <- n_eff_R <- NA_integer_
  eff_H <- NA_real_
}

message(sprintf(paste0(
  "\nFor effector neurons (motor + visceral/circulatory), %d neurons were ",
  "matched to a MANC cell type by morphology (manc_png_match). ",
  "Connectivity matching suggested a MANC match for %d/%d effector neurons. ",
  "Of these, %d confirmed the existing morphology match, %d were type ",
  "change suggestions, and %d were new matches (no prior morphology match)."),
  n_eff_png, n_eff_manc_match, n_eff_total,
  n_eff_unchanged, n_eff_changed, n_eff_new))

if (!is.na(n_eff_X)) {
  message(sprintf(paste0(
    "\nOf the %d effector type change suggestions, %d were confirmed as ",
    "better, %d as worse (old type preferred), and %d as poor matches in ",
    "both cases. Agreement rate: %.1f%%."),
    n_eff_changed, n_eff_X, n_eff_Y, n_eff_R, eff_H))
} else {
  message("  Effector type changes have not yet been reviewed.")
}

message("### banc: VNC type change analysis complete ###")

###############################################################
### Update manc_match / malecns_match from alignment CSVs    ###
###                                                          ###
### Reads ID-to-ID alignment files:                          ###
###   banc_manc_*_alignment_scores.csv (banc_id → manc_id)   ###
###   banc_mcns_*_alignment_scores.csv (banc_id → mcns_id)   ###
### BANC IDs map to root_626 in SeaTable.                    ###
###                                                          ###
### Only pushes if the cell_type of the matched MANC/maleCNS ###
### neuron agrees with manc_cell_type / malecns_cell_type    ###
### in BANC SeaTable. Otherwise the match is skipped.        ###
###                                                          ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################
message("=== Alignment-based manc_match / malecns_match update ===")

# --- Read alignment CSVs (no header) ---
manc_align_file <- "data/codex/id_id_matches/banc_manc_min_syn_1_vnc_with_mn_norm_lr_gdt_node_alignment_scores.csv"
mcns_align_file <- "data/codex/id_id_matches/banc_mcns_min_syn_1_vnc_with_mn_norm_lr_gdt_node_alignment_scores.csv"

if (!file.exists(manc_align_file) || !file.exists(mcns_align_file)) {
  message("  Skipping alignment match update: CSV files not found")
} else {

  manc_align <- readr::read_csv(manc_align_file,
                                 col_names = c("banc_id", "manc_id", "score1", "score2"),
                                 col_types = readr::cols(banc_id = "c", manc_id = "c",
                                                          score1 = "d", score2 = "d"),
                                 show_col_types = FALSE)
  mcns_align <- readr::read_csv(mcns_align_file,
                                 col_names = c("banc_id", "mcns_id", "score1", "score2"),
                                 col_types = readr::cols(banc_id = "c", mcns_id = "c",
                                                          score1 = "d", score2 = "d"),
                                 show_col_types = FALSE)

  message(sprintf("  Alignment rows: MANC=%d, maleCNS=%d",
                  nrow(manc_align), nrow(mcns_align)))

  # --- MANC: look up cell_type for each manc_id via franken_meta ---
  if (!exists("fm_types") || !is.data.frame(fm_types)) {
    fm_types <- franken_meta("SELECT manc_id, cell_type FROM franken_meta") %>%
      dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
      dplyr::mutate(manc_id = as.character(manc_id)) %>%
      dplyr::distinct(manc_id, .keep_all = TRUE)
  }

  manc_validated <- manc_align %>%
    # Join to BANC SeaTable to get _id, manc_cell_type, current manc_match
    dplyr::inner_join(bc %>% dplyr::select(root_626, `_id`, manc_cell_type,
                                            current_manc_match = manc_match),
                      by = c("banc_id" = "root_626")) %>%
    dplyr::filter(!is.na(`_id`), `_id` != "",
                  !is.na(manc_cell_type), manc_cell_type != "") %>%
    # Join to franken_meta to get the matched neuron's cell_type
    dplyr::left_join(fm_types %>% dplyr::rename(match_cell_type = cell_type),
                     by = c("manc_id" = "manc_id")) %>%
    # Only keep where cell_type matches manc_cell_type
    dplyr::filter(!is.na(match_cell_type), match_cell_type == manc_cell_type) %>%
    dplyr::distinct(banc_id, .keep_all = TRUE)

  message(sprintf("  MANC alignment matches with cell_type agreement: %d/%d",
                  nrow(manc_validated), nrow(manc_align)))

  # Only push rows where the match has actually changed
  manc_push <- manc_validated %>%
    dplyr::filter(is.na(current_manc_match) | current_manc_match == "" |
                    current_manc_match != manc_id) %>%
    dplyr::transmute(`_id`, root_626=banc_id, manc_match = manc_id) %>%
    as.data.frame()

  message(sprintf("  MANC matches to update (changed from current): %d", nrow(manc_push)))

  # --- maleCNS: look up cell_type for each mcns_id via maleCNS SeaTable ---
  mcns_types <- tryCatch({
    banctable_query("SELECT malecns_09_id, cell_type FROM malecns",
                    base = "cns_meta") %>%
      dplyr::filter(!is.na(malecns_09_id), !is.na(cell_type), cell_type != "") %>%
      dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
      dplyr::distinct(malecns_09_id, .keep_all = TRUE)
  }, error = function(e) {
    message("  Could not query maleCNS SeaTable: ", e$message)
    # Fallback: try maleCNS GCS meta
    if (exists("malecns_meta") && is.data.frame(malecns_meta)) {
      malecns_meta %>%
        dplyr::filter(!is.na(cell_type), cell_type != "") %>%
        dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
        dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
        dplyr::select(malecns_09_id, cell_type)
    } else {
      NULL
    }
  })

  if (!is.null(mcns_types) && nrow(mcns_types) > 0) {
    mcns_validated <- mcns_align %>%
      # Join to BANC SeaTable to get _id, malecns_cell_type, current malecns_match
      dplyr::inner_join(bc %>% dplyr::select(root_626, `_id`, malecns_cell_type,
                                              current_malecns_match = malecns_match),
                        by = c("banc_id" = "root_626")) %>%
      dplyr::filter(!is.na(`_id`), `_id` != "",
                    !is.na(malecns_cell_type), malecns_cell_type != "") %>%
      # Join to maleCNS types to get the matched neuron's cell_type
      dplyr::left_join(mcns_types %>% dplyr::rename(match_cell_type = cell_type),
                       by = c("mcns_id" = "malecns_09_id")) %>%
      # Only keep where cell_type matches malecns_cell_type
      dplyr::filter(!is.na(match_cell_type), match_cell_type == malecns_cell_type) %>%
      dplyr::distinct(banc_id, .keep_all = TRUE)

    message(sprintf("  maleCNS alignment matches with cell_type agreement: %d/%d",
                    nrow(mcns_validated), nrow(mcns_align)))

    # Only push rows where the match has actually changed
    mcns_push <- mcns_validated %>%
      dplyr::filter(is.na(current_malecns_match) | current_malecns_match == "" |
                      current_malecns_match != mcns_id) %>%
      dplyr::transmute(`_id`, root_626=banc_id, malecns_match = mcns_id) %>%
      as.data.frame()

    message(sprintf("  maleCNS matches to update (changed from current): %d", nrow(mcns_push)))
  } else {
    mcns_push <- data.frame(`_id` = character(0), malecns_match = character(0),
                             check.names = FALSE)
    message("  maleCNS types unavailable; skipping malecns_match update")
  }

  # --- Push to SeaTable ---
  if (nrow(manc_push) > 0) {
    message(sprintf("  Ready to push %d manc_match updates", nrow(manc_push)))
    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = manc_push,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  manc_match push complete")
  }

  if (nrow(mcns_push) > 0) {
    message(sprintf("  Ready to push %d malecns_match updates", nrow(mcns_push)))
    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = mcns_push,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  malecns_match push complete")
  }

  rm(manc_align, mcns_align, manc_validated, manc_push,
     mcns_validated, mcns_push, mcns_types); gc()
}

})


