#' banc-update-celltypes — Derive cell types from cross-dataset matches and push to SeaTable.
#'
#' For each matched dataset (FAFB, MANC, maleCNS, hemibrain, FANC), look
#' up the matched neuron's cell_type and set the dataset-specific column.
#' Also sets `cell_type` itself per the cascade (FAFB > MANC > hemibrain > FANC)
#' if not already set. Cross-references maleCNS via shared cell types.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - `franken_meta()`
#'   - maleCNS SeaTable
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `fafb_cell_type`, `manc_cell_type`,
#'     `malecns_cell_type`, `hemibrain_cell_type`, `fanc_cell_type`,
#'     `cell_type`, `cell_type_source`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`

###########################################################
### Derive cell types from cross-species matches and
### push to BANC seatable
###
### For each matched dataset (FAFB, MANC, maleCNS,
### hemibrain, FANC):
###   - Look up matched neuron's cell_type
###   - Set dataset-specific columns (fafb_cell_type, etc.)
###   - Set cell_type if not already set
###
### Cell type cascade: FAFB > MANC > hemibrain > FANC
### Also: maleCNS cross-referencing via shared cell types
###
### Pushes cell_type columns to seatable independently
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: updating cell types from matches ###")
t_start <- Sys.time()

###########################
### Read current state  ###
###########################

bc <- banctable_query()
bc <- bc %>%
  dplyr::filter(!is.na(root_id), root_id != "0") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Save original values for anti-join (only push rows that actually changed)
ct_orig_cols <- c("_id", "cell_type", "hemilineage", "super_class",
                  "fafb_cell_type", "manc_cell_type", "hemibrain_cell_type",
                  "fanc_cell_type", "malecns_cell_type")
ct_orig_cols <- intersect(ct_orig_cols, colnames(bc))
bc.ct.orig <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(dplyr::all_of(ct_orig_cols)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~dplyr::coalesce(as.character(.), "")))

# Read reference metadata
franken.meta <- franken_meta()
fw.meta <- franken.meta %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::rename(root_783 = fafb_id)
mc.meta <- franken.meta %>%
  dplyr::filter(!is.na(manc_id)) %>%
  dplyr::rename(bodyid = manc_id)
hb.meta <- suppressWarnings(readr::read_csv(file.path(banc.meta.save.path, "hemibrain_meta.csv"),
                                             col_types = hemibrainr:::sql_col_types))
fc.meta <- suppressWarnings(readr::read_csv(file.path(banc.meta.save.path, "fanc_meta.csv"),
                                             col_types = hemibrainr:::sql_col_types))
malecns.meta <- readr::read_csv(file.path(banc.meta.save.path, "malecns_09_meta.csv"),
                                 col_types = banc.col.types)
mcns.meta <- malecns.meta %>%
  dplyr::filter(!is.na(malecns_09_id)) %>%
  dplyr::rename(bodyid = malecns_09_id)

# Read local NBLAST-derived columns from banc_meta.csv
banc.meta <- readr::read_csv(file.path(banc.meta.save.path, "banc_meta.csv"),
                              col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

message("  Deriving FAFB cell types...")
###########################
### FAFB cell types     ###
###########################

fafb.cols <- colnames(banc.meta)[grepl("^fafb", colnames(banc.meta))]
fafb.cols <- setdiff(fafb.cols, c("fafb_nblast_match", "fafb_nblast"))
banc.meta.fafb <- banc.meta[, c("root_id", fafb.cols)]
banc.meta.fafb <- banc.meta.fafb[!duplicated(banc.meta.fafb$root_id), ]

bc <- bc %>%
  dplyr::left_join(banc.meta.fafb, by = "root_id", suffix = c("", ".meta")) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(fafb_match) ~
        fw.meta$cell_type[match(fafb_match, fw.meta$root_783)],
      TRUE ~ cell_type
    ),
    fafb_cell_type = dplyr::case_when(
      !is.na(fafb_match) & fafb_match != "" ~ fw.meta$cell_type[match(fafb_match, fw.meta$root_783)],
      is.na(fafb_cell_type) & "fafb_auto_cell_type" %in% colnames(banc.meta.fafb) ~ fafb_auto_cell_type,
      TRUE ~ fafb_cell_type
    ),
    top_nt = dplyr::case_when(
      (is.na(cell_type) | grepl("auto", cell_type)) & !is.na(fafb_match) ~
        fw.meta$top_nt[match(fafb_match, fw.meta$root_783)],
      TRUE ~ top_nt
    )
  ) %>%
  dplyr::select(-dplyr::ends_with(".meta"))

message("  Deriving MANC cell types...")
###########################
### MANC cell types     ###
###########################

manc.cols <- colnames(banc.meta)[grepl("^manc", colnames(banc.meta))]
manc.cols <- setdiff(manc.cols, c("manc_nblast_match", "manc_nblast"))
banc.meta.manc <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::all_of(manc.cols))

bc <- bc %>%
  dplyr::left_join(banc.meta.manc, by = "root_id", suffix = c("", ".meta")) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(manc_match) ~
        mc.meta$cell_type[match(manc_match, mc.meta$bodyid)],
      TRUE ~ cell_type
    ),
    manc_cell_type = dplyr::case_when(
      !is.na(manc_match) & manc_match != "" ~ mc.meta$cell_type[match(manc_match, mc.meta$bodyid)],
      is.na(manc_cell_type) & "manc_auto_cell_type" %in% colnames(banc.meta.manc) ~ manc_auto_cell_type,
      TRUE ~ manc_cell_type
    ),
    top_nt = dplyr::case_when(
      !grepl("auto", top_nt) & !is.na(top_nt) ~ top_nt,
      (is.na(cell_type) | grepl("auto", cell_type)) & !is.na(manc_match) ~
        mc.meta$top_nt[match(manc_match, mc.meta$bodyid)],
      TRUE ~ top_nt
    )
  ) %>%
  dplyr::select(-dplyr::ends_with(".meta"))

message("  Deriving maleCNS cell types...")
###########################
### maleCNS cell types  ###
###########################

malecns.cols <- colnames(banc.meta)[grepl("^malecns", colnames(banc.meta))]
malecns.cols <- setdiff(malecns.cols, c("malecns_nblast_match", "malecns_nblast"))
banc.meta.malecns <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::all_of(malecns.cols))

# Safe top_nt lookup — maleCNS metadata may lack this column
.mcns_top_nt <- if ("top_nt" %in% colnames(mcns.meta)) mcns.meta$top_nt else rep(NA_character_, nrow(mcns.meta))

bc <- bc %>%
  dplyr::left_join(banc.meta.malecns, by = "root_id", suffix = c("", ".meta")) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(malecns_match) ~
        mcns.meta$cell_type[match(malecns_match, mcns.meta$bodyid)],
      TRUE ~ cell_type
    ),
    malecns_cell_type = dplyr::case_when(
      !is.na(malecns_match) & malecns_match != "" ~ mcns.meta$cell_type[match(malecns_match, mcns.meta$bodyid)],
      is.na(malecns_cell_type) & "malecns_auto_cell_type" %in% colnames(banc.meta.malecns) ~ malecns_auto_cell_type,
      TRUE ~ malecns_cell_type
    ),
    top_nt = dplyr::case_when(
      !grepl("auto", top_nt) & !is.na(top_nt) ~ top_nt,
      (is.na(cell_type) | grepl("auto", cell_type)) & !is.na(malecns_match) ~
        .mcns_top_nt[match(malecns_match, mcns.meta$bodyid)],
      TRUE ~ top_nt
    )
  ) %>%
  dplyr::select(-dplyr::ends_with(".meta"))

message("  Deriving hemibrain cell types...")
###########################
### hemibrain cell types ###
###########################

hemibrain.cols <- colnames(banc.meta)[grepl("^hemibrain", colnames(banc.meta))]
hemibrain.cols <- setdiff(hemibrain.cols, c("hemibrain_nblast_match", "hemibrain_nblast"))
banc.meta.hemibrain <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::all_of(hemibrain.cols)) %>%
  dplyr::filter(!is.na(hemibrain_auto_cell_type))

bc <- bc %>%
  dplyr::left_join(banc.meta.hemibrain, by = "root_id", suffix = c("", ".meta")) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(hemibrain_match) ~
        hb.meta$cell_type[match(hemibrain_match, hb.meta$bodyid)],
      TRUE ~ cell_type
    ),
    hemibrain_cell_type = dplyr::case_when(
      !is.na(hemibrain_match) & hemibrain_match != "" ~
        hb.meta$cell_type[match(hemibrain_match, hb.meta$bodyid)],
      is.na(hemibrain_cell_type) & "hemibrain_auto_cell_type" %in% colnames(banc.meta.hemibrain) ~ hemibrain_auto_cell_type,
      TRUE ~ hemibrain_cell_type
    )
  ) %>%
  dplyr::select(-dplyr::ends_with(".meta"))

message("  Deriving FANC cell types...")
###########################
### FANC cell types     ###
###########################

fanc.cols <- colnames(banc.meta)[grepl("^fanc", colnames(banc.meta))]
fanc.cols <- setdiff(fanc.cols, c("fanc_nblast_match", "fanc_nblast"))
banc.meta.fanc <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::all_of(fanc.cols))

bc <- bc %>%
  dplyr::left_join(banc.meta.fanc, by = "root_id", suffix = c("", ".meta")) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(fanc_match) ~
        fc.meta$cell_type[match(fanc_match, fc.meta$cell_id)],
      TRUE ~ cell_type
    ),
    fanc_cell_type = dplyr::case_when(
      !is.na(fanc_match) & fanc_match != "" ~
        fc.meta$cell_type[match(fanc_match, fc.meta$cell_id)],
      is.na(fanc_cell_type) & "fanc_auto_cell_type" %in% colnames(banc.meta.fanc) ~ fanc_auto_cell_type,
      TRUE ~ fanc_cell_type
    )
  ) %>%
  dplyr::select(-dplyr::ends_with(".meta"))

message("  Applying cell type cascade (FAFB > MANC > hemibrain > FANC)...")
###########################
### Cell type cascade   ###
###########################

# Region-sanity guard for the cell_type cascade.
# Background: the cascade below inherits cell_type from fafb/manc/hemibrain/
# fanc_cell_type whenever a non-auto candidate is available. Without a
# region check this happily injects geometrically impossible labels — e.g.
# a Kenyon-cell or ORN cell_type onto a banc neuron whose region is
# optic_lobe, or a motor-neuron cell_type onto a brain neuron. Cases like
# this surfaced 2026-05-15 (27 newly-injected rows since the previous run).
# Fix: only allow a cross-source cell_type to be inherited when its
# franken_meta-recorded region(s) are compatible with banc's region. The
# check is intentionally permissive: when franken has no record for the
# cell_type, or when banc.region is unset, inheritance is allowed.
.expected_regions_by_ct <- franken.meta %>%
  dplyr::filter(!is.na(cell_type) & cell_type != "" &
                !is.na(region) & region != "") %>%
  dplyr::distinct(cell_type, region) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(.rl = paste0("|", paste(unique(region), collapse = "|"), "|"),
                   .groups = "drop")
.expected_regions <- setNames(.expected_regions_by_ct$.rl,
                              .expected_regions_by_ct$cell_type)

# Vectorised: TRUE iff banc.region is compatible with the inherited
# cell_type's franken-recorded region(s), or we have no info to contradict.
.region_compatible <- function(inherited_ct, banc_region) {
  bracketed <- .expected_regions[as.character(inherited_ct)]
  mapply(function(br, eb) {
    if (is.na(eb) || is.null(eb)) return(TRUE)
    if (is.na(br) || br == "") return(TRUE)
    grepl(paste0("|", br, "|"), eb, fixed = TRUE)
  }, banc_region, bracketed, USE.NAMES = FALSE)
}

# Final cascade: if cell_type still missing, try dataset-specific types.
# Every branch that inherits cell_type from a cross-source column is now
# additionally guarded by .region_compatible() so we never copy a
# geometrically wrong label onto a banc row.
bc <- bc %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(manc_cell_type) & grepl("neck", region) & grepl("ascending|sensory", super_class) &
      !grepl("auto", manc_cell_type) ~ manc_cell_type,
    !is.na(fafb_cell_type) & grepl("neck", region) & grepl("descending", super_class) &
      !grepl("auto", fafb_cell_type) ~ fafb_cell_type,
    !is.na(fafb_cell_type) & !grepl("auto", fafb_cell_type) &
      .region_compatible(fafb_cell_type, region) ~ fafb_cell_type,
    !is.na(manc_cell_type) & !grepl("auto", manc_cell_type) &
      .region_compatible(manc_cell_type, region) ~ manc_cell_type,
    !is.na(hemibrain_cell_type) & !grepl("auto", hemibrain_cell_type) &
      .region_compatible(hemibrain_cell_type, region) ~ hemibrain_cell_type,
    !is.na(fanc_cell_type) & !grepl("auto", fanc_cell_type) &
      .region_compatible(fanc_cell_type, region) ~ fanc_cell_type,
    !is.na(fafb_cell_type) &
      .region_compatible(fafb_cell_type, region) ~ fafb_cell_type,
    !is.na(manc_cell_type) &
      .region_compatible(manc_cell_type, region) ~ manc_cell_type,
    !is.na(hemibrain_cell_type) &
      .region_compatible(hemibrain_cell_type, region) ~ hemibrain_cell_type,
    !is.na(fanc_cell_type) &
      .region_compatible(fanc_cell_type, region) ~ fanc_cell_type,
    TRUE ~ cell_type
  ))

#########################################
### maleCNS cross-reference cell type ###
#########################################

# malecns.meta.ct <- malecns.meta %>%
#   dplyr::select(malecns_09_id, cell_type,
#                 mcns_fafb_cell_type = fafb_cell_type,
#                 mcns_hemibrain_cell_type = hemibrain_cell_type,
#                 mcns_manc_cell_type = manc_cell_type,
#                 mcns_manc_match = manc_match) %>%
#   dplyr::mutate(cell_type = dplyr::coalesce(cell_type, mcns_fafb_cell_type,
#                                              mcns_hemibrain_cell_type, mcns_manc_cell_type)) %>%
#   dplyr::filter(!is.na(cell_type))
# 
# bc$malecns_cell_type <- NA
# bc <- bc %>%
#   dplyr::left_join(malecns.meta.ct %>%
#                      dplyr::filter(!is.na(mcns_fafb_cell_type)) %>%
#                      dplyr::distinct(mcns_fafb_cell_type, .keep_all = TRUE) %>%
#                      dplyr::select(malecns_cell_type_from_fafb = cell_type, mcns_fafb_cell_type),
#                    by = c("fafb_cell_type" = "mcns_fafb_cell_type")) %>%
#   dplyr::left_join(malecns.meta.ct %>%
#                      dplyr::filter(!is.na(mcns_hemibrain_cell_type)) %>%
#                      dplyr::distinct(mcns_hemibrain_cell_type, .keep_all = TRUE) %>%
#                      dplyr::select(malecns_cell_type_from_hemibrain = cell_type, mcns_hemibrain_cell_type),
#                    by = c("hemibrain_cell_type" = "mcns_hemibrain_cell_type")) %>%
#   dplyr::left_join(malecns.meta.ct %>%
#                      dplyr::filter(!is.na(mcns_manc_cell_type)) %>%
#                      dplyr::distinct(mcns_manc_cell_type, .keep_all = TRUE) %>%
#                      dplyr::select(malecns_cell_type_from_manc = cell_type, mcns_manc_cell_type),
#                    by = c("manc_cell_type" = "mcns_manc_cell_type")) %>%
#   dplyr::left_join(malecns.meta.ct %>%
#                      dplyr::filter(!is.na(mcns_manc_match)) %>%
#                      dplyr::distinct(mcns_manc_match, .keep_all = TRUE) %>%
#                      dplyr::select(malecns_cell_type_from_manc_id = cell_type, mcns_manc_match),
#                    by = c("manc_match" = "mcns_manc_match")) %>%
#   dplyr::mutate(
#     malecns_cell_type = dplyr::coalesce(malecns_cell_type_from_fafb, malecns_cell_type_from_manc,
#                                          malecns_cell_type_from_manc_id, malecns_cell_type_from_hemibrain),
#     malecns_cell_type = dplyr::case_when(
#       cell_type %in% !!malecns.meta.ct$cell_type ~ cell_type,
#       TRUE ~ malecns_cell_type
#     )
#   ) %>%
#   dplyr::select(-dplyr::starts_with("malecns_cell_type_from_"))

# Clean auto cell types
bc <- bc %>%
  dplyr::mutate(cell_type = ifelse(grepl("auto", cell_type), NA, cell_type))

###############################################################################
### Cell type cascade re-derivation (2026-05-15)                            ###
###############################################################################
# If the current cell_type came from a specific dataset's prior match (i.e.
# cell_type equals the ORIGINAL *_cell_type for that dataset), and the new
# match in this run produced a different *_cell_type, propagate the new
# value to cell_type — UNLESS status carries a manual-override flag.
#
# Manual-override flags (cell_type stays untouched):
#   - HAS_MANUAL_ANNOTATION       : user typed cell_type by hand
#   - <DATASET>_MATCH_PREFERRED   : user manually pinned cell_type to a
#                                    specific dataset's match. We do NOT
#                                    auto-update even from the pinned
#                                    dataset's new match — once pinned,
#                                    the value is sacrosanct.
#
# This makes the pipeline robust to new PNG additions: the latest reviewed
# match replaces the cell_type derived from the older one for un-overridden
# rows, without manual intervention. Manually-overridden cell_types stay.
.has_manual_override <- function(status) {
  s <- ifelse(is.na(status), "", status)
  grepl("HAS_MANUAL_ANNOTATION|_MATCH_PREFERRED", s)
}

# Original values from before this run's *_cell_type overwrites
.orig_idx <- match(bc$`_id`, bc.ct.orig$`_id`)
.orig_cell_type <- bc.ct.orig$cell_type[.orig_idx]
.orig_fafb_ct   <- bc.ct.orig$fafb_cell_type[.orig_idx]
.orig_manc_ct   <- bc.ct.orig$manc_cell_type[.orig_idx]
.orig_hemi_ct   <- bc.ct.orig$hemibrain_cell_type[.orig_idx]
.orig_fanc_ct   <- bc.ct.orig$fanc_cell_type[.orig_idx]
.orig_mcns_ct   <- bc.ct.orig$malecns_cell_type[.orig_idx]

# Treat "" the same as NA for the equality test
.eq <- function(a, b) {
  ifelse(is.na(a) | a == "" | is.na(b) | b == "", FALSE, a == b)
}
.diff <- function(a, b) {
  ifelse(is.na(a) | a == "" | is.na(b) | b == "", FALSE, a != b)
}

.manual_override <- .has_manual_override(bc$status)

bc <- bc %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    .manual_override ~ cell_type,
    .eq(.orig_cell_type, .orig_fafb_ct) & .diff(fafb_cell_type, .orig_fafb_ct) ~ fafb_cell_type,
    .eq(.orig_cell_type, .orig_manc_ct) & .diff(manc_cell_type, .orig_manc_ct) ~ manc_cell_type,
    .eq(.orig_cell_type, .orig_hemi_ct) & .diff(hemibrain_cell_type, .orig_hemi_ct) ~ hemibrain_cell_type,
    .eq(.orig_cell_type, .orig_fanc_ct) & .diff(fanc_cell_type, .orig_fanc_ct) ~ fanc_cell_type,
    .eq(.orig_cell_type, .orig_mcns_ct) & .diff(malecns_cell_type, .orig_mcns_ct) ~ malecns_cell_type,
    TRUE ~ cell_type
  ))

.n_cascade <- sum(.orig_cell_type != bc$cell_type & !is.na(bc$cell_type), na.rm = TRUE)
.n_manual <- sum(.manual_override, na.rm = TRUE)
message(sprintf("  Cell type cascade: %d cell_types propagated from new matches (%d rows skipped due to manual-override status)",
                .n_cascade, .n_manual))

###############################################################################
### super_class inheritance from unverified high-NBLAST matches (2026-05-15) ##
###############################################################################
# For neurons that are proofread (or roughly_proofread) but have NO super_class
# assigned in SeaTable, adopt the super_class of their best NBLAST match —
# even if that match has not yet been manually reviewed via PNG. Only
# super_class is inherited; cell_type / cell_class / cell_sub_class are left
# untouched (the unverified match may be morphologically similar but not the
# same type — super_class is the broadest, safest level to inherit).
#
# Cascade order matches the cell_type cascade: FAFB > MANC > hemibrain > FANC
# > maleCNS. NBLAST score threshold: >= 0.5 (matches the historical
# "interesting match" floor used elsewhere in the pipeline).
.sc_threshold <- 0.5

# Per-dataset best (highest score) NBLAST hit per pt_root_id, regardless of
# validation status. banc.meta.<x>.nb feathers were loaded in
# banc-nblast-compile.R style above — but this script doesn't load them.
# Read them here.
.read_nb <- function(fname) {
  fp <- file.path(banc.meta.save.path, fname)
  if (!file.exists(fp)) return(NULL)
  tryCatch(arrow::read_feather(fp), error = function(e) NULL)
}
.fafb.nb   <- .read_nb("banc_fafb_783_nblast.feather")
.manc.nb   <- .read_nb("banc_manc_v1.2.1_nblast.feather")
.hemi.nb   <- .read_nb("banc_hemibrain_v1.2.1_nblast.feather")
.fanc.nb   <- .read_nb("banc_fanc_1116_nblast.feather")
.mcns.nb   <- .read_nb("banc_malecns_v0.9_nblast.feather")

.best_match <- function(nb_df) {
  if (is.null(nb_df) || !nrow(nb_df)) return(NULL)
  nb_df %>%
    dplyr::filter(!is.na(score), score >= .sc_threshold,
                  !is.na(match_id), match_id != "") %>%
    dplyr::group_by(pt_root_id) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(pt_root_id, match_id, score)
}

.fafb.best <- .best_match(.fafb.nb)
.manc.best <- .best_match(.manc.nb)
.hemi.best <- .best_match(.hemi.nb)
.fanc.best <- .best_match(.fanc.nb)
.mcns.best <- .best_match(.mcns.nb)

# Look up super_class from each source's meta. Some sources may lack the
# column — guard with %||% NA.
.lookup_sc <- function(best_df, source_meta, source_id_col) {
  if (is.null(best_df) || !nrow(best_df) || !("super_class" %in% colnames(source_meta))) {
    return(setNames(rep(NA_character_, 0), character(0)))
  }
  sc <- source_meta$super_class[match(best_df$match_id, source_meta[[source_id_col]])]
  setNames(sc, best_df$pt_root_id)
}

.fafb_sc <- .lookup_sc(.fafb.best, fw.meta, "root_783")
.manc_sc <- .lookup_sc(.manc.best, mc.meta, "bodyid")
.hemi_sc <- .lookup_sc(.hemi.best, hb.meta, "bodyid")
.fanc_sc <- .lookup_sc(.fanc.best, fc.meta, "cell_id")
.mcns_sc <- .lookup_sc(.mcns.best, mcns.meta, "bodyid")

# Eligible: proofread OR roughly_proofread, super_class empty, AND no
# HAS_MANUAL_ANNOTATION on status (don't auto-populate identity columns
# on neurons the user has flagged as manually annotated, even if a
# specific column happens to be empty).
.proofread_flag <- function(x) !is.na(x) & x %in% c("TRUE", "T", "t", "true", TRUE)
.eligible <- (.proofread_flag(bc$proofread) | .proofread_flag(bc$roughly_proofread)) &
              (is.na(bc$super_class) | bc$super_class == "") &
              !grepl("HAS_MANUAL_ANNOTATION",
                     ifelse(is.na(bc$status), "", bc$status))

# Cascade lookup for each eligible neuron
.new_sc <- rep(NA_character_, nrow(bc))
.new_sc <- ifelse(is.na(.new_sc), .fafb_sc[bc$root_id], .new_sc)
.new_sc <- ifelse(is.na(.new_sc), .manc_sc[bc$root_id], .new_sc)
.new_sc <- ifelse(is.na(.new_sc), .hemi_sc[bc$root_id], .new_sc)
.new_sc <- ifelse(is.na(.new_sc), .fanc_sc[bc$root_id], .new_sc)
.new_sc <- ifelse(is.na(.new_sc), .mcns_sc[bc$root_id], .new_sc)

.apply <- .eligible & !is.na(.new_sc) & .new_sc != ""
if (any(.apply)) {
  bc$super_class[.apply] <- .new_sc[.apply]
  message(sprintf("  super_class inheritance: assigned %d neurons from unverified NBLAST (threshold score >= %.2f)",
                  sum(.apply), .sc_threshold))
} else {
  message("  super_class inheritance: no eligible neurons found")
}

message("  Deriving hemilineage from cross-species matches...")
###########################
### Hemilineage cascade ###
###########################

# Empty hemilineage values
hemi_empty <- c(NA, "", "TBD")

# Strip "auto:" prefix for clean lookup values
strip_auto <- function(x) gsub("^auto:", "", x)

# Per-dataset hemilineage columns are in banc_meta.csv (from banc-meta.R):
#   fafb_ito_lee_hemilineage, manc_hemilineage, malecns_hemilineage,
#   hemibrain_hemilineage, fanc_ito_lee_hemilineage
# These have "auto:" prefix — strip for the unified hemilineage column.

# Also look up directly from match metadata for neurons with a match but
# missing the auto: column (e.g. match added after banc-meta.R ran)
# Safe hemilineage vectors — some reference datasets lack this column
.fw_hemi <- if ("hemilineage" %in% colnames(fw.meta)) fw.meta$hemilineage else rep(NA_character_, nrow(fw.meta))
.mcns_hemi <- if ("hemilineage" %in% colnames(mcns.meta)) mcns.meta$hemilineage else rep(NA_character_, nrow(mcns.meta))
.mc_hemi <- if ("hemilineage" %in% colnames(mc.meta)) mc.meta$hemilineage else rep(NA_character_, nrow(mc.meta))
.hb_hemi <- if ("hemilineage" %in% colnames(hb.meta)) hb.meta$hemilineage else rep(NA_character_, nrow(hb.meta))
.fc_hemi <- if ("hemilineage" %in% colnames(fc.meta)) fc.meta$hemilineage else rep(NA_character_, nrow(fc.meta))

bc <- bc %>%
  dplyr::mutate(
    # Direct lookup from match IDs where dataset hemilineage column is missing
    fafb_hemi_direct = dplyr::case_when(
      !is.na(fafb_match) & fafb_match != "" ~
        .fw_hemi[match(fafb_match, fw.meta$root_783)],
      TRUE ~ NA_character_
    ),
    malecns_hemi_direct = dplyr::case_when(
      !is.na(malecns_match) & malecns_match != "" ~
        .mcns_hemi[match(malecns_match, mcns.meta$bodyid)],
      TRUE ~ NA_character_
    ),
    manc_hemi_direct = dplyr::case_when(
      !is.na(manc_match) & manc_match != "" ~
        .mc_hemi[match(manc_match, mc.meta$bodyid)],
      TRUE ~ NA_character_
    ),
    hemibrain_hemi_direct = dplyr::case_when(
      !is.na(hemibrain_match) & hemibrain_match != "" ~
        .hb_hemi[match(hemibrain_match, hb.meta$bodyid)],
      TRUE ~ NA_character_
    ),
    fanc_hemi_direct = dplyr::case_when(
      !is.na(fanc_match) & fanc_match != "" ~
        .fc_hemi[match(fanc_match, fc.meta$cell_id)],
      TRUE ~ NA_character_
    )
  )

# Resolve per-dataset: prefer banc_meta.csv auto: column, fall back to direct match lookup
# Ensure hemilineage columns exist — some may not be in banc_meta after left_joins.
# case_when eagerly evaluates all formula LHS expressions, so column must exist even if
# the %in% colnames check would be FALSE.
for (.hcol in c("fafb_ito_lee_hemilineage", "malecns_hemilineage", "manc_hemilineage",
                "hemibrain_hemilineage", "fanc_ito_lee_hemilineage")) {
  if (!.hcol %in% colnames(bc)) bc[[.hcol]] <- NA_character_
}

bc <- bc %>%
  dplyr::mutate(
    fafb_hemi = dplyr::case_when(
      !is.na(fafb_ito_lee_hemilineage) &
        !fafb_ito_lee_hemilineage %in% hemi_empty ~ strip_auto(fafb_ito_lee_hemilineage),
      !is.na(fafb_hemi_direct) & !fafb_hemi_direct %in% hemi_empty ~ fafb_hemi_direct,
      TRUE ~ NA_character_
    ),
    malecns_hemi = dplyr::case_when(
      !is.na(malecns_hemilineage) &
        !malecns_hemilineage %in% hemi_empty ~ strip_auto(malecns_hemilineage),
      !is.na(malecns_hemi_direct) & !malecns_hemi_direct %in% hemi_empty ~ malecns_hemi_direct,
      TRUE ~ NA_character_
    ),
    manc_hemi = dplyr::case_when(
      !is.na(manc_hemilineage) &
        !manc_hemilineage %in% hemi_empty ~ strip_auto(manc_hemilineage),
      !is.na(manc_hemi_direct) & !manc_hemi_direct %in% hemi_empty ~ manc_hemi_direct,
      TRUE ~ NA_character_
    ),
    hemibrain_hemi = dplyr::case_when(
      !is.na(hemibrain_hemilineage) &
        !hemibrain_hemilineage %in% hemi_empty ~ strip_auto(hemibrain_hemilineage),
      !is.na(hemibrain_hemi_direct) & !hemibrain_hemi_direct %in% hemi_empty ~ hemibrain_hemi_direct,
      TRUE ~ NA_character_
    ),
    fanc_hemi = dplyr::case_when(
      !is.na(fanc_ito_lee_hemilineage) &
        !fanc_ito_lee_hemilineage %in% hemi_empty ~ strip_auto(fanc_ito_lee_hemilineage),
      !is.na(fanc_hemi_direct) & !fanc_hemi_direct %in% hemi_empty ~ fanc_hemi_direct,
      TRUE ~ NA_character_
    )
  )

# Cascade: FAFB > maleCNS > MANC > hemibrain > FANC
# Only update if current hemilineage is empty
bc <- bc %>%
  dplyr::mutate(
    hemilineage = dplyr::case_when(
      !is.na(hemilineage) & !hemilineage %in% hemi_empty ~ hemilineage,
      !is.na(fafb_hemi) ~ fafb_hemi,
      !is.na(malecns_hemi) ~ malecns_hemi,
      !is.na(manc_hemi) ~ manc_hemi,
      !is.na(hemibrain_hemi) ~ hemibrain_hemi,
      !is.na(fanc_hemi) ~ fanc_hemi,
      TRUE ~ NA_character_   # was: ~ hemilineage; clears TBD / placeholders
                              # that survive when no source has a real value
    )
  ) %>%
  dplyr::select(-dplyr::matches("_hemi_direct$|_hemi$"))

# Hemilineage hygiene guard (per the user 2026-05-15):
#  - Sensory super_classes (sensory / sensory_ascending / sensory_descending)
#    should not carry a hemilineage value — peripheral sensory neurons
#    derive from sensory organs, not central neuroblast lineages. Force
#    NA on the cascade output regardless of which source provided it.
#  - Support-cell super_classes (glia, trachea) also do not carry a
#    hemilineage — they are not neuroblast-lineage cells.
#  - The literal "TBD" is a placeholder, not a hemilineage. Force NA.
# These rules run AFTER the cascade so they override any source value.
bc <- bc %>%
  dplyr::mutate(
    hemilineage = dplyr::case_when(
      !is.na(super_class) & super_class %in%
        c("sensory", "sensory_ascending", "sensory_descending") ~ NA_character_,
      !is.na(super_class) & super_class %in%
        c("glia", "trachea") ~ NA_character_,
      !is.na(hemilineage) & hemilineage == "TBD" ~ NA_character_,
      TRUE ~ hemilineage
    )
  )

# Flow hygiene guard (per the user 2026-05-15):
#  - Support-cell super_classes (glia, trachea) should carry
#    flow = "intrinsic" (the bancpipeline taxonomy historically left
#    these as flow == NA but the project convention is now intrinsic).
# Only fires on these super_classes; any other row's flow is untouched
# by this step (set elsewhere from the per-dataset cascade / SeaTable).
if ("flow" %in% colnames(bc)) {
  bc <- bc %>%
    dplyr::mutate(flow = dplyr::case_when(
      !is.na(super_class) & super_class %in% c("glia", "trachea") ~ "intrinsic",
      TRUE ~ flow
    ))
}

message(sprintf("  Hemilineage assigned for %d neurons",
                sum(!is.na(bc$hemilineage) & !bc$hemilineage %in% hemi_empty)))

###########################
### Push to seatable    ###
###########################

ct_push_cols <- c("_id", "cell_type", "hemilineage", "super_class",
                  "fafb_cell_type", "manc_cell_type", "hemibrain_cell_type",
                  "fanc_cell_type", "malecns_cell_type")
ct_push_cols <- intersect(ct_push_cols, colnames(bc))

bc.ct.update <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(dplyr::all_of(ct_push_cols)) %>%
  as.data.frame()
bc.ct.update[is.na(bc.ct.update)] <- ""

# Anti-join: only push rows where cell type columns actually changed
bc.ct.update <- dplyr::anti_join(bc.ct.update, bc.ct.orig, by = intersect(colnames(bc.ct.update), colnames(bc.ct.orig)))
message(sprintf("Pushing cell type data for %d changed neurons to seatable (skipped %d unchanged)",
                nrow(bc.ct.update), nrow(bc.ct.orig) - nrow(bc.ct.update)))
if (nrow(bc.ct.update)) {
  banctable_update_rows(base = 'banc_meta',
                        table = 'banc_meta',
                        df = bc.ct.update,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

message(sprintf("### banc: cell type update complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
