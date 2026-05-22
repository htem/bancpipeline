#' banc-vnc-sexually-dimorphic — Propagate sexually_dimorphic labels across BANC, MANC, FAFB, maleCNS.
#'
#' Cross-dataset propagation: if a cell_type is `dimorphic` in any dataset
#' and matches a cell_type in another dataset (direct or via bridge columns),
#' the label is propagated. Hard-codes a small set of SNch types whose
#' dimorphism cannot be inferred from cross-dataset evidence.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`: cols `cell_type`, `sexually_dimorphic`,
#'     `manc_cell_type`, `malecns_cell_type`, `fafb_cell_type`
#'   - `franken_meta()` (MANC + FAFB)
#'   - maleCNS SeaTable (`manc_cell_type` bridge)
#'   - `data/codex/banc_vnc_cell_type_and_dimorphism.csv`
#'   - `data/codex/manc_vnc_cell_type_and_dimorphism.csv`
#'   - `data/codex/banc_effector_dimorphism.csv`
#'   - `data/codex/mcns_type_dimorph_scores_wrt_banc.csv`
#'
#' @section Writes:
#'   - SeaTable `banc_meta` / `franken_meta` / maleCNS:
#'     col `sexually_dimorphic` ∈ {dimorphic, isomorphic, female-specific, male-specific}
#'     (push calls COMMENTED OUT — uncomment to apply)
#'
#' @section Notes:
#'   - Hardcoded male-specific: SNch08, SNch13, SNch15, SNch16, SNch03.
#'   - Hardcoded dimorphic: SNch09, SNch07 (BANC: SNch09f / SNch07f).

###########################################################
### Update sexually_dimorphic labels across four datasets
###
### Data sources:
###   SeaTable cell_type + sexually_dimorphic columns from:
###     banc_meta (BANC, female) — with bridge columns:
###       manc_cell_type, malecns_cell_type, fafb_cell_type
###     franken_meta (MANC, male; FAFB, female)
###     malecns (maleCNS, male) — with manc_cell_type bridge
###   CSV feeder files in data/codex/:
###     banc_vnc_cell_type_and_dimorphism.csv (BANC types)
###     manc_vnc_cell_type_and_dimorphism.csv (MANC types)
###     banc_effector_dimorphism.csv (BANC effector neurons)
###     mcns_type_dimorph_scores_wrt_banc.csv (maleCNS types)
###
### Cross-dataset dimorphism propagation:
###   If a cell_type is "dimorphic" in any dataset and matches
###   cell_type in another dataset, the type gets "dimorphic"
###   in all datasets. Matching uses direct cell_type and BANC
###   bridge columns (manc_cell_type, malecns_cell_type,
###   fafb_cell_type) plus maleCNS manc_cell_type bridge.
###
### Hardcoded: SNch08/13/15/16/03 male-specific,
###            SNch09/07 dimorphic (BANC: SNch09f/07f)
###
### Allowed SeaTable values: dimorphic, isomorphic,
###   female-specific, male-specific
###
### Targets (SeaTable pushes COMMENTED OUT):
###   1. banc_meta.sexually_dimorphic   (base: banc_meta)
###   2. franken_meta.sexually_dimorphic (base: cns_meta)
###   3. malecns.sexually_dimorphic     (base: cns_meta)
###
### Protected labels: "female-specific" and "dimorphic" in
### banc_meta are never overwritten by this script.
###
### Usage: Rscript banc/matching/banc-vnc-sexually-dimorphic.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: updating sexually_dimorphic ###")
t_start <- Sys.time()

library(ggplot2)
library(dplyr)

# Output directory for plots and review CSVs
plot.dir <- file.path("inst", "images")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

# Labels that must not be overwritten in banc_meta
protected <- c("female-specific", "dimorphic")

###########################
### Query all datasets ###
###########################

message("\n=== Querying datasets from SeaTable ===")

# 1. BANC (banc_meta) — female; bridge columns map to other datasets
bc <- banctable_query_cached(paste0(
  "SELECT _id, root_id, root_626, region, cell_type, manc_cell_type, malecns_cell_type, ",
  "fafb_cell_type, super_class, cell_class, cell_sub_class, sexually_dimorphic ",
  "FROM banc_meta"))
if (!"root_id" %in% colnames(bc) && "banc_746_id" %in% colnames(bc)) {
  bc$root_id <- as.character(bc$banc_746_id)
  if (!"root_626" %in% colnames(bc)) bc$root_626 <- bc$root_id
  if (!"_id" %in% colnames(bc)) bc$`_id` <- paste0("gcs_", seq_len(nrow(bc)))
  if (!"sexually_dimorphic" %in% colnames(bc)) bc$sexually_dimorphic <- NA_character_
  message("  Using GCS feather fallback (banc_746_id -> root_id)")
}
bc <- bc %>%
  dplyr::filter(!is.na(root_id), root_id != "0") %>%
  dplyr::mutate(root_626 = as.character(root_626)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("  BANC: %d neurons", nrow(bc)))

# 2. MANC (franken_meta, dataset=MANC) — male
fm <- tryCatch({
  res <- banctable_query_cached(paste0(
    "SELECT _id, manc_id, cell_type, super_class, cell_class, cell_sub_class, ",
    "sexually_dimorphic, dataset FROM franken_meta"),
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  if (!"manc_id" %in% colnames(res)) {
    id_col <- intersect(c("manc_121_id", "bodyid"), colnames(res))[1]
    if (!is.na(id_col)) {
      res$manc_id <- as.character(res[[id_col]])
      message(sprintf("  Using GCS feather fallback (%s -> manc_id)", id_col))
    }
  }
  if (!"dataset" %in% colnames(res)) res$dataset <- "MANC"
  if (!"_id" %in% colnames(res)) res$`_id` <- paste0("gcs_", seq_len(nrow(res)))
  if (!"sexually_dimorphic" %in% colnames(res)) res$sexually_dimorphic <- NA_character_
  res %>%
    dplyr::filter(grepl("MANC", dataset, ignore.case = TRUE), !is.na(manc_id)) %>%
    dplyr::mutate(manc_id = as.character(manc_id)) %>%
    dplyr::distinct(manc_id, .keep_all = TRUE)
}, error = function(e) {
  warning("Could not query franken_meta: ", e$message)
  NULL
})
if (!is.null(fm)) message(sprintf("  MANC: %d neurons", nrow(fm)))

# 3. FAFB (franken_meta, dataset=FAFB/FlyWire) — female
fafb_df <- tryCatch({
  res <- banctable_query_cached(
    "SELECT cell_type, super_class, sexually_dimorphic, dataset FROM franken_meta",
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  res %>%
    dplyr::filter(grepl("FAFB|FlyWire", dataset, ignore.case = TRUE),
                  !is.na(cell_type), cell_type != "")
}, error = function(e) {
  warning("Could not load FAFB data: ", e$message)
  NULL
})
if (!is.null(fafb_df)) message(sprintf("  FAFB: %d neurons", nrow(fafb_df)))

# 4. maleCNS — male
mc <- tryCatch({
  res <- banctable_query_cached(paste0(
    "SELECT _id, cell_type, manc_cell_type, super_class, cell_class, ",
    "cell_sub_class, sexually_dimorphic FROM malecns"),
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  if (!"_id" %in% colnames(res)) res$`_id` <- paste0("gcs_", seq_len(nrow(res)))
  if (!"sexually_dimorphic" %in% colnames(res)) res$sexually_dimorphic <- NA_character_
  res %>%
    dplyr::filter(!is.na(cell_type), cell_type != "") %>%
    dplyr::distinct(cell_type, .keep_all = TRUE)
}, error = function(e) {
  warning("Could not query malecns: ", e$message)
  NULL
})
if (!is.null(mc)) message(sprintf("  maleCNS: %d neurons", nrow(mc)))

###############################################
### Read CSV feeder files                   ###
###############################################

message("\n=== Reading dimorphism from CSV feeder files ===")

.valid_labels <- c("dimorphic", "female-specific", "male-specific")
.dim_rank <- c("dimorphic" = 3, "female-specific" = 2, "male-specific" = 2)

# Normalize CSV labels to SeaTable convention (hyphens)
.normalize_dim <- function(x) {
  dplyr::case_when(
    x %in% c("known_dimorphic", "dimorphic") ~ "dimorphic",
    x %in% c("known_female_specific", "female_specific", "female-specific") ~ "female-specific",
    x %in% c("known_male_specific", "male_specific", "male-specific",
             "sex_specific") ~ "male-specific",
    TRUE ~ NA_character_
  )
}

# 1. BANC VNC types + dimorphism
banc_csv_path <- "data/codex/banc_vnc_cell_type_and_dimorphism.csv"
banc_csv_raw <- if (file.exists(banc_csv_path)) {
  readr::read_csv(banc_csv_path, col_types = banc.col.types,
                   show_col_types = FALSE) %>%
    dplyr::transmute(root_626 = ID, cell_type = Type,
                     dim_label = .normalize_dim(Dimorphism))
} else NULL
# Type-level (for unified_type_dim)
banc_csv <- if (!is.null(banc_csv_raw)) {
  banc_csv_raw %>%
    dplyr::select(cell_type, dim_label) %>%
    dplyr::filter(dim_label %in% .valid_labels)
} else NULL
# Neuron-level (for direct join to bc by root_626)
# Keeps ALL rows (not just non-isomorphic) so csv_cell_type can be used for type lookup
banc_csv_neurons <- if (!is.null(banc_csv_raw)) {
  banc_csv_raw %>%
    dplyr::transmute(root_626, csv_cell_type = cell_type,
                     dim_csv = dplyr::if_else(dim_label %in% .valid_labels,
                                              dim_label, NA_character_)) %>%
    dplyr::distinct(root_626, .keep_all = TRUE)
} else NULL
if (!is.null(banc_csv)) {
  n_with_dim <- sum(!is.na(banc_csv_neurons$dim_csv))
  message(sprintf("  BANC CSV: %d neurons total (%d unique types, %d with explicit dimorphism)",
    nrow(banc_csv_neurons), dplyr::n_distinct(banc_csv$cell_type), n_with_dim))
}

# 2. MANC VNC types + dimorphism
manc_csv_path <- "data/codex/manc_vnc_cell_type_and_dimorphism.csv"
manc_csv <- if (file.exists(manc_csv_path)) {
  readr::read_csv(manc_csv_path, col_types = banc.col.types,
                   show_col_types = FALSE) %>%
    dplyr::transmute(cell_type = Type,
                     dim_label = .normalize_dim(Dimorphism)) %>%
    dplyr::filter(dim_label %in% .valid_labels)
} else NULL
if (!is.null(manc_csv)) message(sprintf("  MANC CSV: %d non-isomorphic neurons", nrow(manc_csv)))

# # 3. BANC effector dimorphism (neuron-level, applied later)
# eff_csv_path <- "data/codex/banc_effector_dimorphism.csv"
# banc_eff_csv <- if (file.exists(eff_csv_path)) {
#   readr::read_csv(eff_csv_path, col_types = readr::cols(root_id = "c"),
#                    show_col_types = FALSE) %>%
#     dplyr::filter(sexually_dimorphic %in% .valid_labels)
# } else NULL
# if (!is.null(banc_eff_csv)) message(sprintf("  BANC effector CSV: %d neurons", nrow(banc_eff_csv)))

# 4. maleCNS type-level dimorphism scores
mcns_csv_path <- "data/codex/mcns_type_dimorph_scores_wrt_banc.csv"
mcns_csv <- if (file.exists(mcns_csv_path)) {
  readr::read_csv(mcns_csv_path, show_col_types = FALSE, col_types = banc.col.types) %>%
    dplyr::transmute(cell_type = type,
                     dim_label = .normalize_dim(dimorphism)) %>%
    dplyr::filter(dim_label %in% .valid_labels)
} else NULL
if (!is.null(mcns_csv)) message(sprintf("  maleCNS CSV: %d non-isomorphic types", nrow(mcns_csv)))

###############################################
### Build unified type-level dimorphism     ###
###############################################

message("\n=== Building type-level dimorphism from SeaTable + CSVs ===")

# Per-dataset type-level dimorphism (strongest label per cell_type)
.type_dim <- function(df, label_col = "sexually_dimorphic") {
  df %>%
    dplyr::filter(!is.na(cell_type), cell_type != "",
                  .data[[label_col]] %in% .valid_labels) %>%
    dplyr::mutate(.rank = .dim_rank[.data[[label_col]]]) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(cell_type, dim_label = .data[[label_col]])
}

# Helper: merge CSV type dims into SeaTable type dims (strongest label wins)
.merge_type_dims <- function(st_dim, csv_dim) {
  if (is.null(csv_dim) || nrow(csv_dim) == 0) return(st_dim)
  csv_td <- csv_dim %>%
    dplyr::mutate(.rank = .dim_rank[dim_label]) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell_type, dim_label)
  dplyr::bind_rows(st_dim, csv_td) %>%
    dplyr::mutate(.rank = .dim_rank[dim_label]) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell_type, dim_label)
}

# SeaTable type dims
banc_type_dim  <- .type_dim(bc)
manc_type_dim  <- if (!is.null(fm)) .type_dim(fm) else
                    data.frame(cell_type = character(0), dim_label = character(0))
fafb_type_dim  <- if (!is.null(fafb_df)) .type_dim(fafb_df) else
                    data.frame(cell_type = character(0), dim_label = character(0))
mcns_type_dim  <- if (!is.null(mc)) .type_dim(mc, "sexually_dimorphic") else
                    data.frame(cell_type = character(0), dim_label = character(0))

message(sprintf("  SeaTable type dims: BANC=%d, MANC=%d, FAFB=%d, maleCNS=%d",
                nrow(banc_type_dim), nrow(manc_type_dim),
                nrow(fafb_type_dim), nrow(mcns_type_dim)))

# Merge CSV type dims into SeaTable type dims
banc_type_dim <- .merge_type_dims(banc_type_dim, banc_csv)
manc_type_dim <- .merge_type_dims(manc_type_dim, manc_csv)
mcns_type_dim <- .merge_type_dims(mcns_type_dim, mcns_csv)

message(sprintf("  After CSV merge: BANC=%d, MANC=%d, FAFB=%d, maleCNS=%d",
                nrow(banc_type_dim), nrow(manc_type_dim),
                nrow(fafb_type_dim), nrow(mcns_type_dim)))

# Collect types marked "dimorphic" in any dataset
dimorphic_pool <- unique(c(
  banc_type_dim$cell_type[banc_type_dim$dim_label == "dimorphic"],
  manc_type_dim$cell_type[manc_type_dim$dim_label == "dimorphic"],
  fafb_type_dim$cell_type[fafb_type_dim$dim_label == "dimorphic"],
  mcns_type_dim$cell_type[mcns_type_dim$dim_label == "dimorphic"]
))

# Cross-dataset propagation using BANC bridge columns:
# If a BANC cell_type's linked type in another dataset is dimorphic → BANC type is dimorphic.
# If a BANC cell_type is dimorphic → all linked types in other datasets are dimorphic.
banc_bridges <- bc %>%
  dplyr::filter(!is.na(cell_type), cell_type != "") %>%
  dplyr::distinct(cell_type, manc_cell_type, malecns_cell_type, fafb_cell_type)

# Forward: linked type is dimorphic → BANC type is dimorphic
banc_linked <- banc_bridges %>%
  dplyr::filter(
    manc_cell_type %in% dimorphic_pool |
    malecns_cell_type %in% dimorphic_pool |
    fafb_cell_type %in% dimorphic_pool
  ) %>%
  dplyr::pull(cell_type)
dimorphic_pool <- unique(c(dimorphic_pool, banc_linked))

# Backward: BANC type is dimorphic → linked types are dimorphic
linked_back <- banc_bridges %>%
  dplyr::filter(cell_type %in% dimorphic_pool)
dimorphic_pool <- unique(c(dimorphic_pool,
  linked_back$manc_cell_type[!is.na(linked_back$manc_cell_type)],
  linked_back$malecns_cell_type[!is.na(linked_back$malecns_cell_type)],
  linked_back$fafb_cell_type[!is.na(linked_back$fafb_cell_type)]
))

# Also: maleCNS manc_cell_type bridge
if (!is.null(mc) && "manc_cell_type" %in% colnames(mc)) {
  mcns_dimorphic <- mc$cell_type[!is.na(mc$sexually_dimorphic) & mc$sexually_dimorphic == "dimorphic"]
  mcns_manc_linked <- mc$manc_cell_type[mc$cell_type %in% mcns_dimorphic &
                                          !is.na(mc$manc_cell_type)]
  dimorphic_pool <- unique(c(dimorphic_pool, mcns_manc_linked))
  # And back: manc types in pool → linked maleCNS types
  mcns_from_pool <- mc$cell_type[mc$manc_cell_type %in% dimorphic_pool &
                                   !is.na(mc$manc_cell_type)]
  dimorphic_pool <- unique(c(dimorphic_pool, mcns_from_pool))
}

# Build unified type-level lookup: dimorphic types from pool, others keep strongest label
all_type_labels <- dplyr::bind_rows(banc_type_dim, manc_type_dim, fafb_type_dim, mcns_type_dim)
unified_type_dim <- all_type_labels %>%
  dplyr::filter(!cell_type %in% dimorphic_pool) %>%
  dplyr::mutate(.rank = .dim_rank[dim_label]) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(cell_type, dimorphism = dim_label)
unified_type_dim <- dplyr::bind_rows(
  unified_type_dim,
  data.frame(cell_type = dimorphic_pool, dimorphism = "dimorphic",
             stringsAsFactors = FALSE)
) %>%
  dplyr::filter(!grepl("unknown|nerve",cell_type),
                !is.na(cell_type),
                cell_type!="")

message(sprintf("  Unified: %d types (%s)",
                nrow(unified_type_dim),
                paste(names(table(unified_type_dim$dimorphism)),
                      table(unified_type_dim$dimorphism), sep = "=", collapse = ", ")))
message(sprintf("  %d types propagated as dimorphic via cross-dataset matching",
                length(dimorphic_pool) - sum(all_type_labels$dim_label == "dimorphic")))

#############################################
### Hardcoded cell type dimorphism        ###
#############################################

# Manually curated dimorphism for specific sensory neuron types:
#   Male-specific: SNch08, SNch13, SNch15, SNch16, SNch03
#   Dimorphic: SNch09, SNch07 (in BANC with "f" suffix: SNch09f, SNch07f)
hardcoded_dim <- data.frame(
  cell_type = c("SNch08", "SNch13", "SNch15", "SNch16", "SNch03", "SNch09", "SNch07"),
  dimorphism = c("male-specific", "male-specific", "male-specific", "male-specific",
                 "male-specific", "dimorphic", "dimorphic"),
  stringsAsFactors = FALSE
)
unified_type_dim <- dplyr::bind_rows(
  unified_type_dim %>% dplyr::filter(!cell_type %in% hardcoded_dim$cell_type),
  hardcoded_dim
)
message(sprintf("  After hardcoded: %d types total", nrow(unified_type_dim)))

###########################
### fix_super_class     ###
###########################

# Fix super_class based on cell_type naming conventions:
#   AN* = ascending, DN* = descending, SN* = sensory, SA[letter] = sensory ascending
fix_super_class <- function(df, ct_col = "cell_type") {
  if (!ct_col %in% colnames(df)) return(df)
  ct <- df[[ct_col]]
  fixed <- dplyr::case_when(
    grepl("^AN", ct) ~ "ascending",
    grepl("^DN", ct) ~ "descending",
    grepl("^SN", ct) ~ "sensory",
    grepl("^SA[A-Za-z]", ct) ~ "ascending",
    TRUE ~ df$super_class
  )
  n_changed <- sum(fixed != df$super_class, na.rm = TRUE)
  if (n_changed > 0) {
    message(sprintf("  Fixed %d super_class assignments by cell_type prefix", n_changed))
  }
  df$super_class <- fixed
  df
}

###########################
### Apply fix_super_class ##
###########################

bc <- fix_super_class(bc, "cell_type")
if (!is.null(fm)) fm <- fix_super_class(fm, "cell_type")
if (!is.null(mc)) {
  mcns_ct_col <- if ("manc_cell_type" %in% colnames(mc)) "manc_cell_type" else "cell_type"
  mc <- fix_super_class(mc, mcns_ct_col)
}

###########################
### 1. BANC seatable   ###
###########################

message("\n=== 1. BANC banc_meta ===")

# Join unified type-level dimorphism via cell_type, bridge columns, and "f" suffix
bc.merged <- bc %>%
  dplyr::left_join(unified_type_dim %>% dplyr::rename(dim_ct = dimorphism),
                   by = "cell_type")

# Also match via manc_cell_type, malecns_cell_type, fafb_cell_type bridge columns
for (.bridge in c("manc_cell_type", "malecns_cell_type", "fafb_cell_type")) {
  .col <- paste0("dim_", sub("_cell_type$", "", .bridge))
  if (.bridge %in% colnames(bc.merged)) {
    bc.merged <- bc.merged %>%
      dplyr::left_join(unified_type_dim %>% dplyr::rename(!!.col := dimorphism),
                       by = stats::setNames("cell_type", .bridge))
  } else {
    bc.merged[[.col]] <- NA_character_
  }
}

# BANC cell_types may have "f" suffix (female variant, e.g. SNch09f -> SNch09)
bc.merged <- bc.merged %>%
  dplyr::mutate(.ct_stripped = sub("f$", "", cell_type)) %>%
  dplyr::left_join(unified_type_dim %>% dplyr::rename(dim_ct_stripped = dimorphism),
                   by = c(".ct_stripped" = "cell_type")) %>%
  dplyr::select(-.ct_stripped)

# BANC VNC CSV: neuron-level data (cell_type + explicit dimorphism, by root_626)
if (!is.null(banc_csv_neurons)) {
  bc.merged <- bc.merged %>%
    dplyr::left_join(banc_csv_neurons, by = "root_626")
  n_csv_match <- sum(!is.na(bc.merged$csv_cell_type))
  n_csv_dim <- sum(!is.na(bc.merged$dim_csv))
  message(sprintf("  %d neurons matched BANC VNC CSV by root_626 (%d with explicit dimorphism)",
                  n_csv_match, n_csv_dim))
} else {
  bc.merged$csv_cell_type <- NA_character_
  bc.merged$dim_csv <- NA_character_
}

# Look up unified_type_dim using CSV cell_type (catches neurons missing cell_type in SeaTable)
bc.merged <- bc.merged %>%
  dplyr::left_join(unified_type_dim %>% dplyr::rename(dim_csv_type = dimorphism),
                   by = c("csv_cell_type" = "cell_type"))

# Priority: CSV explicit dimorphism > unified type (SeaTable cell_type > CSV cell_type >
#           bridge > stripped) > terminalia sensory > existing metadata > isomorphic
bc.merged <- bc.merged %>%
  dplyr::mutate(
    is_terminalia_sensory = super_class == "sensory" &
      (grepl("terminalia", cell_class, ignore.case = TRUE) |
       grepl("terminalia", cell_sub_class, ignore.case = TRUE)),
    .raw_label = dplyr::case_when(
      !is.na(dim_csv) ~ dim_csv,
      !is.na(dim_ct) ~ dim_ct,
      !is.na(dim_csv_type) ~ dim_csv_type,
      !is.na(dim_manc) ~ dim_manc,
      !is.na(dim_malecns) ~ dim_malecns,
      !is.na(dim_fafb) ~ dim_fafb,
      !is.na(dim_ct_stripped) ~ dim_ct_stripped,
      is_terminalia_sensory ~ "dimorphic",
      sexually_dimorphic %in% .valid_labels ~ sexually_dimorphic,
      !is.na(cell_type) & cell_type != "" ~ "isomorphic",
      !is.na(csv_cell_type) & csv_cell_type != "" ~ "isomorphic",
      TRUE ~ NA_character_
    ),
    # BANC is female: a type present here can't be male-specific → dimorphic
    sexually_dimorphic_new = dplyr::if_else(
      .raw_label == "male-specific", "dimorphic", .raw_label)
  ) %>%
  dplyr::select(-dplyr::any_of(c("dim_csv", "csv_cell_type", "dim_csv_type",
                                  "dim_ct", "dim_manc", "dim_malecns", "dim_fafb",
                                  "dim_ct_stripped", "is_terminalia_sensory", ".raw_label")))

n_from_type <- sum(!is.na(bc.merged$sexually_dimorphic_new) &
                    bc.merged$sexually_dimorphic_new != "isomorphic", na.rm = TRUE)
n_terminalia <- sum(bc.merged$super_class == "sensory" &
                     (grepl("terminalia", bc.merged$cell_class, ignore.case = TRUE) |
                      grepl("terminalia", bc.merged$cell_sub_class, ignore.case = TRUE)),
                    na.rm = TRUE)
message(sprintf("  %d neurons with non-isomorphic labels from type-level matching", n_from_type))
message(sprintf("  %d terminalia sensory neurons -> dimorphic", n_terminalia))

# Respect protected labels: never overwrite "female-specific" or "dimorphic"
bc.changed <- bc.merged %>%
  dplyr::filter(
    !is.na(sexually_dimorphic_new),
    !sexually_dimorphic %in% protected
  )

message(sprintf("  %d neurons need updating (respecting %d protected labels)",
                nrow(bc.changed),
                sum(bc.merged$sexually_dimorphic %in% protected, na.rm = TRUE)))

# Save proposed changes for review
if (nrow(bc.changed) > 0) {
  review_csv <- file.path(plot.dir, "sexually_dimorphic_proposed_changes.csv")
  readr::write_csv(
    bc.changed %>% dplyr::select(`_id`, root_id, sexually_dimorphic, sexually_dimorphic_new),
    review_csv
  )
  message(sprintf("  Saved: %s (%d rows)", review_csv, nrow(bc.changed)))
}

### Push to banc_meta — COMMENTED OUT
if (nrow(bc.changed) > 0) {
  push.df <- bc.changed %>%
    dplyr::filter(super_class != "sensory",
                  super_class != "motor" & region != "central_brain",
                  super_class != "visceral_circulatory" & region != "central_brain",
                  super_class != "glia",
                  super_class != "trachea") %>%
    dplyr::transmute(`_id`, root_id, sexually_dimorphic = sexually_dimorphic_new) %>%
    as.data.frame()
  message(sprintf("  Pushing sexually_dimorphic for %d neurons to banc_meta", nrow(push.df)))
  banctable_update_rows(base = 'banc_meta',
                        table = "banc_meta",
                        df = push.df,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

###########################
### 2. MANC franken_meta ##
###########################

message("\n=== 2. MANC franken_meta ===")

if (is.null(fm)) {
  message("  Skipping MANC section (SeaTable unavailable)")
  fm.merged <- NULL
} else {
message(sprintf("  %d MANC neurons in franken_meta", nrow(fm)))

# Join unified type-level dimorphism by cell_type
fm.merged <- fm %>%
  dplyr::left_join(unified_type_dim %>% dplyr::rename(unified_label = dimorphism),
                   by = "cell_type")

# Priority: unified type > terminalia sensory > existing metadata > isomorphic
fm.merged <- fm.merged %>%
  dplyr::mutate(
    is_terminalia_sensory = super_class == "sensory" &
      (grepl("terminalia", cell_class, ignore.case = TRUE) |
       grepl("terminalia", cell_sub_class, ignore.case = TRUE)),
    .raw_label = dplyr::case_when(
      !is.na(unified_label) ~ unified_label,
      is_terminalia_sensory ~ "dimorphic",
      sexually_dimorphic %in% .valid_labels ~ sexually_dimorphic,
      !is.na(cell_type) & cell_type != "" ~ "isomorphic",
      TRUE ~ NA_character_
    ),
    # MANC is male: a type present here can't be female-specific → dimorphic
    sexually_dimorphic_new = dplyr::if_else(
      .raw_label == "female-specific", "dimorphic", .raw_label)
  ) %>%
  dplyr::select(-dplyr::any_of(c("unified_label", "is_terminalia_sensory", ".raw_label")))

fm.changed <- fm.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new),
                is.na(sexually_dimorphic) | sexually_dimorphic == "" |
                  sexually_dimorphic != sexually_dimorphic_new)

n_from_unified <- sum(!is.na(fm.merged$sexually_dimorphic_new) &
                       fm.merged$sexually_dimorphic_new != "isomorphic", na.rm = TRUE)
message(sprintf("  %d MANC neurons with non-isomorphic labels from unified type transfer; %d need updating",
                n_from_unified, nrow(fm.changed)))

# Save MANC proposed changes
if (nrow(fm.changed) > 0) {
  manc_review <- file.path(plot.dir, "manc_dimorphism_proposed_changes.csv")
  readr::write_csv(
    fm.changed %>% dplyr::select(`_id`, manc_id, sexually_dimorphic, sexually_dimorphic_new),
    manc_review
  )
  message(sprintf("  Saved: %s (%d rows)", manc_review, nrow(fm.changed)))
}

### Push to franken_meta — COMMENTED OUT
if (nrow(fm.changed) > 0) {
  push.fm <- fm.changed %>%
    dplyr::transmute(`_id`, manc_id, sexually_dimorphic = sexually_dimorphic_new) %>%
    as.data.frame()
  message(sprintf("  Pushing sexually_dimorphic for %d MANC neurons to franken_meta", nrow(push.fm)))
  # banctable_update_rows(base = 'cns_meta',
  #                       table = "franken_meta",
  #                       df = push.fm,
  #                       append_allowed = FALSE,
  #                       chunksize = 1000)
 }
} # end if fm not null

###########################
### 3. maleCNS         ###
###########################

message("\n=== 3. maleCNS ===")

if (is.null(mc)) {
  message("  Skipping maleCNS section (SeaTable unavailable)")
  mc.merged <- NULL
} else {
message(sprintf("  %d typed maleCNS neurons in malecns table", nrow(mc)))

# Join unified type-level dimorphism by cell_type
mc.merged <- mc %>%
  dplyr::left_join(unified_type_dim %>% dplyr::rename(dim_ct = dimorphism),
                   by = "cell_type")

# Also try matching via manc_cell_type bridge
if ("manc_cell_type" %in% colnames(mc.merged)) {
  mc.merged <- mc.merged %>%
    dplyr::left_join(unified_type_dim %>% dplyr::rename(dim_manc = dimorphism),
                     by = c("manc_cell_type" = "cell_type"))
} else {
  mc.merged$dim_manc <- NA_character_
}

# Priority: unified (cell_type) > unified (manc_cell_type) > terminalia sensory > existing > isomorphic
mc.merged <- mc.merged %>%
  dplyr::mutate(
    is_terminalia_sensory = super_class == "sensory" &
      (grepl("terminalia", cell_class, ignore.case = TRUE) |
       grepl("terminalia", cell_sub_class, ignore.case = TRUE)),
    .raw_label = dplyr::case_when(
      !is.na(dim_ct) ~ dim_ct,
      !is.na(dim_manc) ~ dim_manc,
      is_terminalia_sensory ~ "dimorphic",
      sexually_dimorphic %in% .valid_labels ~ sexually_dimorphic,
      !is.na(cell_type) & cell_type != "" ~ "isomorphic",
      TRUE ~ NA_character_
    ),
    # maleCNS is male: a type present here can't be female-specific → dimorphic
    sexually_dimorphic_new = dplyr::if_else(
      .raw_label == "female-specific", "dimorphic", .raw_label)
  ) %>%
  dplyr::select(-dplyr::any_of(c("dim_ct", "dim_manc", "is_terminalia_sensory", ".raw_label")))

mc.changed <- mc.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new),
                is.na(sexually_dimorphic) | sexually_dimorphic == "" |
                  sexually_dimorphic != sexually_dimorphic_new)

n_from_unified <- sum(!is.na(mc.merged$sexually_dimorphic_new) &
                       mc.merged$sexually_dimorphic_new != "isomorphic", na.rm = TRUE)
message(sprintf("  %d maleCNS neurons with non-isomorphic labels; %d need updating",
                n_from_unified, nrow(mc.changed)))

# Save maleCNS proposed changes
if (nrow(mc.changed) > 0) {
  malecns_review <- file.path(plot.dir, "malecns_dimorphism_proposed_changes.csv")
  readr::write_csv(
    mc.changed %>% dplyr::select(`_id`, cell_type, sexually_dimorphic, sexually_dimorphic_new),
    malecns_review
  )
  message(sprintf("  Saved: %s (%d rows)", malecns_review, nrow(mc.changed)))
}

### Push to malecns — COMMENTED OUT
if (nrow(mc.changed) > 0) {
  push.mc <- mc.changed %>%
    dplyr::filter(super_class != "sensory",#|cell_type %in% c("SNxx17","SNxx12","SNxx10","SNxx09","SNxx08","SNxx07","SNxx04","SNta03","SNpp52","SNch01","LgLG8","LgLG5"),
                  super_class != "glia",
                  super_class != "trachea") %>%
    dplyr::transmute(`_id`, cell_type, sexually_dimorphic = sexually_dimorphic_new) %>%
    as.data.frame()
  message(sprintf("  Pushing sexually_dimorphic for %d maleCNS neurons to malecns table", nrow(push.mc)))
  banctable_update_rows(base = 'cns_meta',
                        table = "malecns",
                        df = push.mc,
                        append_allowed = FALSE,
                        chunksize = 1000)
}
} # end if mc not null

###########################
### Save dimorphism CSVs ###
###########################

message("\n=== Saving dimorphism assignment CSVs ===")

# Save BANC dimorphism mapping (for density script and eventual SeaTable push)
banc_dimorphism_out <- bc.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
  dplyr::select(dplyr::any_of(c("_id", "root_id", "root_626", "cell_type",
                                  "manc_cell_type", "super_class",
                                  "sexually_dimorphic", "sexually_dimorphic_new")))
banc_dim_file <- file.path(plot.dir, "banc_dimorphism_assignments.csv")
readr::write_csv(banc_dimorphism_out, banc_dim_file)
message(sprintf("  Saved: %s (%d rows)", banc_dim_file, nrow(banc_dimorphism_out)))

# Save MANC dimorphism mapping
if (!is.null(fm.merged)) {
  manc_dimorphism_out <- fm.merged %>%
    dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
    dplyr::select(dplyr::any_of(c("_id", "manc_id", "cell_type", "super_class",
                                    "sexually_dimorphic", "sexually_dimorphic_new")))
  manc_dim_file <- file.path(plot.dir, "manc_dimorphism_assignments.csv")
  readr::write_csv(manc_dimorphism_out, manc_dim_file)
  message(sprintf("  Saved: %s (%d rows)", manc_dim_file, nrow(manc_dimorphism_out)))
}

# Save maleCNS dimorphism mapping
if (!is.null(mc.merged)) {
  mcns_dimorphism_out <- mc.merged %>%
    dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
    dplyr::select(dplyr::any_of(c("_id", "cell_type", "manc_cell_type", "super_class",
                                    "sexually_dimorphic", "sexually_dimorphic_new")))
  mcns_dim_file <- file.path(plot.dir, "malecns_dimorphism_assignments.csv")
  readr::write_csv(mcns_dimorphism_out, mcns_dim_file)
  message(sprintf("  Saved: %s (%d rows)", mcns_dim_file, nrow(mcns_dimorphism_out)))
}

###########################
### Breakdown Tables   ###
###########################

message("\n=== Dimorphism Breakdown Tables ===")

all_labels <- c("isomorphic", "dimorphic", "female-specific", "male-specific")

### BANC ---
banc_bd <- bc.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
  dplyr::mutate(label = sexually_dimorphic_new,
                super_class = ifelse(is.na(super_class) | super_class == "", "unknown", super_class))

banc_by_neuron <- banc_bd %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

banc_by_type <- banc_bd %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

# Add totals row
banc_by_neuron <- dplyr::bind_rows(
  banc_by_neuron,
  data.frame(super_class = "TOTAL",
             as.list(colSums(banc_by_neuron[, -1, drop = FALSE])))
)
banc_by_type <- dplyr::bind_rows(
  banc_by_type,
  data.frame(super_class = "TOTAL",
             as.list(colSums(banc_by_type[, -1, drop = FALSE])))
)

message("\n--- BANC: Neuron count by super_class x dimorphism ---")
print(knitr::kable(banc_by_neuron))
message("\n--- BANC: Cell type count by super_class x dimorphism ---")
print(knitr::kable(banc_by_type))

### MANC ---
if (!is.null(fm.merged)) {
manc_bd <- fm.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
  dplyr::mutate(label = sexually_dimorphic_new,
                super_class = ifelse(is.na(super_class) | super_class == "", "unknown", super_class))

manc_by_neuron <- manc_bd %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

manc_by_type <- manc_bd %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

manc_by_neuron <- dplyr::bind_rows(
  manc_by_neuron,
  data.frame(super_class = "TOTAL",
             as.list(colSums(manc_by_neuron[, -1, drop = FALSE])))
)
manc_by_type <- dplyr::bind_rows(
  manc_by_type,
  data.frame(super_class = "TOTAL",
             as.list(colSums(manc_by_type[, -1, drop = FALSE])))
)

message("\n--- MANC: Neuron count by super_class x dimorphism ---")
print(knitr::kable(manc_by_neuron))
message("\n--- MANC: Cell type count by super_class x dimorphism ---")
print(knitr::kable(manc_by_type))
} else {
  message("\n  MANC: skipped (data unavailable)")
}

### maleCNS ---
if (!is.null(mc.merged)) {
mcns_bd <- mc.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
  dplyr::mutate(label = sexually_dimorphic_new,
                super_class = ifelse(is.na(super_class) | super_class == "", "unknown", super_class))

mcns_by_neuron <- mcns_bd %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

mcns_by_type <- mcns_bd %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::count(super_class, label, name = "n") %>%
  tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  dplyr::arrange(super_class)

mcns_by_neuron <- dplyr::bind_rows(
  mcns_by_neuron,
  data.frame(super_class = "TOTAL",
             as.list(colSums(mcns_by_neuron[, -1, drop = FALSE])))
)
mcns_by_type <- dplyr::bind_rows(
  mcns_by_type,
  data.frame(super_class = "TOTAL",
             as.list(colSums(mcns_by_type[, -1, drop = FALSE])))
)

message("\n--- maleCNS: Neuron count by super_class x dimorphism ---")
print(knitr::kable(mcns_by_neuron))
message("\n--- maleCNS: Cell type count by super_class x dimorphism ---")
print(knitr::kable(mcns_by_type))
} else {
  message("\n  maleCNS: skipped (data unavailable)")
}

###########################
### Plots              ###
###########################

message("\n=== Plots ===")

###########################
### Plot 1: Heatmap     ###
###########################

# Confusion-matrix: current seatable label (y) vs proposed label (x) for BANC
label_levels <- c("unlabelled", "isomorphic", "dimorphic", "female-specific", "male-specific")

bc.sourced <- bc.merged %>%
  dplyr::filter(!is.na(sexually_dimorphic_new)) %>%
  dplyr::mutate(
    old_label = dplyr::if_else(is.na(sexually_dimorphic) | sexually_dimorphic == "",
                                "unlabelled", sexually_dimorphic)
  )

heatmap_data <- bc.sourced %>%
  dplyr::mutate(
    old_label = factor(old_label, levels = label_levels),
    sexually_dimorphic_new = factor(sexually_dimorphic_new, levels = label_levels)
  ) %>%
  dplyr::count(old_label, sexually_dimorphic_new, .drop = FALSE) %>%
  dplyr::filter(n > 0)

plot_dimorphism_heatmap <- ggplot(heatmap_data,
             aes(x = sexually_dimorphic_new, y = old_label, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 3.5,
            colour = ifelse(heatmap_data$n > max(heatmap_data$n) * 0.5, "white", "black")) +
  scale_fill_gradient(low = "grey90", high = paper.cols["dimorphic"], name = "Neurons",
                      trans = "log1p") +
  labs(title = "BANC: current vs proposed sexually_dimorphic labels",
       x = "Proposed label (from unified type-level transfer)", y = "Current seatable label") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10))

ggsave(file.path(plot.dir, "sexually_dimorphic_current_vs_proposed.pdf"),
       plot = plot_dimorphism_heatmap, width = 7, height = 5, dpi = 300, bg = "white")
message("  Saved: sexually_dimorphic_current_vs_proposed.pdf")

##############################################
### Plot 2: % dimorphic by super_class     ###
##############################################

# For typed BANC neurons, show proportion that are dimorphic/female-specific
# by super_class. Uses CSV label where available, else existing seatable label.
non_neuronal <- c("glia", "trachea", "")

bc.final <- bc.merged %>%
  dplyr::filter(!is.na(cell_type), cell_type != "",
                !is.na(super_class), !super_class %in% non_neuronal) %>%
  dplyr::mutate(
    final_label = dplyr::case_when(
      sexually_dimorphic %in% protected ~ sexually_dimorphic,
      !is.na(sexually_dimorphic_new) ~ sexually_dimorphic_new,
      !is.na(sexually_dimorphic) & sexually_dimorphic != "" ~ sexually_dimorphic,
      TRUE ~ "isomorphic"
    )
  )

dimorphic_by_sc <- bc.final %>%
  dplyr::mutate(is_dimorphic = final_label %in% c("dimorphic", "female-specific")) %>%
  dplyr::group_by(super_class) %>%
  dplyr::summarise(
    n_dimorphic = sum(is_dimorphic),
    n_total = dplyr::n(),
    pct_dimorphic = n_dimorphic / n_total * 100,
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(pct_dimorphic))

# Assign region: brain, VNC, or projection (ascending/descending/sensory span both)
assign_region <- function(sc) {
  dplyr::case_when(
    sc %in% c("central_brain_intrinsic", "optic_lobe_intrinsic",
              "visual_projection", "visual_centrifugal") ~ "brain",
    sc %in% c("ventral_nerve_cord_intrinsic", "motor", "visceral_circulatory",
              "ascending_visceral_circulatory", "sensory_descending",
              "sensory_ascending", "unknown_visceral_circulatory") ~ "VNC",
    sc %in% c("ascending", "descending", "sensory") ~ "projection",
    TRUE ~ "other"
  )
}

label_by_sc <- bc.final %>%
  dplyr::mutate(region = assign_region(super_class)) %>%
  dplyr::count(super_class, region, final_label) %>%
  dplyr::group_by(super_class) %>%
  dplyr::mutate(prop = n / sum(n), total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(final_label %in% c("dimorphic", "female-specific"))

sc_totals <- label_by_sc %>%
  dplyr::group_by(super_class, region) %>%
  dplyr::summarise(prop_sum = sum(prop), bar_total = sum(n), .groups = "drop")

sc_order <- dimorphic_by_sc$super_class
label_by_sc$super_class <- factor(label_by_sc$super_class, levels = sc_order)
sc_totals$super_class <- factor(sc_totals$super_class, levels = sc_order)
label_by_sc$region <- factor(label_by_sc$region, levels = c("brain", "projection", "VNC"))
sc_totals$region <- factor(sc_totals$region, levels = c("brain", "projection", "VNC"))

plot_dimorphism_by_superclass <- ggplot(label_by_sc, aes(x = super_class, y = prop, fill = final_label)) +
  geom_col(color = NA) +
  geom_text(data = sc_totals,
            aes(x = super_class, y = prop_sum + 0.02,
                label = paste0("n=", bar_total), fill = NULL),
            size = 2.8, hjust = 0.5) +
  scale_fill_manual(values = paper.cols) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ region, scales = "free_x") +
  labs(title = "Sexually dimorphic neurons by super_class",
       subtitle = "% dimorphic + female-specific (typed BANC neurons)",
       x = NULL, y = "Proportion", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1))

ggsave(file.path(plot.dir, "sexually_dimorphic_by_superclass.pdf"),
       plot = plot_dimorphism_by_superclass, width = 14, height = 5, dpi = 300, bg = "white")
message("  Saved: sexually_dimorphic_by_superclass.pdf")

##############################################
### Plot 2b-d: % dimorphic by super_class  ###
### for MANC, maleCNS, FAFB               ###
##############################################

# Helper: build the dimorphic-by-superclass plot for any dataset
make_dimorphic_by_sc_plot <- function(df, label_col, dataset_name, sex) {
  non_neuronal <- c("glia", "trachea", "")
  df_final <- df %>%
    dplyr::filter(!is.na(cell_type), cell_type != "",
                  !is.na(super_class), !super_class %in% non_neuronal) %>%
    dplyr::mutate(final_label = !!rlang::sym(label_col))
  if (nrow(df_final) == 0) return(NULL)

  # For male datasets, dimorphic labels are "dimorphic" and "male-specific"
  # For female datasets, "dimorphic" and "female-specific"
  if (sex == "male") {
    dim_labels <- c("dimorphic", "male-specific")
  } else {
    dim_labels <- c("dimorphic", "female-specific")
  }

  df_final <- df_final %>%
    dplyr::mutate(region = assign_region(super_class))

  lbl <- df_final %>%
    dplyr::count(super_class, region, final_label) %>%
    dplyr::group_by(super_class) %>%
    dplyr::mutate(prop = n / sum(n), total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(final_label %in% dim_labels)
  if (nrow(lbl) == 0) return(NULL)

  tots <- lbl %>%
    dplyr::group_by(super_class, region) %>%
    dplyr::summarise(prop_sum = sum(prop), bar_total = sum(n), .groups = "drop")

  # Order by % dimorphic
  sc_pct <- df_final %>%
    dplyr::mutate(is_dim = final_label %in% dim_labels) %>%
    dplyr::group_by(super_class) %>%
    dplyr::summarise(pct = sum(is_dim) / dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(pct))

  lbl$super_class <- factor(lbl$super_class, levels = sc_pct$super_class)
  tots$super_class <- factor(tots$super_class, levels = sc_pct$super_class)
  lbl$region <- factor(lbl$region, levels = c("brain", "projection", "VNC"))
  tots$region <- factor(tots$region, levels = c("brain", "projection", "VNC"))

  ggplot(lbl, aes(x = super_class, y = prop, fill = final_label)) +
    geom_col(color = NA) +
    geom_text(data = tots,
              aes(x = super_class, y = prop_sum + 0.02,
                  label = paste0("n=", bar_total), fill = NULL),
              size = 2.8, hjust = 0.5) +
    scale_fill_manual(values = paper.cols) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ region, scales = "free_x") +
    labs(title = sprintf("Sexually dimorphic neurons by super_class (%s)", dataset_name),
         subtitle = sprintf("%% dimorphic + %s (typed neurons)",
                            if (sex == "male") "male-specific" else "female-specific"),
         x = NULL, y = "Proportion", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1))
}

# MANC
if (!is.null(fm.merged)) {
  plot_dimorphism_by_superclass_manc <- make_dimorphic_by_sc_plot(fm.merged, "sexually_dimorphic_new", "MANC", "male")
  if (!is.null(plot_dimorphism_by_superclass_manc)) {
    ggsave(file.path(plot.dir, "sexually_dimorphic_by_superclass_manc.pdf"),
           plot = plot_dimorphism_by_superclass_manc, width = 14, height = 5, dpi = 300, bg = "white")
    message("  Saved: sexually_dimorphic_by_superclass_manc.pdf")
  }
}

# maleCNS
if (!is.null(mc.merged)) {
  plot_dimorphism_by_superclass_malecns <- make_dimorphic_by_sc_plot(mc.merged, "sexually_dimorphic_new", "maleCNS", "male")
  if (!is.null(plot_dimorphism_by_superclass_malecns)) {
    ggsave(file.path(plot.dir, "sexually_dimorphic_by_superclass_malecns.pdf"),
           plot = plot_dimorphism_by_superclass_malecns, width = 14, height = 5, dpi = 300, bg = "white")
    message("  Saved: sexually_dimorphic_by_superclass_malecns.pdf")
  }
}

# FAFB — use already-queried fafb_df
fafb.merged <- if (!is.null(fafb_df) && nrow(fafb_df) > 0) {
  fafb_neurons <- fafb_df %>%
    dplyr::left_join(unified_type_dim %>% dplyr::rename(unified_dim = dimorphism),
                     by = "cell_type") %>%
    dplyr::mutate(
      .raw_label = dplyr::case_when(
        !is.na(unified_dim) ~ unified_dim,
        sexually_dimorphic %in% .valid_labels ~ sexually_dimorphic,
        !is.na(cell_type) & cell_type != "" ~ "isomorphic",
        TRUE ~ NA_character_
      ),
      # FAFB is female: a type present here can't be male-specific → dimorphic
      sexually_dimorphic_new = dplyr::if_else(
        .raw_label == "male-specific", "dimorphic", .raw_label)
    ) %>%
    dplyr::select(-unified_dim, -.raw_label)
  fix_super_class(fafb_neurons, "cell_type")
} else {
  NULL
}

if (!is.null(fafb.merged)) {
  plot_dimorphism_by_superclass_fafb <- make_dimorphic_by_sc_plot(fafb.merged, "sexually_dimorphic_new", "FAFB", "female")
  if (!is.null(plot_dimorphism_by_superclass_fafb)) {
    ggsave(file.path(plot.dir, "sexually_dimorphic_by_superclass_fafb.pdf"),
           plot = plot_dimorphism_by_superclass_fafb, width = 14, height = 5, dpi = 300, bg = "white")
    message("  Saved: sexually_dimorphic_by_superclass_fafb.pdf")
  }
}

###################################
### Plot 3: Final label summary ###
###################################

label_summary <- bc.final %>%
  dplyr::count(final_label) %>%
  dplyr::mutate(final_label = factor(final_label,
                                      levels = c("isomorphic", "dimorphic",
                                                 "female-specific", "unlabelled")))

plot_dimorphism_label_distribution <- ggplot(label_summary, aes(x = final_label, y = n, fill = final_label)) +
  geom_col(color = NA, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = paper.cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Final sexually_dimorphic label distribution",
       subtitle = "Typed BANC neurons",
       x = NULL, y = "Neurons") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(plot.dir, "sexually_dimorphic_label_distribution.pdf"),
       plot = plot_dimorphism_label_distribution, width = 7, height = 5, dpi = 300, bg = "white")
message("  Saved: sexually_dimorphic_label_distribution.pdf")

message(sprintf("\n### banc: sexually_dimorphic update complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
