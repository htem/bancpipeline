#' banc-cx-fix — Fix fafb_/hemibrain_ matches for central-complex neurons.
#'
#' CX neurons often have good `cell_type` but poor `fafb_cell_type` /
#' `fafb_match` and stale `hemibrain_cell_type` / `hemibrain_match`. This
#' script picks the best matching FAFB id and hemibrain bodyid per
#' cell_type using NBLAST scores from GCS.
#'
#' @section Reads:
#'   - SeaTable `banc_meta` (CX neurons via `cell_class` / `cns_network`)
#'   - `banc_fafb_783_nblast.feather`, `banc_hemibrain_v1.2.1_nblast.feather` (GCS-cached)
#'   - `franken_meta()`
#'
#' @section Writes:
#'   - `data/codex/cx_fafb_fixes.csv` — proposed fafb/hemibrain match swaps

###########################################################
### Fix fafb_cell_type, fafb_match, hemibrain_cell_type,
### and hemibrain_match for CX neurons
###
### Central complex neurons (cell_class or cns_network
### contains "central_complex") have good cell_type but
### often poor fafb_cell_type/fafb_match and sometimes
### stale hemibrain_cell_type/hemibrain_match.
###
### Uses NBLAST scores from GCS to pick the best matching
### fafb_id and hemibrain bodyid for each cell_type.
###
### Output: data/codex/cx_fafb_fixes.csv
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: CX fafb + hemibrain fix ###")

###########################
### Read data           ###
###########################

bc <- banctable_query(
  "SELECT _id, root_888, root_id, cell_type, cell_class, cns_network, cell_type_source, fafb_cell_type, fafb_match, hemibrain_cell_type, hemibrain_match, side FROM banc_meta"
) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

# CX neurons: cell_class OR cns_network contains "central_complex"
cx <- bc %>%
  dplyr::filter(
    grepl("central_complex", cell_class, ignore.case = TRUE) |
    grepl("central_complex", cns_network, ignore.case = TRUE)
  ) %>%
  dplyr::filter(!is.na(cell_type), cell_type != "")
message(sprintf("  CX neurons with cell_type: %d", nrow(cx)))

###########################
### Load reference data ###
###########################

# franken_meta for FAFB lookups
message("  Loading franken_meta...")
fm <- franken_meta()
fafb_pool <- fm %>%
  dplyr::filter(!is.na(hemibrain_type), hemibrain_type != "",
                !is.na(cell_type), cell_type != "",
                !is.na(fafb_id), fafb_id != "") %>%
  dplyr::select(fafb_id, fm_hemibrain_type = hemibrain_type,
                fm_cell_type = cell_type, fm_side = side) %>%
  dplyr::mutate(fafb_id = as.character(fafb_id))

fm_cell_types <- unique(fafb_pool$fm_cell_type)
fm_hb_types <- unique(fafb_pool$fm_hemibrain_type)

# Hemibrain meta
message("  Loading hemibrain meta...")
hb_meta_path <- "/Users/papers/BANC-project/data/meta/hemibrain_meta.csv"
if (!file.exists(hb_meta_path)) {
  hb_meta_path <- file.path(banc.meta.save.path, "hemibrain_meta.csv")
}
hb_meta <- readr::read_csv(hb_meta_path, col_types = readr::cols(.default = "c"),
                             show_col_types = FALSE) %>%
  dplyr::filter(!is.na(bodyid), bodyid != "",
                !is.na(cell_type), cell_type != "")
message(sprintf("  Hemibrain meta: %d neurons with cell_type", nrow(hb_meta)))

# NBLAST scores from GCS (cached)
message("  Loading NBLAST feathers...")
nblast_cache <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache, showWarnings = FALSE)

load_nblast_feather <- function(gcs_path) {
  local_file <- file.path(nblast_cache, basename(gcs_path))
  if (!file.exists(local_file)) {
    # Also check /tmp/nblast_cache from previous runs
    alt_file <- file.path("/tmp/nblast_cache", basename(gcs_path))
    if (file.exists(alt_file)) {
      local_file <- alt_file
    } else {
      message(sprintf("    Downloading %s...", basename(gcs_path)))
      system2("gsutil", c("cp", gcs_path, local_file), stdout = FALSE, stderr = FALSE)
    }
  } else {
    message(sprintf("    Using cached %s", basename(local_file)))
  }
  arrow::read_feather(local_file)
}

fafb_nblast <- load_nblast_feather(
  "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_fafb_783_nblast.feather"
)
hb_nblast <- load_nblast_feather(
  "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_hemibrain_v1.2.1_nblast.feather"
)
message(sprintf("  FAFB NBLAST: %d rows, Hemibrain NBLAST: %d rows",
                nrow(fafb_nblast), nrow(hb_nblast)))

###########################
### Helper functions    ###
###########################

# Find fafb_cell_type via hemibrain_type lookup (1:1 only)
find_via_hemibrain <- function(hb_type, pool) {
  matches <- pool %>%
    dplyr::filter(fm_hemibrain_type == hb_type) %>%
    dplyr::distinct(fm_cell_type)
  if (nrow(matches) == 1) return(matches$fm_cell_type[1])
  NA_character_
}

# Resolve the fafb_cell_type for a given BANC cell_type + hemibrain_cell_type
resolve_fafb_cell_type <- function(ct, hb) {
  # Strategy 1: cell_type exists directly in franken_meta
  if (ct %in% fm_cell_types) return(list(type = ct, strategy = "direct_cell_type"))

  # Strategy 2: hemibrain_cell_type matches franken_meta hemibrain_type
  if (!is.na(hb) && hb != "" && hb %in% fm_hb_types) {
    fm_ct <- find_via_hemibrain(hb, fafb_pool)
    if (!is.na(fm_ct)) return(list(type = fm_ct, strategy = "hemibrain_type_direct"))
  }

  # Strategy 3: strip suffix (e.g. ER2_c -> ER2)
  if (!is.na(hb) && hb != "") {
    hb_base <- sub("_[a-z]$", "", hb)
    if (hb_base != hb && hb_base %in% fm_hb_types) {
      fm_ct <- find_via_hemibrain(hb_base, fafb_pool)
      if (!is.na(fm_ct)) return(list(type = fm_ct, strategy = "hemibrain_base_type"))
    }
  }

  list(type = NA_character_, strategy = "no_match")
}

# Pre-index NBLAST data for fast lookup
message("  Pre-indexing NBLAST data...")
fafb_nblast_idx <- split(
  data.frame(match_id = as.character(fafb_nblast$match_id),
             score = fafb_nblast$score, stringsAsFactors = FALSE),
  fafb_nblast$root_888
)
hb_nblast_idx <- split(
  data.frame(match_id = as.character(hb_nblast$match_id),
             score = hb_nblast$score, stringsAsFactors = FALSE),
  hb_nblast$root_888
)

# Pre-index type -> IDs
fafb_type_ids <- split(fafb_pool$fafb_id, fafb_pool$fm_cell_type)
hb_type_ids <- split(hb_meta$bodyid, hb_meta$cell_type)

message(sprintf("  FAFB index: %d neurons, HB index: %d neurons",
                length(fafb_nblast_idx), length(hb_nblast_idx)))

# Find best FAFB match by NBLAST score
find_best_fafb_by_nblast <- function(root_888_id, target_type, side) {
  type_ids <- fafb_type_ids[[target_type]]
  if (is.null(type_ids) || length(type_ids) == 0) return(NA_character_)

  hits <- fafb_nblast_idx[[root_888_id]]
  if (!is.null(hits) && nrow(hits) > 0) {
    type_mask <- hits$match_id %in% type_ids
    if (any(type_mask)) {
      type_hits <- hits[type_mask, ]
      return(type_hits$match_id[which.max(type_hits$score)])
    }
  }

  # Fallback: any FAFB neuron of that type, prefer same side
  candidates <- fafb_pool[fafb_pool$fm_cell_type == target_type, ]
  if (!is.na(side) && side != "") {
    same_side <- candidates[candidates$fm_side == side, ]
    if (nrow(same_side) > 0) candidates <- same_side
  }
  as.character(candidates$fafb_id[1])
}

# Find best hemibrain match by NBLAST score
find_best_hb_by_nblast <- function(root_888_id, target_type) {
  type_ids <- hb_type_ids[[target_type]]
  if (is.null(type_ids) || length(type_ids) == 0)
    return(list(id = NA_character_, type = NA_character_))

  hits <- hb_nblast_idx[[root_888_id]]
  if (!is.null(hits) && nrow(hits) > 0) {
    type_mask <- hits$match_id %in% type_ids
    if (any(type_mask)) {
      type_hits <- hits[type_mask, ]
      return(list(id = type_hits$match_id[which.max(type_hits$score)],
                  type = target_type))
    }
  }

  list(id = as.character(type_ids[1]), type = target_type)
}

###########################
### FAFB fixes          ###
###########################

message("\n  === FAFB fixes ===")
needs_fafb <- cx %>%
  dplyr::filter(is.na(fafb_cell_type) | fafb_cell_type == "" | fafb_cell_type != cell_type)
message(sprintf("  Needing fafb_cell_type fix: %d", nrow(needs_fafb)))

fafb_results <- needs_fafb %>%
  dplyr::select(root_888, `_id`, root_id, cell_type, hemibrain_cell_type,
                old_fafb_cell_type = fafb_cell_type, old_fafb_match = fafb_match,
                cell_type_source, side) %>%
  dplyr::mutate(new_fafb_cell_type = NA_character_,
                new_fafb_match = NA_character_,
                fafb_strategy = NA_character_)

for (i in seq_len(nrow(fafb_results))) {
  ct <- fafb_results$cell_type[i]
  hb <- fafb_results$hemibrain_cell_type[i]
  rid <- fafb_results$root_888[i]
  sid <- fafb_results$side[i]

  res <- resolve_fafb_cell_type(ct, hb)
  fafb_results$new_fafb_cell_type[i] <- res$type
  fafb_results$fafb_strategy[i] <- res$strategy

  if (!is.na(res$type)) {
    fafb_results$new_fafb_match[i] <- find_best_fafb_by_nblast(rid, res$type, sid)
  }
}

tab <- table(fafb_results$fafb_strategy)
for (s in names(tab)) message(sprintf("    %s: %d", s, tab[s]))

###########################
### Hemibrain fixes     ###
###########################

message("\n  === Hemibrain fixes ===")
needs_hb <- cx %>%
  dplyr::filter(!is.na(hemibrain_cell_type), hemibrain_cell_type != "",
                cell_type != hemibrain_cell_type)
message(sprintf("  cell_type != hemibrain_cell_type: %d", nrow(needs_hb)))

# Resolve the correct hemibrain_cell_type: use cell_type if it exists in hb_meta,
# else try stripping suffix, else use cell_type directly as hemibrain type
hb_types <- unique(hb_meta$cell_type)

hb_results <- needs_hb %>%
  dplyr::select(root_888, `_id`, root_id, cell_type, cell_type_source, side,
                old_hemibrain_cell_type = hemibrain_cell_type,
                old_hemibrain_match = hemibrain_match) %>%
  dplyr::mutate(new_hemibrain_cell_type = NA_character_,
                new_hemibrain_match = NA_character_,
                hb_strategy = NA_character_)

for (i in seq_len(nrow(hb_results))) {
  ct <- hb_results$cell_type[i]
  rid <- hb_results$root_888[i]

  # Strategy 1: cell_type exists directly in hemibrain meta
  if (ct %in% hb_types) {
    best <- find_best_hb_by_nblast(rid, ct)
    hb_results$new_hemibrain_cell_type[i] <- best$type
    hb_results$new_hemibrain_match[i] <- best$id
    hb_results$hb_strategy[i] <- "direct"
    next
  }

  # Strategy 2: strip suffix (e.g. ER2_c -> ER2)
  ct_base <- sub("_[a-z]$", "", ct)
  if (ct_base != ct && ct_base %in% hb_types) {
    best <- find_best_hb_by_nblast(rid, ct_base)
    hb_results$new_hemibrain_cell_type[i] <- best$type
    hb_results$new_hemibrain_match[i] <- best$id
    hb_results$hb_strategy[i] <- "base_type"
    next
  }

  hb_results$hb_strategy[i] <- "no_match"
}

tab <- table(hb_results$hb_strategy)
for (s in names(tab)) message(sprintf("    %s: %d", s, tab[s]))

###########################
### Merge and report    ###
###########################

# Combine FAFB and hemibrain fixes into a single results table
# Some neurons may appear in both
all_ids <- union(fafb_results$root_888, hb_results$root_888)

combined <- cx %>%
  dplyr::filter(root_888 %in% all_ids) %>%
  dplyr::select(root_888, `_id`, root_id, cell_type, cell_type_source, side,
                old_fafb_cell_type = fafb_cell_type, old_fafb_match = fafb_match,
                old_hemibrain_cell_type = hemibrain_cell_type, old_hemibrain_match = hemibrain_match)

# Join FAFB fixes
combined <- combined %>%
  dplyr::left_join(
    fafb_results %>% dplyr::select(root_888, new_fafb_cell_type, new_fafb_match, fafb_strategy),
    by = "root_888"
  )

# Join hemibrain fixes
combined <- combined %>%
  dplyr::left_join(
    hb_results %>% dplyr::select(root_888, new_hemibrain_cell_type, new_hemibrain_match, hb_strategy),
    by = "root_888"
  )

message(sprintf("\n  === Combined ==="))
message(sprintf("  Total neurons with any fix: %d", nrow(combined)))
message(sprintf("  With FAFB fix: %d", sum(!is.na(combined$new_fafb_cell_type))))
message(sprintf("  With hemibrain fix: %d", sum(!is.na(combined$new_hemibrain_cell_type))))
message(sprintf("  With both: %d",
                sum(!is.na(combined$new_fafb_cell_type) & !is.na(combined$new_hemibrain_cell_type))))

###########################
### Save CSV            ###
###########################

out_file <- "data/codex/cx_fafb_fixes.csv"
readr::write_csv(combined, out_file)
message(sprintf("\n  Saved %d rows to %s", nrow(combined), out_file))

###########################
### SeaTable update     ###
###########################

# Build push dataframe: only rows with actual changes
push_df <- combined %>%
  dplyr::filter(!is.na(new_fafb_cell_type) | !is.na(new_hemibrain_cell_type)) %>%
  dplyr::transmute(
    `_id`,
    fafb_cell_type = dplyr::if_else(!is.na(new_fafb_cell_type), new_fafb_cell_type, old_fafb_cell_type),
    fafb_match = dplyr::if_else(!is.na(new_fafb_match), new_fafb_match, old_fafb_match),
    hemibrain_cell_type = dplyr::if_else(!is.na(new_hemibrain_cell_type),
                                          new_hemibrain_cell_type, old_hemibrain_cell_type),
    hemibrain_match = dplyr::if_else(!is.na(new_hemibrain_match),
                                      new_hemibrain_match, old_hemibrain_match),
    cell_type_source = append_status(cell_type_source, "bates")
  ) %>%
  as.data.frame()

# Clean NAs for seatable
for (col in c("fafb_cell_type", "fafb_match", "hemibrain_cell_type", "hemibrain_match")) {
  push_df[[col]][is.na(push_df[[col]])] <- ""
}

message(sprintf("\n  SeaTable push: %d rows", nrow(push_df)))
message(sprintf("    fafb_cell_type: %d non-empty", sum(push_df$fafb_cell_type != "")))
message(sprintf("    fafb_match: %d non-empty", sum(push_df$fafb_match != "")))
message(sprintf("    hemibrain_cell_type: %d non-empty", sum(push_df$hemibrain_cell_type != "")))
message(sprintf("    hemibrain_match: %d non-empty", sum(push_df$hemibrain_match != "")))

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_df,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  SeaTable update complete")

})
