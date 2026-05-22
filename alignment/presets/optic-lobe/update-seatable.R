#' update-seatable (optic-lobe preset) — Push optic-lobe alignment predictions to SeaTable.
#'
#' Loaded by `alignment/banc-alignment-update-seatable.R --region optic-lobe`.
#' Uploads connectivity-alignment predictions joining on `root_888`. Writes
#' only to the `fafb_alignment_*` / `fafb_ntac_*` column family.
#'
#' @section Reads:
#'   - `[results_file]` positional (default `data/optic_lobe/banc_optic_both_alignment_optic_lobe_full.csv`)
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `fafb_alignment_cell_type`,
#'     `fafb_alignment_match`, `fafb_alignment_confidence`,
#'     `fafb_alignment_runner_up`, `fafb_ntac_cell_type`
#'
#' @section CLI (via dispatcher):
#'   --region optic-lobe [results_file]
#'   --dry-run

###########################################################
### Optic Lobe Alignment: SeaTable Update
###
### Uploads connectivity alignment predictions to BANC
### SeaTable, joining on root_888.
###
### Columns updated:
###   - fafb_alignment_cell_type: assigned cell type
###   - fafb_alignment_match: best FAFB neuron match (by NBLAST within assigned type)
###   - fafb_alignment_confidence: confidence (best_score - runner_up_score)
###   - fafb_alignment_runner_up: runner-up type
###   - fafb_ntac_cell_type: NTAC within-dataset prediction (if NTAC results available)
###
### Usage:
###   Rscript alignment/banc-alignment-update-seatable.R --region optic-lobe [results_file]
###   Rscript alignment/banc-alignment-update-seatable.R --region optic-lobe  # default: optic_lobe_full
###   Rscript alignment/banc-alignment-update-seatable.R --region optic-lobe --dry-run
###########################################################
source("banc/banc-startup.R")

local({

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
args <- args[args != "--dry-run"]

# Default to the bilateral full run
results_file <- if (length(args) > 0) {
  args[1]
} else {
  "data/optic_lobe/banc_optic_both_alignment_optic_lobe_full.csv"
}

if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}

message(sprintf("=== Preparing SeaTable update from %s ===", results_file))
if (dry_run) message("  *** DRY RUN — no changes will be pushed ***")

# Load alignment results
alignment <- readr::read_csv(results_file,
  col_types = readr::cols(root_888 = "c", best_target_match = "c"),
  show_col_types = FALSE)
message(sprintf("  Alignment results: %d neurons (%d with cell_type, %d with FAFB match)",
                nrow(alignment),
                sum(alignment$assigned_cell_type != "", na.rm = TRUE),
                sum(!is.na(alignment$best_target_match) & alignment$best_target_match != "")))

# Query SeaTable for _id and root_888 mapping
message("  Querying SeaTable for _id mapping...")
bc <- banctable_query(
  "SELECT _id, root_id, root_888, cell_type, fafb_cell_type FROM banc_meta"
) %>%
  dplyr::filter(!is.na(root_888), root_888 != "") %>%
  dplyr::mutate(root_888 = as.character(root_888)) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)
message(sprintf("  SeaTable: %d neurons with root_888", nrow(bc)))

# Build push dataframe — only for neurons WITHOUT existing GT cell_type
# (anchored neurons already have correct types; we only push predictions for untyped neurons)
has_gt <- bc %>%
  dplyr::filter(!is.na(cell_type) & cell_type != "" & !grepl("^auto:", cell_type)) %>%
  dplyr::pull(root_888)

push_df <- alignment %>%
  dplyr::filter(assigned_cell_type != "" | !is.na(best_target_match)) %>%
  dplyr::filter(!root_888 %in% has_gt) %>%
  dplyr::inner_join(bc %>% dplyr::select(root_888, `_id`), by = "root_888") %>%
  dplyr::transmute(
    `_id` = `_id`,
    fafb_alignment_cell_type = ifelse(is.na(assigned_cell_type), "", assigned_cell_type),
    fafb_alignment_match = ifelse(is.na(best_target_match), "", best_target_match),
    fafb_alignment_confidence = round(confidence, 4),
    fafb_alignment_runner_up = ifelse(is.na(runner_up_type), "", runner_up_type)
  ) %>%
  as.data.frame()

message(sprintf("  Skipping %d neurons with existing GT cell_type", length(has_gt)))

# Load NTAC results (if available alongside alignment results)
ntac_file <- sub("_alignment_", "_ntac_", sub("alignment", "ntac", results_file))
# Also try standard naming: same dir, banc_{prefix}_{side}_ntac_full.csv
if (!file.exists(ntac_file)) {
  ntac_file <- file.path(dirname(results_file),
    sub("alignment.*", "ntac_full.csv", basename(results_file)))
}
if (file.exists(ntac_file)) {
  message(sprintf("  Loading NTAC results: %s", ntac_file))
  ntac <- readr::read_csv(ntac_file,
    col_types = readr::cols(root_888 = "c"), show_col_types = FALSE)
  push_df <- push_df %>%
    dplyr::left_join(
      ntac %>%
        dplyr::transmute(root_888 = as.character(root_888),
                         fafb_ntac_cell_type = ifelse(is.na(ntac_cell_type) | ntac_cell_type == "?",
                                                      "", as.character(ntac_cell_type))) %>%
        dplyr::inner_join(bc %>% dplyr::select(root_888, `_id`), by = "root_888") %>%
        dplyr::select(`_id`, fafb_ntac_cell_type),
      by = "_id"
    )
  push_df$fafb_ntac_cell_type[is.na(push_df$fafb_ntac_cell_type)] <- ""
  n_ntac <- sum(push_df$fafb_ntac_cell_type != "")
  message(sprintf("  NTAC: %d non-empty predictions", n_ntac))
} else {
  message(sprintf("  NTAC results not found at %s, skipping", ntac_file))
}

# Summary
n_typed <- sum(push_df$fafb_alignment_cell_type != "")
n_matched <- sum(push_df$fafb_alignment_match != "")
n_unique_types <- length(unique(push_df$fafb_alignment_cell_type[push_df$fafb_alignment_cell_type != ""]))

message(sprintf("  Ready to push: %d rows", nrow(push_df)))
message(sprintf("    fafb_alignment_cell_type: %d non-empty (%d unique types)", n_typed, n_unique_types))
message(sprintf("    fafb_alignment_match: %d non-empty", n_matched))
conf_numeric <- as.numeric(push_df$fafb_alignment_confidence)
message(sprintf("    fafb_alignment_confidence: mean=%.4f, median=%.4f",
                mean(conf_numeric, na.rm = TRUE),
                median(conf_numeric, na.rm = TRUE)))

# Check for conflicts: neurons where assigned type differs from existing cell_type
conflict_check <- alignment %>%
  dplyr::filter(assigned_cell_type != "") %>%
  dplyr::inner_join(bc %>% dplyr::select(root_888, cell_type, fafb_cell_type), by = "root_888") %>%
  dplyr::filter(
    !is.na(cell_type) & cell_type != "" & !grepl("^auto:", cell_type) &
      cell_type != assigned_cell_type
  )
if (nrow(conflict_check) > 0) {
  message(sprintf("  WARNING: %d neurons where alignment differs from existing cell_type", nrow(conflict_check)))
  top_conflicts <- conflict_check %>%
    dplyr::count(cell_type, assigned_cell_type, sort = TRUE) %>%
    head(10)
  for (i in seq_len(nrow(top_conflicts))) {
    r <- top_conflicts[i, ]
    message(sprintf("    %s -> %s: %d", r$cell_type, r$assigned_cell_type, r$n))
  }
}

# Push to SeaTable
if (!dry_run) {
  message("  Pushing to SeaTable...")
  banctable_update_rows(base = "banc_meta",
                        table = "banc_meta",
                        df = push_df,
                        append_allowed = FALSE,
                        chunksize = 200)
  message("  SeaTable update complete!")
} else {
  message("  Dry run complete. Remove --dry-run to push.")
}

message("\n=== Done ===")

})
