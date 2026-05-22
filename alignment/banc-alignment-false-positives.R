#' banc-alignment-false-positives — Pull reviewed alignment rejections into the forbidden-matches CSV.
#'
#' Reads `fafb_alignment_decision == F` rows from SeaTable and writes a
#' two-column (`root_id`, `cell_type`) blacklist that the Python aligner
#' consumes via `--forbidden-matches`. Written to BOTH the optic-lobe and
#' whole-brain preset paths so both runs share one source of truth.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`: cols `_id`, `root_id`, `root_888`,
#'     `fafb_alignment_cell_type`, `fafb_alignment_decision`,
#'     `cell_type`, `fafb_cell_type`
#'
#' @section Writes:
#'   - `alignment/presets/optic-lobe/forbidden-matches.csv`
#'   - `alignment/presets/whole-brain/forbidden-matches.csv`
#'
#' @section Invoked by:
#'   `o2/alignment/o2_banc_alignment.sh`,
#'   `o2/alignment/o2_banc_optic_prep.sh`,
#'   `o2/oneshots/o2_banc_optic_sweep.sh`,
#'   `o2/oneshots/o2_banc_wb_probe.sh`
#'
#' @section Notes:
#'   - See `banc-alignment-false-positives-append.R` for the append-mode variant
#'     that preserves on-disk blacklist rows.

###########################################################
### Alignment: forbidden (root_id, cell_type) pairs
###
### Pulls reviewed alignment predictions from BANC SeaTable
### and writes a CSV listing every (root_id, cell_type) pair
### that has been manually rejected (fafb_alignment_decision
### == FALSE/F). The downstream Python aligner can be passed
### `--forbidden-matches PATH` to veto these pairs during
### iterative assignment.
###
### Source columns from banc_meta:
###   - _id (SeaTable row id, kept only for sanity checks)
###   - root_id, root_888 (BANC neuron identifiers)
###   - fafb_alignment_cell_type (the rejected type)
###   - fafb_alignment_decision (T/F/empty review label)
###   - cell_type, fafb_cell_type (used to flag conflicts
###     where a forbidden pair contradicts existing GT)
###
### Output: two-column CSV (`root_id`, `cell_type`), unique
### on the (root_id, cell_type) pair, written to BOTH
###   alignment/presets/optic-lobe/forbidden-matches.csv
###   alignment/presets/whole-brain/forbidden-matches.csv
### so the optic and whole-brain alignment runs can both
### consume the same source of truth.
###
### Usage:
###   Rscript alignment/banc-alignment-false-positives.R
###########################################################
source("banc/banc-startup.R")

local({

message("=== Querying SeaTable for reviewed alignment matches ===")
bc <- banctable_query(
  "SELECT _id, root_id, root_888, fafb_alignment_cell_type, fafb_alignment_decision, cell_type, fafb_cell_type FROM banc_meta"
)
message(sprintf("  Pulled %s rows from banc_meta",
                format(nrow(bc), big.mark = ",")))

# Normalize string columns (banctable returns mixed types)
bc <- bc %>%
  dplyr::mutate(
    root_id = ifelse(is.na(root_id), NA_character_, as.character(root_id)),
    root_888 = ifelse(is.na(root_888), NA_character_, as.character(root_888)),
    fafb_alignment_cell_type = ifelse(is.na(fafb_alignment_cell_type),
                                      NA_character_,
                                      as.character(fafb_alignment_cell_type)),
    fafb_alignment_decision = ifelse(is.na(fafb_alignment_decision),
                                     NA_character_,
                                     as.character(fafb_alignment_decision)),
    cell_type = ifelse(is.na(cell_type), NA_character_, as.character(cell_type)),
    fafb_cell_type = ifelse(is.na(fafb_cell_type), NA_character_,
                            as.character(fafb_cell_type))
  )

# Sanity check: align.py keys neurons by `root_id` from the meta CSV (which is
# loaded straight from banctable_query()). The forbidden CSV must use the same
# value. We verify here that root_id and root_888 agree where both exist; if
# they ever diverge, we prefer root_888 (the current proofread root) and warn.
both_present <- !is.na(bc$root_id) & nzchar(bc$root_id) &
                !is.na(bc$root_888) & nzchar(bc$root_888)
n_disagree <- sum(both_present & bc$root_id != bc$root_888)
if (n_disagree > 0) {
  message(sprintf("  WARNING: %d rows where root_id != root_888; using root_888 for those",
                  n_disagree))
}
# Build canonical id used by align.py: prefer root_888, fall back to root_id
bc$banc_id <- dplyr::if_else(!is.na(bc$root_888) & nzchar(bc$root_888),
                             bc$root_888, bc$root_id)

# Identify rejected rows. The seatable column may come back as text; accept
# any common spelling of false. NA / empty / anything else is treated as
# "not reviewed" and ignored.
decision_norm <- toupper(trimws(bc$fafb_alignment_decision))
is_false <- decision_norm %in% c("F", "FALSE")
is_true  <- decision_norm %in% c("T", "TRUE")

message(sprintf("  Reviewed: %d FALSE, %d TRUE, %d unreviewed",
                sum(is_false), sum(is_true),
                nrow(bc) - sum(is_false) - sum(is_true)))

forbidden <- bc %>%
  dplyr::filter(is_false,
                !is.na(fafb_alignment_cell_type),
                nzchar(fafb_alignment_cell_type),
                !is.na(banc_id),
                nzchar(banc_id))

if (nrow(forbidden) == 0) {
  message("  No FALSE rows with both root_id and fafb_alignment_cell_type; nothing to write.")
  return(invisible(NULL))
}

# Conflict warning: a forbidden (id, type) shouldn't match an existing GT
# cell_type for that neuron. The downstream aligner will let the anchor win,
# but we surface them here so they can be cleaned up in SeaTable.
gt_for_id <- dplyr::case_when(
  !is.na(forbidden$fafb_cell_type) & nzchar(forbidden$fafb_cell_type) &
    !grepl("^auto:", forbidden$fafb_cell_type) ~ forbidden$fafb_cell_type,
  !is.na(forbidden$cell_type) & nzchar(forbidden$cell_type) &
    !grepl("^auto:", forbidden$cell_type) ~ forbidden$cell_type,
  TRUE ~ NA_character_
)
conflict <- !is.na(gt_for_id) & gt_for_id == forbidden$fafb_alignment_cell_type
if (any(conflict)) {
  message(sprintf("  WARNING: %d forbidden pairs match an existing GT cell_type for the same neuron",
                  sum(conflict)))
  warn_df <- forbidden[conflict, c("banc_id", "fafb_alignment_cell_type")] %>%
    dplyr::distinct() %>%
    head(10)
  for (i in seq_len(nrow(warn_df))) {
    message(sprintf("    %s -> %s", warn_df$banc_id[i],
                    warn_df$fafb_alignment_cell_type[i]))
  }
  message("    (anchor wins in align.py; clean these up in SeaTable if intended)")
}

out_df <- forbidden %>%
  dplyr::transmute(root_id = banc_id, cell_type = fafb_alignment_cell_type) %>%
  dplyr::distinct(root_id, cell_type)

n_unique_ids <- dplyr::n_distinct(out_df$root_id)
n_unique_types <- dplyr::n_distinct(out_df$cell_type)
message(sprintf("  Forbidden pairs: %d (%d unique neurons, %d unique types)",
                nrow(out_df), n_unique_ids, n_unique_types))

top_types <- out_df %>%
  dplyr::count(cell_type, sort = TRUE) %>%
  head(10)
if (nrow(top_types) > 0) {
  message("  Top rejected types:")
  for (i in seq_len(nrow(top_types))) {
    message(sprintf("    %-20s %d", top_types$cell_type[i], top_types$n[i]))
  }
}

# Write to both alignment data dirs so optic and whole-brain runs share the
# same source of truth.
out_paths <- c("alignment/presets/optic-lobe/forbidden-matches.csv",
               "alignment/presets/whole-brain/forbidden-matches.csv")
for (p in out_paths) {
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(out_df, p)
  message(sprintf("  Wrote %s (%d rows)", p, nrow(out_df)))
}

message("\n=== Done ===")

})
