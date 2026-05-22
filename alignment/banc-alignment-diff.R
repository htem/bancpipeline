#' banc-alignment-diff — Diff two alignment results CSVs.
#'
#' Compare two `banc_*_alignment_results.csv` runs and summarise per-neuron
#' changes in `assigned_cell_type` / `best_target_match` / `confidence` /
#' `runner_up_type`. Writes a per-row diff CSV (only changed rows) and a
#' human-readable summary on stdout.
#'
#' @section Reads:
#'   - `<NEW_CSV>`, `<OLD_CSV>` — two alignment results CSVs
#'
#' @section Writes:
#'   - `[out.csv]` if a third positional argument is given
#'
#' @section CLI:
#'   <NEW_CSV>   required; the new results CSV
#'   <OLD_CSV>   required; the prior results CSV
#'   [out.csv]   optional; output diff CSV path
#'
#' @section Invoked by:
#'   `o2/oneshots/o2_banc_wb_dryrun_diff.sh`

###########################################################
### Diff a new whole-brain alignment results CSV against
### the previous run's CSV, summarising per-neuron changes
### in assigned_cell_type / best_target_match / confidence
### / runner_up_type.
###
### Usage:
###   Rscript alignment/banc-alignment-diff.R NEW_CSV OLD_CSV [out.csv]
###
### Writes a per-row diff CSV (only rows that differ) plus a
### human-readable summary on stdout.
###########################################################
source("banc/banc-startup.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: banc-alignment-diff.R NEW_CSV OLD_CSV [out.csv]")
}
new_path <- args[1]
old_path <- args[2]
out_path <- if (length(args) >= 3) args[3] else NA_character_

stopifnot(file.exists(new_path), file.exists(old_path))

key_cols <- c("root_888")
cmp_cols <- c("assigned_cell_type", "best_target_match",
              "confidence", "runner_up_type")

read_align <- function(p) {
  readr::read_csv(p,
    col_types = readr::cols(.default = readr::col_character(),
                            confidence = readr::col_double()),
    show_col_types = FALSE)
}

message(sprintf("=== Diff %s vs %s ===", basename(new_path), basename(old_path)))
new_df <- read_align(new_path)
old_df <- read_align(old_path)

# Normalise blanks for character columns (NA == "" for diff purposes).
norm_char <- function(x) ifelse(is.na(x) | x == "", NA_character_, x)
for (c in setdiff(cmp_cols, "confidence")) {
  new_df[[c]] <- norm_char(new_df[[c]])
  old_df[[c]] <- norm_char(old_df[[c]])
}

new_df <- dplyr::distinct(new_df, root_888, .keep_all = TRUE)
old_df <- dplyr::distinct(old_df, root_888, .keep_all = TRUE)

message(sprintf("  new: %d rows; old: %d rows", nrow(new_df), nrow(old_df)))

joined <- dplyr::full_join(
  new_df[, c(key_cols, cmp_cols)],
  old_df[, c(key_cols, cmp_cols)],
  by = key_cols, suffix = c(".new", ".old")
)

# Per-column change flags
eq_char <- function(a, b) ifelse(is.na(a) & is.na(b), TRUE,
                          ifelse(is.na(a) | is.na(b), FALSE, a == b))
eq_num  <- function(a, b, tol = 1e-4) ifelse(is.na(a) & is.na(b), TRUE,
                                       ifelse(is.na(a) | is.na(b), FALSE,
                                              abs(a - b) <= tol))

joined$cell_type_changed <- !eq_char(joined$assigned_cell_type.new,
                                     joined$assigned_cell_type.old)
joined$match_changed     <- !eq_char(joined$best_target_match.new,
                                     joined$best_target_match.old)
joined$runner_up_changed <- !eq_char(joined$runner_up_type.new,
                                     joined$runner_up_type.old)
joined$confidence_changed <- !eq_num(joined$confidence.new,
                                     joined$confidence.old)

any_change <- joined$cell_type_changed | joined$match_changed |
              joined$runner_up_changed | joined$confidence_changed

n_only_new <- sum(!joined$root_888 %in% old_df$root_888, na.rm = TRUE)
n_only_old <- sum(!joined$root_888 %in% new_df$root_888, na.rm = TRUE)

message("")
message(sprintf("  neurons only in new run: %d", n_only_new))
message(sprintf("  neurons only in old run: %d", n_only_old))
message(sprintf("  neurons in both:         %d",
                sum(joined$root_888 %in% new_df$root_888 &
                    joined$root_888 %in% old_df$root_888)))
message("")
message(sprintf("  rows with ANY change:           %d", sum(any_change)))
message(sprintf("    assigned_cell_type changed:  %d",
                sum(joined$cell_type_changed, na.rm = TRUE)))
message(sprintf("    best_target_match changed:     %d",
                sum(joined$match_changed, na.rm = TRUE)))
message(sprintf("    runner_up_type changed:      %d",
                sum(joined$runner_up_changed, na.rm = TRUE)))
message(sprintf("    confidence changed (>1e-4):  %d",
                sum(joined$confidence_changed, na.rm = TRUE)))

# Breakdown of cell_type transitions
ct_changes <- joined %>%
  dplyr::filter(cell_type_changed) %>%
  dplyr::mutate(
    transition = paste0(
      ifelse(is.na(assigned_cell_type.old), "(none)", assigned_cell_type.old),
      " -> ",
      ifelse(is.na(assigned_cell_type.new), "(none)", assigned_cell_type.new)
    )
  ) %>%
  dplyr::count(transition, sort = TRUE)

if (nrow(ct_changes) > 0) {
  message("\n  Top assigned_cell_type transitions (new <- old):")
  top <- head(ct_changes, 20)
  for (i in seq_len(nrow(top))) {
    message(sprintf("    %-40s %d", top$transition[i], top$n[i]))
  }
}

# Persist the row-level diff
if (!is.na(out_path)) {
  diff_out <- joined[any_change, ]
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(diff_out, out_path)
  message(sprintf("\n  Wrote per-row diff: %s (%d rows)",
                  out_path, nrow(diff_out)))
}

message("\n=== Done ===")
