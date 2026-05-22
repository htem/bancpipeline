#' update-seatable (whole-brain preset) — Push whole-brain alignment + NTAC predictions to SeaTable.
#'
#' Loaded by `alignment/banc-alignment-update-seatable.R --region whole-brain`.
#' Uploads connectivity-alignment predictions and the most-anchored NTAC
#' sweep result, joining on `root_888` and `root_850` respectively. Writes
#' only to the `fafb_alignment_*` / `fafb_ntac_*` column family — never to
#' GT cell_type / super_class.
#'
#' @section Reads:
#'   - `[alignment_csv]` positional (default the latest WB production CSV)
#'   - env var `BANC_WB_NTAC_CSV` — override NTAC source CSV
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `fafb_alignment_cell_type`,
#'     `fafb_alignment_match`, `fafb_alignment_confidence`,
#'     `fafb_alignment_runner_up`, `fafb_alignment_super_class`,
#'     `fafb_alignment_decision` (cleared iff cell_type changes),
#'     `fafb_ntac_cell_type`
#'   - `<banc.meta.save.path>/snapshots/<ts>_banc_seatable.csv` (pre-push snapshot)
#'
#' @section CLI (via dispatcher):
#'   --region whole-brain [alignment_csv]
#'   --dry-run
#'
#' @section Invoked by:
#'   `o2/alignment/o2_banc_wb_push.sh`,
#'   `o2/oneshots/o2_banc_wb_dryrun_diff.sh`
#'
#' @section Notes:
#'   - Preserve-on-blank rule: blank push values never overwrite populated
#'     SeaTable cells (guards against the earlier pivot_wider wipe bug).

###########################################################
### Whole-Brain Alignment: SeaTable Update
###
### Uploads connectivity-alignment predictions from the WB-production run
### (and the most-anchored WB NTAC sweep) to BANC SeaTable, joining on
### root_888 for alignment results and root_850 for NTAC results.
###
### We only ever write to the `fafb_alignment_*` and `fafb_ntac_*` column
### family — never to GT columns (cell_type, super_class, etc.). So no
### GT skip-list is applied; every alignment row with a prediction is
### pushed.
###
### Preserve-on-blank rule
### ----------------------
### If our new value for a target column is blank and SeaTable already
### has a populated value there, we KEEP the existing value. This
### protects against an earlier class of bug where pivot_wider() / NA
### gaps in a push frame silently wiped curated cells.
###
### Special case for fafb_alignment_decision
### ----------------------------------------
### A prior T/F review applies to the PRIOR fafb_alignment_cell_type
### call. If we are overwriting that cell_type with a different value,
### the review is stale → we explicitly blank fafb_alignment_decision.
### If our new cell_type equals the existing one, the review is still
### valid → preserve.
###
### Snapshot before push
### --------------------
### A full-table CSV snapshot is written to
###   <banc.meta.save.path>/snapshots/<timestamp>_banc_seatable.csv
### so the push can be reverted if needed.
###
### Columns updated:
###   - fafb_alignment_cell_type
###   - fafb_alignment_match
###   - fafb_alignment_confidence
###   - fafb_alignment_runner_up
###   - fafb_alignment_super_class
###   - fafb_alignment_decision       (wiped iff cell_type changes)
###   - fafb_ntac_cell_type           (from a100_ho000 v850 NTAC,
###                                    joined on root_850)
###
### Usage:
###   Rscript alignment/banc-alignment-update-seatable.R --region whole-brain [alignment_csv]
###   Rscript alignment/banc-alignment-update-seatable.R --region whole-brain --dry-run
###
### Optional env vars:
###   BANC_WB_NTAC_CSV — override path to the NTAC source CSV.
###     Default: data/whole_brain_alignment/
###              banc_brain_both_ntac_sweep_ntac_a100_ho000.csv
###########################################################
source("banc/banc-startup.R")

local({

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
args <- args[args != "--dry-run"]

results_file <- if (length(args) > 0) {
  args[1]
} else {
  "data/whole_brain_alignment_v888v2/banc_brain_both_alignment_production_v888v2_all_tier1.csv"
}
if (!file.exists(results_file)) stop("Results file not found: ", results_file)

ntac_file <- Sys.getenv(
  "BANC_WB_NTAC_CSV",
  unset = "data/whole_brain_alignment/banc_brain_both_ntac_sweep_ntac_a100_ho000.csv"
)

message(sprintf("=== Preparing SeaTable update from %s ===", results_file))
if (dry_run) message("  *** DRY RUN — no changes will be pushed ***")

###########################################################
### Helper: NA-safe character coalesce
###########################################################
.blank_na <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

###########################################################
### Helper: SeaTable retry wrapper (copied from banc-update-seatable.R)
### Survives transient ReadTimeouts on cloud.seatable.io.
###########################################################
.seatable_retry <- function(expr, label, max_tries = 6L, base_wait = 10) {
  .e <- substitute(expr)
  for (.try in seq_len(max_tries)) {
    ok <- tryCatch({ eval(.e, envir = parent.frame()); TRUE },
      error = function(e) {
        msg <- conditionMessage(e)
        message(sprintf("  [%s] attempt %d/%d failed: %s",
                        label, .try, max_tries, msg))
        FALSE
      })
    if (ok) {
      if (.try > 1L) message(sprintf("  [%s] succeeded on retry %d", label, .try))
      return(invisible(TRUE))
    }
    if (.try < max_tries) {
      wait <- min(300, base_wait * 2^(.try - 1))
      message(sprintf("  [%s] sleeping %ds before retry", label, wait))
      Sys.sleep(wait)
    }
  }
  stop(sprintf("[%s] exhausted %d retries", label, max_tries), call. = FALSE)
}

###########################################################
### Load alignment results (root_888 keyed)
###########################################################
alignment <- readr::read_csv(results_file,
  col_types = readr::cols(root_888 = "c", best_target_match = "c"),
  show_col_types = FALSE)
message(sprintf("  Alignment results: %d neurons (%d with cell_type, %d with FAFB match)",
                nrow(alignment),
                sum(alignment$assigned_cell_type != "", na.rm = TRUE),
                sum(!is.na(alignment$best_target_match) & alignment$best_target_match != "")))

###########################################################
### Query SeaTable — fetch _id, key mapping, AND current
### values of target columns (for preserve-on-blank).
###########################################################
message("  Querying SeaTable for _id mapping + current target values...")
.seatable_retry(
  bc <- banctable_query(
    paste(
      "SELECT _id, root_id, root_850, root_888,",
      "fafb_alignment_cell_type, fafb_alignment_match,",
      "fafb_alignment_confidence, fafb_alignment_runner_up,",
      "fafb_alignment_super_class, fafb_alignment_decision,",
      "fafb_ntac_cell_type FROM banc_meta"
    )
  ),
  label = "banctable_query(bc)"
)
bc <- bc %>%
  dplyr::filter(!is.na(root_888), root_888 != "") %>%
  dplyr::mutate(root_888 = as.character(root_888),
                root_850 = as.character(root_850)) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)
message(sprintf("  SeaTable: %d neurons with root_888 (%d also with root_850)",
                nrow(bc), sum(!is.na(bc$root_850) & bc$root_850 != "")))

###########################################################
### Snapshot — full table dump before any write
###########################################################
if (!dry_run) {
  message("  Taking full-table snapshot before push...")
  .seatable_retry(bc.orig <- banctable_query(),
                  label = "banctable_query(bc.orig)")
  dir.create(file.path(banc.meta.save.path, "snapshots"),
             showWarnings = FALSE, recursive = TRUE)
  datetime_string <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  snapshot_file <- file.path(banc.meta.save.path, "snapshots",
                             paste0(datetime_string,
                                    "_banc_seatable_pre_wb_alignment.csv"))
  readr::write_csv(bc.orig, snapshot_file)
  message(sprintf("  Snapshot saved: %s (%d rows)",
                  basename(snapshot_file), nrow(bc.orig)))
  rm(bc.orig)
}

###########################################################
### FAFB super_class lookup (modal super_class per cell_type)
###########################################################
message("  Loading FAFB super_class lookup from franken_meta...")
fafb.meta <- franken_meta(
  "SELECT cell_type, super_class FROM franken_meta"
) %>%
  dplyr::filter(!is.na(cell_type), cell_type != "",
                !is.na(super_class), super_class != "") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(super_class = names(which.max(table(super_class))),
                   .groups = "drop")
message(sprintf("  FAFB type→super_class lookup: %d unique cell_types",
                nrow(fafb.meta)))

###########################################################
### Build push frame — NEW values only at this stage
###########################################################
push_df <- alignment %>%
  dplyr::filter(assigned_cell_type != "" | !is.na(best_target_match)) %>%
  dplyr::inner_join(bc, by = "root_888", suffix = c("", ".bc")) %>%
  dplyr::left_join(fafb.meta, by = c("assigned_cell_type" = "cell_type"),
                   suffix = c("", ".predicted")) %>%
  dplyr::mutate(
    new_cell_type   = .blank_na(assigned_cell_type),
    new_match       = .blank_na(best_target_match),
    # confidence stays numeric — SeaTable column is Number type, not Text.
    # Use NA for "no new value" (preserve-on-blank handled with is.na below).
    new_confidence  = round(confidence, 4),
    new_runner_up   = .blank_na(runner_up_type),
    # super_class is brought in by the fafb.meta left_join. The bc query
    # does NOT select the GT super_class column, so the FAFB column keeps
    # its bare name `super_class` with no collision.
    new_super_class = .blank_na(super_class),
    # NTAC join handled separately below.
    new_ntac        = ""
  ) %>%
  as.data.frame()

message(sprintf("  Initial push frame: %d rows (no GT guard)", nrow(push_df)))

###########################################################
### NTAC join (root_850 keyed) — fill new_ntac
###########################################################
if (file.exists(ntac_file)) {
  message(sprintf("  Loading NTAC results: %s", ntac_file))
  ntac_hdr <- readr::read_csv(ntac_file, n_max = 0, show_col_types = FALSE)
  ntac_key <- if ("root_888" %in% colnames(ntac_hdr)) {
    "root_888"
  } else if ("root_850" %in% colnames(ntac_hdr)) {
    "root_850"
  } else {
    stop("NTAC CSV missing expected root_888 or root_850 column: ", ntac_file)
  }
  message(sprintf("  NTAC key column: %s", ntac_key))
  ntac <- readr::read_csv(ntac_file,
    col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  ntac <- ntac %>%
    dplyr::transmute(
      !!ntac_key := as.character(.data[[ntac_key]]),
      ntac_pred = ifelse(is.na(ntac_cell_type) | ntac_cell_type == "?",
                         "", as.character(ntac_cell_type))
    ) %>%
    dplyr::filter(!is.na(.data[[ntac_key]]), .data[[ntac_key]] != "") %>%
    dplyr::distinct(.data[[ntac_key]], .keep_all = TRUE)
  push_df <- push_df %>%
    dplyr::left_join(ntac, by = ntac_key)
  push_df$new_ntac <- .blank_na(push_df$ntac_pred)
  push_df$ntac_pred <- NULL
  n_ntac_raw <- sum(push_df$new_ntac != "")
  message(sprintf("  NTAC joined on %s: %d non-empty (of %d rows)",
                  ntac_key, n_ntac_raw, nrow(push_df)))
} else {
  message(sprintf("  NTAC results not found at %s — leaving NTAC empty",
                  ntac_file))
}

###########################################################
### Preserve-on-blank merge
###
### For each target column: if NEW value is "" but the EXISTING SeaTable
### value is populated, keep the existing value. Otherwise write NEW.
###
### Special case for fafb_alignment_decision:
###   - if cell_type is changing (new_cell_type != existing): WIPE (set "")
###   - else: preserve existing
###########################################################
existing <- list(
  cell_type   = .blank_na(push_df$fafb_alignment_cell_type),
  match       = .blank_na(push_df$fafb_alignment_match),
  confidence  = suppressWarnings(as.numeric(push_df$fafb_alignment_confidence)),
  runner_up   = .blank_na(push_df$fafb_alignment_runner_up),
  super_class = .blank_na(push_df$fafb_alignment_super_class),
  decision    = .blank_na(push_df$fafb_alignment_decision),
  ntac        = .blank_na(push_df$fafb_ntac_cell_type)
)
new_vals <- list(
  cell_type   = .blank_na(push_df$new_cell_type),
  match       = .blank_na(push_df$new_match),
  confidence  = as.numeric(push_df$new_confidence),  # numeric
  runner_up   = .blank_na(push_df$new_runner_up),
  super_class = .blank_na(push_df$new_super_class),
  ntac        = .blank_na(push_df$new_ntac)
)
.preserve <- function(new, old) ifelse(new == "" & old != "", old, new)
.preserve_num <- function(new, old) ifelse(is.na(new) & !is.na(old), old, new)

final_cell_type   <- .preserve(new_vals$cell_type,   existing$cell_type)
final_match       <- .preserve(new_vals$match,       existing$match)
final_confidence  <- .preserve_num(new_vals$confidence, existing$confidence)
final_runner_up   <- .preserve(new_vals$runner_up,   existing$runner_up)
final_super_class <- .preserve(new_vals$super_class, existing$super_class)
final_ntac        <- .preserve(new_vals$ntac,        existing$ntac)

# decision: wipe only when cell_type actually changes
ct_changes <- final_cell_type != "" & final_cell_type != existing$cell_type
final_decision <- ifelse(ct_changes, "", existing$decision)

push_final <- data.frame(
  `_id`                       = push_df$`_id`,
  fafb_alignment_cell_type    = final_cell_type,
  fafb_alignment_match        = final_match,
  fafb_alignment_confidence   = final_confidence,
  fafb_alignment_runner_up    = final_runner_up,
  fafb_alignment_super_class  = final_super_class,
  fafb_alignment_decision     = final_decision,
  fafb_ntac_cell_type         = final_ntac,
  stringsAsFactors            = FALSE,
  check.names                 = FALSE
)

###########################################################
### Final NA scrub — defense-in-depth against the wipe bug.
###
### Character columns: any residual NA → "" so SeaTable never sees a
### null in a text cell (its API otherwise reads NULL as "clear field").
### Numeric columns (confidence): NA is the natural empty-number sentinel
### and is left as-is — SeaTable Number type accepts NULL for empty.
###########################################################
.char_cols <- vapply(push_final, is.character, logical(1))
.na_count_char <- sum(vapply(push_final[, .char_cols, drop = FALSE],
                             function(x) sum(is.na(x)), integer(1)))
if (.na_count_char > 0) {
  message(sprintf("  WARNING: %d NA cells found in character columns; scrubbing to \"\"",
                  .na_count_char))
  for (.cn in colnames(push_final)[.char_cols]) {
    push_final[[.cn]][is.na(push_final[[.cn]])] <- ""
  }
}

###########################################################
### Summary
###########################################################
n_typed       <- sum(push_final$fafb_alignment_cell_type != "")
n_matched     <- sum(push_final$fafb_alignment_match != "")
n_super       <- sum(push_final$fafb_alignment_super_class != "")
n_decision    <- sum(push_final$fafb_alignment_decision != "")
n_ntac        <- sum(push_final$fafb_ntac_cell_type != "")
n_unique_types <- length(unique(
  push_final$fafb_alignment_cell_type[push_final$fafb_alignment_cell_type != ""]))
n_decision_wiped <- sum(ct_changes)
n_preserved_blank <- sum(
  (new_vals$cell_type == "" & existing$cell_type != "") |
  (new_vals$match == "" & existing$match != "") |
  (new_vals$super_class == "" & existing$super_class != "") |
  (new_vals$ntac == "" & existing$ntac != "")
)

message(sprintf("\n  Final push frame: %d rows", nrow(push_final)))
message(sprintf("    fafb_alignment_cell_type:   %d non-empty (%d unique types)",
                n_typed, n_unique_types))
message(sprintf("    fafb_alignment_match:       %d non-empty", n_matched))
message(sprintf("    fafb_alignment_super_class: %d non-empty", n_super))
message(sprintf("    fafb_alignment_decision:    %d non-empty (wiped %d rows where cell_type changed)",
                n_decision, n_decision_wiped))
message(sprintf("    fafb_ntac_cell_type:        %d non-empty", n_ntac))
message(sprintf("    preserve-on-blank kept %d existing values across cell_type/match/super_class/ntac",
                n_preserved_blank))

###########################################################
### Push
###########################################################
if (!dry_run) {
  # Push in batches with retry. Smaller chunksize (100) keeps each
  # SeaTable API call below the 30s read-timeout; outer batching of
  # 10k rows lets the retry wrapper restart from a known offset if
  # cloud.seatable.io flakes mid-push (the inner banctable_update_rows
  # pblapply does not itself retry).
  message("\n  Pushing to SeaTable in 10k-row batches (chunksize=100)...")
  .batch_size <- 10000L
  n_total <- nrow(push_final)
  n_batches <- ceiling(n_total / .batch_size)
  for (.b in seq_len(n_batches)) {
    .from <- (.b - 1L) * .batch_size + 1L
    .to   <- min(.b * .batch_size, n_total)
    .chunk <- push_final[.from:.to, , drop = FALSE]
    message(sprintf("  Batch %d/%d: rows %d-%d (%d rows)",
                    .b, n_batches, .from, .to, nrow(.chunk)))
    .seatable_retry(
      banctable_update_rows(base = "banc_meta",
                            table = "banc_meta",
                            df = .chunk,
                            append_allowed = FALSE,
                            chunksize = 100L),
      label = sprintf("batch_%02d", .b)
    )
  }
  message("  SeaTable update complete!")
} else {
  message("\n  Dry run complete. Remove --dry-run to push.")
}

message("\n=== Done ===")

})
