#' banc-alignment-fill-cell-type — Fill empty cell_type from NTAC/alignment consensus.
#'
#' For BANC neurons where alignment and NTAC agree on a cell_type and
#' (a) the existing cell_type is blank or `auto:*`, (b) the super_class is
#' populated and matches the alignment-predicted super_class, and (c) the
#' super_class is not in the manual-veto list, push the consensus cell_type
#' back to SeaTable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`: cols `_id`, `root_id`, `cell_type`,
#'     `super_class`, `status`, `fafb_alignment_cell_type`,
#'     `fafb_alignment_super_class`, `fafb_alignment_decision`,
#'     `ntac_cell_type`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `cell_type`, `super_class`, `status`
#'     (status appended with `ALIGNMENT_CELL_TYPE` via `append_status`)
#'   - `<banc.meta.save.path>/snapshots/<ts>_banc_seatable_pre_ct_fill.csv`
#'
#' @section CLI:
#'   --dry-run   do not push to SeaTable
#'
#' @section Notes:
#'   - Veto list: ascending / descending / motor / visceral_circulatory.
#'   - Cell-type veto list: Tm1, Mi9, T5d (NTAC optic-lobe attractor types).

###########################################################
### Fill empty cell_type (or overwrite auto:*) from NTAC/alignment
### consensus
###
### Filters:
###   - alignment_cell_type == ntac_cell_type, both populated
###   - existing cell_type is empty OR begins with "auto:"
###   - fafb_alignment_decision != "F"
###   - existing super_class is populated AND matches predicted
###     fafb_alignment_super_class
###   - super_class NOT IN {ascending, descending, motor,
###                         visceral_circulatory} (manual-veto list)
###   - cell_type NOT IN {Tm1, Mi9, T5d} (NTAC OL-attractor types)
###
### Negative-confidence rows are retained per current request.
###
### Writes:
###   - cell_type           = predicted alignment+NTAC consensus type
###   - super_class         = predicted super_class (idempotent — already
###                            matches existing per the filter)
###   - status              = existing | "ALIGNMENT_CELL_TYPE"
###                            (via banc-functions.R append_status)
###
### Snapshot saved before push to
### <banc.meta.save.path>/snapshots/<timestamp>_banc_seatable_pre_ct_fill.csv
###
### Usage:
###   Rscript alignment/banc-alignment-fill-cell-type.R --dry-run
###   Rscript alignment/banc-alignment-fill-cell-type.R
###########################################################
source("banc/banc-startup.R")

local({

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

message("=== Fill cell_type from NTAC/alignment consensus ===")
if (dry_run) message("  *** DRY RUN — no changes will be pushed ***")

.blank_na <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""; x
}

.seatable_retry <- function(expr, label, max_tries = 6L, base_wait = 10) {
  .e <- substitute(expr)
  for (.try in seq_len(max_tries)) {
    ok <- tryCatch({ eval(.e, envir = parent.frame()); TRUE },
      error = function(e) {
        message(sprintf("  [%s] attempt %d/%d failed: %s",
                        label, .try, max_tries, conditionMessage(e))); FALSE
      })
    if (ok) {
      if (.try > 1L) message(sprintf("  [%s] succeeded on retry %d", label, .try))
      return(invisible(TRUE))
    }
    if (.try < max_tries) {
      wait <- min(300, base_wait * 2^(.try - 1))
      message(sprintf("  [%s] sleeping %ds", label, wait))
      Sys.sleep(wait)
    }
  }
  stop(sprintf("[%s] exhausted %d retries", label, max_tries), call. = FALSE)
}

# Veto lists
BANNED_SUPER_CLASSES <- c("ascending", "descending", "motor",
                          "visceral_circulatory")
NTAC_ATTRACTORS <- c("Tm1", "Mi9", "T5d")
LAMINA_TYPES <- c("R7", "R8", "L1", "L2", "L3", "L4", "L5")

###########################################################
### Re-pull SeaTable
###########################################################
message("  Re-pulling banc_meta...")
.seatable_retry(
  bc <- banctable_query(
    paste("SELECT _id, root_888, cell_type, super_class, status,",
          "fafb_alignment_cell_type, fafb_alignment_super_class,",
          "fafb_alignment_decision, target_ntac_cell_type FROM banc_meta")
  ),
  label = "banctable_query(bc)"
)
bc <- bc %>%
  dplyr::mutate(
    root_888 = .blank_na(root_888),
    cell_type = .blank_na(cell_type),
    super_class = .blank_na(super_class),
    status = .blank_na(status),
    fafb_alignment_cell_type = .blank_na(fafb_alignment_cell_type),
    fafb_alignment_super_class = .blank_na(fafb_alignment_super_class),
    fafb_alignment_decision = .blank_na(fafb_alignment_decision),
    target_ntac_cell_type = .blank_na(target_ntac_cell_type)
  ) %>%
  dplyr::filter(root_888 != "") %>%
  dplyr::distinct(root_888, .keep_all = TRUE)
message(sprintf("  SeaTable: %d rows", nrow(bc)))

###########################################################
### Snapshot before any write
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
                                    "_banc_seatable_pre_ct_fill.csv"))
  readr::write_csv(bc.orig, snapshot_file)
  message(sprintf("  Snapshot saved: %s (%d rows)",
                  basename(snapshot_file), nrow(bc.orig)))
  rm(bc.orig)
}

###########################################################
### Filter
###########################################################
eligible <- bc %>%
  dplyr::filter(
    (cell_type == "" | grepl("^auto:", cell_type)),
    fafb_alignment_cell_type != "",
    target_ntac_cell_type != "",
    fafb_alignment_cell_type == target_ntac_cell_type,
    fafb_alignment_decision != "F",
    # super_class either empty (we'll fill it) OR matches prediction
    (super_class == "" | super_class == fafb_alignment_super_class),
    !super_class %in% BANNED_SUPER_CLASSES,
    !fafb_alignment_super_class %in% BANNED_SUPER_CLASSES,
    !fafb_alignment_cell_type %in% NTAC_ATTRACTORS,
    !fafb_alignment_cell_type %in% LAMINA_TYPES
  )
message(sprintf("\n  Eligible rows: %d", nrow(eligible)))

# Funnel diagnostics
base <- bc %>%
  dplyr::filter((cell_type == "" | grepl("^auto:", cell_type)),
                fafb_alignment_cell_type != "",
                target_ntac_cell_type != "",
                fafb_alignment_cell_type == target_ntac_cell_type)
message(sprintf("  Base agreement pool: %d", nrow(base)))
n_decF <- sum(base$fafb_alignment_decision == "F")
n_no_sc <- sum(base$fafb_alignment_decision != "F" & base$super_class == "")
n_sc_mis <- sum(base$fafb_alignment_decision != "F" &
                base$super_class != "" &
                base$super_class != base$fafb_alignment_super_class)
n_banned_sc <- sum(base$super_class %in% BANNED_SUPER_CLASSES)
n_attractor <- sum(base$fafb_alignment_cell_type %in% NTAC_ATTRACTORS)
n_sc_mis2 <- sum(base$fafb_alignment_decision != "F" &
                 base$super_class != "" &
                 base$super_class != base$fafb_alignment_super_class)
n_lamina <- sum(base$fafb_alignment_cell_type %in% LAMINA_TYPES)
message(sprintf("    excluded by decision == F:                  %d", n_decF))
message(sprintf("    INCLUDED (was-excluded) no existing super_class: %d", n_no_sc))
message(sprintf("    excluded by super_class mismatch:           %d", n_sc_mis2))
message(sprintf("    excluded by banned super_class:             %d", n_banned_sc))
message(sprintf("    excluded by NTAC attractor cell_type:       %d", n_attractor))
message(sprintf("    excluded by lamina cell_type:               %d", n_lamina))

if (nrow(eligible) == 0) {
  message("  Nothing to push — exiting.")
  return(invisible(NULL))
}

# Distribution
message("\n  Super_class distribution:")
dist_sc <- eligible %>% dplyr::count(super_class, sort = TRUE)
for (i in seq_len(nrow(dist_sc))) {
  r <- dist_sc[i, ]
  message(sprintf("    %-30s  %6d", r$super_class, r$n))
}
message("\n  Top 15 cell_types proposed:")
dist_ct <- eligible %>%
  dplyr::count(fafb_alignment_cell_type, sort = TRUE) %>% head(15)
for (i in seq_len(nrow(dist_ct))) {
  r <- dist_ct[i, ]
  message(sprintf("    %-30s  %6d", r$fafb_alignment_cell_type, r$n))
}
n_auto_overwrite <- sum(grepl("^auto:", eligible$cell_type))
message(sprintf("\n  Rows whose existing cell_type carries auto: prefix (will be overwritten): %d",
                n_auto_overwrite))

###########################################################
### Build push frame
###########################################################
push_df <- eligible %>%
  dplyr::mutate(
    new_status = append_status(status, "ALIGNMENT_CELL_TYPE")
  ) %>%
  dplyr::transmute(
    `_id`       = `_id`,
    cell_type   = fafb_alignment_cell_type,
    super_class = fafb_alignment_super_class,
    status      = new_status
  ) %>%
  as.data.frame(stringsAsFactors = FALSE, check.names = FALSE)

# Final NA scrub (defense in depth)
for (.cn in colnames(push_df)) {
  if (is.character(push_df[[.cn]])) {
    push_df[[.cn]][is.na(push_df[[.cn]])] <- ""
  }
}

stopifnot(all(push_df$cell_type != ""))
stopifnot(!any(push_df$super_class %in% BANNED_SUPER_CLASSES))
stopifnot(!any(push_df$cell_type %in% NTAC_ATTRACTORS))
stopifnot(!any(push_df$cell_type %in% LAMINA_TYPES))
stopifnot(!any(is.na(push_df)))

# A few example rows to eyeball
message("\n  First 5 push rows:")
print(utils::head(push_df, 5))

###########################################################
### Push in 10k batches with retry
###########################################################
if (!dry_run) {
  message("\n  Pushing to SeaTable in 10k-row batches (chunksize=100)...")
  .batch_size <- 10000L
  n_total <- nrow(push_df)
  n_batches <- ceiling(n_total / .batch_size)
  for (.b in seq_len(n_batches)) {
    .from <- (.b - 1L) * .batch_size + 1L
    .to   <- min(.b * .batch_size, n_total)
    .chunk <- push_df[.from:.to, , drop = FALSE]
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
