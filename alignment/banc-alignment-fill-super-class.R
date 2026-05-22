#' banc-alignment-fill-super-class — Fill empty `super_class` from NTAC/alignment consensus.
#'
#' For BANC neurons where SeaTable `super_class` is blank, the alignment-implied
#' and NTAC-implied super_classes agree, and the predicted super_class lies in
#' a safe whitelist (optic_lobe_intrinsic / central_brain_intrinsic /
#' visual_projection / visual_centrifugal), push the consensus class plus the
#' modal region / flow back to SeaTable. Preserves existing non-blank values.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`: cols `_id`, `root_id`, `super_class`, `region`,
#'     `flow`, `fafb_alignment_super_class`, `fafb_alignment_decision`,
#'     `ntac_super_class`, etc.
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `super_class`, `region`, `flow`
#'   - `<banc.meta.save.path>/snapshots/<ts>_banc_seatable_pre_sc_fill.csv` (pre-push snapshot)
#'
#' @section CLI:
#'   --dry-run   do not push to SeaTable; report what would change
#'
#' @section Notes:
#'   - Excludes ascending / descending / sensory / motor from auto-fill —
#'     consensus is unreliable there.
#'   - Final NA scrub dodges the previous wipe bug where SeaTable received
#'     empty strings over existing values.

###########################################################
### Fill empty `super_class` from NTAC/alignment consensus
###
### For rows where:
###   - GT `super_class` is empty
###   - `fafb_alignment_decision` is NOT "F"
###   - alignment-implied super_class matches NTAC-implied super_class
###   - consensus super_class ∈ {optic_lobe_intrinsic,
###                              central_brain_intrinsic,
###                              visual_projection,
###                              visual_centrifugal}
###     (ascending / descending / sensory / motor are manually-verified
###     to be wrong from this consensus, so excluded)
###
### Writes super_class + correlated region + flow per the BANC modal
### conventions probed from banc_meta.
###
### Snapshot saved before push to <banc.meta.save.path>/snapshots/.
### Includes preserve-on-blank logic — never writes "" over an existing
### value — and final NA scrub to dodge the previous wipe bug.
###
### Usage:
###   Rscript alignment/banc-alignment-fill-super-class.R --dry-run
###   Rscript alignment/banc-alignment-fill-super-class.R
###########################################################
source("banc/banc-startup.R")

local({

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

message("=== Fill empty super_class from NTAC/alignment consensus ===")
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

# Allowed consensus super_classes + their correlated (region, flow).
# Suffixed `.new` to avoid collision with bc's existing region / flow cols
# during the left_join below.
sc_map <- data.frame(
  super_class = c("optic_lobe_intrinsic",
                  "central_brain_intrinsic",
                  "visual_projection",
                  "visual_centrifugal"),
  region.new  = c("optic_lobe", "central_brain", "optic_lobe", "central_brain"),
  flow.new    = c("intrinsic",  "intrinsic",     "intrinsic",  "intrinsic"),
  stringsAsFactors = FALSE
)

###########################################################
### Re-pull SeaTable (user just made manual changes)
###########################################################
message("  Re-pulling banc_meta (user reports recent manual edits)...")
.seatable_retry(
  bc <- banctable_query(
    paste("SELECT _id, root_888, cell_type, super_class, region, flow,",
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
    region = .blank_na(region),
    flow = .blank_na(flow),
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
                                    "_banc_seatable_pre_sc_fill.csv"))
  readr::write_csv(bc.orig, snapshot_file)
  message(sprintf("  Snapshot saved: %s (%d rows)",
                  basename(snapshot_file), nrow(bc.orig)))
  rm(bc.orig)
}

###########################################################
### NTAC-implied super_class via franken_meta lookup
###########################################################
message("  Loading FAFB super_class lookup for NTAC type translation...")
target.meta <- franken_meta("SELECT cell_type, super_class FROM franken_meta") %>%
  dplyr::filter(!is.na(cell_type), cell_type != "",
                !is.na(super_class), super_class != "") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(super_class = names(which.max(table(super_class))),
                   .groups = "drop") %>%
  dplyr::rename(ntac_sc = super_class)

bc <- bc %>%
  dplyr::left_join(target.meta, by = c("target_ntac_cell_type" = "cell_type")) %>%
  dplyr::mutate(ntac_sc = .blank_na(ntac_sc))

###########################################################
### Eligible rows
###########################################################
LAMINA_TYPES <- c(paste0("L", 1:5), paste0("R", 1:8))
.is_lamina <- function(x) toupper(trimws(x)) %in% LAMINA_TYPES

eligible <- bc %>%
  dplyr::filter(
    super_class == "",                          # empty GT
    fafb_alignment_decision != "F",             # not user-rejected
    fafb_alignment_super_class != "",
    ntac_sc != "",
    fafb_alignment_super_class == ntac_sc,      # consensus
    fafb_alignment_super_class %in% sc_map$super_class,  # allowed only
    !.is_lamina(fafb_alignment_cell_type),      # drop L1-5/R1-8 calls
    !.is_lamina(target_ntac_cell_type)
  ) %>%
  dplyr::left_join(sc_map, by = c("fafb_alignment_super_class" = "super_class"))
message(sprintf("\n  Eligible (consensus, decision!=F, allowed sc, no lamina): %d rows",
                nrow(eligible)))

# Region-conflict skip: if existing region is populated AND disagrees with the
# super_class's implied region, the consensus is likely wrong — drop it.
.region_conflict <- eligible$region != "" & eligible$region != eligible$region.new
n_region_conflict <- sum(.region_conflict)
if (n_region_conflict > 0) {
  message(sprintf("  Skipping %d rows where existing region conflicts with proposed region",
                  n_region_conflict))
  eligible <- eligible[!.region_conflict, , drop = FALSE]
}
message(sprintf("  Final eligible after region-conflict skip: %d rows", nrow(eligible)))

if (nrow(eligible) == 0) {
  message("  Nothing to push — exiting.")
  return(invisible(NULL))
}

###########################################################
### Build push frame — sc + correlated region + flow
### region/flow are PRESERVE-ON-BLANK: only fill if empty in SeaTable.
### (super_class is always written — eligibility already requires it empty.)
###########################################################
push_df <- eligible %>%
  dplyr::transmute(
    `_id`       = `_id`,
    super_class = fafb_alignment_super_class,
    region      = ifelse(region == "", region.new, region),
    flow        = ifelse(flow   == "", flow.new,   flow)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE, check.names = FALSE)

# Final NA scrub (character cols only — defense in depth)
for (.cn in colnames(push_df)) {
  if (is.character(push_df[[.cn]])) {
    push_df[[.cn]][is.na(push_df[[.cn]])] <- ""
  }
}

# Distribution + sanity
message("\n  Push distribution by super_class:")
dist <- push_df %>% dplyr::count(super_class, sort = TRUE)
for (i in seq_len(nrow(dist))) {
  r <- dist[i, ]
  message(sprintf("    %-30s  %6d", r$super_class, r$n))
}

stopifnot(all(push_df$super_class %in% sc_map$super_class))
stopifnot(all(push_df$region != ""))
stopifnot(all(push_df$flow != ""))
stopifnot(!any(is.na(push_df)))

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
