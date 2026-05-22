#' banc-alignment-false-positives-append — Append-mode refresh of alignment blacklist.
#'
#' Mirrors `banc-alignment-false-positives.R` but UNIONS the fresh SeaTable
#' pull with whatever is already on disk so rows already in the blacklist
#' survive even when they no longer carry `fafb_alignment_decision == F`
#' in SeaTable (e.g. cell deleted, decision reverted, manual addition).
#'
#' @section Reads:
#'   - SeaTable `banc_meta`: cols `_id`, `root_id`, `root_888`,
#'     `fafb_alignment_cell_type`, `fafb_alignment_decision`,
#'     `cell_type`, `fafb_cell_type`
#'   - `alignment/presets/optic-lobe/forbidden-matches.csv` (existing)
#'   - `alignment/presets/whole-brain/forbidden-matches.csv` (existing)
#'
#' @section Writes:
#'   - `alignment/presets/optic-lobe/forbidden-matches.csv`
#'   - `alignment/presets/whole-brain/forbidden-matches.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_refresh_blacklist.sh`

###########################################################
### Append-mode refresh of forbidden-matches.csv blacklists
###
### Mirrors banc-alignment-false-positives.R but UNIONS the
### fresh SeaTable pull with whatever is already on disk at
###   alignment/presets/optic-lobe/forbidden-matches.csv
###   alignment/presets/whole-brain/forbidden-matches.csv
### so rows already in the blacklist are preserved even if
### they no longer carry fafb_alignment_decision == F in
### SeaTable (e.g. cell removed, decision reverted, manual
### addition).
###
### Usage:
###   Rscript alignment/banc-alignment-false-positives-append.R
###########################################################
source("banc/banc-startup.R")

local({

out_paths <- c("alignment/presets/optic-lobe/forbidden-matches.csv",
               "alignment/presets/whole-brain/forbidden-matches.csv")

# 1. Existing on-disk blacklists ------------------------------------------
existing <- lapply(out_paths, function(p) {
  if (file.exists(p)) {
    df <- readr::read_csv(p, col_types = readr::cols(.default = readr::col_character()))
    message(sprintf("  Existing %s: %d rows", p, nrow(df)))
    df
  } else {
    message(sprintf("  Existing %s: missing -- treating as empty", p))
    tibble::tibble(root_id = character(), cell_type = character())
  }
})

# 2. Fresh SeaTable pull of fafb_alignment_decision == F -----------------
message("=== Querying SeaTable for reviewed alignment matches ===")
bc <- banctable_query(
  "SELECT _id, root_id, root_888, fafb_alignment_cell_type, fafb_alignment_decision, cell_type, fafb_cell_type FROM banc_meta"
)
message(sprintf("  Pulled %s rows from banc_meta",
                format(nrow(bc), big.mark = ",")))

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

both_present <- !is.na(bc$root_id) & nzchar(bc$root_id) &
                !is.na(bc$root_888) & nzchar(bc$root_888)
n_disagree <- sum(both_present & bc$root_id != bc$root_888)
if (n_disagree > 0) {
  message(sprintf("  WARNING: %d rows where root_id != root_888; using root_888 for those",
                  n_disagree))
}
bc$banc_id <- dplyr::if_else(!is.na(bc$root_888) & nzchar(bc$root_888),
                             bc$root_888, bc$root_id)

decision_norm <- toupper(trimws(bc$fafb_alignment_decision))
is_false <- decision_norm %in% c("F", "FALSE")
is_true  <- decision_norm %in% c("T", "TRUE")
message(sprintf("  Reviewed: %d FALSE, %d TRUE, %d unreviewed",
                sum(is_false), sum(is_true),
                nrow(bc) - sum(is_false) - sum(is_true)))

fresh <- bc %>%
  dplyr::filter(is_false,
                !is.na(fafb_alignment_cell_type),
                nzchar(fafb_alignment_cell_type),
                !is.na(banc_id),
                nzchar(banc_id)) %>%
  dplyr::transmute(root_id = banc_id, cell_type = fafb_alignment_cell_type) %>%
  dplyr::distinct(root_id, cell_type)

message(sprintf("  Fresh F pairs from SeaTable: %d (%d unique neurons, %d unique types)",
                nrow(fresh),
                dplyr::n_distinct(fresh$root_id),
                dplyr::n_distinct(fresh$cell_type)))

# 3. Union + dedupe per output path ---------------------------------------
for (i in seq_along(out_paths)) {
  p <- out_paths[[i]]
  prev <- existing[[i]]
  prev$root_id <- as.character(prev$root_id)
  prev$cell_type <- as.character(prev$cell_type)

  merged <- dplyr::bind_rows(prev[, c("root_id", "cell_type")],
                             fresh[, c("root_id", "cell_type")]) %>%
    dplyr::filter(!is.na(root_id), nzchar(root_id),
                  !is.na(cell_type), nzchar(cell_type)) %>%
    dplyr::distinct(root_id, cell_type) %>%
    dplyr::arrange(root_id, cell_type)

  added <- nrow(merged) - nrow(prev)
  message(sprintf("  %s: %d -> %d rows (+%d new)",
                  p, nrow(prev), nrow(merged), added))

  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(merged, p)
  message(sprintf("  Wrote %s", p))
}

message("\n=== Done ===")

})
