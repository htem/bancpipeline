#' banc-ids888 — Standalone root_888 backfill for SeaTable banc_meta.
#'
#' Lifted from `banc-ids.R` so the v888 backfill can run independently of
#' the rest (which can crash early in `banc_proofreading_notes()`). Reads
#' SeaTable, finds rows with missing / empty / zero `root_888`, resolves
#' their `supervoxel_id` at v888 via CAVE `banc_rootid()`, and writes back
#' via `banctable_update_rows()`. Idempotent.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - CAVE `banc_rootid()` at v888
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `root_id`, `root_888`, `supervoxel_id`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_v888_rebuild.sh`

###############################################################################
### Standalone root_888 backfill for SeaTable banc_meta.
###
### Lifted from banc/update/banc-ids.R lines 297–312 (the bc888 block) so that
### it can run independently of the rest of banc-ids.R (which currently crashes
### early in banc_proofreading_notes()).
###
### Reads SeaTable banc_meta, finds rows with missing/empty/zero root_888,
### resolves their supervoxel_id at v888 via CAVE banc_rootid(), and writes
### root_id/root_888/supervoxel_id back via banctable_update_rows().
###
### Idempotent: re-running on a fully populated table is a no-op.
###############################################################################
source("banc/banc-startup.R")

bancr::choose_banc()

message(sprintf("[%s] querying banc_meta for rows missing root_888...",
                Sys.time()))

bc.all <- banctable_query(
  "SELECT _id, status, root_id, root_888, supervoxel_id, position from banc_meta"
)
message(sprintf("  banc_meta total rows: %s",
                format(nrow(bc.all), big.mark = ",")))

bc888 <- bc.all %>%
  dplyr::select(`_id`, root_888, root_id, supervoxel_id, position) %>%
  dplyr::filter(
    is.na(root_888) | root_888 == "" | root_888 == "0" | root_888 == 0 |
    is.na(supervoxel_id) | supervoxel_id == "0" | supervoxel_id == 0
  )

message(sprintf("  rows needing root_888 backfill: %s",
                format(nrow(bc888), big.mark = ",")))

if (nrow(bc888) == 0L) {
  message("nothing to update -- root_888 already populated for all rows")
} else {
  # 3-pass resolution. Previous version unconditionally overwrote
  # supervoxel_id via banc_xyz2id(position), clobbering valid svids whenever
  # `position` was empty or outside segmentation — leaving root_888 stuck NA.
  # Now: only fill svid where it's missing; final fallback uses root_id via
  # banc_latestid() for rows that still can't be resolved.
  empty_like <- function(x) {
    is.na(x) | x == "" | x == "0" | x == 0
  }

  # ---- Pass 1: existing valid svid -> root_888 via CAVE
  svid_ok <- !empty_like(bc888$supervoxel_id)
  message(sprintf("[%s] pass 1: resolving root_888 from %s existing supervoxel_id values...",
                  Sys.time(), format(sum(svid_ok), big.mark = ",")))
  if (any(svid_ok)) {
    bc888$root_888[svid_ok] <- bancr::banc_rootid(bc888$supervoxel_id[svid_ok],
                                                  version = "888")
  }

  # ---- Pass 2: rows still missing root_888 AND missing svid -> recompute
  # svid from position, then root_888.
  needs_pass2 <- empty_like(bc888$root_888) & empty_like(bc888$supervoxel_id)
  has_position <- !empty_like(bc888$position)
  pass2_idx <- which(needs_pass2 & has_position)
  message(sprintf("[%s] pass 2: resolving %s svids from position...",
                  Sys.time(), format(length(pass2_idx), big.mark = ",")))
  if (length(pass2_idx)) {
    bc888$supervoxel_id[pass2_idx] <- bancr::banc_xyz2id(
      bc888$position[pass2_idx], root = FALSE, rawcoords = TRUE)
    refresh <- pass2_idx[!empty_like(bc888$supervoxel_id[pass2_idx])]
    if (length(refresh)) {
      bc888$root_888[refresh] <- bancr::banc_rootid(
        bc888$supervoxel_id[refresh], version = "888")
    }
  }

  # ---- Pass 3: last-ditch fallback via root_id -> banc_latestid(version=888).
  # Catches rows with no position AND no svid (or where position-based lookup
  # returned 0 because the segment isn't in v888). banc_latestid samples leaves
  # off the current root_id and walks them back to the v888 root.
  still_missing <- empty_like(bc888$root_888)
  has_root_id   <- !empty_like(bc888$root_id)
  pass3_idx <- which(still_missing & has_root_id)
  message(sprintf("[%s] pass 3: %s rows still missing root_888, attempting root_id fallback via banc_latestid()...",
                  Sys.time(), format(length(pass3_idx), big.mark = ",")))
  if (length(pass3_idx)) {
    bc888$root_888[pass3_idx] <- tryCatch(
      bancr::banc_latestid(bc888$root_id[pass3_idx], version = "888"),
      error = function(e) {
        warning(sprintf("banc_latestid fallback failed: %s", e$message))
        rep(NA_character_, length(pass3_idx))
      })
  }

  bc888[is.na(bc888)] <- ""

  n_resolved <- sum(nzchar(bc888$root_888) & bc888$root_888 != "0")
  n_stuck    <- nrow(bc888) - n_resolved
  message(sprintf("  resolved root_888 for %s/%s rows (%s still unresolved)",
                  format(n_resolved,  big.mark = ","),
                  format(nrow(bc888), big.mark = ","),
                  format(n_stuck,     big.mark = ",")))

  message(sprintf("[%s] writing back to SeaTable in chunks of 1000...",
                  Sys.time()))
  banctable_update_rows(
    base           = "banc_meta",
    table          = "banc_meta",
    df             = bc888[, c("_id", "root_id", "root_888", "supervoxel_id")],
    append_allowed = FALSE,
    chunksize      = 1000
  )
  message(sprintf("[%s] done.", Sys.time()))
}
