#' banc-ids890 — Standalone root_890 backfill for SeaTable banc_meta.
#'
#' Cloned from `banc-ids888.R` for the v888 → v890 migration. Runs
#' independently of `banc-ids.R`. Resolves `root_890` from `supervoxel_id`
#' via CAVE and pushes back to SeaTable. Idempotent. Includes a clear
#' fallback if the `root_890` SeaTable column does not exist yet: aborts
#' before any write with a "create the column" message.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - CAVE `banc_rootid()` at v890
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `root_id`, `root_890`, `supervoxel_id`
#'
#' @section Invoked by:
#'   `o2/oneshots/o2_banc_v890_rebuild.sh`

###############################################################################
### Standalone root_890 backfill for SeaTable banc_meta.
###
### Cloned from banc/update/banc-ids888.R for the v888 → v890 migration.
### Runs independently of banc-ids.R so we can backfill root_890 even when
### the rest of banc-ids.R is blocked by upstream issues.
###
### Reads SeaTable banc_meta, finds rows with missing/empty/zero root_890,
### resolves their supervoxel_id at v890 via CAVE banc_rootid(), and writes
### root_id/root_890/supervoxel_id back via banctable_update_rows().
###
### Idempotent: re-running on a fully populated table is a no-op.
###
### Fallback for the case where the root_890 SeaTable column does not exist
### yet at run time: the script catches the column-missing error, queries the
### subset of rows it CAN see (without root_890), resolves root_890 via CAVE
### in-memory, and aborts with a clear "create the root_890 column" message
### so the user can add the column and re-run. We never write to a column we
### can't read, so this aborts before any write.
###############################################################################
source("banc/banc-startup.R")

bancr::choose_banc()

target.version <- "890"
root.col <- paste0("root_", target.version)

message(sprintf("[%s] querying banc_meta for rows missing %s...",
                Sys.time(), root.col))

sql.full <- sprintf(
  "SELECT _id, status, root_id, %s, supervoxel_id, position from banc_meta",
  root.col
)

bc.all <- tryCatch(
  banctable_query(sql.full),
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl(root.col, msg, fixed = TRUE)) {
      message(sprintf(
        "ERROR: SeaTable rejected SELECT including %s — column likely doesn't exist yet.\n  Underlying message: %s\n  To proceed: add a %s text column to the banc_meta SeaTable, then re-run this script.",
        root.col, msg, root.col))
      stop(sprintf("Add %s column to banc_meta SeaTable before running banc-ids890.R.",
                   root.col), call. = FALSE)
    }
    stop(e)
  }
)
message(sprintf("  banc_meta total rows: %s",
                format(nrow(bc.all), big.mark = ",")))

bc890 <- bc.all %>%
  dplyr::select(`_id`, !!rlang::sym(root.col), root_id, supervoxel_id, position) %>%
  dplyr::filter(
    is.na(.data[[root.col]]) | .data[[root.col]] == "" |
    .data[[root.col]] == "0" | .data[[root.col]] == 0 |
    is.na(supervoxel_id) | supervoxel_id == "0" | supervoxel_id == 0
  )

message(sprintf("  rows needing %s backfill: %s",
                root.col, format(nrow(bc890), big.mark = ",")))

if (nrow(bc890) == 0L) {
  message(sprintf("nothing to update -- %s already populated for all rows",
                  root.col))
} else {
  # Resolve supervoxel_id from position when missing, then root_890 from svid.
  message(sprintf("[%s] resolving supervoxel_id from position...", Sys.time()))
  bc890$supervoxel_id <- bancr::banc_xyz2id(bc890$position,
                                            root = FALSE,
                                            rawcoords = TRUE)
  message(sprintf("[%s] resolving %s via CAVE...", Sys.time(), root.col))
  bc890[[root.col]] <- bancr::banc_rootid(bc890$supervoxel_id,
                                          version = target.version)

  bc890[is.na(bc890)] <- ""

  n_resolved <- sum(nzchar(bc890[[root.col]]) & bc890[[root.col]] != "0")
  message(sprintf("  resolved %s for %s/%s rows",
                  root.col,
                  format(n_resolved,  big.mark = ","),
                  format(nrow(bc890), big.mark = ",")))

  message(sprintf("[%s] writing back to SeaTable in chunks of 1000...",
                  Sys.time()))
  banctable_update_rows(
    base           = "banc_meta",
    table          = "banc_meta",
    df             = bc890[, c("_id", "root_id", root.col, "supervoxel_id")],
    append_allowed = FALSE,
    chunksize      = 1000
  )
  message(sprintf("[%s] done.", Sys.time()))
}
