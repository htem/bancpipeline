#' banc-update-matches — Push NBLAST + PNG match columns to SeaTable.
#'
#' Reads PNG review CSVs (mirror, FAFB, MANC, hemibrain, FANC, maleCNS) and
#' NBLAST feather files, combines matches, detects conflicts, and pushes
#' the per-dataset `*_match` / `*_png_match` / `*_nblast_match` / `*_nblast`
#' columns to SeaTable independently.
#'
#' @section Reads:
#'   - `<banc.nblast.{dataset}.save.path>/correct/*.csv` (PNG matches per dataset)
#'   - `<banc.save.path>/banc_<target>_<ver>_nblast.feather` (per-dataset NBLAST)
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: per-dataset `*_match`, `*_png_match`,
#'     `*_nblast_match`, `*_nblast` columns
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`

###########################################################
### Update NBLAST and PNG match columns in BANC seatable
###
### Reads PNG review CSVs (mirror, FAFB, MANC, hemibrain,
### FANC, maleCNS) and NBLAST feather files.
### Combines matches, detects conflicts.
###
### Pushes match columns to seatable independently:
###   *_match, *_png_match, *_nblast_match, *_nblast
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: updating match columns ###")
t_start <- Sys.time()

###########################
### Read current state  ###
###########################

bc <- banctable_query()
bc <- bc %>%
  dplyr::filter(!is.na(root_id), root_id != "0") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Save original match values for anti-join (only push rows that actually changed)
match_orig_cols <- c("_id",
                     "banc_match", "banc_match_supervoxel_id",
                     "banc_nblast_match", "banc_nblast_match_supervoxel_id", "banc_nblast",
                     "fafb_match", "fafb_nblast_match", "fafb_nblast",
                     "manc_match", "manc_nblast_match", "manc_nblast",
                     "hemibrain_match", "hemibrain_nblast_match", "hemibrain_nblast",
                     "fanc_match", "fanc_nblast_match", "fanc_nblast",
                     "malecns_match", "malecns_nblast_match", "malecns_nblast")
match_orig_cols <- intersect(match_orig_cols, colnames(bc))
bc.match.orig <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(dplyr::all_of(match_orig_cols))

###########################
### Read match data     ###
###########################

# PNG review CSVs
mirror.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_mirror_reviewed_matches.csv"),
                                      col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id,
                banc_png_match = match_root_id,
                banc_png_match_supervoxel_id = match_supervoxel_id)

fafb.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_fafb_reviewed_matches.csv"),
                                    col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id, fafb_png_match = match_id)

manc.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_manc_reviewed_matches.csv"),
                                    col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id, manc_png_match = match_id)

hemibrain.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_hemibrain_reviewed_matches.csv"),
                                         col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id, hemibrain_png_match = match_id)

fanc.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_fanc_reviewed_matches.csv"),
                                    col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id, fanc_png_match = match_id)

malecns.matches.df <- readr::read_csv(file.path(banc.meta.save.path, "banc_malecns_reviewed_matches.csv"),
                                       col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id), valid == 't') %>%
  dplyr::rename(root_id = pt_root_id, malecns_png_match = match_id)

# NBLAST feather files — version strings for cross-dataset NBLAST results
nblast.versions <- list(
  fafb = "783",
  manc = "v1.2.1",
  hemibrain = "v1.2.1",
  fanc = "1116",
  malecns = "v0.9"
)

banc.meta.fafb.nb <- arrow::read_feather(
  file.path(banc.meta.save.path, sprintf("banc_fafb_%s_nblast.feather", nblast.versions$fafb))) %>%
  dplyr::filter(validation == TRUE, !grepl("\\_", match_id)) %>%
  dplyr::rename(root_id = pt_root_id, fafb_match = match_id)

banc.meta.manc.nb <- arrow::read_feather(
  file.path(banc.meta.save.path, sprintf("banc_manc_%s_nblast.feather", nblast.versions$manc))) %>%
  dplyr::filter(validation == TRUE, !grepl("\\_", match_id)) %>%
  dplyr::rename(root_id = pt_root_id, manc_match = match_id)

banc.meta.hemibrain.nb <- arrow::read_feather(
  file.path(banc.meta.save.path, sprintf("banc_hemibrain_%s_nblast.feather", nblast.versions$hemibrain))) %>%
  dplyr::filter(validation == TRUE, !grepl("\\_", match_id)) %>%
  dplyr::rename(root_id = pt_root_id, hemibrain_match = match_id)

banc.meta.fanc.nb <- arrow::read_feather(
  file.path(banc.meta.save.path, sprintf("banc_fanc_%s_nblast.feather", nblast.versions$fanc))) %>%
  dplyr::filter(validation == TRUE, !grepl("\\_", match_id)) %>%
  dplyr::rename(root_id = pt_root_id, fanc_match = match_id)

banc.meta.malecns.nb <- arrow::read_feather(
  file.path(banc.meta.save.path, sprintf("banc_malecns_%s_nblast.feather", nblast.versions$malecns))) %>%
  dplyr::filter(validation == TRUE, !grepl("\\_", match_id)) %>%
  dplyr::rename(root_id = pt_root_id, malecns_match = match_id)

banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_mirror_nblast.feather")) %>%
  dplyr::filter(valid == 't') %>%
  dplyr::rename(root_id = pt_root_id,
                banc_match = match_root_id,
                banc_match_supervoxel_id = match_supervoxel_id)

# Read NBLAST scores from banc_meta.csv
banc.meta <- readr::read_csv(file.path(banc.meta.save.path, "banc_meta.csv"),
                              col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Combine all match data
matches.df <- plyr::rbind.fill(mirror.matches.df, fafb.matches.df, manc.matches.df,
                                fanc.matches.df, hemibrain.matches.df, malecns.matches.df,
                                banc.meta.fafb.nb, banc.meta.manc.nb, banc.meta.fanc.nb,
                                banc.meta.hemibrain.nb, banc.meta.malecns.nb,
                                banc.meta.mirror.nb)

###########################
### Update match cols   ###
###########################

update.cols <- c("banc_match", "banc_match_supervoxel_id",
                 "fafb_match", "manc_match", "hemibrain_match", "fanc_match", "malecns_match",
                 "banc_png_match", "banc_png_match_supervoxel_id",
                 "fafb_png_match", "manc_png_match", "hemibrain_png_match",
                 "fanc_png_match", "malecns_png_match")

# Initialize PNG match columns
for (col in c("banc_png_match", "banc_png_match_supervoxel_id",
              "fafb_png_match", "manc_png_match", "hemibrain_png_match",
              "fanc_png_match", "malecns_png_match")) {
  bc[[col]] <- ""
}

# NBLAST score columns from banc_meta
nblast_cols <- c("banc_nblast_match", "banc_nblast_match_supervoxel_id", "banc_nblast",
                 "fafb_nblast_match", "fafb_nblast",
                 "manc_nblast_match", "manc_nblast",
                 "fanc_nblast_match", "fanc_nblast",
                 "hemibrain_nblast_match", "hemibrain_nblast",
                 "malecns_nblast_match", "malecns_nblast")
nblast_cols_available <- intersect(nblast_cols, colnames(banc.meta))
if (length(nblast_cols_available) > 0) {
  bc <- bc %>%
    dplyr::select(-dplyr::any_of(nblast_cols_available)) %>%
    dplyr::left_join(banc.meta[, c("root_id", nblast_cols_available)], by = "root_id")
}

message(sprintf("Updating %d match columns for %d neurons...", length(update.cols), nrow(bc)))
empty_values <- c(NA, "NA", "", " ", "unknown", "issue", "ISSUE", "ISSUE!")
conflicts <- data.frame()
for (update.col in update.cols) {
  if (!update.col %in% colnames(matches.df)) next
  message(sprintf("  Processing column: %s", update.col))

  # Get best new value per root_id from matches data
  new_vals <- matches.df %>%
    dplyr::filter(!is.na(.data[[update.col]]), .data[[update.col]] != "") %>%
    dplyr::distinct(root_id, .keep_all = TRUE) %>%
    dplyr::select(root_id, new_value = dplyr::all_of(update.col))

  # Join new values to bc
  merged <- bc %>%
    dplyr::select(root_id, `_id`, old_value = dplyr::all_of(update.col)) %>%
    dplyr::left_join(new_vals, by = "root_id")

  # Identify rows where we can update (old is empty, new is not)
  can_update <- !is.na(merged$new_value) & merged$new_value != "" &
                (is.na(merged$old_value) | merged$old_value %in% empty_values)
  bc[[update.col]][can_update] <- merged$new_value[can_update]

  # Identify conflicts (old has a real value, new differs)
  has_conflict <- !is.na(merged$new_value) & merged$new_value != "" &
                  !is.na(merged$old_value) & !merged$old_value %in% empty_values &
                  merged$old_value != merged$new_value
  if (any(has_conflict)) {
    conflict_rows <- merged[has_conflict, ]
    for (j in seq_len(nrow(conflict_rows))) {
      message(sprintf("conflict %s: current %s vs suggested %s for %s",
                      update.col, conflict_rows$old_value[j],
                      conflict_rows$new_value[j], conflict_rows$root_id[j]))
    }
    conflicts <- rbind(conflicts, data.frame(
      `_id` = conflict_rows$`_id`, root_id = conflict_rows$root_id,
      type = update.col, current = conflict_rows$old_value,
      suggested = conflict_rows$new_value))
  }
}

# Merge PNG matches into *_match where *_match is empty
bc <- bc %>%
  dplyr::mutate(
    fafb_match = dplyr::case_when(
      is.na(fafb_match) & !fafb_png_match %in% c("no_match", "", "0") ~ fafb_png_match,
      TRUE ~ fafb_match),
    manc_match = dplyr::case_when(
      is.na(manc_match) & !manc_png_match %in% c("no_match", "", "0") ~ manc_png_match,
      TRUE ~ manc_match),
    hemibrain_match = dplyr::case_when(
      is.na(hemibrain_match) & !hemibrain_png_match %in% c("no_match", "", "0") ~ hemibrain_png_match,
      TRUE ~ hemibrain_match),
    fanc_match = dplyr::case_when(
      is.na(fanc_match) & !fanc_png_match %in% c("no_match", "", "0") ~ fanc_png_match,
      TRUE ~ fanc_match),
    banc_match = dplyr::case_when(
      is.na(banc_match) & !banc_png_match %in% c("no_match", "", "0") ~ banc_png_match,
      TRUE ~ banc_match),
    banc_match_supervoxel_id = dplyr::case_when(
      is.na(banc_match_supervoxel_id) & !banc_png_match_supervoxel_id %in% c("no_match", "", "0") ~ banc_png_match_supervoxel_id,
      TRUE ~ banc_match_supervoxel_id)
  )

###########################
### Push to seatable    ###
###########################

# Select only match columns + _id for push
match_push_cols <- c("_id",
                     "banc_match", "banc_match_supervoxel_id",
                     "banc_png_match", "banc_png_match_supervoxel_id",
                     "banc_nblast_match", "banc_nblast_match_supervoxel_id", "banc_nblast",
                     "fafb_match", "fafb_png_match", "fafb_nblast_match", "fafb_nblast",
                     "manc_match", "manc_png_match", "manc_nblast_match", "manc_nblast",
                     "hemibrain_match", "hemibrain_png_match", "hemibrain_nblast_match", "hemibrain_nblast",
                     "fanc_match", "fanc_png_match", "fanc_nblast_match", "fanc_nblast",
                     "malecns_match", "malecns_png_match", "malecns_nblast_match", "malecns_nblast")
match_push_cols <- intersect(match_push_cols, colnames(bc))

bc.match.update <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(dplyr::all_of(match_push_cols)) %>%
  as.data.frame()

# Convert integer64 columns to character — SeaTable expects text for ID columns,
# and integer64 values lose precision when passed through reticulate as doubles
int64_cols <- names(bc.match.update)[sapply(bc.match.update, bit64::is.integer64)]
if (length(int64_cols) > 0) {
  message(sprintf("  Converting %d integer64 columns to character: %s",
                  length(int64_cols), paste(int64_cols, collapse = ", ")))
  for (col in int64_cols) {
    bc.match.update[[col]] <- as.character(bc.match.update[[col]])
  }
}

# Look up SeaTable column types to handle numeric vs text columns correctly
col_info <- tryCatch(
  bancr:::banctable_columns(table = "banc_meta"),
  error = function(e) {
    warning("Could not fetch column types from SeaTable: ", e$message,
            " — falling back to regex-based detection")
    NULL
  }
)
if (!is.null(col_info)) {
  numeric_cols <- col_info$name[col_info$type == "number"]
  numeric_cols <- intersect(numeric_cols, colnames(bc.match.update))
} else {
  numeric_cols <- grep("_nblast$", colnames(bc.match.update), value = TRUE)
}

# Sanitize numeric columns: ensure finite values, round scores
for (col in numeric_cols) {
  vals <- as.numeric(bc.match.update[[col]])
  vals[!is.finite(vals)] <- NA
  bc.match.update[[col]] <- round(vals, 3)
}

# Replace text NAs with empty string (clears SeaTable text cells)
for (col in colnames(bc.match.update)) {
  if (col %in% numeric_cols) next
  bc.match.update[[col]][is.na(bc.match.update[[col]])] <- ""
}

# Diagnose column keys for debugging SeaTable errors
if (!is.null(col_info)) {
  message("  Number-type columns being pushed: ",
          paste(numeric_cols, collapse = ", "))
  # Show SeaTable types for all columns being pushed
  pushed_types <- col_info[col_info$name %in% colnames(bc.match.update), c("name", "type", "key")]
  for (ri in seq_len(nrow(pushed_types))) {
    rtype <- class(bc.match.update[[pushed_types$name[ri]]])[1]
    message(sprintf("  Column '%s' (key=%s): SeaTable=%s, R=%s",
                    pushed_types$name[ri], pushed_types$key[ri],
                    pushed_types$type[ri], rtype))
  }
  # Map column keys to names for error diagnosis
  key_map <- setNames(col_info$name, col_info$key)
  problem_keys <- c("gudM", "1uB4", "Dm7e", "BjFL")
  for (pk in problem_keys) {
    if (pk %in% names(key_map)) {
      cname <- key_map[[pk]]
      message(sprintf("  Column key '%s' = '%s' (type: %s)",
                      pk, cname, col_info$type[col_info$key == pk]))
      if (cname %in% colnames(bc.match.update)) {
        vals <- bc.match.update[[cname]]
        bad <- vals[!vals %in% c("") & suppressWarnings(is.na(as.numeric(vals)))]
        if (length(bad))
          message(sprintf("    Non-numeric values in '%s': %s",
                          cname, paste(head(unique(bad), 10), collapse = ", ")))
      }
    }
  }
}

# Anti-join: only push rows where match columns actually changed
# For anti-join, compare text columns as character and round numeric columns to match push format
bc.match.orig.cmp <- bc.match.orig
for (col in intersect(numeric_cols, colnames(bc.match.orig.cmp))) {
  vals <- as.numeric(bc.match.orig.cmp[[col]])
  vals[!is.finite(vals)] <- NA
  bc.match.orig.cmp[[col]] <- round(vals, 3)
}
for (col in intersect(setdiff(colnames(bc.match.orig.cmp), numeric_cols), colnames(bc.match.update))) {
  bc.match.orig.cmp[[col]][is.na(bc.match.orig.cmp[[col]])] <- ""
  bc.match.orig.cmp[[col]] <- as.character(bc.match.orig.cmp[[col]])
  bc.match.update[[col]] <- as.character(bc.match.update[[col]])
}
shared_cmp_cols <- intersect(colnames(bc.match.update), colnames(bc.match.orig.cmp))
n_before <- nrow(bc.match.update)
bc.match.update <- dplyr::anti_join(bc.match.update, bc.match.orig.cmp, by = shared_cmp_cols)
message(sprintf("Pushing match data for %d changed neurons to seatable (skipped %d unchanged)",
                nrow(bc.match.update), n_before - nrow(bc.match.update)))
if (nrow(bc.match.update)) {
  # Split push: text columns for all rows, numeric columns only for non-NA rows.
  # This avoids numeric NA serialization issues (older bancr may serialize NA as
  # "NA" string which SeaTable rejects for number columns).
  #
  # IMPORTANT: never push _id + only 1 data column — banc_df2updatepayload drops
  # the column name when as.list() receives a scalar, producing a JSON array
  # instead of object, which SeaTable rejects as "failed to parse json".
  text_cols <- setdiff(colnames(bc.match.update), numeric_cols)

  # Retry wrapper for transient SeaTable errors (timeouts, connection resets)
  retry_push <- function(df, max_retries = 3, wait_sec = 15, ...) {
    for (attempt in seq_len(max_retries)) {
      ok <- tryCatch({
        banctable_update_rows(df = df, ...)
        TRUE
      }, error = function(e) {
        if (grepl("Timeout|timed out|ConnectionError|Connection reset|BrokenPipe",
                  e$message, ignore.case = TRUE) && attempt < max_retries) {
          message(sprintf("    Timeout on attempt %d/%d, retrying in %ds...",
                          attempt, max_retries, wait_sec))
          Sys.sleep(wait_sec)
          FALSE
        } else {
          stop(e)
        }
      })
      if (ok) return(invisible(TRUE))
    }
  }

  # Phase 1: push text columns for all rows
  if (length(text_cols) > 1) {
    text_df <- bc.match.update[, text_cols, drop = FALSE]
    message(sprintf("  Phase 1: pushing %d text columns for %d rows",
                    length(text_cols) - 1, nrow(text_df)))
    retry_push(text_df, base = 'banc_meta', table = 'banc_meta',
               append_allowed = FALSE, chunksize = 1000)
  }

  # Phase 2: push numeric columns — batch complete rows, then handle partials
  all_complete <- complete.cases(bc.match.update[, numeric_cols, drop = FALSE])
  if (any(all_complete)) {
    complete_df <- bc.match.update[all_complete, c("_id", numeric_cols), drop = FALSE]
    message(sprintf("  Phase 2a: pushing %d numeric columns for %d complete rows",
                    length(numeric_cols), nrow(complete_df)))
    retry_push(complete_df, base = 'banc_meta', table = 'banc_meta',
               append_allowed = FALSE, chunksize = 1000)
  }

  # Partial rows: group by which numeric columns are non-NA and push together.
  partial <- bc.match.update[!all_complete, , drop = FALSE]
  if (nrow(partial) > 0) {
    na_key <- apply(is.na(partial[, numeric_cols, drop = FALSE]), 1,
                    function(r) paste(which(!r), collapse = "-"))
    for (key in unique(na_key)) {
      if (key == "") next
      rows <- which(na_key == key)
      nonna_idx <- as.integer(strsplit(key, "-")[[1]])
      nonna_cols <- numeric_cols[nonna_idx]
      push_cols <- c("_id", nonna_cols)
      # If only 1 numeric column, include all text columns to avoid
      # single-data-column serialization bug in banc_df2updatepayload
      if (length(nonna_cols) == 1) push_cols <- c(text_cols, nonna_cols)
      group_df <- partial[rows, push_cols, drop = FALSE]
      message(sprintf("  Phase 2b: pushing %s for %d rows",
                      paste(nonna_cols, collapse = "+"), nrow(group_df)))
      retry_push(group_df, base = 'banc_meta', table = 'banc_meta',
                 append_allowed = FALSE, chunksize = 1000)
    }
  }

  message("  Push complete")
}

message(sprintf("### banc: match columns update complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
