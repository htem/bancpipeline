#' banc-update-seatable — Master SeaTable integration orchestrator.
#'
#' Top-level driver: (1) update root_ids in SeaTable, (2) save a snapshot,
#' (3) merge identity data from `banc_ids.csv`, (4) update
#' region/cell_class/nerve/nucleus from `banc_meta.csv`, (5) source
#' update-matches / update-celltypes / update-status, (6) append new rows
#' to SeaTable, (7) detect concurrent changes.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/banc_ids.csv`, `banc_meta.csv`
#'   - SeaTable `banc_meta`
#'   - sourced child scripts: `banc-update-matches.R`,
#'     `banc-update-celltypes.R`, `banc-update-status.R`
#'
#' @section Writes:
#'   - `<banc.meta.save.path>/snapshots/<ts>_banc_seatable.csv` — pre-push snapshot
#'   - `logs/banc-update-seatable_<date>.txt`
#'   - SeaTable `banc_meta` — many columns (cascaded through children)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`

###########################################################
### Master seatable integration: update IDs, snapshot,
### merge identity data, and append new rows
###
### This script handles:
###   1. Update root IDs in seatable
###   2. Save snapshot
###   3. Read + merge identity data from banc_ids.csv
###   4. Update region/cell_class/nerve/nucleus from
###      banc_meta.csv
###   5. Source update-matches, update-celltypes,
###      update-status
###   6. Append new rows to seatable
###   7. Concurrent-change detection
###########################################################
source("banc/banc-startup.R")

# Set up logging
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
.log_file <- file.path(log_dir, sprintf("banc-update-seatable_%s.txt", format(Sys.time(), "%Y-%m-%d")))
.log_con <- file(.log_file, open = "wt")
sink(.log_con, type = "output", split = TRUE)
sink(.log_con, type = "message")
message(sprintf("Log started: %s", Sys.time()))

local({

message("### banc: seatable integration ###")
t_start <- Sys.time()

# SeaTable cloud API has a hardcoded 30s read timeout per chunk; partial
# bulk writes can stall and crash the whole pipeline. Wrap each SeaTable-
# mutating call in retry-with-exponential-backoff. Idempotent operations
# (banctable_updateids re-resolves from CAVE; banctable_update_rows skips
# unchanged rows; banctable_append_rows guards on supervoxel_id+root_id in
# our pre-filter) so re-running converges.
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
      wait <- min(300, base_wait * 2^(.try - 1))   # 10,20,40,80,160,300 cap
      message(sprintf("  [%s] sleeping %ds before retry", label, wait))
      Sys.sleep(wait)
    }
  }
  stop(sprintf("[%s] exhausted %d retries", label, max_tries), call. = FALSE)
}

###########################
### Update IDs          ###
###########################

.seatable_retry(bancr:::banctable_updateids(), label = "banctable_updateids")

###########################
### Save snapshot       ###
###########################

.seatable_retry(bc.orig <- banctable_query(), label = "banctable_query(bc.orig)")
dir.create(file.path(banc.meta.save.path, "snapshots"), showWarnings = FALSE)
datetime_string <- format(Sys.time(), "%Y-%m-%d_%H-%M")
snapshot_file <- file.path(banc.meta.save.path, "snapshots",
                           paste0(datetime_string, "_banc_seatable.csv"))
readr::write_csv(bc.orig, snapshot_file)
message(sprintf("Snapshot saved: %s (%d rows)", basename(snapshot_file), nrow(bc.orig)))

###########################
### Read identity data  ###
###########################

drop.cols <- c("_locked", "_locked_by", "_archived", "_creator",
               "_ctime", "_last_modifier", "_mtime")

.seatable_retry(bc1 <- banctable_query(), label = "banctable_query(bc1)")
bc1[, drop.cols] <- NULL

# Read banc_ids.csv (from banc-ids.R)
banc.ids <- readr::read_csv(file.path(banc.meta.save.path, "banc_ids.csv"),
                            col_types = banc.col.types, 
                            show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Merge identity data
# Pull nucleus columns from banc.ids too, so newly-appended rows get nucleus
# info (banc.ids carries it from CAVE banc_nuclei via banc-ids.R). Without
# this, new rows landed in seatable with blank nucleus_id and had to wait
# for the NEXT cycle to converge — see drift_reports/
# banc_meta_appended_missing_nucleus_*.csv for affected rows.
.nucleus_cols_from_banc_ids <- intersect(
  c("nucleus_id", "nucleus_supervoxel_id",
    "nucleus_position", "nucleus_position_nm"),
  colnames(banc.ids))
bc <- bc1 %>%
  dplyr::rename(root_id.x = root_id, position.x = position) %>%
  dplyr::full_join(banc.ids %>%
                     dplyr::select(dplyr::all_of(c("root_id", "supervoxel_id", "position",
                                                    .nucleus_cols_from_banc_ids))) %>%
                     dplyr::rename_with(~ paste0(.x, ".bi"),
                                        .cols = dplyr::all_of(.nucleus_cols_from_banc_ids)) %>%
                     dplyr::rename(root_id.y = root_id, position.y = position),
                   by = "supervoxel_id") %>%
  dplyr::mutate(
    root_id = dplyr::coalesce(root_id.y, root_id.x),
    position = dplyr::coalesce(position.y, position.x)
  )
# For each nucleus column: prefer banctable_query() value (existing rows have
# what they have), fall back to banc.ids value (fills the new rows). Don't
# overwrite curated SeaTable values with stale banc.ids ones.
for (.nc in .nucleus_cols_from_banc_ids) {
  .bi <- paste0(.nc, ".bi")
  if (.nc %in% colnames(bc) && .bi %in% colnames(bc)) {
    bc[[.nc]] <- ifelse(is.na(bc[[.nc]]) | bc[[.nc]] == "",
                        bc[[.bi]], bc[[.nc]])
    bc[[.bi]] <- NULL
  } else if (.bi %in% colnames(bc)) {
    bc[[.nc]] <- bc[[.bi]]
    bc[[.bi]] <- NULL
  }
}
bc <- bc %>%
  dplyr::select(-dplyr::any_of(c("root_id.x", "root_id.y", "position.x", "position.y"))) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
  dplyr::distinct(position, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(root_id), !is.na(supervoxel_id))

###########################
### Merge banc_meta     ###
###########################

# Read BANC meta local (computed by banc-meta.R)
banc.meta <- readr::read_csv(file.path(banc.meta.save.path, "banc_meta.csv"),
                              col_types = banc.col.types, show_col_types = FALSE) %>%
  dplyr::arrange(cell_type, dplyr::desc(fafb_nblast), dplyr::desc(manc_nblast),
                 dplyr::desc(malecns_nblast), dplyr::desc(banc_nblast)) %>%
  dplyr::filter(!is.na(root_id), !is.na(supervoxel_id), root_id != "0")
banc.meta <- banc_updateids(banc.meta)
readr::write_csv(banc.meta, file.path(banc.meta.save.path, "banc_meta.csv"))
banc.meta <- banc.meta[!duplicated(banc.meta$root_id), ]

# Metric columns to auto-update from banc_meta.
# `proofread` and `roughly_proofread` are intentionally absent — banc-ids.R
# is the sole writer for those flags (it has direct access to the CAVE
# backbone_proofread + proofreading_notes tables and applies a preserve-on-
# uncertainty policy). Including them here previously caused the
# wipe-and-rejoin at L167-168 to lose curator-asserted TRUEs whenever a
# root_id was missing from banc_meta (e.g. no L2 metrics CSV yet).
metric.cols <- c("l2_nodes", "l2_cable_length_um",
                 "input_connections", "output_connections",
                 "input_side_index", "output_side_index",
                 "banc_nblast_match", "banc_nblast_match_supervoxel_id", "banc_nblast",
                 "fafb_nblast_match", "fafb_nblast",
                 "manc_nblast_match", "manc_nblast",
                 "fanc_nblast_match", "fanc_nblast",
                 "hemibrain_nblast_match", "hemibrain_nblast",
                 "malecns_nblast_match", "malecns_nblast")
metric.cols <- intersect(metric.cols, colnames(banc.meta))

# Wipe auto-update columns and re-join from banc_meta
bc[, colnames(bc) %in% metric.cols] <- NULL
bc <- dplyr::full_join(bc, banc.meta[, c("root_id", metric.cols)], by = "root_id")
bc[bc == "NA"] <- NA
bc[bc == ""] <- NA
bc[bc == " "] <- NA
bc[] <- lapply(bc, function(x) {
  if (is.character(x)) x[grepl("^auto", x, ignore.case = TRUE)] <- NA
  x
})
bc$status[bc$status == "NA"] <- ""

# Update region/cell_class/nerve from banc_meta (vectorized)
message("Updating identity columns from banc_meta...")
empty_values <- c(NA, "NA", "", " ", "unknown", "issue", "ISSUE", "ISSUE!")
update.cols <- c("region", "cell_class", "nerve")
for (update.col in update.cols) {
  if (!update.col %in% colnames(banc.meta)) next
  message(sprintf("  Processing column: %s", update.col))

  meta_lookup <- banc.meta %>%
    dplyr::filter(!is.na(.data[[update.col]]) & .data[[update.col]] != "") %>%
    dplyr::distinct(root_id, .keep_all = TRUE) %>%
    dplyr::select(root_id, .meta_value = dplyr::all_of(update.col))

  bc <- bc %>%
    dplyr::left_join(meta_lookup, by = "root_id") %>%
    dplyr::mutate(
      !!update.col := dplyr::if_else(
        is.na(.data[[update.col]]) | .data[[update.col]] %in% empty_values,
        .meta_value,
        .data[[update.col]]
      )
    ) %>%
    dplyr::select(-.meta_value)
}

###########################
### Source update scripts ###
###########################

message("### banc: sourcing update-matches.R ###")
source("banc/update/banc-update-matches.R")

message("### banc: sourcing update-celltypes.R ###")
source("banc/update/banc-update-celltypes.R")

message("### banc: sourcing update-status.R ###")
source("banc/update/banc-update-status.R")

###########################
### Backfill super_class from CAVE cell_info ###
###########################

# For neurons with blank super_class in seatable, fill with glia/trachea
# if the CAVE cell_info table classifies them as such (via banc-ids.R cave_super_class)
if ("cave_super_class" %in% colnames(banc.ids)) {
  cave_sc_lookup <- banc.ids %>%
    dplyr::filter(!is.na(cave_super_class), cave_super_class != "",
                  grepl("glia|trachea", cave_super_class, ignore.case = TRUE)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE) %>%
    dplyr::select(root_id, cave_super_class)

  sc_blank <- is.na(bc$super_class) | bc$super_class == ""
  sc_fill <- dplyr::left_join(
    bc[sc_blank, c("root_id"), drop = FALSE],
    cave_sc_lookup, by = "root_id"
  )
  filled <- !is.na(sc_fill$cave_super_class)
  if (sum(filled) > 0) {
    bc$super_class[sc_blank][filled] <- sc_fill$cave_super_class[filled]
    message(sprintf("  Backfilled super_class from CAVE cell_info: %d neurons (glia/trachea)", sum(filled)))
  }
}

###########################
### Push identity data  ###
###########################

# Build update dataframe
bc.update <- bc %>%
  dplyr::filter(root_id != "0") %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.character), ~tidyr::replace_na(., "")),
                dplyr::across(dplyr::where(is.factor), ~as.factor(tidyr::replace_na(as.character(.), ""))),
                dplyr::across(dplyr::where(is.logical), ~.)) %>%
  dplyr::select(
    `_id`, root_id, position, supervoxel_id,
    nucleus_id, nucleus_position, nucleus_position_nm, nucleus_supervoxel_id,
    super_class
  )
# Dedupe extant rows by `_id`; dedupe new rows (no `_id`) by root_id so they
# don't all collapse on NA `_id`.
.has_id <- !is.na(bc.update$`_id`) & bc.update$`_id` != ""
bc.update <- dplyr::bind_rows(
  bc.update[.has_id, ] %>% dplyr::distinct(`_id`, .keep_all = TRUE),
  bc.update[!.has_id, ] %>% dplyr::distinct(root_id, .keep_all = TRUE)
)

###########################
### Concurrent changes  ###
###########################

# Check if anything changed in the meantime
.seatable_retry(bc2 <- banctable_query(), label = "banctable_query(bc2)")
bc2[, drop.cols] <- NULL
updated.cols <- intersect(colnames(bc2), colnames(bc.update))
bc2.check <- bc2[, updated.cols]
bc2.check <- bc2.check[, intersect(colnames(bc2.check), colnames(bc1))]
bckeep <- bc2.check %>%
  dplyr::semi_join(bc1, by = names(bc2.check))
bckeep <- bckeep$`_id`

# Update rows that haven't changed since we started
bc.update.present <- bc.update %>%
  dplyr::filter(`_id` %in% bckeep, `_id` != "") %>%
  as.data.frame()
# Replace NAs per column type: text → "", numeric → keep NA (serialized as JSON null)
for (col in colnames(bc.update.present)) {
  if (is.numeric(bc.update.present[[col]])) next
  bc.update.present[[col]][is.na(bc.update.present[[col]])] <- ""
}
bc.update.present <- bc.update.present %>%
  dplyr::anti_join(bc2, by = updated.cols)

if (nrow(bc.update.present)) {
  message(sprintf("Modifying %d extant neurons in banctable", nrow(bc.update.present)))
  .seatable_retry(
    banctable_update_rows(base = 'banc_meta',
                          table = "banc_meta",
                          df = bc.update.present,
                          append_allowed = FALSE,
                          chunksize = 1000),
    label = "banctable_update_rows(bc.update.present)")
}

###########################
### Append new rows     ###
###########################

new.rows <- bc.update %>%
  dplyr::filter(is.na(`_id`) | `_id` == "",
                !root_id %in% bc.update.present$root_id,
                !supervoxel_id %in% bc.update.present$supervoxel_id,
                !root_id %in% bc2$root_id,
                !supervoxel_id %in% bc2$supervoxel_id) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  as.data.frame()
new.rows$`_id` <- NULL
new.rows[is.na(new.rows)] <- ""
new.rows$status <- NULL

if (nrow(new.rows)) {
  message(sprintf("Adding %d new neurons to banctable", nrow(new.rows)))
  .seatable_retry(
    banctable_append_rows(base = 'banc_meta',
                          table = "banc_meta",
                          df = new.rows,
                          chunksize = 1000),
    label = "banctable_append_rows(new.rows)")
}

###########################
## DEFENSIVE CHECK: any newly-appended row that lacks nucleus_id but for
## which CAVE has a nucleus assigned to that root_id. Logs a "second-cycle-
## needed" warning and writes a CSV — typically caused by banc-meta.R missing
## a row (no L2 metrics yet) or banc-update-seatable.R dropping nucleus
## columns from banc.ids during the full_join. Report-only.
###########################
try({
  blank_nuc_new <- new.rows %>%
    dplyr::filter(is.na(nucleus_id) | nucleus_id == "")
  if (nrow(blank_nuc_new) > 0L) {
    cave_nuc_check <- bancr::banc_nuclei(rawcoords = TRUE, table = "both") %>%
      dplyr::transmute(cave_nucleus_id = as.character(id),
                       root_id         = as.character(pt_root_id)) %>%
      dplyr::filter(!is.na(cave_nucleus_id), cave_nucleus_id != "",
                    !is.na(root_id),         root_id        != "0")
    needs_second_cycle <- blank_nuc_new %>%
      dplyr::inner_join(cave_nuc_check, by = "root_id")
    if (nrow(needs_second_cycle) > 0L) {
      out_dir <- file.path(banc.meta.save.path, "drift_reports")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      stamp <- format(Sys.time(), "%Y-%m-%d")
      out <- file.path(out_dir,
                       sprintf("banc_meta_appended_missing_nucleus_%s.csv", stamp))
      readr::write_csv(needs_second_cycle, out)
      message(sprintf("  WARN: %d newly-appended rows lack nucleus_id but CAVE has one assigned. Will resolve next cycle. -> %s",
                      nrow(needs_second_cycle), out))
    } else {
      message(sprintf("  OK: %d new rows had blank nucleus_id but CAVE confirmed none for those root_ids",
                      nrow(blank_nuc_new)))
    }
  }
}, silent = FALSE)

# Save fi
readr::write_csv(bc.update, file.path(banc.meta.save.path, "banc_seatable.csv"))
message(sprintf("### banc: seatable integration complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})

# Close logging
message(sprintf("Log ended: %s", Sys.time()))
sink(type = "message")
sink(type = "output")
close(.log_con)
message(sprintf("Log saved: %s", .log_file))
