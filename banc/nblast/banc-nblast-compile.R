#' banc-nblast-compile — Consolidate per-query NBLAST CSVs into per-dataset feathers.
#'
#' Updates compiled feathers from reviewed CSVs + scanned PNG review
#' folders. Run AFTER `banc-nblast-wrong-matches.R`.
#'
#' @section Reads:
#'   - per-query NBLAST CSVs, reviewed CSVs, PNG review folders, existing feathers
#'
#' @section Writes:
#'   - `banc_<target>_<ver>_nblast.feather` (one per dataset, plus mirror + native)
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   BANC-project/R/figures/panels_proofread_matching.R (loads `banc_{fafb,manc,hemibrain,fanc,malecns}_nblast.feather`);
#'   `banc/share/banc-nblast-share-gcs.R` (publishes to GCS).
#'
#' @section Schema:
#'   banc_fafb_783_nblast.md.
#'
#' @section Paper:
#'   Methods §"NBLAST cross-dataset matching".

###########################################################
### Compile NBLAST results from per-query CSVs into
### combined feather files, incorporating manual PNG
### review results.
###
### For each dataset (maleCNS, mirror, FANC, FAFB,
### MANC, hemibrain):
###   1. Read existing compiled feather + reviewed CSV
###   2. Update validated matches from reviewed CSVs
###   3. Scan PNG review folders for new reviews
###   4. Read per-query NBLAST result CSVs
###   5. Combine into updated feather files
###
### Depends on: banc-nblast-wrong-matches.R having run
### first to clean out invalidated matches.
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: compiling NBLAST results ###")
t_start <- Sys.time()

version <- banc.nblast.version

# Register cores
cl <- setup_parallel()
on.exit(stop_parallel(cl), add = TRUE)

###########################
### Helper functions    ###
###########################

# Convert CAVE-format feather columns to internal names used throughout this script
cave_to_internal <- function(df) {
  if ("query_root_id" %in% names(df)) df <- dplyr::rename(df, query_id = query_root_id)
  if ("validation" %in% names(df)) {
    df$valid <- ifelse(df$validation, 't', 'f')
    df$validation <- NULL
  }
  df
}

# Convert internal column names to CAVE cell_match format for feather output
internal_to_cave <- function(df) {
  if ("query_id" %in% names(df)) df <- dplyr::rename(df, query_root_id = query_id)
  if ("valid" %in% names(df)) {
    df$validation <- df$valid == 't'
    df$valid <- NULL
  }
  df
}

# Safe CSV reader with fallback for oversized files
safe_read_nblast_csv <- function(path, col_types, max_rows = 50000) {
  tryCatch(
    readr::read_csv(path, col_types = col_types, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Error reading %s: %s\n  Retrying with row limit (%d rows)...",
                      basename(path), e$message, max_rows))
      tryCatch(
        readr::read_csv(path, col_types = col_types, show_col_types = FALSE,
                        n_max = max_rows),
        error = function(e2) {
          warning(sprintf("  Failed again for %s: %s — skipping", basename(path), e2$message))
          NULL
        }
      )
    }
  )
}

###########################
### Read existing data  ###
###########################

# Existing compiled feather files
banc.meta.fafb.nb <- cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather")))
banc.meta.manc.nb <- cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather")))
banc.meta.fanc.nb <- cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_fanc_1116_nblast.feather")))
banc.meta.hemibrain.nb <- cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather")))
banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_mirror_nblast.feather"))
banc.meta.malecns.nb <- cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_malecns_v0.9_nblast.feather")))

# Native within-BANC NBLAST (no mirroring) — separate feather. First-run
# bootstrap: if the file doesn't exist yet, start from an empty frame with
# the same schema as banc.meta.mirror.nb.
native_feather_path <- file.path(banc.meta.save.path, "banc_native_nblast.feather")
if (file.exists(native_feather_path)) {
  banc.meta.native.nb <- arrow::read_feather(native_feather_path)
} else {
  banc.meta.native.nb <- banc.meta.mirror.nb[0, ]
}

# Existing reviewed match CSVs
mirror.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_mirror_reviewed_matches.csv"),
                                           col_types = banc.col.types,
                                           show_col_types = FALSE)
fafb.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_fafb_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id))
manc.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_manc_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id))
hemibrain.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_hemibrain_reviewed_matches.csv"),
                                              col_types = banc.col.types,
                                              show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id))
fanc.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_fanc_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id))
malecns.matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_malecns_reviewed_matches.csv"),
                                            col_types = banc.col.types,
                                            show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_", match_id))

###########################
### Update verified     ###
###########################

banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(fafb.matches.df.valid[, c("pt_root_id", "match_id", "valid")],
                   by = c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

banc.meta.manc.nb <- banc.meta.manc.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(manc.matches.df.valid[, c("pt_root_id", "match_id", "valid")],
                   by = c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(hemibrain.matches.df.valid[, c("pt_root_id", "match_id", "valid")],
                   by = c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

banc.meta.fanc.nb <- banc.meta.fanc.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(fanc.matches.df.valid[, c("pt_root_id", "match_id", "valid")],
                   by = c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(malecns.matches.df.valid[, c("pt_root_id", "match_id", "valid")],
                   by = c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

###########################
### Get meta data       ###
###########################

banc.meta <- tryCatch({
  message("  Querying SeaTable for current neuron metadata...")
  bm <- banctable_query()
  bm <- banc_updateids(bm)
  bm
}, error = function(e) {
  warning(sprintf("  banctable_query failed: %s — falling back to banc_ids.csv", e$message))
  bm <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_ids.csv"),
                         col_types = banc.col.types,
                         show_col_types = FALSE)
  bm <- banc_updateids(bm)
  bm
})

# Fast ID update using banc.meta lookup, falling back to banc_updateids()
# for IDs not in banc.meta
fast_updateids <- function(ids, max_retries = 2, delay = 30) {
  ids <- as.character(ids)
  idx <- match(ids, banc.meta$root_id)
  found <- !is.na(idx)
  if (all(found)) return(ids)
  missing <- unique(ids[!found])
  if (length(missing)) {
    updated <- NULL
    for (.attempt in seq_len(max_retries + 1)) {
      updated <- tryCatch(banc_updateids(missing), error = function(e) e)
      if (!inherits(updated, "error")) break
      if (.attempt <= max_retries) {
        warning(sprintf("  fast_updateids: attempt %d/%d failed: %s — retrying in %ds...",
                        .attempt, max_retries + 1, updated$message, delay))
        Sys.sleep(delay)
      } else {
        stop(sprintf("fast_updateids failed after %d attempts: %s",
                     max_retries + 1, updated$message))
      }
    }
    ids[!found] <- updated[match(ids[!found], missing)]
  }
  ids
}

fast_updateids_df <- function(df,
                              root.column = "pt_root_id",
                              supervoxel.column = "pt_supervoxel_id",
                              position.column = "pt_position",
                              max_retries = 2, delay = 30) {
  ids <- as.character(df[[root.column]])
  idx <- match(ids, banc.meta$root_id)
  found <- !is.na(idx)
  # Update found rows from banc.meta (regardless of whether all are found)
  if (any(found)) {
    if (supervoxel.column %in% names(df))
      df[[supervoxel.column]][found] <- banc.meta$supervoxel_id[idx[found]]
    if (position.column %in% names(df))
      df[[position.column]][found] <- banc.meta$position[idx[found]]
  }
  if (all(found)) return(df)
  # Only call banc_updateids for unfound rows (not the entire dataframe)
  unfound_idx <- which(!found)
  message(sprintf("  fast_updateids_df: %d/%d IDs not in banc.meta, updating via CAVE...",
                  length(unfound_idx), nrow(df)))
  unfound_df <- df[unfound_idx, , drop = FALSE]
  for (.attempt in seq_len(max_retries + 1)) {
    result <- tryCatch(
      banc_updateids(unfound_df, root.column = root.column,
                     supervoxel.column = supervoxel.column,
                     position.column = position.column),
      error = function(e) e)
    if (!inherits(result, "error")) {
      df[unfound_idx, ] <- result
      return(df)
    }
    if (.attempt <= max_retries) {
      warning(sprintf("  fast_updateids_df: attempt %d/%d failed: %s — retrying in %ds...",
                      .attempt, max_retries + 1, result$message, delay))
      Sys.sleep(delay)
    } else {
      warning(sprintf("  fast_updateids_df: failed after %d attempts: %s — returning %d rows with stale IDs",
                      max_retries + 1, result$message, length(unfound_idx)))
      return(df)
    }
  }
}

# After fast_updateids_df migrates pt_root_id forward, old NBLAST rows become
# stale: pt_root_id points at the CURRENT neuron but query_id still holds the
# pre-edit root_id, so the score represents similarity between the OLD
# morphology and match_id — not the current morphology. This helper drops
# stale rows (pt_root_id != query_id) ONLY for pt_root_ids that also have
# fresh rows (pt_root_id == query_id). Stale rows are preserved for neurons
# that have never been re-NBLASTed after an edit, so we never lose
# "best we have" data.
#
# Expects the INTERNAL column name query_id (cave_to_internal has been applied
# on read). On-disk feather uses query_root_id.
prune_stale_nblast <- function(df, pt_col = "pt_root_id", query_col = "query_id") {
  if (!all(c(pt_col, query_col) %in% names(df))) {
    warning(sprintf("prune_stale_nblast: missing columns (%s/%s), skipping",
                    pt_col, query_col))
    return(df)
  }
  n0 <- nrow(df)
  is_fresh <- df[[pt_col]] == df[[query_col]]
  # pt_root_ids that have at least one fresh row
  fresh_pts <- unique(df[[pt_col]][is_fresh])
  # Drop stale rows whose pt_root_id is in fresh_pts
  drop_mask <- !is_fresh & (df[[pt_col]] %in% fresh_pts)
  if (any(drop_mask)) {
    df <- df[!drop_mask, , drop = FALSE]
    message(sprintf("  prune_stale_nblast: dropped %d stale rows (%d fresh rows preserved, %d stale-only rows kept as best-we-have)",
                    sum(drop_mask), sum(is_fresh),
                    n0 - sum(is_fresh) - sum(drop_mask)))
  } else {
    message(sprintf("  prune_stale_nblast: no redundant stale rows found (%d fresh, %d stale-only)",
                    sum(is_fresh), n0 - sum(is_fresh)))
  }
  df
}

# Fast supervoxel → root_id lookup: banc.meta first, CAVE fallback for missing.
# Deduplicates supervoxels for efficient batch CAVE calls.
fast_rootid_from_sv <- function(svids, ...) {
  svids_char <- as.character(svids)
  unique_svids <- unique(na.omit(svids_char))
  idx <- match(unique_svids, banc.meta$supervoxel_id)
  found <- !is.na(idx)
  result_map <- setNames(rep(NA_character_, length(unique_svids)), unique_svids)
  result_map[found] <- banc.meta$root_id[idx[found]]
  if (any(!found)) {
    message(sprintf("  fast_rootid_from_sv: %d/%d unique supervoxels not in banc.meta, calling CAVE...",
                    sum(!found), length(unique_svids)))
    cave_result <- banc_rootid(unique_svids[!found], ...)
    result_map[!found] <- cave_result
  }
  unname(result_map[svids_char])
}

fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "flywire_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "manc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "hemibrain_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
hb.ids <- unique(hb.meta$bodyid)
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "fanc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(fc.meta$cell_id)
mcns.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "malecns_09_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
mcns.ids <- unique(mcns.meta$malecns_09_id)

# NBLAST folder paths
fafb.nblast.folder <- file.path(banc.nblast.fafb.save.path, "results", version)
mirror.nblast.folder <- file.path(banc.nblast.mirror.save.path, "results")
native.nblast.folder <- file.path(banc.nblast.native.save.path, "results")
manc.nblast.folder <- file.path(banc.nblast.manc.save.path, "results", version)
fanc.nblast.folder <- file.path(banc.nblast.fanc.save.path, "results", version)
hemibrain.nblast.folder <- file.path(banc.nblast.hemibrain.save.path, "results_with_mirrored", version)
malecns.nblast.folder <- file.path(banc.nblast.malecns.save.path, "results", banc.nblast.malecns.version)

#######################################
### Collect manual matching results ###
#######################################

# Helper: collect PNG matches from review folders.
# Returns a character vector of basenames (back-compat) with `mtime` and
# `confidence` attributes attached as a parallel data.frame, plus a `path`
# attribute carrying full paths. Downstream code that just iterates basenames
# keeps working; the recency-aware rank_png_matches() consumes the attributes.
collect_png_matches <- function(correct_match_path, confidence_levels = 1:5) {
  rows <- list()
  for (level in confidence_levels) {
    folder_name <- c("1_perfect", "2_confident", "3_good", "4_likely", "5_possible")[level]
    folder <- file.path(correct_match_path, folder_name)
    if (!dir.exists(folder)) next
    files <- list.files(folder, pattern = "\\.png$", full.names = TRUE)
    if (!length(files)) next
    rows[[length(rows) + 1L]] <- data.frame(
      path = files,
      file = basename(files),
      mtime = as.numeric(file.info(files)$mtime),
      confidence = level,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    out <- character(0)
    attr(out, "meta") <- data.frame(file = character(), path = character(),
                                    mtime = numeric(), confidence = integer())
    return(out)
  }
  df <- do.call(rbind, rows)
  # Newer mtime / lower confidence number (more confident bucket) takes priority
  # in the dedupe, so the same basename appearing under multiple confidence
  # folders keeps the "best, most recent" entry.
  df <- df[order(-df$mtime, df$confidence), , drop = FALSE]
  df <- df[!duplicated(df$file), , drop = FALSE]
  out <- df$file
  attr(out, "meta") <- df
  out
}

# Helper: pick one row per pt_root_id from a parsed PNG-match data.frame.
# Ranking (descending priority):
#   1. is_current   — query_id matches pt_root_id (PNG generated at the
#                     neuron's current root_id, i.e. most recent NBLAST run).
#   2. mtime        — PNG file mtime; newer wins for older-root collisions.
#   3. score        — NBLAST score; tie-break per user spec ("if two matches
#                     have the same root_id but a different hit, choose the
#                     one with the highest NBLAST score").
# Inputs:
#   df         — parsed match df with at minimum pt_root_id, match_id,
#                query_id columns. May also carry `file` (basename).
#   nb_df      — NBLAST feather (rename pt_root_id, match_id, score).
#   png_meta   — attr(matches, "meta") from collect_png_matches; supplies
#                file → mtime mapping. NULL is OK (treat all mtimes as 0).
rank_png_matches <- function(df, nb_df = NULL, png_meta = NULL,
                              key_col = "pt_root_id") {
  if (!nrow(df)) return(df)
  df$.is_current <- as.integer(!is.na(df$query_id) & df$query_id == df[[key_col]])
  if (!is.null(png_meta) && "file" %in% colnames(df)) {
    df$.mtime <- png_meta$mtime[match(df$file, png_meta$file)]
    df$.mtime[is.na(df$.mtime)] <- 0
  } else {
    df$.mtime <- 0
  }
  if (!is.null(nb_df) && nrow(nb_df) &&
      all(c(key_col, "match_id", "score") %in% colnames(nb_df))) {
    join_keys <- c(key_col, "match_id")
    score_idx <- match(paste(df[[key_col]], df$match_id, sep = "\r"),
                       paste(nb_df[[key_col]], nb_df$match_id, sep = "\r"))
    df$.score <- nb_df$score[score_idx]
    df$.score[is.na(df$.score)] <- 0
  } else {
    df$.score <- 0
  }
  df <- df[order(-df$.is_current, -df$.mtime, -df$.score), , drop = FALSE]
  df <- df[!duplicated(df[[key_col]]), , drop = FALSE]
  df$.is_current <- NULL; df$.mtime <- NULL; df$.score <- NULL
  df
}

### Mirror ###
message("Working on BANC-mirror reviewed matches ...")
mirror.matches <- collect_png_matches(banc.mirror.correct.match.path)
mirror.matches.df <- data.frame()
for (file in mirror.matches) {
  mdf <- data.frame(root_id = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*", "", file),
                    supervoxel_id = ifelse(grepl("supervoxel_id", file), gsub(".*_supervoxel_id_|_hit.*|_query.*", "", file), NA),
                    banc_match = gsub(".*_hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*", "", file))
  mirror.matches.df <- rbind(mirror.matches.df, mdf)
}
mirror.matches.df <- mirror.matches.df %>%
  dplyr::anti_join(mirror.matches.df.valid %>% dplyr::filter(valid == 't'),
                   by = c("banc_match" = "match_id", "root_id" = "query_id"))
if (nrow(mirror.matches.df)) {
  # Add reverse matches
  mirror.matches.df.missing <- mirror.matches.df %>%
    dplyr::filter(!banc_match %in% root_id) %>%
    dplyr::rename(root_id = banc_match, banc_match = root_id)
  mirror.matches.df <- rbind(mirror.matches.df, mirror.matches.df.missing) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)

  # Update IDs
  mirror.matches.df$query <- mirror.matches.df$root_id
  mirror.matches.df$match_id <- mirror.matches.df$banc_match
  mirror.matches.df$root_id <- fast_updateids(mirror.matches.df$root_id)
  mirror.matches.df$banc_match <- fast_updateids(mirror.matches.df$banc_match)

  # Assign positions
  mirror.matches.df$position <- banc.meta$position[match(mirror.matches.df$root_id, banc.meta$root_id)]
  mirror.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(mirror.matches.df$root_id, banc.meta$root_id)]
  mirror.matches.df$banc_match_position <- banc.meta$position[match(mirror.matches.df$banc_match, banc.meta$root_id)]
  mirror.matches.df$banc_match_supervoxel_id <- banc.meta$supervoxel_id[match(mirror.matches.df$banc_match, banc.meta$root_id)]

  # Organise
  mirror.matches.df <- mirror.matches.df %>%
    dplyr::rename(pt_root_id = root_id,
                  pt_supervoxel_id = supervoxel_id,
                  pt_position = position,
                  match_root_id = banc_match,
                  match_supervoxel_id = banc_match_supervoxel_id,
                  match_position = banc_match_position,
                  query_id = query) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                  match_root_id, match_supervoxel_id, match_position,
                  query_id, match_id) %>%
    dplyr::bind_rows(mirror.matches.df.valid) %>%
    dplyr::mutate(valid = 't') %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::filter(pt_root_id != match_root_id) %>%
    dplyr::distinct(pt_root_id, match_root_id, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(match_supervoxel_id))

  # Save
  mirror.matches.df <- fast_updateids_df(mirror.matches.df,
                                         root.column = "pt_root_id",
                                         position.column = "pt_position",
                                         supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(mirror.matches.df, file.path(banc.meta.save.path, "banc_mirror_reviewed_matches.csv"))
}

### FAFB ###
message("Working on BANC-FAFB reviewed matches ...")
fafb.matches <- collect_png_matches(banc.fafb.correct.match.path)
fafb.png.meta <- attr(fafb.matches, "meta")
fafb.matches.df <- data.frame()
for (file in fafb.matches) {
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*", "", file),
                    supervoxel_ids = ifelse(grepl("supervoxel_id", file), gsub(".*_supervoxel_id_|_hit.*|_query.*", "", file), NA),
                    fafb_match = gsub(".*root_783_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png", "", file),
                    file = file,
                    stringsAsFactors = FALSE)
  fafb.matches.df <- rbind(fafb.matches.df, mdf)
}
fafb.matches.df <- fafb.matches.df %>%
  dplyr::filter(!grepl("\\_", fafb_match)) %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
fafb.matches.df <- fafb.matches.df %>%
  dplyr::anti_join(fafb.matches.df.valid %>% dplyr::filter(valid == 't'),
                   by = c("fafb_match" = "match_id", "root_id" = "query_id"))
if (nrow(fafb.matches.df)) {
  fafb.matches.df$query <- fafb.matches.df$root_id
  fafb.matches.df <- fast_updateids_df(fafb.matches.df, root.column = "root_id", supervoxel.column = "supervoxel_id")
  fafb.matches.df$position <- banc.meta$position[match(fafb.matches.df$root_id, banc.meta$root_id)]
  fafb.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(fafb.matches.df$root_id, banc.meta$root_id)]

  fafb.matches.df <- fafb.matches.df %>%
    dplyr::left_join(fw.meta[, c("root_783", "cell_type")], by = c("fafb_match" = "root_783")) %>%
    dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                  pt_position = position, query_id = query,
                  match_id = fafb_match, match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, file) %>%
    dplyr::mutate(valid = 't') %>%
    plyr::rbind.fill(fafb.matches.df.valid)
  # Recency-aware dedupe: one match per pt_root_id, preferring rows generated
  # against the CURRENT root_id (query_id == pt_root_id), then by PNG mtime,
  # then by NBLAST score (banc.meta.fafb.nb).
  fafb.matches.df <- rank_png_matches(fafb.matches.df,
                                       nb_df = banc.meta.fafb.nb,
                                       png_meta = fafb.png.meta)
  if ("file" %in% colnames(fafb.matches.df)) fafb.matches.df$file <- NULL

  fafb.matches.df <- fast_updateids_df(fafb.matches.df, root.column = "pt_root_id",
                                       position.column = "pt_position", supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(fafb.matches.df, file.path(banc.meta.save.path, "banc_fafb_reviewed_matches.csv"))
} else {
  fafb.matches.df <- fafb.matches.df.valid
}

### MANC ###
message("Working on BANC-MANC reviewed matches ...")
manc.matches <- collect_png_matches(banc.manc.correct.match.path)
manc.png.meta <- attr(manc.matches, "meta")
manc.matches.df <- data.frame()
for (file in manc.matches) {
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*", "", file),
                    supervoxel_ids = ifelse(grepl("supervoxel_id", file), gsub(".*_supervoxel_id_|_hit.*|_query.*", "", file), NA),
                    manc_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png", "", file),
                    file = file,
                    stringsAsFactors = FALSE)
  manc.matches.df <- rbind(manc.matches.df, mdf)
}
manc.matches.df <- manc.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
manc.matches.df <- manc.matches.df %>%
  dplyr::anti_join(manc.matches.df.valid %>% dplyr::filter(valid == 't'),
                   by = c("manc_match" = "match_id", "root_id" = "query_id"))
if (nrow(manc.matches.df)) {
  manc.matches.df$query <- manc.matches.df$root_id
  manc.matches.df <- fast_updateids_df(manc.matches.df, root.column = "root_id", supervoxel.column = "supervoxel_id")
  manc.matches.df$position <- banc.meta$position[match(manc.matches.df$root_id, banc.meta$root_id)]
  manc.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(manc.matches.df$root_id, banc.meta$root_id)]

  manc.matches.df <- manc.matches.df %>%
    dplyr::left_join(mc.meta[, c("bodyid", "cell_type")], by = c("manc_match" = "bodyid")) %>%
    dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                  pt_position = position, query_id = query,
                  match_id = manc_match, match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, file) %>%
    dplyr::mutate(valid = 't') %>%
    plyr::rbind.fill(manc.matches.df.valid)
  # Recency-aware dedupe — see rank_png_matches() definition above.
  manc.matches.df <- rank_png_matches(manc.matches.df,
                                       nb_df = banc.meta.manc.nb,
                                       png_meta = manc.png.meta)
  if ("file" %in% colnames(manc.matches.df)) manc.matches.df$file <- NULL

  manc.matches.df <- fast_updateids_df(manc.matches.df, root.column = "pt_root_id",
                                       position.column = "pt_position", supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(manc.matches.df, file.path(banc.meta.save.path, "banc_manc_reviewed_matches.csv"))
} else {
  manc.matches.df <- manc.matches.df.valid
}

### hemibrain ###
message("Working on BANC-hemibrain reviewed matches ...")
hemibrain.matches <- collect_png_matches(banc.hemibrain.correct.match.path, confidence_levels = 1:4)
hemibrain.png.meta <- attr(hemibrain.matches, "meta")
hemibrain.matches.df <- data.frame()
for (file in hemibrain.matches) {
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id_.*", "", file),
                    supervoxel_ids = ifelse(grepl("supervoxel_id", file), gsub(".*_supervoxel_id_|_hit.*|_query.*", "", file), NA),
                    hemibrain_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*|\\.png", "", file),
                    file = file,
                    stringsAsFactors = FALSE)
  hemibrain.matches.df <- rbind(hemibrain.matches.df, mdf)
}
hemibrain.matches.df <- hemibrain.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
hemibrain.matches.df <- hemibrain.matches.df %>%
  dplyr::anti_join(hemibrain.matches.df.valid %>% dplyr::filter(valid == 't'),
                   by = c("hemibrain_match" = "match_id", "root_id" = "query_id"))
if (nrow(hemibrain.matches.df)) {
  hemibrain.matches.df$query <- hemibrain.matches.df$root_id
  hemibrain.matches.df <- fast_updateids_df(hemibrain.matches.df, root.column = "root_id", supervoxel.column = "supervoxel_id")
  hemibrain.matches.df$position <- banc.meta$position[match(hemibrain.matches.df$root_id, banc.meta$root_id)]
  hemibrain.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(hemibrain.matches.df$root_id, banc.meta$root_id)]

  hemibrain.matches.df <- hemibrain.matches.df %>%
    dplyr::left_join(hb.meta[, c("bodyid", "cell_type")], by = c("hemibrain_match" = "bodyid")) %>%
    dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                  pt_position = position, query_id = query,
                  match_id = hemibrain_match, match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, file) %>%
    dplyr::mutate(valid = 't') %>%
    plyr::rbind.fill(hemibrain.matches.df.valid)
  hemibrain.matches.df <- rank_png_matches(hemibrain.matches.df,
                                            nb_df = banc.meta.hemibrain.nb,
                                            png_meta = hemibrain.png.meta)
  if ("file" %in% colnames(hemibrain.matches.df)) hemibrain.matches.df$file <- NULL

  hemibrain.matches.df <- fast_updateids_df(hemibrain.matches.df, root.column = "pt_root_id",
                                            position.column = "pt_position", supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(hemibrain.matches.df, file.path(banc.meta.save.path, "banc_hemibrain_reviewed_matches.csv"))
} else {
  hemibrain.matches.df <- hemibrain.matches.df.valid
}

### FANC ###
message("Working on BANC-fanc reviewed matches ...")
fanc.matches <- collect_png_matches(banc.fanc.correct.match.path, confidence_levels = 1:4)
fanc.png.meta <- attr(fanc.matches, "meta")
fanc.matches.df <- data.frame()
for (file in fanc.matches) {
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id_.*", "", file),
                    supervoxel_ids = gsub(".*supervoxel_id_([0-9]+).*", "\\1", file),
                    fanc_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*|\\.png", "", file),
                    file = file,
                    stringsAsFactors = FALSE)
  fanc.matches.df <- rbind(fanc.matches.df, mdf)
}
fanc.matches.df <- fanc.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
fanc.matches.df <- fanc.matches.df %>%
  dplyr::anti_join(fanc.matches.df.valid %>% dplyr::filter(valid == 't'),
                   by = c("fanc_match" = "match_id", "root_id" = "query_id"))
if (nrow(fanc.matches.df)) {
  fanc.matches.df$query <- fanc.matches.df$root_id
  fanc.matches.df <- fast_updateids_df(fanc.matches.df, root.column = "root_id", supervoxel.column = "supervoxel_id")
  fanc.matches.df$position <- banc.meta$position[match(fanc.matches.df$root_id, banc.meta$root_id)]
  fanc.matches.df$fanc_match <- fc.meta$cell_id[match(fanc.matches.df$fanc_match, fc.meta$root_id)]

  fanc.matches.df <- fanc.matches.df %>%
    dplyr::left_join(fc.meta[, c("cell_id", "cell_type")], by = c("fanc_match" = "cell_id")) %>%
    dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                  pt_position = position, query_id = query,
                  match_id = fanc_match, match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, file) %>%
    dplyr::mutate(valid = 't') %>%
    plyr::rbind.fill(fanc.matches.df.valid)
  fanc.matches.df <- rank_png_matches(fanc.matches.df,
                                       nb_df = banc.meta.fanc.nb,
                                       png_meta = fanc.png.meta)
  if ("file" %in% colnames(fanc.matches.df)) fanc.matches.df$file <- NULL

  fanc.matches.df <- fast_updateids_df(fanc.matches.df, root.column = "pt_root_id",
                                       position.column = "pt_position", supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(fanc.matches.df, file.path(banc.meta.save.path, "banc_fanc_reviewed_matches.csv"))
} else {
  fanc.matches.df <- fanc.matches.df.valid
}

### maleCNS ###
message("Working on BANC-maleCNS reviewed matches ...")
malecns.matches <- collect_png_matches(banc.malecns.correct.match.path)
malecns.png.meta <- attr(malecns.matches, "meta")
malecns.matches.df <- data.frame()
for (file in malecns.matches) {
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*", "", file),
                    supervoxel_ids = ifelse(grepl("supervoxel_id", file), stringr::str_match(file, "supervoxel_id_([0-9]+)_")[, 2], NA),
                    malecns_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png", "", file),
                    file = file,
                    stringsAsFactors = FALSE)
  malecns.matches.df <- rbind(malecns.matches.df, mdf)
}
if (nrow(malecns.matches.df)) {
  malecns.matches.df <- malecns.matches.df %>%
    tidyr::separate_rows(root_ids, sep = "_") %>%
    dplyr::rename(root_id = root_ids) %>%
    tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
    dplyr::rename(supervoxel_id = supervoxel_ids)
  malecns.matches.df <- malecns.matches.df %>%
    dplyr::anti_join(malecns.matches.df.valid %>% dplyr::filter(valid == 't'),
                     by = c("malecns_match" = "match_id", "root_id" = "query_id"))
}
if (nrow(malecns.matches.df)) {
  malecns.matches.df$query <- malecns.matches.df$root_id
  malecns.matches.df <- fast_updateids_df(malecns.matches.df, root.column = "root_id", supervoxel.column = "supervoxel_id")
  malecns.matches.df$position <- banc.meta$position[match(malecns.matches.df$root_id, banc.meta$root_id)]
  malecns.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(malecns.matches.df$root_id, banc.meta$root_id)]

  malecns.matches.df <- malecns.matches.df %>%
    dplyr::left_join(mcns.meta[, c("malecns_09_id", "cell_type")], by = c("malecns_match" = "malecns_09_id")) %>%
    dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                  pt_position = position, query_id = query,
                  match_id = malecns_match, match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, file) %>%
    dplyr::mutate(valid = 't') %>%
    plyr::rbind.fill(malecns.matches.df.valid)
  malecns.matches.df <- rank_png_matches(malecns.matches.df,
                                          nb_df = banc.meta.malecns.nb,
                                          png_meta = malecns.png.meta)
  if ("file" %in% colnames(malecns.matches.df)) malecns.matches.df$file <- NULL

  malecns.matches.df <- fast_updateids_df(malecns.matches.df, root.column = "pt_root_id",
                                          position.column = "pt_position", supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(malecns.matches.df, file.path(banc.meta.save.path, "banc_malecns_reviewed_matches.csv"))
} else {
  malecns.matches.df <- malecns.matches.df.valid
}

############################
### maleCNS NBLAST list  ###
############################

message("Working on BANC-maleCNS NBLAST matches ...")
malecns.nblast.files <- list.files(malecns.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.malecns.nb$query_id)
malecns.nblast.files <- malecns.nblast.files[!gsub(".*_|.csv", "", basename(malecns.nblast.files)) %in% done]
batch_size <- 250
total_files <- length(malecns.nblast.files)
if (total_files) {
  num_batches <- ceiling(total_files / batch_size)
  for (batch_idx in seq_len(num_batches)) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, total_files)
    current_batch <- malecns.nblast.files[start_idx:end_idx]
    message(sprintf("Processing maleCNS batch %d/%d (files %d to %d)",
                    batch_idx, num_batches, start_idx, end_idx))
    if (length(current_batch)) {
      by.query.malecns <- foreach::foreach(mfile = current_batch, .errorhandling = 'pass') %do% {
        id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
        sv_id <- stringr::str_match(basename(mfile), "supervoxel_id_([0-9]+)_root_id")[, 2]
        mdf <- safe_read_nblast_csv(mfile, col_types = banc.col.types)
        if (is.null(mdf)) return(NULL)
        numeric_columns <- sapply(mdf, is.numeric)
        mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
        mdf$query <- id
        mdf$supervoxel_id <- sv_id
        mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
          dplyr::group_by(query) %>%
          dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
          dplyr::ungroup()
        mdf
      }
      if (!length(by.query.malecns)) break
      by.query.malecns <- by.query.malecns[unlist(lapply(by.query.malecns, is.data.frame))]

      tryCatch({
      malecns.nblast.scores.all <- dplyr::bind_rows(by.query.malecns)
      malecns.nblast.scores.all$root_location <- NULL
      malecns.nblast.scores.all <- malecns.nblast.scores.all %>%
        dplyr::mutate(root_id = query) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(nb)) %>%
        dplyr::filter(nb <= 1)

      # Update validated matches
      if (nrow(malecns.matches.df)) {
        banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
          dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
          dplyr::select(-valid) %>%
          dplyr::left_join(malecns.matches.df[, c("query_id", "match_id", "valid")],
                           by = c("query_id", "match_id"),
                           relationship = "many-to-many") %>%
          dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))
      }

      banc.meta.malecns.scores.todo <- malecns.nblast.scores.all %>%
        dplyr::anti_join(banc.meta.malecns.nb, by = c("bodyid" = "match_id", "query" = "query_id"))
      if (nrow(banc.meta.malecns.scores.todo)) {
        # Use supervoxel_id (from filename) for primary root_id lookup
        svids.all <- banc.meta.malecns.scores.todo$supervoxel_id
        svids <- unique(na.omit(svids.all))
        root_ids.all <- rep(NA_character_, nrow(banc.meta.malecns.scores.todo))
        if (length(svids)) {
          root_ids_sv <- banc_rootid(svids)
          root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
        }
        # Fall back to fast_updateids for any that failed (0, NA, or no supervoxel)
        bad <- is.na(root_ids.all) | root_ids.all == "0"
        if (any(bad)) {
          fallback <- fast_updateids(banc.meta.malecns.scores.todo$root_id[bad])
          root_ids.all[bad] <- fallback
        }
        banc.meta.malecns.scores.todo$root_id <- root_ids.all
        banc.meta.malecns.scores.todo$position <- banc.meta$position[match(banc.meta.malecns.scores.todo$root_id, banc.meta$root_id)]
        banc.meta.malecns.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.malecns.scores.todo$root_id, banc.meta$root_id)]
        banc.meta.malecns.scores.todo <- subset(banc.meta.malecns.scores.todo, !is.na(position))

        banc.meta.malecns.nb <- banc.meta.malecns.scores.todo %>%
          dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                        pt_position = position, query_id = query,
                        match_id = bodyid, match_cell_type = cell_type, score = nb) %>%
          dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                        query_id, match_id, match_cell_type, score) %>%
          dplyr::arrange(dplyr::desc(score)) %>%
          dplyr::bind_rows(banc.meta.malecns.nb) %>%
          dplyr::select(-valid) %>%
          dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
          dplyr::left_join(malecns.matches.df[, c("pt_root_id", "match_id", "valid")],
                           by = c("pt_root_id", "match_id"),
                           relationship = "many-to-many") %>%
          dplyr::arrange(dplyr::desc(score)) %>%
          dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid)) %>%
          dplyr::filter(!is.na(pt_supervoxel_id), !is.na(pt_root_id))

        malecns.nblast.scores <- fast_updateids_df(banc.meta.malecns.nb, root.column = 'pt_root_id',
                                                   supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
        malecns.nblast.scores <- prune_stale_nblast(malecns.nblast.scores)
        malecns.nblast.scores$root_626 <- banc_rootid(malecns.nblast.scores$pt_supervoxel_id, version = "626")
        malecns.nblast.scores[[paste0("root_", banc.version)]] <- banc_rootid(malecns.nblast.scores$pt_supervoxel_id, version = banc.version)
        arrow::write_feather(internal_to_cave(malecns.nblast.scores), file.path(banc.meta.save.path, "banc_malecns_v0.9_nblast.feather"))
      }
      rm(by.query.malecns, malecns.nblast.scores.all); gc()
      }, error = function(e) {
        warning(sprintf("  maleCNS batch %d/%d failed: %s — skipping", batch_idx, num_batches, e$message))
      })
    }
  }
}

##########################
### Mirror NBLAST list ###
##########################

message("Working on BANC mirror NBLAST matches ...")
mirror.nblast.files <- list.files(mirror.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
mirror.nblast.files <- mirror.nblast.files[grepl("supervoxel_id", mirror.nblast.files)]
done <- unique(banc.meta.mirror.nb$query_id)
mirror.nblast.files <- mirror.nblast.files[!gsub(".*_root_id_|.csv", "", basename(mirror.nblast.files)) %in% done]

if (length(mirror.nblast.files)) {
  message(sprintf("  Reading %d mirror NBLAST files...", length(mirror.nblast.files)))
  n_mirror <- length(mirror.nblast.files)
  by.query.mirror <- foreach::foreach(mfile = mirror.nblast.files, i = seq_along(mirror.nblast.files)) %do% {
    if (i == 1 || i %% 200 == 0 || i == n_mirror)
      message(sprintf("    Mirror: %d/%d (%.0f%%)", i, n_mirror, 100 * i / n_mirror))
    tryCatch({
      id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
        dplyr::group_by(query) %>%
        dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
        dplyr::ungroup()
      mdf
    }, error = function(e) {
      message(sprintf("    Error reading %s: %s", basename(mfile), e$message))
      NULL
    })
  }
  by.query.mirror <- by.query.mirror[unlist(lapply(by.query.mirror, is.data.frame))]

  tryCatch({
  banc.meta.mirror.scores.all <- dplyr::bind_rows(by.query.mirror) %>%
    dplyr::mutate(banc_nblast_match = root_id, root_id = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(nb)) %>%
    dplyr::filter(nb <= 1)

  banc.meta.mirror.scores.todo <- banc.meta.mirror.scores.all %>%
    dplyr::anti_join(banc.meta.mirror.nb, by = c("banc_nblast_match" = "match_id", "query" = "query_id"))
  if (nrow(banc.meta.mirror.scores.todo)) {
    # Update query ids via supervoxel_id → banc_rootid, with fast_updateids fallback
    queries.all <- banc.meta.mirror.scores.todo$root_id
    svids.all <- banc.meta$supervoxel_id[match(queries.all, banc.meta$root_id)]
    svids <- unique(na.omit(svids.all))
    root_ids.all <- rep(NA_character_, length(queries.all))
    if (length(svids)) {
      root_ids_sv <- banc_rootid(svids)
      root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
    }
    bad <- is.na(root_ids.all) | root_ids.all == "0"
    if (any(bad)) {
      fallback <- fast_updateids(queries.all[bad])
      root_ids.all[bad] <- fallback
    }
    banc.meta.mirror.scores.todo$root_id <- root_ids.all
    banc.meta.mirror.scores.todo$position <- banc.meta$position[match(banc.meta.mirror.scores.todo$root_id, banc.meta$root_id)]
    banc.meta.mirror.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.mirror.scores.todo$root_id, banc.meta$root_id)]

    # Update match ids via supervoxel_id → banc_rootid, with fast_updateids fallback
    targets.all <- banc.meta.mirror.scores.todo$banc_nblast_match
    tgt_svids.all <- banc.meta$supervoxel_id[match(targets.all, banc.meta$root_id)]
    tgt_svids <- unique(na.omit(tgt_svids.all))
    banc_nblast_matches <- rep(NA_character_, length(targets.all))
    if (length(tgt_svids)) {
      tgt_rids <- banc_rootid(tgt_svids)
      banc_nblast_matches[!is.na(tgt_svids.all)] <- tgt_rids[match(tgt_svids.all[!is.na(tgt_svids.all)], tgt_svids)]
    }
    tgt_bad <- is.na(banc_nblast_matches) | banc_nblast_matches == "0"
    if (any(tgt_bad)) {
      tgt_fallback <- fast_updateids(targets.all[tgt_bad])
      banc_nblast_matches[tgt_bad] <- tgt_fallback
    }
    banc_nblast_matches.all <- banc_nblast_matches
    banc.meta.mirror.scores.todo$banc_nblast_match_latest <- banc_nblast_matches.all
    banc.meta.mirror.scores.todo$banc_match_position <- banc.meta$position[match(banc.meta.mirror.scores.todo$banc_nblast_match_latest, banc.meta$root_id)]
    banc.meta.mirror.scores.todo$banc_match_supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.mirror.scores.todo$banc_nblast_match_latest, banc.meta$root_id)]

    # Organise
    banc.meta.mirror.scores <- banc.meta.mirror.scores.todo %>%
      dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                    pt_position = position, match_root_id = banc_nblast_match_latest,
                    match_supervoxel_id = banc_match_supervoxel_id,
                    match_position = banc_match_position,
                    query_id = query, match_id = banc_nblast_match, score = nb) %>%
      dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                    match_root_id, match_supervoxel_id, match_position,
                    query_id, match_id, score) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::bind_rows(banc.meta.mirror.nb) %>%
      dplyr::select(-valid) %>%
      dplyr::filter(pt_root_id != match_root_id) %>%
      dplyr::distinct(pt_root_id, match_root_id, .keep_all = TRUE) %>%
      dplyr::left_join(mirror.matches.df[, c("pt_root_id", "match_id", "valid")],
                       by = c("pt_root_id", "match_id")) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

    # Save
    banc.meta.mirror.scores <- fast_updateids_df(banc.meta.mirror.scores)
    # Use deduplicated supervoxel lookups for versioned root IDs
    unique_pt_sv <- unique(na.omit(as.character(banc.meta.mirror.scores$pt_supervoxel_id)))
    message(sprintf("  Mirror compile: looking up %d unique pt supervoxels for versions 626/850...", length(unique_pt_sv)))
    pt_v626 <- banc_rootid(unique_pt_sv, version = "626")
    pt_v850 <- banc_rootid(unique_pt_sv, version = banc.version)
    pt_v626_map <- setNames(pt_v626, unique_pt_sv)
    pt_v850_map <- setNames(pt_v850, unique_pt_sv)
    banc.meta.mirror.scores$root_626 <- unname(pt_v626_map[as.character(banc.meta.mirror.scores$pt_supervoxel_id)])
    banc.meta.mirror.scores[[paste0("root_", banc.version)]] <- unname(pt_v850_map[as.character(banc.meta.mirror.scores$pt_supervoxel_id)])
    if ("match_supervoxel_id" %in% names(banc.meta.mirror.scores)) {
      unique_match_sv <- unique(na.omit(as.character(banc.meta.mirror.scores$match_supervoxel_id)))
      message(sprintf("  Mirror compile: looking up %d unique match supervoxels for versions 626/850...", length(unique_match_sv)))
      match_v626 <- banc_rootid(unique_match_sv, version = "626")
      match_v850 <- banc_rootid(unique_match_sv, version = banc.version)
      match_v626_map <- setNames(match_v626, unique_match_sv)
      match_v850_map <- setNames(match_v850, unique_match_sv)
      banc.meta.mirror.scores$match_root_626 <- unname(match_v626_map[as.character(banc.meta.mirror.scores$match_supervoxel_id)])
      banc.meta.mirror.scores[[paste0("match_root_", banc.version)]] <- unname(match_v850_map[as.character(banc.meta.mirror.scores$match_supervoxel_id)])
    }
    arrow::write_feather(banc.meta.mirror.scores, file.path(banc.meta.save.path, "banc_mirror_nblast.feather"))
    # Update banc.meta.mirror.nb so the refresh section doesn't overwrite with stale data
    banc.meta.mirror.nb <- banc.meta.mirror.scores
  }
  }, error = function(e) {
    warning(sprintf("  Mirror NBLAST processing failed: %s — skipping", e$message))
  })
}

##########################
### Native (BANC self) ###
##########################
###
### Within-BANC NBLAST without the L-R mirror step. Same schema as the mirror
### feather (so it can be read by the same downstream code), but distinct
### output file: banc_native_nblast.feather. There is no per-match manual
### review yet, so the `valid` column is initialised to 'f' for every row.
message("Working on BANC native NBLAST matches ...")
native.nblast.files <- list.files(native.nblast.folder, pattern = "\\.csv",
                                  full.names = TRUE, recursive = FALSE)
native.nblast.files <- native.nblast.files[grepl("supervoxel_id", native.nblast.files)]
done.native <- unique(banc.meta.native.nb$query_id)
native.nblast.files <- native.nblast.files[
  !gsub(".*_root_id_|.csv", "", basename(native.nblast.files)) %in% done.native
]

if (length(native.nblast.files)) {
  message(sprintf("  Reading %d native NBLAST files...", length(native.nblast.files)))
  n_native <- length(native.nblast.files)
  by.query.native <- foreach::foreach(mfile = native.nblast.files,
                                      i = seq_along(native.nblast.files)) %do% {
    if (i == 1 || i %% 200 == 0 || i == n_native)
      message(sprintf("    Native: %d/%d (%.0f%%)", i, n_native, 100 * i / n_native))
    tryCatch({
      id  <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
        dplyr::group_by(query) %>%
        dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
        dplyr::ungroup()
      mdf
    }, error = function(e) {
      message(sprintf("    Error reading %s: %s", basename(mfile), e$message))
      NULL
    })
  }
  by.query.native <- by.query.native[unlist(lapply(by.query.native, is.data.frame))]

  tryCatch({
    banc.meta.native.scores.all <- dplyr::bind_rows(by.query.native) %>%
      dplyr::mutate(banc_nblast_match = root_id, root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::filter(nb <= 1)

    banc.meta.native.scores.todo <- banc.meta.native.scores.all %>%
      dplyr::anti_join(banc.meta.native.nb,
                       by = c("banc_nblast_match" = "match_id", "query" = "query_id"))
    if (nrow(banc.meta.native.scores.todo)) {
      # Update query ids via supervoxel_id -> banc_rootid, fast_updateids fallback
      queries.all <- banc.meta.native.scores.todo$root_id
      svids.all   <- banc.meta$supervoxel_id[match(queries.all, banc.meta$root_id)]
      svids       <- unique(na.omit(svids.all))
      root_ids.all <- rep(NA_character_, length(queries.all))
      if (length(svids)) {
        root_ids_sv <- banc_rootid(svids)
        root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
      }
      bad <- is.na(root_ids.all) | root_ids.all == "0"
      if (any(bad)) {
        fallback <- fast_updateids(queries.all[bad])
        root_ids.all[bad] <- fallback
      }
      banc.meta.native.scores.todo$root_id       <- root_ids.all
      banc.meta.native.scores.todo$position      <- banc.meta$position[match(banc.meta.native.scores.todo$root_id, banc.meta$root_id)]
      banc.meta.native.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.native.scores.todo$root_id, banc.meta$root_id)]

      # Update match ids the same way
      targets.all     <- banc.meta.native.scores.todo$banc_nblast_match
      tgt_svids.all   <- banc.meta$supervoxel_id[match(targets.all, banc.meta$root_id)]
      tgt_svids       <- unique(na.omit(tgt_svids.all))
      banc_nblast_matches <- rep(NA_character_, length(targets.all))
      if (length(tgt_svids)) {
        tgt_rids <- banc_rootid(tgt_svids)
        banc_nblast_matches[!is.na(tgt_svids.all)] <- tgt_rids[match(tgt_svids.all[!is.na(tgt_svids.all)], tgt_svids)]
      }
      tgt_bad <- is.na(banc_nblast_matches) | banc_nblast_matches == "0"
      if (any(tgt_bad)) {
        tgt_fallback <- fast_updateids(targets.all[tgt_bad])
        banc_nblast_matches[tgt_bad] <- tgt_fallback
      }
      banc.meta.native.scores.todo$banc_nblast_match_latest <- banc_nblast_matches
      banc.meta.native.scores.todo$banc_match_position      <- banc.meta$position[match(banc.meta.native.scores.todo$banc_nblast_match_latest, banc.meta$root_id)]
      banc.meta.native.scores.todo$banc_match_supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.native.scores.todo$banc_nblast_match_latest, banc.meta$root_id)]

      banc.meta.native.scores <- banc.meta.native.scores.todo %>%
        dplyr::rename(pt_root_id          = root_id,
                      pt_supervoxel_id    = supervoxel_id,
                      pt_position         = position,
                      match_root_id       = banc_nblast_match_latest,
                      match_supervoxel_id = banc_match_supervoxel_id,
                      match_position      = banc_match_position,
                      query_id            = query,
                      match_id            = banc_nblast_match,
                      score               = nb) %>%
        dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                      match_root_id, match_supervoxel_id, match_position,
                      query_id, match_id, score) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::bind_rows(banc.meta.native.nb %>%
                          (function(df) if ("valid" %in% names(df)) dplyr::select(df, -valid) else df)) %>%
        dplyr::filter(pt_root_id != match_root_id) %>%
        dplyr::distinct(pt_root_id, match_root_id, .keep_all = TRUE) %>%
        dplyr::mutate(valid = 'f')

      # Versioned root IDs (root_626, root_888) via cached supervoxel lookups
      banc.meta.native.scores <- fast_updateids_df(banc.meta.native.scores)
      unique_pt_sv <- unique(na.omit(as.character(banc.meta.native.scores$pt_supervoxel_id)))
      message(sprintf("  Native compile: looking up %d unique pt supervoxels for versions 626/888...", length(unique_pt_sv)))
      if (length(unique_pt_sv)) {
        pt_v626 <- banc_rootid(unique_pt_sv, version = "626")
        pt_v888 <- banc_rootid(unique_pt_sv, version = banc.version)
        pt_v626_map <- setNames(pt_v626, unique_pt_sv)
        pt_v888_map <- setNames(pt_v888, unique_pt_sv)
        banc.meta.native.scores$root_626 <- unname(pt_v626_map[as.character(banc.meta.native.scores$pt_supervoxel_id)])
        banc.meta.native.scores[[paste0("root_", banc.version)]] <- unname(pt_v888_map[as.character(banc.meta.native.scores$pt_supervoxel_id)])
      }
      if ("match_supervoxel_id" %in% names(banc.meta.native.scores)) {
        unique_match_sv <- unique(na.omit(as.character(banc.meta.native.scores$match_supervoxel_id)))
        message(sprintf("  Native compile: looking up %d unique match supervoxels...", length(unique_match_sv)))
        if (length(unique_match_sv)) {
          match_v626 <- banc_rootid(unique_match_sv, version = "626")
          match_v888 <- banc_rootid(unique_match_sv, version = banc.version)
          match_v626_map <- setNames(match_v626, unique_match_sv)
          match_v888_map <- setNames(match_v888, unique_match_sv)
          banc.meta.native.scores$match_root_626 <- unname(match_v626_map[as.character(banc.meta.native.scores$match_supervoxel_id)])
          banc.meta.native.scores[[paste0("match_root_", banc.version)]] <- unname(match_v888_map[as.character(banc.meta.native.scores$match_supervoxel_id)])
        }
      }
      arrow::write_feather(banc.meta.native.scores, native_feather_path)
      banc.meta.native.nb <- banc.meta.native.scores
    }
  }, error = function(e) {
    warning(sprintf("  Native NBLAST processing failed: %s — skipping", e$message))
  })
}

########################
### FANC NBLAST list ###
########################

message("Working on BANC-FANC NBLAST matches ...")
fanc.nblast.files <- list.files(fanc.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.fanc.nb$query_id)
fanc.nblast.files <- fanc.nblast.files[!gsub(".*_|.csv", "", basename(fanc.nblast.files)) %in% done]
if (length(fanc.nblast.files)) {
  by.query.fanc <- foreach::foreach(mfile = fanc.nblast.files, .errorhandling = 'pass') %do% {
    id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
    mdf <- safe_read_nblast_csv(mfile, col_types = banc.col.types)
    if (is.null(mdf)) return(NULL)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
    mdf$query <- id
    mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
      dplyr::group_by(query) %>%
      dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
      dplyr::ungroup()
    mdf
  }
  if (length(by.query.fanc)) {
    by.query.fanc <- by.query.fanc[unlist(lapply(by.query.fanc, is.data.frame))]

    tryCatch({
    fanc.nblast.scores.all <- dplyr::bind_rows(by.query.fanc) %>%
      dplyr::mutate(root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::distinct(nb, fanc_id, cell_id, cell_type, query, root_id, supervoxel_id) %>%
      dplyr::filter(nb <= 1)

    banc.meta.fanc.scores.todo <- fanc.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.fanc.nb, by = c("cell_id" = "match_id", "query" = "query_id"))
    if (nrow(banc.meta.fanc.scores.todo)) {
      # Look up supervoxel_id from banc.meta using (possibly stale) root_id
      svids.all <- banc.meta$supervoxel_id[match(banc.meta.fanc.scores.todo$query, banc.meta$root_id)]
      svids <- unique(na.omit(svids.all))
      root_ids.all <- rep(NA_character_, length(svids.all))
      if (length(svids)) {
        root_ids_sv <- banc_rootid(svids)
        root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
      }
      # Fall back to fast_updateids for any that failed (0, NA, or no supervoxel)
      bad <- is.na(root_ids.all) | root_ids.all == "0"
      if (any(bad)) {
        fallback <- fast_updateids(banc.meta.fanc.scores.todo$query[bad])
        root_ids.all[bad] <- fallback
      }
      banc.meta.fanc.scores.todo$root_id <- root_ids.all
      banc.meta.fanc.scores.todo$position <- banc.meta$position[match(banc.meta.fanc.scores.todo$root_id, banc.meta$root_id)]
      banc.meta.fanc.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.fanc.scores.todo$root_id, banc.meta$root_id)]

      banc.meta.fanc.nb <- banc.meta.fanc.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                      pt_position = position, query_id = query,
                      match_id = cell_id, match_cell_type = cell_type, score = nb) %>%
        dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, query_id, match_id, match_cell_type, score) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::bind_rows(banc.meta.fanc.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
        dplyr::left_join(fanc.matches.df[, c("pt_root_id", "match_id", "valid")],
                         by = c("pt_root_id", "match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

      fanc.nblast.scores <- fast_updateids_df(banc.meta.fanc.nb, root.column = 'pt_root_id',
                                              supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      fanc.nblast.scores <- prune_stale_nblast(fanc.nblast.scores)
      fanc.nblast.scores$root_626 <- banc_rootid(fanc.nblast.scores$pt_supervoxel_id, version = "626")
      fanc.nblast.scores[[paste0("root_", banc.version)]] <- banc_rootid(fanc.nblast.scores$pt_supervoxel_id, version = banc.version)
      arrow::write_feather(internal_to_cave(fanc.nblast.scores), file.path(banc.meta.save.path, "banc_fanc_1116_nblast.feather"))
    }
    }, error = function(e) {
      warning(sprintf("  FANC NBLAST processing failed: %s — skipping", e$message))
    })
  }
}

########################
### FAFB NBLAST list ###
########################

message("Working on BANC-FAFB NBLAST matches ...")
fafb.nblast.files <- list.files(fafb.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.fafb.nb$query_id)
fafb.nblast.files <- fafb.nblast.files[!gsub(".*_|.csv", "", basename(fafb.nblast.files)) %in% done]
batch_size <- 250
total_files <- length(fafb.nblast.files)
num_batches <- ceiling(total_files / batch_size)
fafb.nblast.files <- sample(fafb.nblast.files)
for (batch_idx in seq_len(num_batches)) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, total_files)
  current_batch <- fafb.nblast.files[start_idx:end_idx]
  message(sprintf("Processing FAFB batch %d/%d (files %d to %d)",
                  batch_idx, num_batches, start_idx, end_idx))
  if (length(current_batch)) {
    by.query.fafb <- NULL
    by.query.fafb <- foreach::foreach(mfile = current_batch, .errorhandling = 'pass') %do% {
      id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
      mdf <- safe_read_nblast_csv(mfile, col_types = banc.col.types)
      if (!is.null(mdf) && length(mdf)) {
        numeric_columns <- sapply(mdf, is.numeric)
        mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
        mdf$query <- id
        mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
          dplyr::group_by(query) %>%
          dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
          dplyr::ungroup()
        mdf
      } else {
        NULL
      }
    }
    if (!length(by.query.fafb)) break
    by.query.fafb <- by.query.fafb[unlist(lapply(by.query.fafb, is.data.frame))]

    tryCatch({
    fafb.nblast.scores.all <- dplyr::bind_rows(by.query.fafb) %>%
      dplyr::mutate(root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::distinct(nb, root_783, nucleus_id, cell_type, query, root_id) %>%
      dplyr::filter(nb <= 1)

    # Update validated matches
    banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
      dplyr::select(-valid) %>%
      dplyr::left_join(fafb.matches.df[, c("query_id", "match_id", "valid")],
                       by = c("query_id", "match_id"),
                       relationship = "many-to-many") %>%
      dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

    banc.meta.fafb.scores.todo <- fafb.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.fafb.nb, by = c("root_783" = "match_id", "query" = "query_id"))
    if (nrow(banc.meta.fafb.scores.todo)) {
      # Look up supervoxel_id from banc.meta using (possibly stale) root_id
      svids.all <- banc.meta$supervoxel_id[match(banc.meta.fafb.scores.todo$root_id, banc.meta$root_id)]
      svids <- unique(na.omit(svids.all))
      root_ids.all <- rep(NA_character_, length(svids.all))
      if (length(svids)) {
        root_ids_sv <- banc_rootid(svids)
        root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
      }
      # Fall back to fast_updateids for any that failed (0, NA, or no supervoxel)
      bad <- is.na(root_ids.all) | root_ids.all == "0"
      if (any(bad)) {
        fallback <- fast_updateids(banc.meta.fafb.scores.todo$root_id[bad])
        root_ids.all[bad] <- fallback
      }
      banc.meta.fafb.scores.todo$root_id <- root_ids.all
      banc.meta.fafb.scores.todo$position <- banc.meta$position[match(banc.meta.fafb.scores.todo$root_id, banc.meta$root_id)]
      banc.meta.fafb.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.fafb.scores.todo$root_id, banc.meta$root_id)]

      banc.meta.fafb.nb <- banc.meta.fafb.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                      pt_position = position, query_id = query,
                      match_id = root_783, match_cell_type = cell_type, score = nb) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                      query_id, match_id, match_cell_type, score) %>%
        dplyr::bind_rows(banc.meta.fafb.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
        dplyr::left_join(fafb.matches.df[, c("pt_root_id", "match_id", "valid")],
                         by = c("pt_root_id", "match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid)) %>%
        dplyr::filter(!is.na(pt_supervoxel_id), !is.na(score), score <= 1)

      fafb.nblast.scores <- fast_updateids_df(banc.meta.fafb.nb, root.column = 'pt_root_id',
                                              supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      fafb.nblast.scores <- prune_stale_nblast(fafb.nblast.scores)
      fafb.nblast.scores$root_626 <- banc_rootid(fafb.nblast.scores$pt_supervoxel_id, version = "626")
      fafb.nblast.scores[[paste0("root_", banc.version)]] <- banc_rootid(fafb.nblast.scores$pt_supervoxel_id, version = banc.version)
      arrow::write_feather(internal_to_cave(fafb.nblast.scores), file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather"))
    }
    }, error = function(e) {
      warning(sprintf("  FAFB batch %d/%d failed: %s — skipping", batch_idx, num_batches, e$message))
    })
  }
}

########################
### MANC NBLAST list ###
########################

message("Working on BANC-MANC NBLAST matches ...")
manc.nblast.files <- list.files(manc.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.manc.nb$query_id)
manc.nblast.files <- manc.nblast.files[!gsub(".*_|.csv", "", basename(manc.nblast.files)) %in% done]
batch_size <- 250
total_files <- length(manc.nblast.files)
num_batches <- ceiling(total_files / batch_size)
for (batch_idx in seq_len(num_batches)) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, total_files)
  current_batch <- manc.nblast.files[start_idx:end_idx]
  message(sprintf("Processing MANC batch %d/%d (files %d to %d)",
                  batch_idx, num_batches, start_idx, end_idx))
  if (length(current_batch)) {
    by.query.manc <- foreach::foreach(mfile = current_batch, .errorhandling = 'pass') %do% {
      id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
      mdf <- safe_read_nblast_csv(mfile, col_types = banc.col.types)
      if (is.null(mdf)) return(NULL)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
        dplyr::group_by(query) %>%
        dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
        dplyr::ungroup()
      mdf
    }
    if (!length(by.query.manc)) break
    by.query.manc <- by.query.manc[unlist(lapply(by.query.manc, is.data.frame))]

    tryCatch({
    manc.nblast.scores.all <- dplyr::bind_rows(by.query.manc)
    manc.nblast.scores.all$root_location <- NULL
    manc.nblast.scores.all <- manc.nblast.scores.all %>%
      dplyr::mutate(root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::filter(nb <= 1)

    # Update validated matches
    banc.meta.manc.nb <- banc.meta.manc.nb %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
      dplyr::select(-valid) %>%
      dplyr::left_join(manc.matches.df[, c("query_id", "match_id", "valid")],
                       by = c("query_id", "match_id"),
                       relationship = "many-to-many") %>%
      dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

    banc.meta.manc.scores.todo <- manc.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.manc.nb, by = c("bodyid" = "match_id", "query" = "query_id"))
    if (nrow(banc.meta.manc.scores.todo)) {
      # Look up supervoxel_id from banc.meta using (possibly stale) root_id
      svids.all <- banc.meta$supervoxel_id[match(banc.meta.manc.scores.todo$root_id, banc.meta$root_id)]
      svids <- unique(na.omit(svids.all))
      root_ids.all <- rep(NA_character_, length(svids.all))
      if (length(svids)) {
        root_ids_sv <- banc_rootid(svids)
        root_ids.all[!is.na(svids.all)] <- root_ids_sv[match(svids.all[!is.na(svids.all)], svids)]
      }
      # Fall back to fast_updateids for any that failed (0, NA, or no supervoxel)
      bad <- is.na(root_ids.all) | root_ids.all == "0"
      if (any(bad)) {
        fallback <- fast_updateids(banc.meta.manc.scores.todo$root_id[bad])
        root_ids.all[bad] <- fallback
      }
      banc.meta.manc.scores.todo$root_id <- root_ids.all
      banc.meta.manc.scores.todo$position <- banc.meta$position[match(banc.meta.manc.scores.todo$root_id, banc.meta$root_id)]
      banc.meta.manc.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.manc.scores.todo$root_id, banc.meta$root_id)]

      banc.meta.manc.nb <- banc.meta.manc.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                      pt_position = position, query_id = query,
                      match_id = bodyid, match_cell_type = cell_type, score = nb) %>%
        dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                      query_id, match_id, match_cell_type, score) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::bind_rows(banc.meta.manc.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
        dplyr::left_join(manc.matches.df[, c("pt_root_id", "match_id", "valid")],
                         by = c("pt_root_id", "match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid)) %>%
        dplyr::filter(!is.na(pt_supervoxel_id), !is.na(pt_root_id))

      manc.nblast.scores <- fast_updateids_df(banc.meta.manc.nb, root.column = 'pt_root_id',
                                              supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      manc.nblast.scores <- prune_stale_nblast(manc.nblast.scores)
      manc.nblast.scores$root_626 <- banc_rootid(manc.nblast.scores$pt_supervoxel_id, version = "626")
      manc.nblast.scores[[paste0("root_", banc.version)]] <- banc_rootid(manc.nblast.scores$pt_supervoxel_id, version = banc.version)
      arrow::write_feather(internal_to_cave(manc.nblast.scores), file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
    }
    }, error = function(e) {
      warning(sprintf("  MANC batch %d/%d failed: %s — skipping", batch_idx, num_batches, e$message))
    })
  }
}

#############################
### Hemibrain NBLAST list ###
#############################

message("Working on BANC-hemibrain NBLAST matches ...")
hemibrain.nblast.files <- list.files(hemibrain.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.hemibrain.nb$query_id)
hemibrain.nblast.files <- hemibrain.nblast.files[!gsub(".*_|.csv", "", basename(hemibrain.nblast.files)) %in% done]
if (length(hemibrain.nblast.files)) {
  by.query.hemibrain <- foreach::foreach(mfile = hemibrain.nblast.files, .errorhandling = 'pass') %do% {
    id <- gsub(".*_root_id_|\\.csv", "", basename(mfile))
    svid <- gsub(".*supervoxel_id_|_root_id_.*|\\.csv", "", basename(mfile))
    mdf <- safe_read_nblast_csv(mfile, col_types = banc.col.types)
    if (is.null(mdf)) return(NULL)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
    mdf$query <- id
    mdf$query_supervoxel_id <- svid
    mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
      dplyr::group_by(query) %>%
      dplyr::filter(nb >= 0.1 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
      dplyr::ungroup()
    mdf
  }
  if (length(by.query.hemibrain)) {
    by.query.hemibrain <- by.query.hemibrain[unlist(lapply(by.query.hemibrain, is.data.frame))]

    tryCatch({
    hemibrain.nblast.scores.all <- dplyr::bind_rows(by.query.hemibrain) %>%
      dplyr::mutate(root_id = query, supervoxel_id = query_supervoxel_id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb), 2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::filter(nb <= 1)

    # Update validated matches
    banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
      dplyr::select(-valid) %>%
      dplyr::left_join(hemibrain.matches.df[, c("query_id", "match_id", "valid")],
                       by = c("query_id", "match_id"),
                       relationship = "many-to-many") %>%
      dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid))

    banc.meta.hemibrain.scores.todo <- hemibrain.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.hemibrain.nb, by = c("bodyid" = "match_id", "query" = "query_id"))
    if (nrow(banc.meta.hemibrain.scores.todo)) {
      svids.all <- banc.meta.hemibrain.scores.todo$supervoxel_id
      svids <- unique(svids.all)
      root_ids <- banc_rootid(svids)
      root_ids.all <- root_ids[match(svids.all, svids)]
      # Fall back to fast_updateids for any that failed (0 or NA)
      bad <- is.na(root_ids.all) | root_ids.all == "0"
      if (any(bad)) {
        fallback <- fast_updateids(banc.meta.hemibrain.scores.todo$root_id[bad])
        root_ids.all[bad] <- fallback
      }
      banc.meta.hemibrain.scores.todo$root_id <- root_ids.all
      banc.meta.hemibrain.scores.todo$position <- banc.meta$position[match(banc.meta.hemibrain.scores.todo$root_id, banc.meta$root_id)]

      banc.meta.hemibrain.nb <- banc.meta.hemibrain.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id,
                      pt_position = position, query_id = query,
                      match_id = bodyid, match_cell_type = cell_type, score = nb) %>%
        dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                      query_id, match_id, match_cell_type, score) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::bind_rows(banc.meta.hemibrain.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
        dplyr::left_join(hemibrain.matches.df[, c("pt_root_id", "match_id", "valid")],
                         by = c("pt_root_id", "match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid), 'f', valid)) %>%
        dplyr::filter(!is.na(pt_supervoxel_id))

      hemibrain.nblast.scores <- fast_updateids_df(banc.meta.hemibrain.nb, root.column = 'pt_root_id',
                                                   supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      hemibrain.nblast.scores <- prune_stale_nblast(hemibrain.nblast.scores)
      hemibrain.nblast.scores$root_626 <- banc_rootid(hemibrain.nblast.scores$pt_supervoxel_id, version = "626")
      hemibrain.nblast.scores[[paste0("root_", banc.version)]] <- banc_rootid(hemibrain.nblast.scores$pt_supervoxel_id, version = banc.version)
      arrow::write_feather(internal_to_cave(hemibrain.nblast.scores), file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))
    }
    }, error = function(e) {
      warning(sprintf("  Hemibrain NBLAST processing failed: %s — skipping", e$message))
    })
  }
}

######################################################
### Refresh pt_root_id, root_626, root_888 in all ###
### feathers to keep IDs current.                  ###
######################################################

message("### Refreshing root IDs in all NBLAST feathers ###")

# Helper: update IDs and rewrite a feather
refresh_feather <- function(nb_df, feather_path, use_cave_format = TRUE) {
  tryCatch({
    message(sprintf("  Refreshing %s (%d rows)...", basename(feather_path), nrow(nb_df)))
    updated <- fast_updateids_df(nb_df, root.column = 'pt_root_id',
                                 supervoxel.column = 'pt_supervoxel_id',
                                 position.column = 'pt_position')
    # Drop redundant stale rows (pt_root_id migrated by fast_updateids_df
    # but query_id left at pre-edit value) when a fresh row exists for the
    # same pt_root_id. Stale-only pt_root_ids are preserved as best-we-have.
    updated <- prune_stale_nblast(updated)
    # Deduplicate supervoxels for versioned lookups
    unique_sv <- unique(na.omit(as.character(updated$pt_supervoxel_id)))
    message(sprintf("  Looking up %d unique supervoxels for versions 626/850...", length(unique_sv)))
    v626 <- banc_rootid(unique_sv, version = "626")
    v850 <- banc_rootid(unique_sv, version = banc.version)
    v626_map <- setNames(v626, unique_sv)
    v850_map <- setNames(v850, unique_sv)
    updated$root_626 <- unname(v626_map[as.character(updated$pt_supervoxel_id)])
    updated[[paste0("root_", banc.version)]] <- unname(v850_map[as.character(updated$pt_supervoxel_id)])
    if (use_cave_format) {
      arrow::write_feather(internal_to_cave(updated), feather_path)
    } else {
      arrow::write_feather(updated, feather_path)
    }
    message(sprintf("  Done: %s", basename(feather_path)))
    updated
  }, error = function(e) {
    warning(sprintf("  Failed to refresh %s: %s", basename(feather_path), e$message))
    nb_df
  })
}

banc.meta.malecns.nb <- refresh_feather(banc.meta.malecns.nb,
  file.path(banc.meta.save.path, "banc_malecns_v0.9_nblast.feather"))
banc.meta.manc.nb <- refresh_feather(banc.meta.manc.nb,
  file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
banc.meta.fafb.nb <- refresh_feather(banc.meta.fafb.nb,
  file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather"))
banc.meta.fanc.nb <- refresh_feather(banc.meta.fanc.nb,
  file.path(banc.meta.save.path, "banc_fanc_1116_nblast.feather"))
banc.meta.hemibrain.nb <- refresh_feather(banc.meta.hemibrain.nb,
  file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))
# Mirror needs special handling: also update match_root_id via match_supervoxel_id
tryCatch({
  message(sprintf("  Refreshing banc_mirror_nblast.feather (%d rows)...", nrow(banc.meta.mirror.nb)))
  banc.meta.mirror.nb <- fast_updateids_df(banc.meta.mirror.nb, root.column = 'pt_root_id',
                                            supervoxel.column = 'pt_supervoxel_id',
                                            position.column = 'pt_position')
  # Drop redundant stale rows on the query side (same logic as cross-dataset
  # feathers). Mirror's query_id is the pt_root_id at NBLAST compute time;
  # fast_updateids_df migrates pt_root_id forward, producing stale rows.
  banc.meta.mirror.nb <- prune_stale_nblast(banc.meta.mirror.nb)
  # Deduplicate pt supervoxels for versioned lookups
  unique_pt_sv <- unique(na.omit(as.character(banc.meta.mirror.nb$pt_supervoxel_id)))
  message(sprintf("  Mirror refresh: looking up %d unique pt supervoxels for versions 626/850...", length(unique_pt_sv)))
  pt_v626 <- banc_rootid(unique_pt_sv, version = "626")
  pt_v850 <- banc_rootid(unique_pt_sv, version = banc.version)
  pt_v626_map <- setNames(pt_v626, unique_pt_sv)
  pt_v850_map <- setNames(pt_v850, unique_pt_sv)
  banc.meta.mirror.nb$root_626 <- unname(pt_v626_map[as.character(banc.meta.mirror.nb$pt_supervoxel_id)])
  banc.meta.mirror.nb[[paste0("root_", banc.version)]] <- unname(pt_v850_map[as.character(banc.meta.mirror.nb$pt_supervoxel_id)])
  if ("match_supervoxel_id" %in% names(banc.meta.mirror.nb)) {
    # Use banc.meta for current root_id (fast), CAVE fallback for missing
    banc.meta.mirror.nb$match_root_id <- fast_rootid_from_sv(banc.meta.mirror.nb$match_supervoxel_id)
    # Deduplicate match supervoxels for versioned lookups
    unique_match_sv <- unique(na.omit(as.character(banc.meta.mirror.nb$match_supervoxel_id)))
    message(sprintf("  Mirror refresh: looking up %d unique match supervoxels for versions 626/850...", length(unique_match_sv)))
    match_v626 <- banc_rootid(unique_match_sv, version = "626")
    match_v850 <- banc_rootid(unique_match_sv, version = banc.version)
    match_v626_map <- setNames(match_v626, unique_match_sv)
    match_v850_map <- setNames(match_v850, unique_match_sv)
    banc.meta.mirror.nb$match_root_626 <- unname(match_v626_map[as.character(banc.meta.mirror.nb$match_supervoxel_id)])
    banc.meta.mirror.nb[[paste0("match_root_", banc.version)]] <- unname(match_v850_map[as.character(banc.meta.mirror.nb$match_supervoxel_id)])
  }
  arrow::write_feather(banc.meta.mirror.nb, file.path(banc.meta.save.path, "banc_mirror_nblast.feather"))
  message("  Done: banc_mirror_nblast.feather")
}, error = function(e) {
  warning(sprintf("  Failed to refresh mirror feather: %s", e$message))
})

message(sprintf("### banc: NBLAST compilation complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
