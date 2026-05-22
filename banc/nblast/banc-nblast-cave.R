#' banc-nblast-cave — Sync compiled NBLAST match feathers to CAVE `cell_match` tables.
#'
#' Runs after `banc-nblast-compile.R`; deletes stale annotations (neurons
#' re-segmented since NBLAST) and uploads fresh ones.
#'
#' @section Reads:
#'   - per-dataset compiled NBLAST feathers, CAVE `cell_match` tables
#'
#' @section Writes:
#'   - CAVE `cell_match` tables (delete + upload)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast.sh`, `o2/production/o2_banc_nblast_compile_chain.sh`.

###########################################################
### Sync compiled NBLAST match feather files to CAVE tables
###
### Runs downstream of banc-nblast-compile.R.
### For each dataset (maleCNS, FAFB, hemibrain, MANC, FANC):
###   1. Read compiled feather (CAVE-format columns)
###   2. Read current CAVE table state
###   3. Delete stale annotations (neuron re-segmented since NBLAST)
###   4. Upload new annotations not yet in CAVE
###
### Uses cell_match schema:
###   pt (BoundSpatialPoint), query_root_id (int64),
###   match_id (string), score (float), validation (bool)
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: syncing NBLAST matches to CAVE ###")
t_start <- Sys.time()

bancr::choose_banc()

###############################################
### One-time table creation (run once only) ###
###############################################
# Uncomment and run interactively to create each table.
# Tables use cell_match schema with PUBLIC read access.
# Scores are bidirectional mean normalised NBLAST (UseAlpha=TRUE, smat_alpha.fcwb).
# Inclusion threshold: score >= 0.3 or top 5 matches per query neuron.
#
# IMPORTANT: voxel_resolution must be c(4, 4, 45) to match the BANC segmentation.
# CAVE resolves supervoxel_id from pt_position using this resolution.
#
# Columns:
#   pt_position      — BANC neuron position in raw voxel coordinates (4x4x45 nm)
#   pt_root_id       — current BANC root_id (updated by CAVE materialisation)
#   pt_supervoxel_id — BANC supervoxel_id (stable across re-segmentation)
#   query_root_id    — BANC root_id at the time the NBLAST was run
#   match_id         — ID in the target dataset (bodyid, root_783, cell_id, etc.)
#   score            — bidirectional mean normalised NBLAST score (0-1)
#   validation       — 't' if confirmed by manual review, 'f' if not
#   valid            — CAVE system column: 't' for current, 'f' for stale/deleted
#                      (set to 'f' by delete_annotation when NBLAST results are superseded)
#
# bancr::banc_cave_new_table("banc_malecns_nblast_v2", "cell_match",
#   voxel_resolution = c(4, 4, 45),
#   description = paste("Bidirectional mean normalised NBLAST matches: BANC -> maleCNS v0.9.",
#     "Scores use UseAlpha=TRUE with smat_alpha.fcwb.",
#     "Threshold: score >= 0.3 or top 5 per query.",
#     "Columns: pt (BoundSpatialPoint), query_root_id (BANC root_id when NBLASTed),",
#     "match_id (maleCNS bodyid), score (0-1),",
#     "validation ('t' if confirmed by manual review, 'f' if not),",
#     "valid (CAVE system: 'f' when superseded by re-NBLAST after re-segmentation)."),
#   write_permission = "GROUP", read_permission = "PUBLIC")
#
# bancr::banc_cave_new_table("banc_fafb_nblast_v2", "cell_match",
#   voxel_resolution = c(4, 4, 45),
#   description = paste("Bidirectional mean normalised NBLAST matches: BANC -> FAFB 783.",
#     "Scores use UseAlpha=TRUE with smat_alpha.fcwb.",
#     "Threshold: score >= 0.3 or top 5 per query.",
#     "Columns: pt (BoundSpatialPoint), query_root_id (BANC root_id when NBLASTed),",
#     "match_id (FAFB root_783), score (0-1),",
#     "validation ('t' if confirmed by manual review, 'f' if not),",
#     "valid (CAVE system: 'f' when superseded by re-NBLAST after re-segmentation)."),
#   write_permission = "GROUP", read_permission = "PUBLIC")
#
# bancr::banc_cave_new_table("banc_hemibrain_nblast_v2", "cell_match",
#   voxel_resolution = c(4, 4, 45),
#   description = paste("Bidirectional mean normalised NBLAST matches: BANC -> hemibrain v1.2.1.",
#     "Scores use UseAlpha=TRUE with smat_alpha.fcwb.",
#     "Threshold: score >= 0.3 or top 5 per query.",
#     "Columns: pt (BoundSpatialPoint), query_root_id (BANC root_id when NBLASTed),",
#     "match_id (hemibrain bodyid), score (0-1),",
#     "validation ('t' if confirmed by manual review, 'f' if not),",
#     "valid (CAVE system: 'f' when superseded by re-NBLAST after re-segmentation)."),
#   write_permission = "GROUP", read_permission = "PUBLIC")
#
# bancr::banc_cave_new_table("banc_manc_nblast_v2", "cell_match",
#   voxel_resolution = c(4, 4, 45),
#   description = paste("Bidirectional mean normalised NBLAST matches: BANC -> MANC v1.2.1.",
#     "Scores use UseAlpha=TRUE with smat_alpha.fcwb.",
#     "Threshold: score >= 0.3 or top 5 per query.",
#     "Columns: pt (BoundSpatialPoint), query_root_id (BANC root_id when NBLASTed),",
#     "match_id (MANC bodyid), score (0-1),",
#     "validation ('t' if confirmed by manual review, 'f' if not),",
#     "valid (CAVE system: 'f' when superseded by re-NBLAST after re-segmentation)."),
#   write_permission = "GROUP", read_permission = "PUBLIC")
#
# bancr::banc_cave_new_table("banc_fanc_nblast_v2", "cell_match",
#   voxel_resolution = c(4, 4, 45),
#   description = paste("Bidirectional mean normalised NBLAST matches: BANC -> FANC 1116.",
#     "Scores use UseAlpha=TRUE with smat_alpha.fcwb.",
#     "Threshold: score >= 0.3 or top 5 per query.",
#     "Columns: pt (BoundSpatialPoint), query_root_id (BANC root_id when NBLASTed),",
#     "match_id (FANC cell_id), score (0-1),",
#     "validation ('t' if confirmed by manual review, 'f' if not),",
#     "valid (CAVE system: 'f' when superseded by re-NBLAST after re-segmentation)."),
#   write_permission = "GROUP", read_permission = "PUBLIC")

###########################
### Helper function     ###
###########################

#' Sync a single NBLAST feather file to its CAVE table
#' @param feather_path Path to compiled feather file (CAVE-format columns)
#' @param table_name CAVE table name
#' @param banc.meta Data frame with current root_id/supervoxel_id mapping
#' @param client CAVE client
sync_nblast_cave_table <- function(feather_path, table_name, banc.meta, client) {
  message(sprintf("  Syncing %s ...", table_name))

  # Read compiled feather (already in CAVE column format)
  compiled <- arrow::read_feather(feather_path)
  message(sprintf("    Compiled feather: %d rows", nrow(compiled)))

  # Handle old-format feathers (query_id/valid) from before CAVE rename
  if ("query_id" %in% names(compiled)) compiled <- dplyr::rename(compiled, query_root_id = query_id)
  if ("valid" %in% names(compiled) && !"validation" %in% names(compiled)) {
    compiled$validation <- compiled$valid == "t"
    compiled$valid <- NULL
  }

  # Ensure key columns are character for consistent matching
  compiled$pt_supervoxel_id <- as.character(compiled$pt_supervoxel_id)
  compiled$query_root_id <- as.character(compiled$query_root_id)
  compiled$match_id <- as.character(compiled$match_id)

  # Read current CAVE table
  cave_data <- tryCatch({
    bancr::banc_cave_query(table_name, live = 2)
  }, error = function(e) {
    message(sprintf("    Could not read CAVE table (may not exist yet): %s", e$message))
    data.frame()
  })
  message(sprintf("    Current CAVE table: %d rows", nrow(cave_data)))

  # --- Delete stale annotations ---
  # Stale = CAVE rows where query_root_id is no longer the current root_id
  # (the BANC neuron was re-segmented since the NBLAST was run)
  n_deleted <- 0
  if (nrow(cave_data) && "query_root_id" %in% names(cave_data)) {
    cave_data$pt_supervoxel_id <- as.character(cave_data$pt_supervoxel_id)
    cave_data$query_root_id <- as.character(cave_data$query_root_id)
    cave_data$match_id <- as.character(cave_data$match_id)

    current_roots <- banc.meta$root_id[match(
      cave_data$pt_supervoxel_id,
      as.character(banc.meta$supervoxel_id)
    )]
    stale_mask <- !is.na(current_roots) &
      cave_data$query_root_id != as.character(current_roots)
    stale_ids <- cave_data$id[stale_mask]

    if (length(stale_ids)) {
      message(sprintf("    Deleting %d stale annotations ...", length(stale_ids)))
      batch_size <- 200
      for (start in seq(1, length(stale_ids), by = batch_size)) {
        end <- min(start + batch_size - 1, length(stale_ids))
        tryCatch(
          client$annotation$delete_annotation(table_name, stale_ids[start:end]),
          error = function(e) warning(sprintf("    Delete batch error: %s", e$message))
        )
      }
      n_deleted <- length(stale_ids)
    }
  }

  # --- Upload new annotations ---
  # Deduplicate: skip rows already in CAVE (by pt_supervoxel_id + match_id)
  if (nrow(cave_data) && "match_id" %in% names(cave_data)) {
    new_rows <- compiled %>%
      dplyr::anti_join(cave_data, by = c("pt_supervoxel_id", "match_id"))
  } else {
    new_rows <- compiled
  }

  # Also deduplicate within the feather itself (same neuron + same match)
  new_rows <- dplyr::distinct(new_rows, pt_supervoxel_id, match_id, .keep_all = TRUE)

  # Drop rows with invalid supervoxel_id (can't anchor to segmentation)
  n_bad_sv <- sum(is.na(new_rows$pt_supervoxel_id) | new_rows$pt_supervoxel_id == "0")
  if (n_bad_sv > 0) {
    message(sprintf("    Dropping %d rows with pt_supervoxel_id=0 or NA", n_bad_sv))
    new_rows <- new_rows %>%
      dplyr::filter(!is.na(pt_supervoxel_id), pt_supervoxel_id != "0")
  }
  message(sprintf("    New rows to upload: %d", nrow(new_rows)))

  n_uploaded <- 0
  n_failed <- 0
  if (nrow(new_rows)) {
    np <- reticulate::import("numpy")
    upload_batch_size <- 500

    # Stage and upload in batches of upload_batch_size
    batch_starts <- seq(1, nrow(new_rows), by = upload_batch_size)
    n_batches <- length(batch_starts)

    for (b in seq_along(batch_starts)) {
      start_i <- batch_starts[b]
      end_i <- min(start_i + upload_batch_size - 1, nrow(new_rows))
      batch <- new_rows[start_i:end_i, ]

      # Stage all rows in this batch
      stage <- client$annotation$stage_annotations(table_name)
      n_staged <- 0
      for (i in seq_len(nrow(batch))) {
        pos <- as.numeric(trimws(strsplit(as.character(batch$pt_position[i]), ",")[[1]]))
        if (length(pos) != 3 || any(is.na(pos))) next

        stage$add(
          pt_position = np$array(as.integer(pos)),
          query_root_id = np$int64(batch$query_root_id[i]),
          match_id = as.character(batch$match_id[i]),
          score = as.numeric(batch$score[i]),
          validation = isTRUE(batch$validation[i])
        )
        n_staged <- n_staged + 1
      }

      if (n_staged == 0) next

      # Upload with retry: on failure, wait 5s and try once more
      ok <- tryCatch({
        client$annotation$upload_staged_annotations(stage)
        TRUE
      }, error = function(e) {
        warning(sprintf("    Batch %d/%d failed: %s — retrying in 5s ...",
                        b, n_batches, e$message))
        Sys.sleep(5)
        tryCatch({
          client$annotation$upload_staged_annotations(stage)
          TRUE
        }, error = function(e2) {
          warning(sprintf("    Batch %d/%d retry failed: %s — skipping",
                          b, n_batches, e2$message))
          FALSE
        })
      })
      stage$clear_annotations()

      if (ok) {
        n_uploaded <- n_uploaded + n_staged
      } else {
        n_failed <- n_failed + n_staged
      }

      if (b %% 10 == 0 || b == n_batches) {
        message(sprintf("    Batch %d/%d — uploaded: %d, failed: %d",
                        b, n_batches, n_uploaded, n_failed))
      }
    }
  }

  message(sprintf("    Done: %d deleted, %d uploaded, %d failed",
                  n_deleted, n_uploaded, n_failed))
}

###########################
### Read metadata       ###
###########################

banc.meta <- readr::read_csv(
  file = file.path(banc.meta.save.path, "banc_ids.csv"),
  col_types = banc.col.types,
  show_col_types = FALSE
)

# CAVE client
client <- bancr::banc_cave_client()

###########################
### Sync each dataset   ###
###########################

# Dataset definitions: table_name -> feather filename
datasets <- list(
  list(table = "banc_malecns_nblast_v2", feather = "banc_malecns_v0.9_nblast.feather"),
  list(table = "banc_fafb_nblast_v2",    feather = "banc_fafb_783_nblast.feather"),
  list(table = "banc_hemibrain_nblast_v2", feather = "banc_hemibrain_v1.2.1_nblast.feather"),
  list(table = "banc_manc_nblast_v2",    feather = "banc_manc_v1.2.1_nblast.feather"),
  list(table = "banc_fanc_nblast_v2",    feather = "banc_fanc_1116_nblast.feather")
)

for (ds in datasets) {
  feather_path <- file.path(banc.meta.save.path, ds$feather)
  if (!file.exists(feather_path)) {
    message(sprintf("  Skipping %s — feather not found: %s", ds$table, feather_path))
    next
  }

  tryCatch(
    sync_nblast_cave_table(feather_path, ds$table, banc.meta, client),
    error = function(e) {
      message(sprintf("  ERROR syncing %s: %s", ds$table, e$message))
    }
  )
}

message(sprintf("### banc: CAVE sync complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
