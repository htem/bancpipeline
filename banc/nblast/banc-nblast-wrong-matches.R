#' banc-nblast-wrong-matches — Process `*_PNG_MATCH_WRONG` SeaTable flags into all artefacts.
#'
#' Invalidates bad matches in compiled feathers + reviewed CSVs + meta CSV,
#' removes PNGs, clears SeaTable flag. Must run BEFORE `banc-nblast-compile.R`.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, compiled NBLAST feathers, reviewed CSVs, PNG folders
#'
#' @section Writes:
#'   - SeaTable `banc_meta` (clears `*_PNG_MATCH_WRONG`)
#'   - drops invalidated entries from feathers / CSVs / PNG folders
#'
#' @section Invoked by:
#'   production v888 rebuild chain.

###########################################################
### Handle wrong PNG match deletion for all datasets
###
### For each dataset (MANC, maleCNS, FAFB, hemibrain, BANC mirror):
###   1. Read *_PNG_MATCH_WRONG status from seatable
###   2. Invalidate in feather + reviewed CSV + banc_meta.csv
###   3. Remove wrong PNG files from review folders
###   4. Sync between storage paths (A ↔ B)
###   5. Update seatable to clear the wrong matches
###
### This must run before banc-nblast-compile.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: handling wrong PNG matches ###")
t_start <- Sys.time()

# Read latest seatable data
bancr:::banctable_updateids()
bc <- banctable_query()
conflicts <- data.frame()

# Storage paths
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'

###########################
### Helper function     ###
###########################

# Generic function to handle wrong matches for a given dataset.
# Pass cell_type_col / auto_cell_type_col = NULL for same-dataset (mirror) matches
# where there is no separate auto-derived cell type column to clear.
handle_wrong_matches <- function(bc, dataset,
                                 status_flag,
                                 match_col, png_match_col,
                                 cell_type_col = NULL,
                                 feather_file, reviewed_csv_file,
                                 auto_cell_type_col = NULL,
                                 correct_match_path,
                                 match_subdir,
                                 conflicts) {

  # Read wrong matches — use exact match to avoid substring collisions
  # (e.g. MANC_PNG_MATCH_WRONG matching inside MALECNS_PNG_MATCH_WRONG)
  has_flag <- function(status, flag) {
    grepl(paste0("(^|,)\\s*", flag, "\\s*(,|$)"), status)
  }

  bc.flagged <- bc %>%
    dplyr::distinct(`_id`, .keep_all = TRUE) %>%
    dplyr::filter(has_flag(status, status_flag),
                  !is.na(root_id))

  # Split: neurons WITH png_match can do full processing;
  # those WITHOUT just need the orphaned status tag removed
  select_cols <- c("_id", "status", "root_id",
                   match_col, png_match_col, "cell_type")
  if (!is.null(cell_type_col)) select_cols <- c(select_cols, cell_type_col)
  select_cols <- c(select_cols, "valid")

  bc.png.wrong <- bc.flagged %>%
    dplyr::filter(!is.na(.data[[png_match_col]]),
                  .data[[png_match_col]] != "") %>%
    dplyr::mutate(valid = 'f') %>%
    dplyr::select(dplyr::all_of(select_cols)) %>%
    plyr::rbind.fill(conflicts) %>%
    dplyr::filter(!is.na(.data[[png_match_col]]),
                  .data[[png_match_col]] != "",
                  !is.na(root_id)) %>%
    as.data.frame()

  bc.orphaned <- bc.flagged %>%
    dplyr::filter(is.na(.data[[png_match_col]]) | .data[[png_match_col]] == "") %>%
    dplyr::select(`_id`, status) %>%
    dplyr::mutate(status = subtract_status(status, status_flag)) %>%
    as.data.frame()

  if (nrow(bc.orphaned)) {
    message(sprintf("  Cleaning %d orphaned %s status tags (no png_match set)",
                    nrow(bc.orphaned), dataset))
    banctable_update_rows(base = 'banc_meta',
                          table = "banc_meta",
                          df = bc.orphaned,
                          append_allowed = FALSE,
                          chunksize = 1000)
  }

  if (!nrow(bc.png.wrong)) {
    message(sprintf("  No wrong %s matches to process", dataset))
    return(invisible())
  }

  message(sprintf("  Processing %d wrong %s matches", nrow(bc.png.wrong), dataset))

  # Invalidate in feather file
  nb.feather <- arrow::read_feather(file.path(banc.meta.save.path, feather_file)) %>%
    dplyr::anti_join(bc.png.wrong,
                     by = stats::setNames(c(png_match_col, "root_id"),
                                          c("match_id", "pt_root_id")))

  # Invalidate in reviewed CSV
  matches.df.valid <- readr::read_csv(file = file.path(banc.meta.save.path, reviewed_csv_file),
                                      col_types = banc.col.types,
                                      show_col_types = FALSE) %>%
    dplyr::left_join(bc.png.wrong %>%
                       dplyr::select(dplyr::all_of(c("root_id", png_match_col, "valid"))),
                     by = stats::setNames(c(png_match_col, "root_id"),
                                          c("match_id", "pt_root_id"))) %>%
    dplyr::mutate(valid = dplyr::case_when(
      !is.na(valid.y) ~ valid.y,
      TRUE ~ valid.x
    )) %>%
    dplyr::select(-valid.x, -valid.y)

  # Invalidate auto cell type in banc_meta.csv (skipped for mirror matches
  # where no separate auto cell type column exists)
  banc.meta <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_meta.csv"),
                                col_types = banc.col.types,
                                show_col_types = FALSE)
  if (!is.null(auto_cell_type_col)) {
    banc.meta <- banc.meta %>%
      dplyr::mutate(!!auto_cell_type_col := ifelse(
        root_id %in% bc.png.wrong$root_id, NA, .data[[auto_cell_type_col]]
      )) %>%
      dplyr::mutate(!!auto_cell_type_col := ifelse(
        supervoxel_id %in% bc.png.wrong$supervoxel_id, NA, .data[[auto_cell_type_col]]
      ))
  }

  # Save updated files
  arrow::write_feather(nb.feather, file.path(banc.meta.save.path, feather_file))
  readr::write_csv(matches.df.valid, file = file.path(banc.meta.save.path, reviewed_csv_file))
  readr::write_csv(banc.meta, file = file.path(banc.meta.save.path, "banc_meta.csv"))

  # Remove saved PNG files
  sync_files(path(B, paste0("matching/", match_subdir, "/correct/")),
             path(A, paste0("matching/", match_subdir, "/correct/")),
             extensions = c("png"), move.old = FALSE)
  inventory <- list.files(correct_match_path, recursive = TRUE, full.names = TRUE)
  for (i in seq_len(nrow(bc.png.wrong))) {
    query <- bc.png.wrong[i, "root_id"]
    match <- bc.png.wrong[i, png_match_col]
    poss <- inventory[grepl(match, inventory)]
    if (!length(poss)) next
    # Keep poss aligned with extracted ids: regmatches drops un-matched files
    # so the original code's `poss[del]` could index the wrong PNGs.
    m <- regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE)
    poss <- poss[m != -1]
    if (!length(poss)) next
    ids <- regmatches(basename(poss), regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE))
    ids <- banc_latestid(ids)
    del <- which(ids == query)
    if (length(del)) {
      file.remove(poss[del])
      message("  deleting: ", paste(basename(poss[del]), collapse = ", "))
    }
  }
  sync_files(path(A, paste0("matching/", match_subdir, "/correct/")),
             path(B, paste0("matching/", match_subdir, "/correct/")),
             extensions = c("png"), move.old = FALSE)

  # Update seatable: clear the wrong match columns
  bc.new <- bc.png.wrong %>%
    dplyr::rowwise() %>%
    {
      if (!is.null(cell_type_col)) {
        dplyr::mutate(.,
          cell_type = dplyr::case_when(
            is.na(.data[[cell_type_col]]) | .data[[cell_type_col]] == "" ~ cell_type,
            is.na(cell_type) | cell_type == "" ~ cell_type,
            is.na(.data[[png_match_col]]) | .data[[png_match_col]] == "" ~ cell_type,
            is.na(.data[[match_col]]) | .data[[match_col]] == "" ~ cell_type,
            .data[[match_col]] != .data[[png_match_col]] ~ cell_type,
            (cell_type == .data[[cell_type_col]]) & (.data[[match_col]] == .data[[png_match_col]]) ~ "",
            TRUE ~ cell_type
          ),
          !!cell_type_col := dplyr::case_when(
            is.na(.data[[cell_type_col]]) | .data[[cell_type_col]] == "" ~ .data[[cell_type_col]],
            is.na(.data[[png_match_col]]) | .data[[png_match_col]] == "" ~ .data[[cell_type_col]],
            is.na(.data[[match_col]]) | .data[[match_col]] == "" ~ "",
            .data[[match_col]] != .data[[png_match_col]] ~ .data[[cell_type_col]],
            (.data[[match_col]] == .data[[png_match_col]]) ~ "",
            TRUE ~ .data[[cell_type_col]]
          )
        )
      } else .
    } %>%
    dplyr::mutate(
      !!match_col := dplyr::case_when(
        is.na(.data[[png_match_col]]) | .data[[png_match_col]] == "" ~ .data[[match_col]],
        is.na(.data[[match_col]]) | .data[[match_col]] == "" ~ .data[[match_col]],
        (.data[[match_col]] == .data[[png_match_col]]) ~ "",
        TRUE ~ .data[[match_col]]
      )
    ) %>%
    dplyr::mutate(!!png_match_col := '') %>%
    dplyr::mutate(status = subtract_status(status, status_flag)) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(
      c("_id", match_col, png_match_col, "cell_type", cell_type_col, "status")
    ))) %>%
    as.data.frame()

  banctable_update_rows(base = 'banc_meta',
                        table = "banc_meta",
                        df = bc.new,
                        append_allowed = FALSE,
                        chunksize = 1000)

  message(sprintf("  Done processing %s wrong matches", dataset))
}

###########################
### MANC                ###
###########################

handle_wrong_matches(
  bc = bc, dataset = "MANC",
  status_flag = "MANC_PNG_MATCH_WRONG",
  match_col = "manc_match", png_match_col = "manc_png_match",
  cell_type_col = "manc_cell_type",
  feather_file = "banc_manc_v1.2.1_nblast.feather",
  reviewed_csv_file = "banc_manc_reviewed_matches.csv",
  auto_cell_type_col = "manc_auto_cell_type",
  correct_match_path = banc.manc.correct.match.path,
  match_subdir = "manc",
  conflicts = conflicts
)

###########################
### maleCNS             ###
###########################

handle_wrong_matches(
  bc = bc, dataset = "maleCNS",
  status_flag = "MALECNS_PNG_MATCH_WRONG",
  match_col = "malecns_match", png_match_col = "malecns_png_match",
  cell_type_col = "malecns_cell_type",
  feather_file = "banc_malecns_v0.9_nblast.feather",
  reviewed_csv_file = "banc_malecns_reviewed_matches.csv",
  auto_cell_type_col = "malecns_auto_cell_type",
  correct_match_path = banc.malecns.correct.match.path,
  match_subdir = "malecns",
  conflicts = conflicts
)

###########################
### FAFB                ###
###########################

handle_wrong_matches(
  bc = bc, dataset = "FAFB",
  status_flag = "FAFB_PNG_MATCH_WRONG",
  match_col = "fafb_match", png_match_col = "fafb_png_match",
  cell_type_col = "fafb_cell_type",
  feather_file = "banc_fafb_783_nblast.feather",
  reviewed_csv_file = "banc_fafb_reviewed_matches.csv",
  auto_cell_type_col = "fafb_auto_cell_type",
  correct_match_path = banc.fafb.correct.match.path,
  match_subdir = "fafb",
  conflicts = conflicts
)

###########################
### hemibrain           ###
###########################

handle_wrong_matches(
  bc = bc, dataset = "hemibrain",
  status_flag = "HEMIBRAIN_PNG_MATCH_WRONG",
  match_col = "hemibrain_match", png_match_col = "hemibrain_png_match",
  cell_type_col = "hemibrain_cell_type",
  feather_file = "banc_hemibrain_v1.2.1_nblast.feather",
  reviewed_csv_file = "banc_hemibrain_reviewed_matches.csv",
  auto_cell_type_col = "hemibrain_auto_cell_type",
  correct_match_path = banc.hemibrain.correct.match.path,
  match_subdir = "hemibrain",
  conflicts = conflicts
)

###########################
### BANC (mirror)       ###
###########################
# Same-dataset left-right mirror match. No separate auto cell type column —
# cell_type is the canonical type and must not be cleared by a wrong mirror match.

handle_wrong_matches(
  bc = bc, dataset = "BANC",
  status_flag = "BANC_PNG_MATCH_WRONG",
  match_col = "banc_match", png_match_col = "banc_png_match",
  cell_type_col = NULL,
  feather_file = "banc_mirror_nblast.feather",
  reviewed_csv_file = "banc_mirror_reviewed_matches.csv",
  auto_cell_type_col = NULL,
  correct_match_path = banc.mirror.correct.match.path,
  match_subdir = "mirror",
  conflicts = conflicts
)

message(sprintf("### banc: wrong PNG match handling complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
