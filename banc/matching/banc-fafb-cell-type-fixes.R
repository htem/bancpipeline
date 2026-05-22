#' banc-fafb-cell-type-fixes — Audit fafb_match quality via NBLAST top-10.
#'
#' For every BANC neuron with a `fafb_match`, check whether that match
#' appears in the top-10 NBLAST hits. Auto-resolves rows where a top hit
#' shares the same `fafb_cell_type` (swap fafb_match, keep type) and flags
#' the rest for manual review with neuroglancer links.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - `banc_fafb_783_nblast.feather`
#'
#' @section Writes:
#'   - `data/codex/fafb_match_fixes.csv` — auto-resolvable swaps
#'   - `data/codex/fafb_match_mismatches.csv` — manual-review queue

###########################################################
### Audit fafb_match quality using NBLAST scores
###
### For every BANC neuron with a fafb_match, check whether
### that match appears in the top 10 NBLAST hits. If not,
### flag it as a potential error.
###
### Resolvable: the neuron has a top NBLAST hit with the
### same fafb_cell_type → swap fafb_match, keep type.
###
### Mismatches: the top NBLAST hit has a different type
### → output for manual review with neuroglancer links.
###
### Output:
###   data/codex/fafb_match_fixes.csv       (auto-resolvable)
###   data/codex/fafb_match_mismatches.csv   (manual review)
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: FAFB match audit via NBLAST ###")

###########################
### Read data           ###
###########################

bc <- banctable_query(
  "SELECT _id, root_888, root_id, cell_type, cell_class, fafb_cell_type, fafb_match, hemibrain_cell_type, side, cell_type_source, status FROM banc_meta"
) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

# Neurons with a fafb_match
has_match <- bc %>%
  dplyr::filter(!is.na(fafb_match), fafb_match != "")
message(sprintf("  BANC neurons with fafb_match: %d", nrow(has_match)))

# Load franken_meta
message("  Loading franken_meta...")
fm <- franken_meta()
fafb_meta <- fm %>%
  dplyr::filter(!is.na(fafb_id), fafb_id != "") %>%
  dplyr::select(fafb_id, fm_cell_type = cell_type, fm_side = side) %>%
  dplyr::mutate(fafb_id = as.character(fafb_id)) %>%
  dplyr::distinct(fafb_id, .keep_all = TRUE)

# Build type -> fafb_ids lookup
fafb_type_ids <- split(fafb_meta$fafb_id, fafb_meta$fm_cell_type)

# Load NBLAST scores
message("  Loading FAFB NBLAST feather...")
nblast_cache <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache, showWarnings = FALSE)

local_file <- file.path(nblast_cache, "banc_fafb_783_nblast.feather")
if (!file.exists(local_file)) {
  alt_file <- file.path("/tmp/nblast_cache", "banc_fafb_783_nblast.feather")
  if (file.exists(alt_file)) {
    local_file <- alt_file
  } else {
    message("    Downloading banc_fafb_783_nblast.feather...")
    system2("gsutil", c("cp",
      "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_fafb_783_nblast.feather",
      local_file), stdout = FALSE, stderr = FALSE)
  }
}
fafb_nblast <- arrow::read_feather(local_file)
message(sprintf("  FAFB NBLAST: %d rows", nrow(fafb_nblast)))

# Pre-index NBLAST by root_888, keeping top 20 per neuron (sorted by score desc)
message("  Pre-indexing NBLAST (top 20 per neuron)...")
nblast_top <- fafb_nblast %>%
  dplyr::mutate(match_id = as.character(match_id), root_888 = as.character(root_888)) %>%
  dplyr::group_by(root_888) %>%
  dplyr::slice_max(order_by = score, n = 20, with_ties = FALSE) %>%
  dplyr::ungroup()
nblast_idx <- split(
  data.frame(match_id = nblast_top$match_id, score = nblast_top$score,
             match_cell_type = as.character(nblast_top$match_cell_type),
             stringsAsFactors = FALSE),
  nblast_top$root_888
)
rm(nblast_top); gc()
message(sprintf("  Index: %d neurons", length(nblast_idx)))

###########################
### Audit               ###
###########################

message("  Auditing fafb_match quality...")

results <- has_match %>%
  dplyr::select(root_888, `_id`, root_id, cell_type, fafb_cell_type,
                fafb_match, hemibrain_cell_type, cell_type_source, side)
results$fafb_match_rank <- NA_integer_
results$fafb_match_score <- NA_real_
results$top1_match_id <- NA_character_
results$top1_cell_type <- NA_character_
results$top1_score <- NA_real_
results$status <- NA_character_

pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)
for (i in seq_len(nrow(results))) {
  setTxtProgressBar(pb, i)
  rid <- results$root_888[i]
  current_match <- results$fafb_match[i]

  hits <- nblast_idx[[rid]]
  if (is.null(hits) || nrow(hits) == 0) {
    results$status[i] <- "no_nblast"
    next
  }

  # Record top hit
  results$top1_match_id[i] <- hits$match_id[1]
  results$top1_cell_type[i] <- hits$match_cell_type[1]
  results$top1_score[i] <- hits$score[1]

  # Find rank of current fafb_match
  rank_pos <- match(current_match, hits$match_id)
  if (!is.na(rank_pos)) {
    results$fafb_match_rank[i] <- rank_pos
    results$fafb_match_score[i] <- hits$score[rank_pos]
    results$status[i] <- if (rank_pos <= 10) "ok" else "low_rank"
  } else {
    results$status[i] <- "not_in_top20"
  }
}
close(pb)

###########################
### Classify issues     ###
###########################

tab <- table(results$status)
message("\n  === Audit results ===")
for (s in names(tab)) message(sprintf("    %s: %d", s, tab[s]))

# Flagged neurons: not in top 10
flagged <- results %>%
  dplyr::filter(status %in% c("low_rank", "not_in_top20"))
message(sprintf("\n  Flagged (fafb_match not in top 10): %d", nrow(flagged)))

if (nrow(flagged) == 0) {
  message("  No issues found. Exiting.")
  return(invisible())
}

# Classify flagged neurons:
# 1. "resolvable" — a top NBLAST hit has the same fafb_cell_type → swap match
# 2. "mismatch" — top NBLAST hit has a different type → manual review

flagged$resolution <- NA_character_
flagged$new_fafb_match <- NA_character_
flagged$new_fafb_cell_type <- NA_character_
flagged$suggested_type <- NA_character_

for (i in seq_len(nrow(flagged))) {
  rid <- flagged$root_888[i]
  current_type <- flagged$fafb_cell_type[i]

  hits <- nblast_idx[[rid]]
  if (is.null(hits) || nrow(hits) == 0) {
    flagged$resolution[i] <- "no_nblast"
    next
  }

  # Can we find a top-10 hit with the same fafb_cell_type?
  if (!is.na(current_type) && current_type != "") {
    same_type_ids <- fafb_type_ids[[current_type]]
    if (!is.null(same_type_ids) && length(same_type_ids) > 0) {
      top10 <- hits[1:min(10, nrow(hits)), ]
      same_type_in_top10 <- top10[top10$match_id %in% same_type_ids, ]
      if (nrow(same_type_in_top10) > 0) {
        flagged$resolution[i] <- "resolvable"
        flagged$new_fafb_match[i] <- same_type_in_top10$match_id[1]
        flagged$new_fafb_cell_type[i] <- current_type
        next
      }
    }
  }

  # Not resolvable — suggest the top NBLAST hit
  # Look up cell_type from franken_meta (not NBLAST match_cell_type, which is hemibrain_type)
  flagged$resolution[i] <- "mismatch"
  top_id <- hits$match_id[1]
  flagged$new_fafb_match[i] <- top_id
  fm_row <- fafb_meta[fafb_meta$fafb_id == top_id, ]
  flagged$suggested_type[i] <- if (nrow(fm_row) > 0) fm_row$fm_cell_type[1] else hits$match_cell_type[1]
}

res_tab <- table(flagged$resolution)
message("\n  === Resolution ===")
for (s in names(res_tab)) message(sprintf("    %s: %d", s, res_tab[s]))

###########################
### Resolvable fixes    ###
###########################

resolvable <- flagged %>% dplyr::filter(resolution == "resolvable")
if (nrow(resolvable) > 0) {
  fix_file <- "data/codex/fafb_match_fixes.csv"
  readr::write_csv(resolvable, fix_file)
  message(sprintf("\n  Saved %d resolvable fixes to %s", nrow(resolvable), fix_file))

  # SeaTable update: swap fafb_match to the better NBLAST hit (same type)
  # Re-read from CSV so this block can be run standalone in RStudio
  resolvable <- readr::read_csv(fix_file, col_types = readr::cols(.default = "c"),
                                 show_col_types = FALSE)
  push_resolvable <- resolvable %>%
    dplyr::transmute(
      `_id`,
      fafb_match = new_fafb_match,
      cell_type_source = append_status(cell_type_source, "bates")
    ) %>%
    as.data.frame()
  push_resolvable$fafb_match[is.na(push_resolvable$fafb_match)] <- ""

  message(sprintf("  SeaTable push (resolvable): %d rows, %d with new fafb_match",
                  nrow(push_resolvable), sum(push_resolvable$fafb_match != "")))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_resolvable,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete (resolvable)")
}

###########################
### Mismatches for review ##
###########################

mismatches <- flagged %>% dplyr::filter(resolution == "mismatch")
if (nrow(mismatches) > 0) {
  message(sprintf("\n  Building neuroglancer links for %d mismatches...", nrow(mismatches)))

  # Decode base neuroglancer scene
  ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4899366325190656"
  ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
  ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
  ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                      return = "text", cache = TRUE)
  ngl_base <- fafbseg::ngl_decode_scene(
    fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

  ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))

  banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)
  if (is.na(banc_layer_idx)) banc_layer_idx <- grep("banc|segmentation", ngl_ls$name, ignore.case = TRUE)[1]

  # FAFB layers: "seatable type" = current match, "connectivity type" = suggested match, "top nblast" = top nblast
  fafb_st_idx <- grep("seatable.type", ngl_ls$name, ignore.case = TRUE)
  fafb_st_idx <- if (length(fafb_st_idx) >= 1) fafb_st_idx[1] else NA_integer_
  fafb_conn_idx <- grep("connectivity.type", ngl_ls$name, ignore.case = TRUE)
  fafb_conn_idx <- if (length(fafb_conn_idx) >= 1) fafb_conn_idx[1] else NA_integer_
  fafb_nbl_idx <- grep("nblast.type", ngl_ls$name, ignore.case = TRUE)
  fafb_nbl_idx <- if (length(fafb_nbl_idx) >= 1) fafb_nbl_idx[1] else NA_integer_

  if (is.na(fafb_st_idx) || is.na(fafb_conn_idx)) {
    fafb_layer_idxs <- grep("fafb|flywire", ngl_ls$name, ignore.case = TRUE)
    if (is.na(fafb_st_idx)) fafb_st_idx <- if (length(fafb_layer_idxs) >= 1) fafb_layer_idxs[1] else NA_integer_
    if (is.na(fafb_conn_idx)) fafb_conn_idx <- if (length(fafb_layer_idxs) >= 2) fafb_layer_idxs[2] else fafb_st_idx
  }

  ngl_urls <- character(nrow(mismatches))
  pb <- txtProgressBar(min = 0, max = nrow(mismatches), style = 3)
  for (i in seq_len(nrow(mismatches))) {
    setTxtProgressBar(pb, i)
    row <- mismatches[i, ]
    tryCatch({
      sc <- ngl_base

      # BANC neuron
      banc_rid <- if (!is.na(row$root_id) && row$root_id != "") row$root_id else row$root_888
      if (!is.na(banc_layer_idx)) {
        sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
        sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
      }

      # Current fafb_match (seatable type layer)
      if (!is.na(fafb_st_idx)) {
        sc[["layers"]][[fafb_st_idx]][["segments"]] <- as.character(row$fafb_match)
        sc[["layers"]][[fafb_st_idx]][["hiddenSegments"]] <- NULL
      }

      # Suggested new match (connectivity type layer)
      if (!is.na(fafb_conn_idx) && !is.na(row$new_fafb_match)) {
        sc[["layers"]][[fafb_conn_idx]][["segments"]] <- as.character(row$new_fafb_match)
        sc[["layers"]][[fafb_conn_idx]][["hiddenSegments"]] <- NULL
      }

      ngl_urls[i] <- as.character(sc)
    }, error = function(e) {
      ngl_urls[i] <<- NA_character_
    })
  }
  close(pb)

  mismatches$neuroglancer_url <- ngl_urls
  mismatches$accept_new <- NA_character_

  out <- mismatches %>%
    dplyr::arrange(dplyr::desc(top1_score)) %>%
    dplyr::select(`_id`, root_888, root_id, cell_type, fafb_cell_type, fafb_match,
                  fafb_match_rank, fafb_match_score,
                  suggested_match = new_fafb_match, suggested_type,
                  top1_score, hemibrain_cell_type, cell_type_source, side,
                  neuroglancer_url, accept_new)

  mismatch_file <- "data/codex/fafb_match_mismatches.csv"
  readr::write_csv(out, mismatch_file)
  message(sprintf("  Saved %d mismatches to %s", nrow(out), mismatch_file))
}

###########################
### Summary             ###
###########################

message("\n  === Final summary ===")
message(sprintf("  Total audited: %d", nrow(results)))
message(sprintf("  OK (fafb_match in top 10): %d", sum(results$status == "ok")))
message(sprintf("  No NBLAST data: %d", sum(results$status == "no_nblast")))
message(sprintf("  Flagged: %d", nrow(flagged)))
message(sprintf("    Resolvable (same type, better match): %d", sum(flagged$resolution == "resolvable", na.rm = TRUE)))
message(sprintf("    Mismatch (different type, needs review): %d", sum(flagged$resolution == "mismatch", na.rm = TRUE)))

###############################################################
### SeaTable update: reviewed mismatches                    ###
### accept_new == "T": accept suggested match + type        ###
### accept_new == "A": mark wrong, clear fafb_match/type    ###
### Both T and A: append FAFB_PNG_MATCH_WRONG to status     ###
###############################################################

reviewed_csv_path <- "data/codex/fafb_match_mismatches_reviewed.csv"
if (file.exists(reviewed_csv_path)) {
  message("\n=== SeaTable update: reviewed FAFB match mismatches ===")

  reviewed <- readr::read_csv(reviewed_csv_path,
                               col_types = readr::cols(.default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new %in% c("T", "A"))

  message(sprintf("  %d reviewed rows (T=%d, A=%d)",
                  nrow(reviewed),
                  sum(reviewed$accept_new == "T"),
                  sum(reviewed$accept_new == "A")))

  if (nrow(reviewed) > 0) {
    # Join to seatable for _id and current values if _id missing
    if (!"_id" %in% names(reviewed) || all(is.na(reviewed$`_id`))) {
      reviewed <- reviewed %>%
        dplyr::select(-dplyr::any_of("_id")) %>%
        dplyr::left_join(
          bc %>% dplyr::select(root_888, `_id`, current_status = status,
                                current_cell_type = cell_type,
                                current_fafb_cell_type = fafb_cell_type,
                                current_cell_type_source = cell_type_source),
          by = "root_888"
        )
    } else {
      reviewed <- reviewed %>%
        dplyr::left_join(
          bc %>% dplyr::select(root_888, current_status = status,
                                current_cell_type = cell_type,
                                current_fafb_cell_type = fafb_cell_type,
                                current_cell_type_source = cell_type_source),
          by = "root_888"
        )
    }
    reviewed <- reviewed %>% dplyr::filter(!is.na(`_id`), `_id` != "")

    # Look up the cell_type of the suggested FAFB match from franken_meta
    reviewed <- reviewed %>%
      dplyr::left_join(
        fafb_meta %>% dplyr::select(fafb_id, new_fafb_cell_type = fm_cell_type),
        by = c("suggested_match" = "fafb_id")
      )

    # Build push dataframe
    push_reviewed <- reviewed %>%
      dplyr::transmute(
        `_id`,
        # T: accept suggested match, set fafb_cell_type from franken_meta
        # A: clear both
        fafb_match = dplyr::case_when(
          accept_new == "T" ~ suggested_match,
          accept_new == "A" ~ "",
          TRUE ~ NA_character_
        ),
        fafb_cell_type = dplyr::case_when(
          accept_new == "T" ~ new_fafb_cell_type,
          accept_new == "A" ~ "",
          TRUE ~ NA_character_
        ),
        # T: if old cell_type == old fafb_cell_type, replace cell_type with new fafb_cell_type
        # A: if old cell_type == old fafb_cell_type, clear cell_type
        cell_type = dplyr::case_when(
          accept_new == "T" & !is.na(current_cell_type) & !is.na(current_fafb_cell_type) &
            current_cell_type == current_fafb_cell_type ~ new_fafb_cell_type,
          accept_new == "A" & !is.na(current_cell_type) & !is.na(current_fafb_cell_type) &
            current_cell_type == current_fafb_cell_type ~ "",
          TRUE ~ NA_character_
        ),
        status = {
          s <- append_status(current_status, "FAFB_PNG_MATCH_WRONG")
          su <- status_update
          su[is.na(su)] <- ""
          ifelse(su != "", append_status(s, su), s)
        },
        cell_type_source = append_status(current_cell_type_source, "bates")
      ) %>%
      as.data.frame()

    # Drop columns that are NA (no change) — keep only columns we actually set
    # For cell_type, only include rows where we're clearing it
    if (all(is.na(push_reviewed$cell_type))) {
      push_reviewed$cell_type <- NULL
    } else {
      # Only set cell_type for rows where we're actually clearing it
      push_reviewed$cell_type[is.na(push_reviewed$cell_type)] <- ""
    }

    push_reviewed$fafb_match[is.na(push_reviewed$fafb_match)] <- ""
    push_reviewed$fafb_cell_type[is.na(push_reviewed$fafb_cell_type)] <- ""

    n_accept <- sum(reviewed$accept_new == "T")
    n_reject <- sum(reviewed$accept_new == "A")
    n_ct_cleared <- if ("cell_type" %in% names(push_reviewed)) sum(push_reviewed$cell_type == "", na.rm = TRUE) else 0

    message(sprintf("  Pushing %d rows:", nrow(push_reviewed)))
    message(sprintf("    T (accept new match): %d", n_accept))
    message(sprintf("    A (clear fafb match): %d", n_reject))
    message(sprintf("    cell_type cleared (was == fafb_cell_type): %d", n_ct_cleared))
    message(sprintf("    status: FAFB_PNG_MATCH_WRONG appended to all"))

    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = push_reviewed,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  SeaTable update complete (reviewed mismatches)")
  }
} else {
  message("\n  Skipping reviewed mismatch update: ", reviewed_csv_path, " not found")
}

})
