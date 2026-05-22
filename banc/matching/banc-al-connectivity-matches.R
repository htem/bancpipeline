#' banc-al-connectivity-matches — Detect BANC↔FAFB type mismatches for AL neurons.
#'
#' Reads Codex BANC↔FAFB connectivity alignment scores for the antennal-lobe
#' cutout, joins to SeaTable to compare against `fafb_cell_type`, emits a
#' review CSV with per-row neuroglancer links. Companion: `banc-al-review-images.R`
#' renders PNGs.
#'
#' @section Reads:
#'   - `data/codex/al/banc_fafb_min_syn_1_cutout_antennal_lobe_..._alignment_scores.csv.gz`
#'   - SeaTable `banc_meta`
#'   - `franken_meta()` (FAFB cell_type lookup)
#'
#' @section Writes:
#'   - `data/codex/al_connectivity_mismatches.csv`

###########################################################
### Compare FAFB connectivity matches vs seatable types
### for antennal lobe (AL) neurons
###
### Reads codex connectivity alignment scores (BANC↔FAFB),
### joins to seatable to get fafb_cell_type, compares
### against the FAFB cell_type of the matched neuron.
### Outputs a CSV of mismatches with neuroglancer links.
###
### Input:  data/codex/al/banc_fafb_min_syn_1_cutout_antennal_lobe_...csv.gz
### Output: data/codex/al_connectivity_mismatches.csv
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: AL connectivity match analysis (BANC↔FAFB) ###")

###########################
### Read data           ###
###########################

# Read connectivity alignment scores (no header)
# Columns: banc_id, fafb_id, score, n_synapses
conn <- readr::read_csv(
  "data/codex/al/banc_fafb_min_syn_1_cutout_antennal_lobe_norm_lr_gdt_node_alignment_scores.csv.gz",
  col_names = c("banc_id", "fafb_id", "score", "n_synapses"),
  col_types = readr::cols(banc_id = "c", fafb_id = "c", score = "d", n_synapses = "d"),
  show_col_types = FALSE
)

# Remove dummy/test rows (IDs < 15 digits are synthetic)
conn <- conn %>% dplyr::filter(nchar(banc_id) > 15, nchar(fafb_id) > 15)
message(sprintf("  Loaded %d connectivity matches", nrow(conn)))

# Keep best match per BANC neuron (highest score)
conn <- conn %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::distinct(banc_id, .keep_all = TRUE)
message(sprintf("  %d unique BANC neurons after dedup", nrow(conn)))

# Query BANC seatable
bc <- banctable_query(
  "SELECT _id, root_id, root_888, supervoxel_id, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, hemibrain_cell_type, fafb_match, malecns_cell_type, malecns_match, side, status FROM banc_meta"
) %>%
  dplyr::filter(!is.na(root_888), root_888 != "") %>%
  dplyr::mutate(
    root_888 = as.character(root_888),
    root_id = as.character(root_id),
    supervoxel_id = as.character(supervoxel_id),
    fafb_match = as.character(fafb_match)
  ) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)
message(sprintf("  SeaTable: %d neurons with root_888", nrow(bc)))

# Read FAFB metadata via franken_meta
message("  Loading FAFB metadata from franken_meta()...")
fm <- franken_meta()
fafb_meta <- fm %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::select(fafb_root_id = fafb_id, fafb_cell_type_meta = cell_type,
                fafb_side = side, fafb_super_class = super_class) %>%
  dplyr::mutate(fafb_root_id = as.character(fafb_root_id))
message(sprintf("  FAFB meta: %d neurons", nrow(fafb_meta)))

###########################
### Join and compare    ###
###########################

# Join connectivity matches to BANC seatable (banc_id = root_888)
df <- conn %>%
  dplyr::inner_join(bc, by = c("banc_id" = "root_888"))
message(sprintf("  %d matches joined to seatable", nrow(df)))

# Join FAFB metadata to get the connectivity-matched FAFB neuron's cell_type
df <- df %>%
  dplyr::left_join(
    fafb_meta %>% dplyr::distinct(fafb_root_id, .keep_all = TRUE),
    by = c("fafb_id" = "fafb_root_id")
  )

# The "connectivity type" is the cell_type of the matched FAFB neuron
df <- df %>%
  dplyr::rename(connectivity_type = fafb_cell_type_meta)

# Determine the "seatable type" — fafb_cell_type from BANC seatable
df <- df %>%
  dplyr::mutate(seatable_type = fafb_cell_type)

# Find mismatches: connectivity_type vs seatable fafb_cell_type
# A mismatch is where both exist and differ, OR connectivity type exists but seatable doesn't
df <- df %>%
  dplyr::mutate(
    has_conn_type = !is.na(connectivity_type) & connectivity_type != "",
    has_st_type = !is.na(seatable_type) & seatable_type != "",
    is_mismatch = dplyr::case_when(
      has_conn_type & has_st_type ~ connectivity_type != seatable_type,
      has_conn_type & !has_st_type ~ TRUE,   # new type from connectivity
      !has_conn_type & has_st_type ~ FALSE,   # connectivity match has no type — not useful
      TRUE ~ FALSE
    )
  )

mismatches <- df %>% dplyr::filter(is_mismatch)
message(sprintf("  Mismatches: %d / %d (%.1f%%)",
                nrow(mismatches), nrow(df), 100 * nrow(mismatches) / nrow(df)))

if (nrow(mismatches) == 0) {
  message("  No mismatches found. Nothing to write.")
  return(invisible())
}

###############################################
### FAFB neuron lookup for neuroglancer     ###
###############################################

# For visualising the seatable type (fafb_cell_type) in neuroglancer:
#   1. If fafb_match exists, use that specific FAFB neuron
#   2. Else, find FAFB neurons with same side + cell_type
#   3. If fafb_cell_type is blank, fall back to hemibrain_cell_type
#      and find any FAFB neuron with that type + same side → mark hemibrain_transfer

# Build FAFB lookup table: cell_type + side → root_ids
fafb_by_type_side <- fafb_meta %>%
  dplyr::filter(!is.na(fafb_cell_type_meta), fafb_cell_type_meta != "") %>%
  dplyr::group_by(fafb_cell_type_meta, fafb_side) %>%
  dplyr::summarise(fafb_ids = list(fafb_root_id), .groups = "drop")

get_fafb_ids_for_seatable_type <- function(row_fafb_match, row_seatable_type,
                                            row_hb_type, row_side) {
  hemibrain_transfer <- FALSE
  type_to_use <- row_seatable_type

  # If fafb_match exists, use it directly
  if (!is.na(row_fafb_match) && row_fafb_match != "") {
    return(list(ids = row_fafb_match, hemibrain_transfer = FALSE))
  }

  # If seatable fafb_cell_type is blank, try hemibrain_cell_type
  if (is.na(type_to_use) || type_to_use == "") {
    if (!is.na(row_hb_type) && row_hb_type != "") {
      type_to_use <- row_hb_type
      hemibrain_transfer <- TRUE
    } else {
      return(list(ids = character(0), hemibrain_transfer = FALSE))
    }
  }

  # Look up FAFB neurons by type + side
  side_to_use <- if (!is.na(row_side) && row_side != "") row_side else NA
  matches <- fafb_by_type_side %>%
    dplyr::filter(fafb_cell_type_meta == type_to_use)
  if (!is.na(side_to_use)) {
    side_matches <- matches %>% dplyr::filter(fafb_side == side_to_use)
    if (nrow(side_matches) > 0) matches <- side_matches
  }
  ids <- unlist(matches$fafb_ids)
  if (is.null(ids)) ids <- character(0)
  list(ids = ids, hemibrain_transfer = hemibrain_transfer)
}

###########################################
### Neuroglancer URLs                   ###
###########################################

message("  Building neuroglancer links...")

# Decode the base neuroglancer scene
ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4899366325190656"
ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                    return = "text", cache = TRUE)
ngl_base <- fafbseg::ngl_decode_scene(
  fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

# Find layer indices
ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))
message(sprintf("  Scene layers: %s", paste(ngl_ls$name, collapse = ", ")))

# BANC segmentation layer
banc_layer_idx <- match("v850 neurons", ngl_ls$name)
if (is.na(banc_layer_idx)) banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)
if (is.na(banc_layer_idx)) banc_layer_idx <- grep("banc|segmentation", ngl_ls$name, ignore.case = TRUE)[1]

# FAFB layers: match by name to avoid positional mix-ups
# "connectivity type" layer shows the connectivity-matched FAFB neuron
# "seatable type" layer shows the existing morphology-based FAFB match
fafb_conn_idx <- grep("connectivity.type", ngl_ls$name, ignore.case = TRUE)
fafb_conn_idx <- if (length(fafb_conn_idx) >= 1) fafb_conn_idx[1] else NA_integer_
fafb_st_idx <- grep("seatable.type", ngl_ls$name, ignore.case = TRUE)
fafb_st_idx <- if (length(fafb_st_idx) >= 1) fafb_st_idx[1] else NA_integer_
# Fallback: positional if name matching fails
if (is.na(fafb_conn_idx) || is.na(fafb_st_idx)) {
  fafb_layer_idxs <- grep("fafb|flywire", ngl_ls$name, ignore.case = TRUE)
  if (is.na(fafb_conn_idx)) fafb_conn_idx <- if (length(fafb_layer_idxs) >= 1) fafb_layer_idxs[1] else NA_integer_
  if (is.na(fafb_st_idx)) fafb_st_idx <- if (length(fafb_layer_idxs) >= 2) fafb_layer_idxs[2] else fafb_conn_idx
}

message(sprintf("  BANC layer: '%s' [%d]",
                if (!is.na(banc_layer_idx)) ngl_ls$name[banc_layer_idx] else "NOT FOUND",
                banc_layer_idx))
message(sprintf("  FAFB layers: conn='%s' [%d], seatable='%s' [%d]",
                if (!is.na(fafb_conn_idx)) ngl_ls$name[fafb_conn_idx] else "NOT FOUND", fafb_conn_idx,
                if (!is.na(fafb_st_idx)) ngl_ls$name[fafb_st_idx] else "NOT FOUND", fafb_st_idx))

first_error <- NULL
ngl_urls <- character(nrow(mismatches))
hb_transfers <- logical(nrow(mismatches))

for (i in seq_len(nrow(mismatches))) {
  row <- mismatches[i, ]

  # Look up FAFB IDs for the seatable type
  st_result <- get_fafb_ids_for_seatable_type(
    row$fafb_match, row$seatable_type, row$hemibrain_cell_type, row$side
  )
  hb_transfers[i] <- st_result$hemibrain_transfer

  tryCatch({
    sc <- ngl_base

    # Set BANC neuron (use root_id for modern segmentation, fall back to banc_id=root_888)
    banc_rid <- if (!is.na(row$root_id) && row$root_id != "") row$root_id else row$banc_id
    if (!is.na(banc_layer_idx)) {
      sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
      sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
    }

    # Set connectivity-matched FAFB neuron
    if (!is.na(fafb_conn_idx)) {
      sc[["layers"]][[fafb_conn_idx]][["segments"]] <- as.character(row$fafb_id)
      sc[["layers"]][[fafb_conn_idx]][["hiddenSegments"]] <- NULL
    }

    # Set seatable-type FAFB neurons
    if (!is.na(fafb_st_idx) && fafb_st_idx != fafb_conn_idx && length(st_result$ids) > 0) {
      sc[["layers"]][[fafb_st_idx]][["segments"]] <- as.character(st_result$ids)
      sc[["layers"]][[fafb_st_idx]][["hiddenSegments"]] <- NULL
    }

    ngl_urls[i] <- as.character(sc)
  }, error = function(e) {
    if (is.null(first_error)) first_error <<- conditionMessage(e)
    ngl_urls[i] <<- NA_character_
  })
}
if (!is.null(first_error)) message(sprintf("  First NGL error: %s", first_error))
message(sprintf("  Generated %d/%d neuroglancer links",
                sum(!is.na(ngl_urls)), nrow(mismatches)))

mismatches$hemibrain_transfer <- hb_transfers
mismatches$neuroglancer_url <- ngl_urls

###########################
### Save CSV            ###
###########################

out <- mismatches %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::transmute(
    root_888 = banc_id,
    root_id,
    fafb_id,
    score,
    n_synapses,
    super_class,
    cell_class,
    cell_sub_class,
    side,
    seatable_type,
    connectivity_type,
    hemibrain_cell_type,
    hemibrain_transfer,
    fafb_match,
    neuroglancer_url
  )

out_file <- "data/codex/al_connectivity_mismatches.csv"
readr::write_csv(out, out_file)
message(sprintf("  Saved %d mismatches to %s", nrow(out), out_file))

# Summary
message("\n  === Summary ===")
message(sprintf("  Total connectivity matches: %d", nrow(df)))
message(sprintf("  Mismatches: %d", nrow(out)))
message(sprintf("  With seatable type: %d", sum(out$seatable_type != "" & !is.na(out$seatable_type))))
message(sprintf("  Using hemibrain_cell_type fallback: %d", sum(out$hemibrain_transfer)))
message(sprintf("  New types (no seatable type): %d",
                sum(is.na(out$seatable_type) | out$seatable_type == "")))

###############################################################
### SeaTable update: reviewed AL type changes                ###
### accept_new == "T": update fafb_match, fafb_cell_type,    ###
### and cell_type from the connectivity-matched FAFB neuron. ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################

reviewed_csv_path <- "data/codex/al_connectivity_mismatches_reviewed.csv"
if (file.exists(reviewed_csv_path)) {
  message("=== SeaTable update: reviewed AL type changes (accept_new == T) ===")

  # Read reviewed CSV and filter to accepted
  accepted <- readr::read_csv(reviewed_csv_path,
                               col_types = readr::cols(.default = "c"),
                               show_col_types = FALSE) %>%
    # Handle old reviewed CSVs that have root_626 instead of root_888
    {if ("root_626" %in% names(.) && !"root_888" %in% names(.)) dplyr::rename(., root_888 = root_626) else .} %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "T")

  message(sprintf("  %d accepted type changes", nrow(accepted)))

  # Join to seatable to get _id and current values
  accepted <- accepted %>%
    dplyr::left_join(bc %>% dplyr::select(root_888, `_id`,
                                           current_cell_type = cell_type,
                                           current_fafb_cell_type = fafb_cell_type,
                                           current_fafb_match = fafb_match,
                                           current_malecns_cell_type = malecns_cell_type),
                     by = "root_888") %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Look up cell_type of the connectivity-matched FAFB neuron from franken_meta
  accepted <- accepted %>%
    dplyr::left_join(
      fafb_meta %>% dplyr::distinct(fafb_root_id, .keep_all = TRUE) %>%
        dplyr::select(fafb_root_id, new_fafb_cell_type = fafb_cell_type_meta),
      by = c("fafb_id" = "fafb_root_id")
    )

  # Resolve the effective new fafb_cell_type for each row
  accepted <- accepted %>%
    dplyr::mutate(
      effective_fafb_cell_type = dplyr::if_else(!is.na(new_fafb_cell_type),
                                                 new_fafb_cell_type, connectivity_type)
    )

  # --- malecns_cell_type from maleCNS SeaTable (fafb_cell_type → cell_type mapping) ---
  mcns_type_map <- tryCatch({
    banctable_query("SELECT malecns_09_id, cell_type, fafb_cell_type FROM malecns",
                    base = "cns_meta") %>%
      dplyr::filter(!is.na(fafb_cell_type), fafb_cell_type != "",
                    !is.na(cell_type), cell_type != "") %>%
      dplyr::distinct(fafb_cell_type, .keep_all = TRUE) %>%
      dplyr::select(fafb_cell_type, malecns_ct = cell_type)
  }, error = function(e) {
    message("  Could not query maleCNS SeaTable: ", e$message)
    NULL
  })

  if (!is.null(mcns_type_map) && nrow(mcns_type_map) > 0) {
    accepted <- accepted %>%
      dplyr::left_join(mcns_type_map,
                       by = c("effective_fafb_cell_type" = "fafb_cell_type"))
    n_malecns <- sum(!is.na(accepted$malecns_ct))
    message(sprintf("  malecns_cell_type resolved for %d/%d neurons via maleCNS SeaTable",
                    n_malecns, nrow(accepted)))
  } else {
    accepted$malecns_ct <- NA_character_
    message("  maleCNS type mapping unavailable; malecns_cell_type will not be updated")
  }

  # Build push data frame
  push_df <- accepted %>%
    dplyr::transmute(
      `_id`,
      fafb_match = fafb_id,
      fafb_cell_type = effective_fafb_cell_type,
      cell_type = dplyr::case_when(
        !is.na(seatable_type) & !is.na(current_cell_type) &
          current_cell_type == seatable_type ~ effective_fafb_cell_type,
        is.na(current_cell_type) | current_cell_type == "" ~ effective_fafb_cell_type,
        TRUE ~ current_cell_type
      ),
      malecns_cell_type = dplyr::case_when(
        !is.na(malecns_ct) &
          (is.na(current_malecns_cell_type) | current_malecns_cell_type == "") ~ malecns_ct,
        !is.na(malecns_ct) ~ malecns_ct,
        TRUE ~ current_malecns_cell_type
      ),
      cell_type_source = "bates"
    ) %>%
    as.data.frame()

  push_df$fafb_match[is.na(push_df$fafb_match)] <- ""
  push_df$fafb_cell_type[is.na(push_df$fafb_cell_type)] <- ""
  push_df$cell_type[is.na(push_df$cell_type)] <- ""
  push_df$malecns_cell_type[is.na(push_df$malecns_cell_type)] <- ""

  n_ct_update <- sum(push_df$cell_type != "" &
    push_df$cell_type != accepted$current_cell_type, na.rm = TRUE)
  n_malecns_update <- sum(push_df$malecns_cell_type != "" &
    (is.na(accepted$current_malecns_cell_type) |
       accepted$current_malecns_cell_type == "" |
       push_df$malecns_cell_type != accepted$current_malecns_cell_type), na.rm = TRUE)

  message(sprintf("  Pushing update for %d neurons:", nrow(push_df)))
  message(sprintf("    fafb_match: %d non-empty", sum(push_df$fafb_match != "")))
  message(sprintf("    fafb_cell_type: %d updated", sum(push_df$fafb_cell_type != "")))
  message(sprintf("    cell_type: %d updated (where old == seatable_type or was blank)",
                  n_ct_update))
  message(sprintf("    malecns_cell_type: %d updated", n_malecns_update))
  message(sprintf("    cell_type_source: set to 'bates' for all %d rows", nrow(push_df)))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_df,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(accepted, push_df); gc()
} else {
  message("  Skipping SeaTable update: reviewed CSV not found at ", reviewed_csv_path)
}

###############################################################
### SeaTable update: non-mismatch fafb_match                ###
### For connectivity matches that are NOT mismatches,        ###
### set fafb_match to the best connectivity match FAFB ID.   ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################

message("=== SeaTable update: non-mismatch fafb_match ===")

non_mismatches <- df %>%
  dplyr::filter(!is_mismatch) %>%
  dplyr::filter(!is.na(`_id`), `_id` != "")

push_nonmismatch <- non_mismatches %>%
  dplyr::transmute(
    `_id`,
    root_id,
    fafb_match = as.character(fafb_id)
  ) %>%
  dplyr::filter(!is.na(fafb_match), fafb_match != "") %>%
  as.data.frame()

message(sprintf("  %d non-mismatch neurons: setting fafb_match", nrow(push_nonmismatch)))

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_nonmismatch,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  SeaTable update complete")

rm(non_mismatches, push_nonmismatch); gc()

})
