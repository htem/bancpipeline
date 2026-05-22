#' banc-update-status — Update status flags + side assignments in SeaTable.
#'
#' Computes status flags (LR_TYPE_CONFLICT, SIDE_CONFLICT, UNROOTED,
#' TOO_SMALL, FAFB_MATCH_WRONG_SIDE, MANC_MATCH_WRONG_SIDE, NO_*_MATCH,
#' TRACING_ISSUE_RESOLVED, CONFLICT_RESOLVED, *_TYPE_CONFLICT) and side
#' from `root_position_nm`. Pushes status + side to SeaTable independently.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `status`, `side`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`
#'
#' @section Notes:
#'   - RETIRED 2026-05-15: `*_ALT_MATCH` flags. The block below strips any
#'     residual ALT_MATCH tags as a one-time migration.

###########################################################
### Update status flags and side assignments in BANC
### seatable
###
### Status flags managed:
###   LR_TYPE_CONFLICT, SIDE_CONFLICT, UNROOTED, TOO_SMALL,
###   FAFB_MATCH_WRONG_SIDE, MANC_MATCH_WRONG_SIDE,
###   NO_*_MATCH, TRACING_ISSUE_RESOLVED,
###   CONFLICT_RESOLVED, *_TYPE_CONFLICT
###
### RETIRED 2026-05-15: *_ALT_MATCH flags. Previously set whenever
### `cell_type != <dataset>_cell_type` to flag "match-and-cell_type
### disagree" collisions. Now that banc-nblast-compile.R picks one
### most-recent match per neuron (rank_png_matches) and banc-update-
### celltypes.R cascades the new cell_type when the existing value
### came from this dataset's prior match (status-guarded), the
### ALT_MATCH flags are redundant. The block below strips any
### residual ALT_MATCH tags from existing status values as a
### one-time migration.
###
### Also computes side from root_position_nm.
###
### Pushes status + side to seatable independently
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: updating status flags and sides ###")
t_start <- Sys.time()

###########################
### Read current state  ###
###########################

bc <- banctable_query()
bc <- bc %>%
  dplyr::filter(!is.na(root_id), root_id != "0") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Save original values for anti-join (only push rows that actually changed)
bc.orig <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(`_id`, status, side) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~dplyr::coalesce(as.character(.), "")))

# Read reference metadata for side comparisons
franken.meta <- franken_meta()
fw.meta <- franken.meta %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::rename(root_783 = fafb_id)
mc.meta <- franken.meta %>%
  dplyr::filter(!is.na(manc_id)) %>%
  dplyr::rename(bodyid = manc_id)

message("  Stripping retired *_ALT_MATCH tags (see header)...")
###########################
### ALT_MATCH retire    ###
###########################
# One-time migration: remove any existing *_ALT_MATCH tags from status.
# Safe to keep this idempotently — subtract_status is a no-op when the
# tag is absent. mcns.meta / fc.meta are still loaded below because
# later blocks (FAFB_MATCH_WRONG_SIDE, MANC_MATCH_WRONG_SIDE, etc.)
# use them.
bc <- bc %>%
  dplyr::mutate(status = subtract_status(status, "FAFB_ALT_MATCH")) %>%
  dplyr::mutate(status = subtract_status(status, "MANC_ALT_MATCH")) %>%
  dplyr::mutate(status = subtract_status(status, "MALECNS_ALT_MATCH")) %>%
  dplyr::mutate(status = subtract_status(status, "FANC_ALT_MATCH")) %>%
  dplyr::mutate(status = subtract_status(status, "HEMIBRAIN_ALT_MATCH"))

# Reference metadata still needed downstream (MATCH_WRONG_SIDE etc.)
mcns.meta <- readr::read_csv(file.path(banc.meta.save.path, "malecns_09_meta.csv"),
                              col_types = banc.col.types) %>%
  dplyr::filter(!is.na(malecns_09_id)) %>%
  dplyr::rename(bodyid = malecns_09_id)
fc.meta <- suppressWarnings(readr::read_csv(file.path(banc.meta.save.path, "fanc_meta.csv"),
                                             col_types = hemibrainr:::sql_col_types))

message("  Computing status flags (UNROOTED, TOO_SMALL, SIDE_CONFLICT, TYPE_CONFLICT, etc.)...")
###########################
### Status wrangling    ###
###########################

bc <- bc %>%
  # Fix root_position_nm from nucleus if missing
  dplyr::mutate(root_position_nm = dplyr::case_when(
    is.na(root_position_nm) | root_position_nm == "" ~ nucleus_position_nm,
    TRUE ~ root_position_nm
  )) %>%
  # LR_TYPE_CONFLICT: mirror match has different cell type (vectorized)
  dplyr::mutate(
    .match_ct = cell_type[match(banc_match, root_id)],
    status = dplyr::case_when(
      is.na(banc_match) | banc_match %in% c("NA", "", " ") ~ status,
      is.na(cell_type) | is.na(.match_ct) ~ status,
      cell_type != .match_ct ~ append_status(status, "LR_TYPE_CONFLICT"),
      cell_type == .match_ct ~ subtract_status(status, "LR_TYPE_CONFLICT"),
      TRUE ~ status
    )
  ) %>%
  dplyr::select(-.match_ct) %>%
  # SIDE_CONFLICT: mirror match on same side (vectorized)
  dplyr::mutate(
    .match_side = side[match(banc_match, root_id)],
    status = dplyr::case_when(
      is.na(banc_match) | banc_match %in% c("NA", "", " ") ~ status,
      side == .match_side ~ append_status(status, "SIDE_CONFLICT"),
      side != .match_side ~ subtract_status(status, "SIDE_CONFLICT"),
      TRUE ~ status
    )
  ) %>%
  dplyr::select(-.match_side) %>%
  # UNROOTED: no root position
  dplyr::mutate(status = dplyr::case_when(
    is.na(root_position_nm) | root_position_nm == "" ~ append_status(status, "UNROOTED"),
    !is.na(root_position_nm) & root_position_nm != "" & grepl("UNROOTED", status) ~
      subtract_status(status, "UNROOTED"),
    TRUE ~ status
  )) %>%
  # Clean up resolved issues
  dplyr::mutate(status = dplyr::case_when(
    grepl("TRACING_ISSUE_RESOLVED", status) ~ subtract_status(status, "TRACING_ISSUE"),
    TRUE ~ status
  )) %>%
  # Remove NO_*_MATCH when match exists
  dplyr::mutate(status = dplyr::case_when(
    !is.na(banc_match) & banc_match != '' ~ subtract_status(status, "NO_MIRROR_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(fafb_match) & fafb_match != '' ~ subtract_status(status, "NO_FAFB_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(fanc_match) & fanc_match != '' ~ subtract_status(status, "NO_FANC_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(manc_match) & manc_match != '' ~ subtract_status(status, "NO_MANC_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(hemibrain_match) & hemibrain_match != '' ~ subtract_status(status, "NO_HEMIBRAIN_MATCH"),
    TRUE ~ status
  )) %>%
  # Clean up resolved conflicts
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED", status) ~ subtract_status(status, "FAFB_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED", status) ~ subtract_status(status, "MANC_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED", status) ~ subtract_status(status, "HEMIBRAIN_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED", status) ~ subtract_status(status, "FANC_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  # TOO_SMALL: cable length <= 50 um
  dplyr::mutate(status = dplyr::case_when(
    grepl("NOT_TOO_SMALL|TOO_SMALL_RESOLVED", status) ~ subtract_status(status, "TOO_SMALL"),
    is.na(l2_cable_length_um) ~ status,
    l2_cable_length_um > 50 ~ subtract_status(status, "TOO_SMALL"),
    l2_cable_length_um <= 50 ~ append_status(status, "TOO_SMALL"),
    TRUE ~ status
  )) %>%
  # FAFB_MATCH_WRONG_SIDE (vectorized)
  dplyr::mutate(
    .fafb_side = fw.meta$side[match(fafb_match, fw.meta$root_783)],
    status = dplyr::case_when(
      is.na(fafb_match) | is.na(side) ~ status,
      grepl("ascending|sensory_ascending|ascending, sensory neuron", cell_class) ~
        subtract_status(status, "FAFB_MATCH_WRONG_SIDE"),
      side == .fafb_side ~ subtract_status(status, "FAFB_MATCH_WRONG_SIDE"),
      side != .fafb_side ~ append_status(status, "FAFB_MATCH_WRONG_SIDE"),
      TRUE ~ status
    )
  ) %>%
  dplyr::select(-.fafb_side) %>%
  # MANC_MATCH_WRONG_SIDE (vectorized)
  dplyr::mutate(
    .manc_side = mc.meta$side[match(manc_match, mc.meta$bodyid)],
    status = dplyr::case_when(
      is.na(manc_match) | is.na(side) ~ status,
      grepl("descending", cell_class) ~ subtract_status(status, "MANC_MATCH_WRONG_SIDE"),
      side == .manc_side ~ subtract_status(status, "MANC_MATCH_WRONG_SIDE"),
      side != .manc_side ~ append_status(status, "MANC_MATCH_WRONG_SIDE"),
      TRUE ~ status
    )
  ) %>%
  dplyr::select(-.manc_side)

message("  Updating side assignments...")
###########################
### Update sides        ###
###########################

# Compute side for neurons that have a root position but no side
root_no_side <- (!is.na(bc$root_position_nm) & bc$root_position_nm != "") &
  (bc$side == "" | is.na(bc$side))
root_no_side[is.na(root_no_side)] <- FALSE
if (sum(root_no_side)) {
  roots <- nat::xyzmatrix(bc$root_position_nm[root_no_side])
  lrdiffs <- bancr:::banc_lr_position(roots, units = "nm")
  sides <- ifelse(lrdiffs > 0, "right", "left")
  bc$side[root_no_side] <- sides
  message(sprintf("  Computed side for %d neurons", sum(root_no_side)))
}

# Safety: preserve existing side values from seatable.  If a neuron had a
# side in the original query but our bc dataframe now has it blank (e.g. due
# to NA propagation), restore the original value so we never erase sides.
bc_side_now_blank <- (is.na(bc$side) | bc$side == "")
bc_side_had_value <- bc.orig$`_id`[bc.orig$side != ""]
restore_mask <- bc_side_now_blank & bc$`_id` %in% bc_side_had_value
if (sum(restore_mask, na.rm = TRUE)) {
  orig_lookup <- setNames(bc.orig$side, bc.orig$`_id`)
  bc$side[restore_mask] <- orig_lookup[bc$`_id`[restore_mask]]
  message(sprintf("  Preserved %d existing side values that would have been blanked",
                  sum(restore_mask, na.rm = TRUE)))
}

###########################
### Clean and push      ###
###########################

bc$status <- gsub("NA\\,|\\,NA", "", bc$status)

status.update <- bc %>%
  dplyr::filter(!is.na(`_id`), `_id` != "") %>%
  dplyr::select(`_id`, status, side) %>%
  as.data.frame()
status.update[is.na(status.update)] <- ""

# Never push a blank side when seatable already has a value
status.update <- status.update %>%
  dplyr::left_join(bc.orig %>% dplyr::select(`_id`, orig_side = side), by = "_id") %>%
  dplyr::mutate(side = dplyr::case_when(
    side == "" & orig_side != "" ~ orig_side,
    TRUE ~ side
  )) %>%
  dplyr::select(-orig_side)

# Anti-join: only push rows where status or side actually changed
status.update <- dplyr::anti_join(status.update, bc.orig, by = c("_id", "status", "side"))
message(sprintf("Pushing status + side for %d changed neurons to seatable (skipped %d unchanged)",
                nrow(status.update), nrow(bc.orig) - nrow(status.update)))
if (nrow(status.update)) {
  banctable_update_rows(base = 'banc_meta',
                        table = 'banc_meta',
                        df = status.update,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

message(sprintf("### banc: status + side update complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
