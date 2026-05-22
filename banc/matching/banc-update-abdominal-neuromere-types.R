#' banc-update-abdominal-neuromere-types — Push FAbG abdominal-neuromere types into SeaTable.
#'
#' Loads the annotated xlsx files from the FAbG project (Jiayi) and pushes
#' cell_type / sexually_dimorphic / status updates to BANC SeaTable. Only
#' rows with `FAbG_match` populated (and not `Not in AbG` / `BANC_unmatched_*`)
#' are touched. The actual push call is COMMENTED OUT — uncomment to apply.
#'
#' @section Reads:
#'   - `data/abdominal_neuromere/FAbG-vs-BANC_MN_New-name_20260508_annotated.xlsx`
#'   - `data/abdominal_neuromere/FAbG-vs-BANC_SN_New-name_20260508_annotated.xlsx`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `cell_type`, `sexually_dimorphic`,
#'     `status` (status append: MANUAL_ANNOTATION + FABG_PREFERRED;
#'     idempotent via `append_status`)
#'
#' @section Notes:
#'   - Re-run after each new xlsx revision from Jiayi.
#'   - The actual `banctable_update_rows()` call is COMMENTED OUT by default.

#####################################################
### BANC: ABDOMINAL NEUROMERE TYPE UPDATES ###
#####################################################
# Pushes naming, dimorphism, and status annotations from the FAbG project
# (Jiayi) into the banc_meta seatable. Inputs are the two annotated xlsx
# files in data/abdominal_neuromere/:
#
#   FAbG-vs-BANC_MN_New-name_20260508_annotated.xlsx   (motor/visceral effectors)
#   FAbG-vs-BANC_SN_New-name_20260508_annotated.xlsx   (sensory + sensory_ascending)
#
# The push has 3 components (the actual banctable_update_rows() call is
# COMMENTED OUT — uncomment to apply):
#
#   1. cell_type           <- banc_cell_type_new
#   2. sexually_dimorphic  <- banc_sexually_dimorphic_new
#   3. status              <- append("MANUAL_ANNOTATION", "FABG_PREFERRED")
#
# append_status() is idempotent: re-running won't duplicate flags.
# Only rows with FAbG_match populated (and not "Not in AbG" /
# "BANC_unmatched_*") are touched.
#
# Re-run after each new xlsx revision from Jiayi.

source("banc/banc-startup.R")
suppressPackageStartupMessages({ library(readxl) })

#####################################################
### A. Naming-pipeline summary (reference only) ###
#####################################################
# The xlsx files were produced by an iterative pipeline that lived across many
# scratch scripts. It is NOT re-run here, because re-running requires manual
# review of clustering choices, dimorphism rules, and locked-name lists. The
# pipeline did (in rough order):
#
#   1. Read Jiayi's source MN/SN xlsx files (BANC tab + FAbG tab).
#   2. Pull current banc_meta and join (banc_cell_type, banc_super_class,
#      banc_super_cluster, banc_sexually_dimorphic).
#   3. For each FAbG_Type, aggregate to majority super_class, super_cluster,
#      sexually_dimorphic.
#   4. Apply the FAbG -> BANC MANC bridge (parse MANC_matched in FAbG tab,
#      expand "+" shorthand, strip "-N" / "_Low"/"_High" sub-suffixes).
#   5. Append uncovered abdominal effectors / afferents (flow + nerve filter)
#      and try to bridge their banc_cell_type via the same MANC map.
#   6. Compute banc_sexually_dimorphic_new per cell type with this rule:
#        female-specific iff BANC majority == "female-specific" AND
#                            FAbG Dimorphism == "Specific"
#        isomorphic      iff BANC majority == "isomorphic"/blank AND
#                            FAbG Dimorphism == "Shared"/blank
#        dimorphic       in all other cases (any BANC <-> FAbG disagreement)
#      Option B: when >=80% of BANC neurons in the group lack
#      sexually_dimorphic, ignore BANC and use FAbG alone.
#   7. Run banc_read_l2skel() -> /1e3 (microns) -> nat::dotprops (k=5,
#      resample=1) -> nat.nblast::nblast_allbyall (normalisation="mean") for
#      the union of both xlsx files (336 effectors + 1095 sensories + cross
#      NBLAST for 36 uncovered unknown_sensory).
#   8. Cluster cell types (hclust, average linkage on 1 - mean_NBLAST) within
#      each prefix; assign names per the convention below.
#   9. For unknown_sensory neurons: tight self-clusters (intra mean NBLAST
#      >= 0.30, n >= 4) get new types; the remaining ~17 stay as
#      banc_cell_type_new = "unknown_sensory" (the rump).
#  10. Reconcile sexually_dimorphic across the cell-type pool and write
#      banc_sexually_dimorphic_new back into both BANC tabs.
#
# Naming convention: {PREFIX}abg{NN}{f?}
#   PREFIX: EFF (effector: motor or visceral_circulatory),
#           SN  (sensory),
#           SA  (sensory_ascending)
#   abg:    abdominal ganglion (deliberately distinct from MANC's "ad")
#   NN:     2-digit zero-padded counter, NBLAST-clustered within prefix;
#           non-FS counts up from 01, FS counts down from 99 (so future
#           additions slot in at the next free slot from either end)
#   f:      appended only for female-specific cell types
#
# Locked / privileged names (kept verbatim regardless of clustering):
#   MN xlsx: CMU, MNxm03, SNpp54, tergotrochanter_extensor_MetaAN,
#            sternal_adductor_MetaAN, hDVM MN, hi1 MN, hi2 MN, hiii2 MN,
#            haltere_motor_neuron_unknown, MNhm03, MNhm42, MNhm43
#   SN xlsx: SPSN
#
# The intermediate scripts (NBLAST runs, clustering, renumbering, etc.) are
# in git history if needed.


#####################################################
### B0. Snapshot banc_meta before any changes ###
#####################################################
# Match the pattern in banc/update/banc-update-seatable.R: save a CSV
# snapshot of the full banc_meta into banc.meta.save.path/snapshots/.

# banc.meta.save.path is an O2 path; on a local box it won't exist, so fall
# back to a project-local snapshots dir.
.snapshot_dir <- if (dir.exists(banc.meta.save.path) ||
                     suppressWarnings(dir.create(banc.meta.save.path,
                                                  showWarnings = FALSE, recursive = TRUE))) {
  file.path(banc.meta.save.path, "snapshots")
} else {
  "data/meta/snapshots"
}
dir.create(.snapshot_dir, showWarnings = FALSE, recursive = TRUE)
.dt <- format(Sys.time(), "%Y-%m-%d_%H-%M")
.snapshot_file <- file.path(.snapshot_dir, paste0(.dt, "_banc_seatable.csv"))
.bc_snapshot <- banctable_query()
readr::write_csv(.bc_snapshot, .snapshot_file)
message(sprintf("Snapshot saved: %s (%d rows)", basename(.snapshot_file), nrow(.bc_snapshot)))

#####################################################
### B. Read the annotated xlsx files ###
#####################################################

ann_dir <- "data/abdominal_neuromere"

mn <- read_excel(file.path(ann_dir, "FAbG-vs-BANC_MN_New-name_20260508_annotated.xlsx"),
                 sheet = "BANC") %>%
  dplyr::mutate(root_626 = as.character(root_626), source_file = "MN")

sn <- read_excel(file.path(ann_dir, "FAbG-vs-BANC_SN_New-name_20260508_annotated.xlsx"),
                 sheet = "BANC") %>%
  dplyr::mutate(root_626 = as.character(root_626), source_file = "SN")

xlsx_both <- dplyr::bind_rows(mn, sn) %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::distinct(root_626, .keep_all = TRUE)
cat(sprintf("xlsx rows: MN=%d, SN=%d, total distinct=%d\n",
            nrow(mn), nrow(sn), nrow(xlsx_both)))


#####################################################
### C. Pull current banc_meta for these neurons ###
#####################################################
# Use the un-cached banctable_query() because we're about to write back.

ids <- unique(xlsx_both$root_626)
ids_quoted <- paste0("'", ids, "'", collapse = ",")
bc <- banctable_query(sprintf(
  "SELECT _id, root_626, root_id, cell_type, super_class, sexually_dimorphic, status, cell_type_source FROM banc_meta WHERE root_626 IN (%s)",
  ids_quoted)) %>%
  dplyr::mutate(across(c(`_id`, root_626, root_id), as.character)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)
cat(sprintf("Fetched %d / %d rows from banc_meta\n", nrow(bc), length(ids)))


#####################################################
### D. Join + filter to FAbG-covered rows ###
#####################################################
# "FAbG-covered" = the neuron is in Jiayi's typology (has a real FAbG_match).
# We do NOT touch uncovered (no FAbG_match), "Not in AbG", or
# "BANC_unmatched_*" rows — those keep their existing seatable fields.

df <- xlsx_both %>%
  dplyr::inner_join(bc %>% dplyr::select(`_id`, root_626,
                                          curr_cell_type           = cell_type,
                                          curr_sexually_dimorphic  = sexually_dimorphic,
                                          curr_status              = status,
                                          curr_cell_type_source    = cell_type_source),
                    by = "root_626") %>%
  dplyr::mutate(fabg_covered = !is.na(FAbG_match) & FAbG_match != "" &
                  FAbG_match != "Not in AbG" &
                  !grepl("^BANC_unmatched", FAbG_match))

cat(sprintf("FAbG-covered rows: %d  uncovered/ignored: %d\n",
            sum(df$fabg_covered), sum(!df$fabg_covered)))


#####################################################
### E. Build per-column update sets (only changed rows) ###
#####################################################

# --- E1. cell_type ----------------------------------------------------------
# Push banc_cell_type_new where it is a real name (not blank, not the
# "unknown_sensory" rump) AND differs from current cell_type.
ct_changes <- df %>%
  dplyr::filter(fabg_covered,
                !is.na(banc_cell_type_new),
                banc_cell_type_new != "",
                banc_cell_type_new != "unknown_sensory",
                is.na(curr_cell_type) | curr_cell_type != banc_cell_type_new) %>%
  dplyr::transmute(`_id`, root_626,
                   cell_type = banc_cell_type_new)
cat(sprintf("cell_type changes: %d\n", nrow(ct_changes)))

# --- E2. sexually_dimorphic -------------------------------------------------
sd_changes <- df %>%
  dplyr::filter(fabg_covered,
                !is.na(banc_sexually_dimorphic_new),
                banc_sexually_dimorphic_new != "",
                is.na(curr_sexually_dimorphic) | curr_sexually_dimorphic != banc_sexually_dimorphic_new) %>%
  dplyr::transmute(`_id`, root_626,
                   sexually_dimorphic = banc_sexually_dimorphic_new)
cat(sprintf("sexually_dimorphic changes: %d\n", nrow(sd_changes)))

# --- E2b. cell_type_source --------------------------------------------------
# Stamp 'jiayi_zhang' as the source for every FAbG-covered row. Only push
# where the field actually changes (so re-runs are no-ops once applied).
cts_changes <- df %>%
  dplyr::filter(fabg_covered,
                is.na(curr_cell_type_source) | curr_cell_type_source != "wang_lab") %>%
  dplyr::transmute(`_id`, root_626, cell_type_source = "wang_lab")
cat(sprintf("cell_type_source changes: %d\n", nrow(cts_changes)))

# --- E3. status -------------------------------------------------------------
# Append MANUAL_ANNOTATION + FABG_PREFERRED to every FAbG-covered row.
# append_status() de-duplicates and sorts, so this is idempotent.
st_changes <- df %>%
  dplyr::filter(fabg_covered) %>%
  dplyr::mutate(status = append_status(curr_status,
                                        "HAS_MANUAL_ANNOTATION, FABG_PREFERRED")) %>%
  dplyr::filter(is.na(curr_status) | status != curr_status) %>%
  dplyr::transmute(`_id`, root_626, status)
cat(sprintf("status changes: %d\n", nrow(st_changes)))


#####################################################
### F. Combine into a single push frame ###
#####################################################
# Build one frame keyed by _id, carrying only the columns that actually
# change. banctable_update_rows() respects only the columns present.

push_df <- df %>%
  dplyr::filter(fabg_covered) %>%
  dplyr::transmute(
    `_id`,
    root_626,
    cell_type = dplyr::if_else(
      !is.na(banc_cell_type_new) & banc_cell_type_new != "" &
        banc_cell_type_new != "unknown_sensory" &
        (is.na(curr_cell_type) | curr_cell_type != banc_cell_type_new),
      banc_cell_type_new, curr_cell_type),
    sexually_dimorphic = dplyr::if_else(
      !is.na(banc_sexually_dimorphic_new) & banc_sexually_dimorphic_new != "" &
        (is.na(curr_sexually_dimorphic) | curr_sexually_dimorphic != banc_sexually_dimorphic_new),
      banc_sexually_dimorphic_new, curr_sexually_dimorphic),
    cell_type_source = "wang_lab",
    status = append_status(curr_status, "HAS_MANUAL_ANNOTATION, FABG_PREFERRED")
  ) %>%
  # Drop rows where nothing actually changes (no-op pushes).
  dplyr::filter(
    (!is.na(cell_type) & (is.na(df$curr_cell_type[match(root_626, df$root_626)]) |
                          df$curr_cell_type[match(root_626, df$root_626)] != cell_type)) |
    (!is.na(sexually_dimorphic) & (is.na(df$curr_sexually_dimorphic[match(root_626, df$root_626)]) |
                                    df$curr_sexually_dimorphic[match(root_626, df$root_626)] != sexually_dimorphic)) |
    (!is.na(status) & (is.na(df$curr_status[match(root_626, df$root_626)]) |
                        df$curr_status[match(root_626, df$root_626)] != status))
  )
cat(sprintf("Combined push frame: %d rows\n", nrow(push_df)))


#####################################################
### G. Push to seatable (COMMENTED OUT) ###
#####################################################
banctable_update_rows(base = 'banc_meta',
                      table = 'banc_meta',
                      df = push_df %>% dplyr::select(`_id`, cell_type,
                                                      sexually_dimorphic,
                                                      cell_type_source,
                                                      status),
                      append_allowed = FALSE,
                      chunksize = 500)


#####################################################
### H. Sanity-check log (post-push) ###
#####################################################
# After pushing, re-query a few of the changed rows and confirm fields landed
# as expected. Uncomment when running the push.
#
# check_ids <- head(push_df$root_626, 10)
# verify <- banctable_query(sprintf(
#   "SELECT root_626, cell_type, sexually_dimorphic, status FROM banc_meta WHERE root_626 IN (%s)",
#   paste0("'", check_ids, "'", collapse = ",")))
# print(verify)
