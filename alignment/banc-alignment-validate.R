#' banc-alignment-validate — Validate optic-lobe alignment results.
#'
#' Read-only validation: (1) holdout accuracy on intrinsic neurons with known
#' types, (2) neurotransmitter consistency vs FAFB type consensus, (3)
#' mismatch detection vs existing SeaTable annotations. Does not push to
#' SeaTable.
#'
#' @section Reads:
#'   - `data/optic_lobe/banc_optic_<side>_alignment_results.csv`
#'   - `data/optic_lobe/banc_optic_<side>_seeds.csv`
#'   - `data/optic_lobe/banc_optic_<side>_meta.csv`
#'   - SeaTable `banc_meta` (for mismatch detection)
#'
#' @section Writes:
#'   - `data/optic_lobe/banc_optic_<side>_holdout_accuracy.csv`
#'   - `data/optic_lobe/banc_optic_<side>_holdout_confusions.csv`
#'   - `data/optic_lobe/banc_optic_<side>_nt_mismatches.csv`
#'   - `data/optic_lobe/banc_optic_<side>_mismatches.csv`
#'
#' @section CLI:
#'   [side]   `right` (default) | `left` | `both`

###########################################################
### Optic Lobe Cross-Dataset Matching: Validation
###
### Validates alignment results against:
###   1. Holdout accuracy (intrinsic neurons with known type)
###   2. Neurotransmitter consistency
###   3. Mismatch detection vs existing SeaTable annotations
###
### NO SeaTable push — validation only.
###########################################################
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

local({

######################
### CONFIGURATION  ###
######################

args <- commandArgs(trailingOnly = TRUE)
side_to_process <- if (length(args) > 0) args[1] else "right"
data_dir <- "data/optic_lobe"
# Target dataset + region — fixed for the optic-lobe paper run. When wiring
# a non-FAFB target or a different region in, parameterise these. Mirror the
# parameter set the prep preset emitted so v2 filenames line up.
target_name    <- "fafb"
region_name    <- "optic-lobe"
nblast_version <- "783"
syn_source     <- Sys.getenv("BANC_SYN_SOURCE",
                              unset = banc.synapse.source.default)

######################
### LOAD DATA      ###
######################

message("=== Loading alignment results ===")

results <- readr::read_csv(
  alignment_path("align-results", query = "banc", target = target_name,
                 region = region_name, side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = data_dir),
  col_types = readr::cols(root_888 = "c", best_target_match = "c"),
  show_col_types = FALSE)

seeds <- readr::read_csv(
  alignment_path("prep-seeds", query = "banc", target = target_name,
                 region = region_name, side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = data_dir),
  col_types = readr::cols(root_888 = "c"),
  show_col_types = FALSE)

banc.meta <- readr::read_csv(
  alignment_path("prep-banc-meta", query = "banc", target = target_name,
                 region = region_name, side = side_to_process,
                 vq = banc.version, syn = syn_source,
                 ext = "csv", dir = data_dir),
  col_types = readr::cols(root_id = "c"),
  show_col_types = FALSE)

target.meta <- readr::read_csv(
  alignment_path("prep-target-meta", query = "banc", target = target_name,
                 region = region_name, side = side_to_process,
                 vq = banc.version, vt = nblast_version,
                 ext = "csv", dir = data_dir),
  col_types = readr::cols(target_id = "c"),
  show_col_types = FALSE)

message(sprintf("  Results: %d neurons, %d typed",
                nrow(results),
                sum(results$assigned_cell_type != "", na.rm = TRUE)))

######################################
### STEP 2.1: HOLDOUT ACCURACY     ###
######################################

message("\n=== Step 2.1: Holdout accuracy ===")

# Holdout = intrinsic neurons with known cell_type (is_holdout == TRUE)
holdout <- seeds %>%
  dplyr::filter(is_holdout == TRUE, !is.na(cell_type), cell_type != "") %>%
  dplyr::select(root_888, true_type = cell_type)

holdout_results <- results %>%
  dplyr::inner_join(holdout, by = "root_888") %>%
  dplyr::filter(assigned_cell_type != "")

if (nrow(holdout_results) > 0) {
  holdout_results$correct <- holdout_results$assigned_cell_type == holdout_results$true_type

  overall_acc <- mean(holdout_results$correct) * 100
  message(sprintf("  Overall holdout accuracy: %.1f%% (%d/%d)",
                  overall_acc,
                  sum(holdout_results$correct),
                  nrow(holdout_results)))

  # Per-type accuracy
  per_type <- holdout_results %>%
    dplyr::group_by(true_type) %>%
    dplyr::summarise(
      n = dplyr::n(),
      correct = sum(correct),
      accuracy = mean(correct) * 100,
      .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))

  message(sprintf("  Per-type accuracy (top 20 by count):"))
  top_types <- head(per_type, 20)
  for (i in seq_len(nrow(top_types))) {
    r <- top_types[i, ]
    message(sprintf("    %-20s %3d neurons, %5.1f%% correct", r$true_type, r$n, r$accuracy))
  }

  # Confusion: most common misassignments
  wrong <- holdout_results %>%
    dplyr::filter(!correct) %>%
    dplyr::count(true_type, assigned_cell_type, sort = TRUE)

  if (nrow(wrong) > 0) {
    message(sprintf("\n  Top 20 confusions:"))
    for (i in seq_len(min(20, nrow(wrong)))) {
      r <- wrong[i, ]
      message(sprintf("    %s -> %s (%d times)", r$true_type, r$assigned_cell_type, r$n))
    }
  }

  # Save
  readr::write_csv(per_type,
    alignment_path("validate-holdout-accuracy", query = "banc", target = target_name,
                   region = region_name, side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir))
  readr::write_csv(wrong,
    alignment_path("validate-holdout-confusions", query = "banc", target = target_name,
                   region = region_name, side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir))
} else {
  message("  No holdout neurons evaluated (none assigned a type)")
}

######################################
### STEP 2.2: NT CONSISTENCY CHECK ###
######################################

message("\n=== Step 2.2: Neurotransmitter consistency ===")

# Get NT for BANC neurons
banc_nt <- banc.meta %>%
  dplyr::select(root_888 = root_id, banc_nt = top_nt, super_class, cell_type) %>%
  dplyr::distinct(root_888, .keep_all = TRUE)

# Override: photoreceptors get histamine
banc_nt$banc_nt[grepl("photoreceptor", banc_nt$cell_type, ignore.case = TRUE) |
                  banc_nt$cell_type %in% c("R1-6", "R7", "R8")] <- "histamine"

# Get consensus NT per target-dataset type
target_type_nt <- target.meta %>%
  dplyr::mutate(
    # Override: photoreceptors get histamine
    target_top_nt = dplyr::if_else(
      grepl("photoreceptor", target_cell_class, ignore.case = TRUE) |
        target_cell_type %in% c("R1-6", "R7", "R8"),
      "histamine", target_top_nt)) %>%
  dplyr::filter(!is.na(target_cell_type), target_cell_type != "",
                !is.na(target_top_nt), target_top_nt != "") %>%
  dplyr::group_by(target_cell_type) %>%
  dplyr::summarise(
    target_nt = names(sort(table(target_top_nt), decreasing = TRUE))[1],
    nt_agreement = max(table(target_top_nt)) / dplyr::n(),
    .groups = "drop")

# Join and check consistency
nt_check <- results %>%
  dplyr::filter(assigned_cell_type != "") %>%
  dplyr::left_join(banc_nt %>% dplyr::select(root_888, banc_nt), by = "root_888") %>%
  dplyr::left_join(target_type_nt, by = c("assigned_cell_type" = "target_cell_type")) %>%
  dplyr::filter(!is.na(banc_nt), banc_nt != "", !is.na(target_nt))

if (nrow(nt_check) > 0) {
  nt_check$nt_match <- tolower(nt_check$banc_nt) == tolower(nt_check$target_nt)
  nt_consistency <- mean(nt_check$nt_match) * 100

  message(sprintf("  NT consistency: %.1f%% (%d/%d neurons with both BANC NT and FAFB type NT)",
                  nt_consistency, sum(nt_check$nt_match), nrow(nt_check)))

  nt_mismatches <- nt_check %>%
    dplyr::filter(!nt_match) %>%
    dplyr::count(banc_nt, target_nt, assigned_cell_type, sort = TRUE)

  if (nrow(nt_mismatches) > 0) {
    message("  Top NT mismatches:")
    for (i in seq_len(min(10, nrow(nt_mismatches)))) {
      r <- nt_mismatches[i, ]
      message(sprintf("    BANC %s != FAFB %s (type %s, %d neurons)",
                      r$banc_nt, r$target_nt, r$assigned_cell_type, r$n))
    }
  }

  # Flag NT-inconsistent assignments in results
  nt_flagged <- nt_check %>%
    dplyr::filter(!nt_match) %>%
    dplyr::select(root_888, banc_nt, target_nt, assigned_cell_type)
  readr::write_csv(nt_flagged,
    alignment_path("validate-nt-mismatches", query = "banc", target = target_name,
                   region = region_name, side = side_to_process,
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir))
} else {
  message("  No neurons with both BANC NT and FAFB type NT — skipping")
}

######################################
### STEP 2.3: EXISTING ANNOTATION  ###
###           MISMATCH DETECTION   ###
######################################

message("\n=== Step 2.3: Mismatch detection ===")

# Compare against existing fafb_cell_type in SeaTable
existing <- banc.meta %>%
  dplyr::select(root_888 = root_id, existing_fafb_type = fafb_cell_type) %>%
  dplyr::filter(!is.na(existing_fafb_type), existing_fafb_type != "") %>%
  dplyr::distinct(root_888, .keep_all = TRUE)

mismatches <- results %>%
  dplyr::filter(assigned_cell_type != "") %>%
  dplyr::inner_join(existing, by = "root_888") %>%
  dplyr::filter(assigned_cell_type != existing_fafb_type)

message(sprintf("  Neurons with existing fafb_cell_type: %d", nrow(existing)))
message(sprintf("  Mismatches (connectivity != existing): %d", nrow(mismatches)))

if (nrow(mismatches) > 0) {
  mismatch_path <- alignment_path("validate-mismatches", query = "banc", target = target_name,
                                   region = region_name, side = side_to_process,
                                   vq = banc.version, vt = nblast_version,
                                   ext = "csv", dir = data_dir)
  readr::write_csv(mismatches, mismatch_path)
  message(sprintf("  Saved: %s", mismatch_path))
}

######################################
### SUMMARY                        ###
######################################

message("\n=== Validation complete ===")
message(sprintf("Side: %s", side_to_process))
message(sprintf("Total neurons: %d", nrow(results)))
message(sprintf("Typed: %d (%.1f%%)",
                sum(results$assigned_cell_type != "", na.rm = TRUE),
                100 * mean(results$assigned_cell_type != "", na.rm = TRUE)))
if (exists("overall_acc")) message(sprintf("Holdout accuracy: %.1f%%", overall_acc))
if (exists("nt_consistency")) message(sprintf("NT consistency: %.1f%%", nt_consistency))
message(sprintf("Mismatches vs existing: %d", nrow(mismatches)))

})
