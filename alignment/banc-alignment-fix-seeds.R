#' banc-alignment-fix-seeds — Regenerate whole-brain alignment seeds in place.
#'
#' One-shot fixer: rebuilds `banc_brain_<side>_seeds.csv` from the existing
#' meta + NBLAST CSVs using the corrected `anchor_super_classes` list
#' (adds compound forms like `sensory_ascending`). Avoids a full re-prep
#' cycle when only the seed policy has changed.
#'
#' @section Reads:
#'   - `<output_dir>/banc_brain_<side>_meta.csv`
#'   - `<output_dir>/fafb_brain_<side>_meta.csv`
#'   - `<output_dir>/banc_fafb_brain_<side>_nblast.csv`
#'
#' @section Writes:
#'   - `<output_dir>/banc_brain_<side>_seeds.csv`
#'
#' @section CLI:
#'   <output_dir>   required; directory holding the prep outputs
#'   [side]         optional; defaults to `both`

###########################################################
### Regenerate banc_brain_<side>_seeds.csv in-place using
### the fixed anchor_super_classes list (adds compound forms
### like sensory_ascending). Avoids a full re-prep cycle.
###
### Usage:
###   Rscript alignment/banc-alignment-fix-seeds.R \
###     <output_dir> [side]
###
### Defaults: side=both. Reads {banc,fafb}_brain_<side>_meta.csv
### and banc_fafb_brain_<side>_nblast.csv from <output_dir>,
### writes banc_brain_<side>_seeds.csv back to the same dir.
###########################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript banc-alignment-fix-seeds.R <output_dir> [side]")
output_dir <- args[1]
side <- if (length(args) >= 2) args[2] else "both"

stopifnot(dir.exists(output_dir))

# Mirror the whole-brain prep parameter set so v2 filenames line up.
target_name    <- "fafb"
region_name    <- "whole-brain"
nblast_version <- "783"
syn_source     <- Sys.getenv("BANC_SYN_SOURCE",
                              unset = banc.synapse.source.default)

# Match the (now-fixed) prep config
anchor_super_classes <- c("visual_centrifugal", "visual_projection",
                          "sensory", "sensory_ascending", "sensory_descending",
                          "motor", "endocrine", "efferent",
                          "ascending", "ascending_visceral_circulatory",
                          "descending")

banc_meta_path <- alignment_path("prep-banc-meta", query = "banc", target = target_name,
                                  region = region_name, side = side,
                                  vq = banc.version, syn = syn_source,
                                  ext = "csv", dir = output_dir)
target_meta_path <- alignment_path("prep-target-meta", query = "banc", target = target_name,
                                    region = region_name, side = side,
                                    vq = banc.version, vt = nblast_version,
                                    ext = "csv", dir = output_dir)
nblast_path    <- alignment_path("prep-nblast", query = "banc", target = target_name,
                                  region = region_name, side = side,
                                  vq = banc.version, vt = nblast_version,
                                  ext = "csv", dir = output_dir)
seeds_path     <- alignment_path("prep-seeds", query = "banc", target = target_name,
                                  region = region_name, side = side,
                                  vq = banc.version, vt = nblast_version,
                                  ext = "csv", dir = output_dir)

stopifnot(file.exists(banc_meta_path), file.exists(target_meta_path), file.exists(nblast_path))

message(sprintf("Regenerating seeds at %s", seeds_path))

# Detect root col from existing seeds (preserve schema)
existing_cols <- if (file.exists(seeds_path)) {
  names(read_csv(seeds_path, n_max = 0, show_col_types = FALSE))
} else character(0)
root_col <- existing_cols[grepl("^root_[0-9]+$", existing_cols)][1]
if (is.na(root_col)) root_col <- "root_888"
message(sprintf("  Using root col: %s", root_col))

banc.pool <- read_csv(banc_meta_path, show_col_types = FALSE,
                      col_types = cols(root_id = col_character(),
                                        supervoxel_id = col_character()))
target.pool <- read_csv(target_meta_path, show_col_types = FALSE,
                      col_types = cols(target_id = col_character()))
nblast.pool <- read_csv(nblast_path, show_col_types = FALSE,
                        col_types = cols(banc_id = col_character(),
                                          target_id = col_character(),
                                          nblast_score = col_double()))

banc_pool_ids <- unique(banc.pool$root_id)

# Tier 1: typed neurons; anchor vs holdout by super_class
typed_banc <- banc.pool %>%
  filter(!is.na(cell_type), cell_type != "") %>%
  select(!!root_col := root_id, cell_type, super_class, region)

anchor_typed <- typed_banc %>%
  filter(super_class %in% anchor_super_classes) %>%
  mutate(tier = 1L, is_holdout = FALSE) %>%
  select(-super_class, -region)

holdout_typed <- typed_banc %>%
  filter(!super_class %in% anchor_super_classes) %>%
  mutate(tier = 1L, is_holdout = TRUE) %>%
  select(-super_class, -region)

tier1 <- bind_rows(anchor_typed, holdout_typed)

message(sprintf("  Tier 1: %d total (%d anchor, %d holdout)",
                nrow(tier1), sum(!tier1$is_holdout), sum(tier1$is_holdout)))

# Anchor breakdown by super_class for sanity
sc_breakdown <- typed_banc %>%
  mutate(anchor = super_class %in% anchor_super_classes) %>%
  count(super_class, anchor) %>%
  arrange(desc(n))
message("  Anchor by super_class:")
for (i in seq_len(nrow(sc_breakdown))) {
  message(sprintf("    %-35s %s (n=%d)",
                  sc_breakdown$super_class[i],
                  if (sc_breakdown$anchor[i]) "ANCHOR" else "holdout",
                  sc_breakdown$n[i]))
}

# Tier 2: NBLAST-initialised
tier1_ids <- tier1[[root_col]]
if (nrow(nblast.pool) > 0) {
  best_nblast <- nblast.pool %>%
    filter(!banc_id %in% tier1_ids) %>%
    group_by(banc_id) %>%
    slice_max(nblast_score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(target.pool %>% select(target_id, target_cell_type) %>%
              distinct(target_id, .keep_all = TRUE), by = "target_id") %>%
    filter(!is.na(fafb_cell_type), fafb_cell_type != "")

  tier2 <- best_nblast %>%
    select(!!root_col := banc_id, cell_type = fafb_cell_type) %>%
    mutate(tier = 2L, is_holdout = FALSE)
} else {
  tier2 <- tibble(!!root_col := character(), cell_type = character(),
                  tier = integer(), is_holdout = logical())
}
message(sprintf("  Tier 2 (NBLAST): %d", nrow(tier2)))

tier3_ids <- setdiff(banc_pool_ids, c(tier1[[root_col]], tier2[[root_col]]))
tier3 <- tibble(!!root_col := tier3_ids, cell_type = NA_character_,
                tier = 3L, is_holdout = FALSE)
message(sprintf("  Tier 3 (unassigned): %d", nrow(tier3)))

seeds <- bind_rows(tier1, tier2, tier3)
write_csv(seeds, seeds_path)
message(sprintf("Wrote %s (%d rows)", seeds_path, nrow(seeds)))
