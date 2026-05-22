#' banc-alignment-sweep-eval — Evaluate alignment sweep runs.
#'
#' Reads all sweep outputs under a given data directory, computes per-run
#' metrics + cross-run agreement, and writes summary CSVs + a comparison PDF.
#' Read-only on SeaTable.
#'
#' @section Reads:
#'   - `data/<dataset_dir>/banc_<prefix>_both_meta.csv`
#'   - `data/<dataset_dir>/fafb_<prefix>_both_meta.csv`
#'   - `data/<dataset_dir>/banc_<prefix>_both_seeds.csv`
#'   - `data/<dataset_dir>/fafb_type_capacity.csv`
#'   - `data/<dataset_dir>/banc_<prefix>_both_alignment_sweep_align_*.csv`
#'   - `data/<dataset_dir>/banc_<prefix>_both_ntac_sweep_ntac_*.csv`
#'
#' @section Writes:
#'   - `data/<dataset_dir>/sweep_eval/sweep_metrics.csv`
#'   - `data/<dataset_dir>/sweep_eval/sweep_agreement.csv`
#'   - `data/<dataset_dir>/sweep_eval/sweep_consensus.csv`
#'   - `data/<dataset_dir>/sweep_eval/sweep_disagreements.csv`
#'   - `data/<dataset_dir>/sweep_eval/sweep_plots.pdf`
#'
#' @section CLI:
#'   <dataset>   one of `whole_brain` (default) | `optic`

###########################################################
### Alignment sweep evaluator
###
### Reads all sweep outputs under a given data directory,
### computes per-run metrics and cross-run agreement, and
### writes summary CSVs + a comparison PDF. Does NOT touch
### SeaTable — this script is intentionally read-only on
### the annotation database.
###
### Usage:
###   Rscript alignment/banc-alignment-sweep-eval.R whole_brain
###   Rscript alignment/banc-alignment-sweep-eval.R optic
###
### Inputs (under data/<dataset_dir>/):
###   banc_<prefix>_both_meta.csv
###   fafb_<prefix>_both_meta.csv
###   banc_<prefix>_both_seeds.csv
###   fafb_type_capacity.csv
###   banc_<prefix>_both_alignment_sweep_align_*.csv    (glob)
###   banc_<prefix>_both_ntac_sweep_ntac_*.csv          (glob)
###
### Outputs (under data/<dataset_dir>/sweep_eval/):
###   sweep_metrics.csv         -- one row per run
###   sweep_agreement.csv       -- pairwise Jaccard agreement between align runs
###   sweep_consensus.csv       -- per-neuron assignment where all 7 align runs agree
###   sweep_disagreements.csv   -- typed neurons where A (0% anchors) and G (100%) disagree
###   sweep_plots.pdf           -- scaling curves, per-super_class accuracy, agreement heatmap
###########################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) > 0) args[1] else "whole_brain"
stopifnot(dataset %in% c("whole_brain", "optic"))

if (dataset == "whole_brain") {
  data_dir   <- "data/whole_brain_alignment"
  file_prefix <- "brain"
} else {
  data_dir   <- "data/optic_lobe"
  file_prefix <- "optic"
}

out_dir <- file.path(data_dir, "sweep_eval")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("=== Sweep eval: %s (%s) ===", dataset, data_dir))

######################################################################
## Load shared inputs
######################################################################
banc_meta <- read_csv(file.path(data_dir, sprintf("banc_%s_both_meta.csv", file_prefix)),
                      show_col_types = FALSE) %>%
  mutate(root_id = as.character(root_id))

fafb_meta <- read_csv(file.path(data_dir, sprintf("fafb_%s_both_meta.csv", file_prefix)),
                      show_col_types = FALSE) %>%
  mutate(fafb_id = as.character(fafb_id))

seeds <- read_csv(file.path(data_dir, sprintf("banc_%s_both_seeds.csv", file_prefix)),
                  show_col_types = FALSE) %>%
  mutate(root_888 = as.character(root_888))

# The "ground truth" cell_type is what the seeds file records for tiers 1-2,
# pulled from SeaTable (tier 1) or the best NBLAST match (tier 2). Tier 3 has
# no label and is excluded from accuracy calculations.
gt <- seeds %>%
  filter(tier == 1L, !is.na(cell_type), cell_type != "") %>%
  select(root_888, true_type = cell_type)

super_class_lookup <- banc_meta %>%
  select(root_888 = root_id, super_class) %>%
  distinct(root_888, .keep_all = TRUE)

fafb_type_counts <- fafb_meta %>%
  filter(!is.na(cell_type), cell_type != "") %>%
  count(cell_type, name = "fafb_count")

######################################################################
## Discover sweep runs
######################################################################
parse_run_name <- function(path) {
  nm <- basename(path)
  # banc_brain_both_alignment_sweep_align_a050_ho050.csv
  # banc_brain_both_ntac_sweep_ntac_a080_ho020.csv
  m <- str_match(nm,
    sprintf("^banc_%s_both_(alignment|ntac)_sweep_(align|ntac)_a(\\d{3})_ho(\\d{3})\\.csv$",
            file_prefix))
  if (is.na(m[1, 1])) return(NULL)
  list(path = path,
       file_tag = m[1, 2],   # "alignment" or "ntac"
       method   = m[1, 3],   # "align" or "ntac"
       anchor_pct  = as.integer(m[1, 4]),
       holdout_pct = as.integer(m[1, 5]),
       run_id   = sprintf("%s_a%s_ho%s", m[1, 3], m[1, 4], m[1, 5]))
}

candidate_files <- list.files(data_dir,
  pattern = sprintf("^banc_%s_both_(alignment|ntac)_sweep_.*\\.csv$", file_prefix),
  full.names = TRUE)

runs <- compact(map(candidate_files, parse_run_name))
if (length(runs) == 0) stop("No sweep outputs found in ", data_dir)

message(sprintf("  Found %d sweep runs", length(runs)))

######################################################################
## Load each run; attach GT and super_class
######################################################################
load_run <- function(run) {
  df <- read_csv(run$path, show_col_types = FALSE)
  # Vintage drift: three NTAC sweep CSVs (a080/a090/a100) were written before
  # the v888 column rename and use `root_850`. Most root_ids are stable across
  # v850→v888 (proofreading-only edits), so we treat root_850 as a substitute
  # join key — rows whose root_id flipped will simply drop on the seeds join.
  # CAVEAT: the alignment code itself has also drifted since those CSVs were
  # written; their accuracy numbers are indicative, not strictly apples-to-
  # apples with the newer runs.
  if (!"root_888" %in% names(df) && "root_850" %in% names(df)) {
    message(sprintf("  NOTE %s: using root_850 as root_888 (cross-vintage; code-drift caveat)",
                    run$run_id))
    df <- df %>% rename(root_888 = root_850)
  }
  if (!"root_888" %in% names(df)) {
    message(sprintf("  SKIP %s: no root_888/root_850 column found", run$run_id))
    return(NULL)
  }
  df <- df %>% mutate(root_888 = as.character(root_888))

  # Align runs: assigned_cell_type + is_anchor (holdout = !is_anchor & typed)
  # NTAC runs:  ntac_cell_type + is_holdout (written directly since we now use
  #             the shared stratified splitter)
  if ("assigned_cell_type" %in% names(df)) {
    df <- df %>% rename(pred_type = assigned_cell_type)
    if (!"is_anchor" %in% names(df)) df$is_anchor <- NA
    df$is_holdout <- ifelse(is.na(df$is_anchor), NA, !df$is_anchor)
  } else if ("ntac_cell_type" %in% names(df)) {
    df <- df %>% rename(pred_type = ntac_cell_type)
    if (!"is_holdout" %in% names(df)) df$is_holdout <- NA
    df$is_anchor <- ifelse(is.na(df$is_holdout), NA, !df$is_holdout)
  } else {
    stop("Run ", run$run_id, " has no recognised prediction column")
  }

  df %>%
    select(root_888, pred_type, is_anchor, is_holdout) %>%
    mutate(pred_type = ifelse(is.na(pred_type), "", pred_type)) %>%
    left_join(gt, by = "root_888") %>%
    left_join(super_class_lookup, by = "root_888")
}

run_data <- set_names(map(runs, load_run), map_chr(runs, "run_id"))
# Drop runs that load_run skipped (vintage mismatch). Keep `runs` index aligned.
keep <- !map_lgl(run_data, is.null)
run_data <- run_data[keep]
runs <- runs[keep]
message(sprintf("  Using %d runs after vintage filter", length(runs)))

######################################################################
## Per-run metrics
######################################################################
compute_metrics <- function(run, df) {
  typed <- df %>% filter(!is.na(true_type))
  # Full acc: on all typed neurons (whether or not they were held out).
  # Holdout acc: only on non-anchor rows (align only; NTAC gets NA).
  full_correct <- typed$pred_type == typed$true_type
  full_acc <- if (nrow(typed) > 0) 100 * mean(full_correct) else NA_real_
  full_acc_macro <- if (nrow(typed) > 0) {
    typed %>%
      group_by(true_type) %>%
      summarise(acc = mean(pred_type == true_type), .groups = "drop") %>%
      pull(acc) %>% mean() * 100
  } else NA_real_

  # Holdout accuracy: use is_holdout for both align (derived from !is_anchor)
  # and NTAC (written directly). At 100% anchors neither has a holdout set.
  ho <- typed %>% filter(!is.na(is_holdout) & is_holdout)
  if (nrow(ho) > 0) {
    holdout_acc <- 100 * mean(ho$pred_type == ho$true_type)
    n_holdout <- nrow(ho)
  } else {
    holdout_acc <- NA_real_
    n_holdout <- 0L
  }

  # Per-super_class accuracy over typed rows
  per_sc <- typed %>%
    filter(!is.na(super_class), super_class != "") %>%
    group_by(super_class) %>%
    summarise(acc = 100 * mean(pred_type == true_type),
              n = n(), .groups = "drop")

  # Type counts vs FAFB
  pred_counts <- df %>%
    filter(pred_type != "") %>%
    count(pred_type, name = "banc_count") %>%
    rename(cell_type = pred_type)

  ct_joined <- full_join(fafb_type_counts, pred_counts, by = "cell_type") %>%
    mutate(fafb_count = replace_na(fafb_count, 0),
           banc_count = replace_na(banc_count, 0))

  type_count_spearman <- if (nrow(ct_joined) > 2) {
    suppressWarnings(cor(ct_joined$fafb_count, ct_joined$banc_count, method = "spearman"))
  } else NA_real_

  type_coverage_pct <- 100 * sum(pred_counts$banc_count > 0) / max(1, nrow(fafb_type_counts))

  list(
    summary = tibble(
      run_id       = run$run_id,
      method       = run$method,
      anchor_pct   = run$anchor_pct,
      holdout_pct  = run$holdout_pct,
      n_neurons    = nrow(df),
      n_assigned   = sum(df$pred_type != ""),
      n_typed_gt   = nrow(typed),
      n_holdout    = n_holdout,
      full_acc     = full_acc,
      full_acc_macro = full_acc_macro,
      holdout_acc  = holdout_acc,
      type_coverage_pct = type_coverage_pct,
      type_count_spearman = type_count_spearman
    ),
    per_sc = per_sc %>% mutate(run_id = run$run_id, .before = 1)
  )
}

metrics_list <- map2(runs, run_data, compute_metrics)
metrics <- bind_rows(map(metrics_list, "summary")) %>%
  arrange(method, anchor_pct)
per_sc_df <- bind_rows(map(metrics_list, "per_sc"))

write_csv(metrics,   file.path(out_dir, "sweep_metrics.csv"))
write_csv(per_sc_df, file.path(out_dir, "sweep_per_super_class.csv"))

######################################################################
## Cross-run agreement (align runs only, pairwise)
######################################################################
align_ids <- map_chr(runs, "run_id")[map_chr(runs, "method") == "align"]
if (length(align_ids) >= 2) {
  align_preds <- map_dfc(align_ids, function(id) {
    set_names(list(run_data[[id]]$pred_type), id)
  })
  align_preds$root_888 <- run_data[[align_ids[1]]]$root_888

  pairs <- expand.grid(a = align_ids, b = align_ids, stringsAsFactors = FALSE) %>%
    filter(a < b)
  agree <- pairs %>%
    rowwise() %>%
    mutate(
      agree_pct = {
        pa <- align_preds[[a]]
        pb <- align_preds[[b]]
        mask <- pa != "" & pb != ""
        if (sum(mask) == 0) NA_real_ else 100 * mean(pa[mask] == pb[mask])
      }
    ) %>%
    ungroup()
  write_csv(agree, file.path(out_dir, "sweep_agreement.csv"))

  # Consensus: neurons where all align runs predict the same non-empty type
  consensus <- align_preds %>%
    rowwise() %>%
    mutate(
      preds = list(c_across(all_of(align_ids))),
      unique_preds = list(unique(preds[preds != ""])),
      is_consensus = length(unique_preds) == 1 && length(preds[preds == ""]) == 0
    ) %>%
    ungroup() %>%
    filter(is_consensus) %>%
    mutate(consensus_type = map_chr(unique_preds, 1)) %>%
    select(root_888, consensus_type)
  write_csv(consensus, file.path(out_dir, "sweep_consensus.csv"))
  message(sprintf("  Consensus across %d align runs: %d neurons",
                  length(align_ids), nrow(consensus)))

  # Disagreements between extremes (0% vs 100% anchors, if both present)
  a_zero <- align_ids[grepl("a000_", align_ids)]
  a_full <- align_ids[grepl("a100_", align_ids)]
  if (length(a_zero) == 1 && length(a_full) == 1) {
    disagree <- tibble(
      root_888  = run_data[[a_zero]]$root_888,
      pred_a000 = run_data[[a_zero]]$pred_type,
      pred_a100 = run_data[[a_full]]$pred_type,
      true_type = run_data[[a_zero]]$true_type,
      super_class = run_data[[a_zero]]$super_class
    ) %>%
      filter(pred_a000 != pred_a100, pred_a000 != "" | pred_a100 != "")
    write_csv(disagree, file.path(out_dir, "sweep_disagreements.csv"))
    message(sprintf("  Disagreements (a000 vs a100): %d neurons", nrow(disagree)))
  }
}

######################################################################
## Plots
######################################################################
pdf(file.path(out_dir, "sweep_plots.pdf"), width = 10, height = 7)

# 1. Scaling curve: full_acc vs anchor_pct
p1 <- ggplot(metrics, aes(x = anchor_pct, y = full_acc, colour = method, group = method)) +
  geom_line(size = 1) + geom_point(size = 3) +
  labs(title = sprintf("%s: full accuracy (all typed neurons) vs anchor fraction", dataset),
       x = "Anchor fraction (%)", y = "Accuracy (%)", colour = "Method") +
  theme_minimal(base_size = 12)
print(p1)

# 2. Holdout acc vs anchor_pct — both methods, same stratified split
ho_metrics <- metrics %>% filter(!is.na(holdout_acc))
if (nrow(ho_metrics) > 0) {
  p2 <- ggplot(ho_metrics, aes(x = anchor_pct, y = holdout_acc, colour = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 3) +
    labs(title = sprintf("%s: holdout accuracy vs anchor fraction", dataset),
         x = "Anchor fraction (%)", y = "Holdout accuracy (%)", colour = "Method") +
    theme_minimal(base_size = 12)
  print(p2)
}

# 3. Type count Spearman vs anchor fraction
p3 <- ggplot(metrics, aes(x = anchor_pct, y = type_count_spearman, colour = method)) +
  geom_line(size = 1) + geom_point(size = 3) +
  labs(title = sprintf("%s: type-count Spearman ρ (BANC vs FAFB)", dataset),
       x = "Anchor fraction (%)", y = "Spearman ρ", colour = "Method") +
  theme_minimal(base_size = 12)
print(p3)

# 4. Per-super_class accuracy across runs
if (nrow(per_sc_df) > 0) {
  p4 <- per_sc_df %>%
    left_join(metrics %>% select(run_id, method, anchor_pct), by = "run_id") %>%
    ggplot(aes(x = anchor_pct, y = acc, colour = method)) +
    geom_line() + geom_point() +
    facet_wrap(~ super_class, scales = "free_y") +
    labs(title = sprintf("%s: per-super_class accuracy", dataset),
         x = "Anchor fraction (%)", y = "Accuracy (%)", colour = "Method") +
    theme_minimal(base_size = 10)
  print(p4)
}

# 5. Agreement heatmap (align only)
if (exists("agree") && nrow(agree) > 0) {
  agree_sym <- bind_rows(
    agree %>% rename(x = a, y = b),
    agree %>% rename(x = b, y = a),
    tibble(x = align_ids, y = align_ids, agree_pct = 100)
  )
  p5 <- ggplot(agree_sym, aes(x = x, y = y, fill = agree_pct)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.0f", agree_pct)), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
    labs(title = sprintf("%s: pairwise agreement (align runs)", dataset),
         x = NULL, y = NULL, fill = "% agree") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p5)
}

dev.off()

message(sprintf("Done. Outputs under %s/", out_dir))
