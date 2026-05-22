#' banc-alignment-plots — Generate publication-quality validation plots for the optic-lobe alignment.
#'
#' Faceted figures comparing BANC alignment assignments to FAFB ground truth
#' for the optic lobe. Falls back to per-side outputs if no bilateral
#' `banc_optic_both_alignment_optic_lobe_full.csv` is found.
#'
#' @section Reads:
#'   - `data/optic_lobe/banc_optic_both_alignment_optic_lobe_full.csv` (preferred)
#'   - `data/optic_lobe/banc_optic_{right,left}_alignment_results.csv` (fallback)
#'   - `data/optic_lobe/banc_optic_{side}_meta.csv`, `..._seeds.csv`
#'
#' @section Writes:
#'   - `inst/images/alignment/optic_lobes/*.png` — per-figure PNGs
#'
#' @section Paper:
#'   Methods §"Optic lobe cell-type alignment".

###########################################################
### Optic Lobe Alignment: Validation Plots
###
### Generates publication-quality validation figures
### comparing BANC alignment results to FAFB ground truth.
### Faceted by side (right | left).
###
### Output: inst/images/alignment/optic_lobes/*.png
###########################################################
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

local({

message("### banc: Optic lobe alignment validation plots ###")

######################
### CONFIGURATION  ###
######################

args <- commandArgs(trailingOnly = TRUE)
# Default to bilateral "both" output; fall back to per-side if not found
data_dir <- "data/optic_lobe"
output_dir <- "inst/images/alignment/optic_lobes"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
target_name    <- "fafb"
region_name    <- "optic-lobe"
nblast_version <- "783"
syn_source     <- Sys.getenv("BANC_SYN_SOURCE",
                              unset = banc.synapse.source.default)
bilateral_file <- alignment_path("align-optic-lobe-full", query = "banc",
                                  target = target_name, region = region_name,
                                  side = "both",
                                  vq = banc.version, vt = nblast_version,
                                  ext = "csv", dir = data_dir)
use_bilateral <- file.exists(bilateral_file)
sides <- if (use_bilateral) "both" else c("right", "left")

# BANC-Project color palette
banc_col   <- "#2a5f75"
fafb_col   <- "#FC1707"
accent_col <- "#EE4244"

# Super_class colors from BANC-Project
sc_colors <- c(
  optic_lobe_intrinsic = "#52aee3",
  visual_centrifugal   = "#121B56",
  visual_projection    = "#3E437E",
  ascending            = "#089c39",
  central              = "#2a5f75",
  descending           = "#F6B83C",
  sensory              = "#EBEF37",
  motor                = "#ee3232",
  unknown              = "#999999"
)

# LaCroix palette for tier colors
tier_colors <- c("1" = "#5BB6E4", "2" = "#4BA747", "3" = "#EE5B32")

######################
### LOAD DATA      ###
######################

all_data <- list()

if (use_bilateral) {
  # Bilateral: load single combined file, split by side from metadata
  message(sprintf("  Loading bilateral results: %s", bilateral_file))
  results <- readr::read_csv(bilateral_file,
    col_types = readr::cols(root_888 = "c", best_target_match = "c"), show_col_types = FALSE)
  banc_meta <- readr::read_csv(
    alignment_path("prep-banc-meta", query = "banc", target = target_name,
                   region = region_name, side = "both",
                   vq = banc.version, syn = syn_source,
                   ext = "csv", dir = data_dir),
    col_types = readr::cols(root_id = "c"), show_col_types = FALSE)
  target_meta <- readr::read_csv(
    alignment_path("prep-target-meta", query = "banc", target = target_name,
                   region = region_name, side = "both",
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir),
    col_types = readr::cols(target_id = "c"), show_col_types = FALSE)
  nblast <- readr::read_csv(
    alignment_path("prep-nblast", query = "banc", target = target_name,
                   region = region_name, side = "both",
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir),
    col_types = readr::cols(banc_id = "c", target_id = "c"), show_col_types = FALSE)
  seeds <- readr::read_csv(
    alignment_path("prep-seeds", query = "banc", target = target_name,
                   region = region_name, side = "both",
                   vq = banc.version, vt = nblast_version,
                   ext = "csv", dir = data_dir),
    col_types = readr::cols(root_888 = "c"), show_col_types = FALSE)

  # Add side to results from metadata
  side_lookup <- banc_meta %>%
    dplyr::select(root_id, side) %>%
    dplyr::mutate(side = tolower(side)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  results <- results %>%
    dplyr::left_join(side_lookup, by = c("root_888" = "root_id"))

  # Split into per-side views
  for (s in c("right", "left")) {
    side_results <- results %>% dplyr::filter(side == s)
    side_target <- target_meta %>% dplyr::filter(tolower(side) == s)
    side_nblast <- nblast %>% dplyr::filter(banc_id %in% side_results$root_888)
    side_seeds <- seeds %>% dplyr::filter(root_888 %in% side_results$root_888)
    all_data[[s]] <- list(results = side_results, nblast = side_nblast,
                          target_meta = side_target, seeds = side_seeds,
                          banc_meta = banc_meta %>% dplyr::filter(tolower(side) == s))
    message(sprintf("  %s: %d results, %d NBLAST scores",
                    s, nrow(side_results), nrow(side_nblast)))
  }
  # Keep full bilateral refs for NT plot
  all_data[["_banc_meta"]] <- banc_meta
  all_data[["_target_meta"]] <- target_meta
  all_data[["_results"]] <- results
} else {
  # Per-side loading (legacy)
  for (side in sides) {
    results_file <- alignment_path("align-results", query = "banc", target = target_name,
                                    region = region_name, side = side,
                                    vq = banc.version, vt = nblast_version,
                                    ext = "csv", dir = data_dir)
    if (!file.exists(results_file)) {
      message(sprintf("  Skipping %s side: results not found", side))
      next
    }
    results <- readr::read_csv(results_file,
      col_types = readr::cols(root_888 = "c", best_target_match = "c"), show_col_types = FALSE)
    nblast <- readr::read_csv(
      alignment_path("prep-nblast", query = "banc", target = target_name,
                     region = region_name, side = side,
                     vq = banc.version, vt = nblast_version,
                     ext = "csv", dir = data_dir),
      col_types = readr::cols(banc_id = "c", target_id = "c"), show_col_types = FALSE)
    target_meta <- readr::read_csv(
      alignment_path("prep-target-meta", query = "banc", target = target_name,
                     region = region_name, side = side,
                     vq = banc.version, vt = nblast_version,
                     ext = "csv", dir = data_dir),
      col_types = readr::cols(target_id = "c"), show_col_types = FALSE)
    seeds <- readr::read_csv(
      alignment_path("prep-seeds", query = "banc", target = target_name,
                     region = region_name, side = side,
                     vq = banc.version, vt = nblast_version,
                     ext = "csv", dir = data_dir),
      col_types = readr::cols(root_888 = "c"), show_col_types = FALSE)
    all_data[[side]] <- list(results = results, nblast = nblast,
                             target_meta = target_meta, seeds = seeds)
    message(sprintf("  Loaded %s: %d results, %d NBLAST scores",
                    side, nrow(results), nrow(nblast)))
  }
}

# Use right/left for all plots regardless of loading mode
sides <- intersect(c("right", "left"), names(all_data))
if (length(sides) == 0) {
  message("  No data found. Exiting.")
  return(invisible())
}

######################
### PLOT A: NBLAST ###
######################

message("\n=== Plot A: NBLAST score distributions ===")

nblast_plot_data <- dplyr::bind_rows(lapply(sides, function(side) {
  d <- all_data[[side]]
  typed <- d$results %>% dplyr::filter(assigned_cell_type != "")

  # Connectivity-assigned match NBLAST score
  assigned_scores <- typed %>%
    dplyr::filter(!is.na(best_target_match), best_target_match != "") %>%
    dplyr::inner_join(d$nblast,
                      by = c("root_888" = "banc_id", "best_target_match" = "target_id")) %>%
    dplyr::transmute(root_888, score = nblast_score, type = "Connectivity match", side = side)

  # Overall top NBLAST score (best match regardless of type)
  top_scores <- d$nblast %>%
    dplyr::group_by(banc_id) %>%
    dplyr::slice_max(nblast_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::filter(banc_id %in% typed$root_888) %>%
    dplyr::transmute(root_888 = banc_id, score = nblast_score,
                     type = "Top overall match", side = side)

  dplyr::bind_rows(assigned_scores, top_scores)
}))

# Count neurons without FAFB match per side
no_match_counts <- sapply(sides, function(side) {
  d <- all_data[[side]]
  typed <- d$results %>% dplyr::filter(assigned_cell_type != "")
  sum(is.na(typed$best_target_match) | typed$best_target_match == "")
})

p_nblast <- ggplot2::ggplot(nblast_plot_data, ggplot2::aes(x = score, fill = type)) +
  ggplot2::geom_density(alpha = 0.6, linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = c(0.2, 0.3, 0.5),
                      linetype = "dashed", color = "grey40", linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = c("Connectivity match" = banc_col,
                                        "Top overall match" = "grey70")) +
  ggplot2::facet_wrap(~side, scales = "free_y") +
  ggplot2::labs(x = "NBLAST score", y = "Density",
                title = "NBLAST score validation",
                subtitle = paste0("Neurons without FAFB match excluded: ",
                                  paste(sprintf("%s=%d", names(no_match_counts),
                                                no_match_counts), collapse = ", ")),
                fill = NULL) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom",
                 strip.text = ggplot2::element_text(size = 14, face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(output_dir, "nblast_validation.png"), p_nblast,
                width = 10, height = 5, dpi = 300, bg = "transparent")
message("  Saved: nblast_validation.png")

######################
### PLOT B: COUNTS ###
######################

message("\n=== Plot B: Type counts (BANC vs FAFB) ===")

count_plot_data <- dplyr::bind_rows(lapply(sides, function(side) {
  d <- all_data[[side]]

  banc_counts <- d$results %>%
    dplyr::filter(assigned_cell_type != "") %>%
    dplyr::count(cell_type = assigned_cell_type, name = "banc_count")

  target_counts <- d$target_meta %>%
    dplyr::filter(!is.na(target_cell_type), target_cell_type != "") %>%
    dplyr::count(cell_type = target_cell_type, name = "target_count")

  # Super_class per type (majority vote from target dataset)
  target_sc <- d$target_meta %>%
    dplyr::filter(!is.na(target_cell_type), target_cell_type != "",
                  !is.na(target_super_class)) %>%
    dplyr::group_by(target_cell_type) %>%
    dplyr::summarise(super_class = names(sort(table(target_super_class),
                                              decreasing = TRUE))[1],
                     .groups = "drop") %>%
    dplyr::rename(cell_type = target_cell_type)

  dplyr::full_join(banc_counts, target_counts, by = "cell_type") %>%
    dplyr::mutate(banc_count = tidyr::replace_na(banc_count, 0),
                  target_count = tidyr::replace_na(target_count, 0)) %>%
    dplyr::left_join(target_sc, by = "cell_type") %>%
    dplyr::mutate(super_class = tidyr::replace_na(super_class, "unknown"),
                  side = side,
                  abs_diff = abs(banc_count - target_count))
}))

# Label top 10 outliers per side
count_plot_data <- count_plot_data %>%
  dplyr::group_by(side) %>%
  dplyr::mutate(label = ifelse(dplyr::row_number(dplyr::desc(abs_diff)) <= 10,
                                cell_type, NA_character_)) %>%
  dplyr::ungroup()

p_counts <- ggplot2::ggplot(count_plot_data,
    ggplot2::aes(x = target_count + 1, y = banc_count + 1, color = super_class)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(alpha = 0.7, size = 1.5) +
  ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 2.5,
                            max.overlaps = 15, segment.color = "grey70",
                            na.rm = TRUE) +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::scale_color_manual(values = sc_colors, na.value = "grey60") +
  ggplot2::facet_wrap(~side) +
  ggplot2::labs(x = "FAFB neuron count (per type)", y = "BANC neuron count (per type)",
                title = "Cell type neuron counts: BANC vs FAFB",
                color = "Super class") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom",
                 legend.direction = "horizontal",
                 strip.text = ggplot2::element_text(size = 14, face = "bold"),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 3)))

ggplot2::ggsave(file.path(output_dir, "type_counts.png"), p_counts,
                width = 12, height = 6, dpi = 300, bg = "transparent")
message("  Saved: type_counts.png")

######################################
### PLOT B2: PER-SUPER_CLASS COUNTS ##
######################################

message("\n=== Plot B2: Per-super_class type counts ===")

# Compute % coverage (BANC / FAFB) per type
count_plot_data <- count_plot_data %>%
  dplyr::mutate(pct_coverage = ifelse(target_count > 0,
                                       100 * banc_count / target_count, 0))

# Filter to super_classes with enough types
sc_with_data <- count_plot_data %>%
  dplyr::filter(target_count > 0) %>%
  dplyr::count(super_class) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::pull(super_class)

count_by_sc <- count_plot_data %>%
  dplyr::filter(super_class %in% sc_with_data, target_count > 0)

# Label highly divergent types (coverage < 50% or > 200%, and FAFB count > 10)
count_by_sc <- count_by_sc %>%
  dplyr::mutate(
    divergent_label = ifelse(
      (pct_coverage < 50 | pct_coverage > 200) & target_count >= 10,
      cell_type, NA_character_
    )
  )

p_coverage <- ggplot2::ggplot(count_by_sc,
    ggplot2::aes(x = target_count, y = pct_coverage)) +
  ggplot2::geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(ggplot2::aes(color = super_class), alpha = 0.6, size = 1.5) +
  ggrepel::geom_text_repel(ggplot2::aes(label = divergent_label), size = 2.2,
                            max.overlaps = 20, segment.color = "grey70",
                            na.rm = TRUE) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_color_manual(values = sc_colors, na.value = "grey60") +
  ggplot2::facet_grid(super_class ~ side, scales = "free_y") +
  ggplot2::labs(x = "FAFB neuron count (per type)", y = "% of FAFB count covered by BANC",
                title = "Type coverage by super_class",
                subtitle = "Labels on types with <50% or >200% coverage (FAFB count >= 10)") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(size = 10),
                 strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
                 panel.grid.minor = ggplot2::element_blank())

n_sc <- length(sc_with_data)
ggplot2::ggsave(file.path(output_dir, "type_coverage_by_super_class.png"), p_coverage,
                width = 10, height = 3 + 2 * n_sc, dpi = 300, bg = "transparent")
message("  Saved: type_coverage_by_super_class.png")

############################
### PLOT C: CONFIDENCE   ###
############################

message("\n=== Plot C: Confidence distribution ===")

conf_plot_data <- dplyr::bind_rows(lapply(sides, function(side) {
  d <- all_data[[side]]
  d$results %>%
    dplyr::filter(assigned_cell_type != "") %>%
    dplyr::mutate(side = side, tier = as.character(tier),
                  has_fafb_match = !is.na(best_target_match) & best_target_match != "")
}))

# Median confidence per side
medians <- conf_plot_data %>%
  dplyr::group_by(side) %>%
  dplyr::summarise(med = median(confidence, na.rm = TRUE), .groups = "drop")

# No-match count per side
no_match_conf <- conf_plot_data %>%
  dplyr::group_by(side) %>%
  dplyr::summarise(n_no_match = sum(!has_fafb_match),
                   pct = 100 * mean(!has_fafb_match), .groups = "drop")

p_conf <- ggplot2::ggplot(conf_plot_data, ggplot2::aes(x = confidence, fill = tier)) +
  ggplot2::geom_histogram(bins = 60, alpha = 0.7, position = "stack") +
  ggplot2::geom_vline(data = medians, ggplot2::aes(xintercept = med),
                      linetype = "dashed", color = "grey30") +
  ggplot2::geom_text(data = medians,
                     ggplot2::aes(x = med, y = Inf,
                                  label = sprintf("median=%.3f", med)),
                     inherit.aes = FALSE, vjust = 1.5, hjust = -0.1, size = 3.5) +
  ggplot2::scale_fill_manual(values = tier_colors,
                             labels = c("1" = "Tier 1 (SeaTable)",
                                        "2" = "Tier 2 (NBLAST)",
                                        "3" = "Tier 3 (unassigned)")) +
  ggplot2::facet_wrap(~side, scales = "free_y") +
  ggplot2::labs(x = "Confidence (best - runner-up score)", y = "Count",
                title = "Assignment confidence distribution",
                subtitle = paste0("Neurons without FAFB match: ",
                                  paste(sprintf("%s=%d (%.1f%%)", no_match_conf$side,
                                                no_match_conf$n_no_match,
                                                no_match_conf$pct), collapse = ", ")),
                fill = "Initialization") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom",
                 strip.text = ggplot2::element_text(size = 14, face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(output_dir, "confidence_distribution.png"), p_conf,
                width = 10, height = 5, dpi = 300, bg = "transparent")
message("  Saved: confidence_distribution.png")

######################################
### PLOT C2: NTAC AGREEMENT        ###
######################################

message("\n=== Plot C2: NTAC vs alignment agreement ===")

# Load NTAC results (combined bilateral or per-side)
ntac_file <- alignment_path("ntac-full", query = "banc", target = target_name,
                             region = region_name, side = "both",
                             vq = banc.version, vt = nblast_version,
                             ext = "csv", dir = data_dir)
if (!file.exists(ntac_file)) {
  ntac_r <- alignment_path("ntac-full", query = "banc", target = target_name,
                            region = region_name, side = "right",
                            vq = banc.version, vt = nblast_version,
                            ext = "csv", dir = data_dir)
  ntac_l <- alignment_path("ntac-full", query = "banc", target = target_name,
                            region = region_name, side = "left",
                            vq = banc.version, vt = nblast_version,
                            ext = "csv", dir = data_dir)
  if (file.exists(ntac_r) && file.exists(ntac_l)) {
    ntac_combined <- dplyr::bind_rows(
      readr::read_csv(ntac_r, col_types = readr::cols(root_888 = "c"), show_col_types = FALSE),
      readr::read_csv(ntac_l, col_types = readr::cols(root_888 = "c"), show_col_types = FALSE)
    )
  } else {
    ntac_combined <- NULL
  }
} else {
  ntac_combined <- readr::read_csv(ntac_file,
    col_types = readr::cols(root_888 = "c"), show_col_types = FALSE)
}

if (!is.null(ntac_combined)) {
  # Get GT neurons to exclude
  banc_meta_all <- if (!is.null(all_data[["_banc_meta"]])) {
    all_data[["_banc_meta"]]
  } else {
    dplyr::bind_rows(lapply(sides, function(s) {
      if (!is.null(all_data[[s]]$banc_meta)) all_data[[s]]$banc_meta else NULL
    }))
  }

  gt_ids <- banc_meta_all %>%
    dplyr::filter(!is.na(cell_type) & cell_type != "" & !grepl("^auto:", cell_type)) %>%
    dplyr::pull(root_id) %>% unique()

  # Build comparison: alignment vs NTAC for non-GT neurons
  align_all <- dplyr::bind_rows(lapply(sides, function(s) {
    all_data[[s]]$results %>%
      dplyr::filter(assigned_cell_type != "") %>%
      dplyr::mutate(side = s)
  }))

  # Add super_class from metadata
  sc_lookup <- banc_meta_all %>%
    dplyr::select(root_id, super_class) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)

  ntac_comparison <- align_all %>%
    dplyr::filter(!root_888 %in% gt_ids) %>%
    dplyr::inner_join(
      ntac_combined %>%
        dplyr::filter(!is.na(ntac_cell_type) & ntac_cell_type != "?") %>%
        dplyr::select(root_888, ntac_cell_type) %>%
        dplyr::mutate(ntac_cell_type = as.character(ntac_cell_type)),
      by = "root_888"
    ) %>%
    dplyr::left_join(sc_lookup, by = c("root_888" = "root_id")) %>%
    dplyr::mutate(
      agree = assigned_cell_type == ntac_cell_type,
      conf_bin = cut(confidence,
                     breaks = c(-Inf, 0.02, 0.05, 0.1, 0.2, Inf),
                     labels = c("<0.02", "0.02-0.05", "0.05-0.1", "0.1-0.2", ">0.2"))
    )

  n_total <- nrow(ntac_comparison)
  n_agree <- sum(ntac_comparison$agree)
  message(sprintf("  NTAC comparison: %d non-GT neurons, %d agree (%.1f%%)",
                  n_total, n_agree, 100 * n_agree / max(n_total, 1)))

  # Plot 1: Agreement by super_class + side
  ntac_by_sc <- ntac_comparison %>%
    dplyr::filter(!is.na(super_class), super_class != "") %>%
    dplyr::group_by(side, super_class) %>%
    dplyr::summarise(n = dplyr::n(),
                     pct_agree = 100 * mean(agree),
                     .groups = "drop") %>%
    dplyr::filter(n >= 10)

  p_ntac_sc <- ggplot2::ggplot(ntac_by_sc,
      ggplot2::aes(x = reorder(super_class, pct_agree), y = pct_agree, fill = super_class)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%\n(n=%s)", pct_agree,
                                                     format(n, big.mark = ","))),
                       hjust = -0.1, size = 3) +
    ggplot2::scale_fill_manual(values = sc_colors, na.value = "grey60") +
    ggplot2::coord_flip(ylim = c(0, 105)) +
    ggplot2::facet_wrap(~side) +
    ggplot2::labs(x = NULL, y = "Agreement (%)",
                  title = "Alignment vs NTAC agreement (non-GT neurons)",
                  subtitle = sprintf("Overall: %.1f%% agree (%s neurons)",
                                     100 * n_agree / max(n_total, 1),
                                     format(n_total, big.mark = ","))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(file.path(output_dir, "ntac_agreement_by_super_class.png"), p_ntac_sc,
                  width = 10, height = 5, dpi = 300, bg = "transparent")
  message("  Saved: ntac_agreement_by_super_class.png")

  # Plot 2: Agreement by confidence bin
  ntac_by_conf <- ntac_comparison %>%
    dplyr::group_by(side, conf_bin) %>%
    dplyr::summarise(n = dplyr::n(),
                     pct_agree = 100 * mean(agree),
                     .groups = "drop")

  p_ntac_conf <- ggplot2::ggplot(ntac_by_conf,
      ggplot2::aes(x = conf_bin, y = pct_agree, fill = pct_agree)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%\n(%s)", pct_agree,
                                                     format(n, big.mark = ","))),
                       vjust = -0.3, size = 3.2) +
    ggplot2::scale_fill_gradient(low = accent_col, high = banc_col, guide = "none") +
    ggplot2::facet_wrap(~side) +
    ggplot2::coord_cartesian(ylim = c(0, 105)) +
    ggplot2::labs(x = "Alignment confidence bin", y = "Agreement with NTAC (%)",
                  title = "Alignment vs NTAC agreement by confidence (non-GT neurons)",
                  subtitle = "Higher confidence assignments agree more with NTAC") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(file.path(output_dir, "ntac_agreement_by_confidence.png"), p_ntac_conf,
                  width = 10, height = 5, dpi = 300, bg = "transparent")
  message("  Saved: ntac_agreement_by_confidence.png")
} else {
  message("  NTAC results not found, skipping")
}

################################
### PLOT D: NT AGREEMENT     ###
################################

message("\n=== Plot D: Neurotransmitter agreement ===")

# Build NT comparison: BANC neuron's predicted NT vs consensus NT of assigned FAFB type
# FAFB consensus NT per type
target_meta_all <- if (!is.null(all_data[["_target_meta"]])) {
  all_data[["_target_meta"]]
} else {
  dplyr::bind_rows(lapply(sides, function(s) all_data[[s]]$target_meta))
}
target_type_nt <- target_meta_all %>%
  dplyr::filter(!is.na(target_cell_type), target_cell_type != "",
                !is.na(target_top_nt), target_top_nt != "") %>%
  dplyr::group_by(target_cell_type) %>%
  dplyr::summarise(target_consensus_nt = names(sort(table(tolower(target_top_nt)),
                                                    decreasing = TRUE))[1],
                   .groups = "drop")

# Build per-neuron NT comparison
nt_data <- dplyr::bind_rows(lapply(sides, function(side) {
  d <- all_data[[side]]
  results <- d$results %>% dplyr::filter(assigned_cell_type != "")

  # Get BANC neuron's NT from metadata, with histamine override for R7/R8
  if (!is.null(d$banc_meta)) {
    banc_nt <- d$banc_meta %>%
      dplyr::select(root_id, banc_nt = top_nt, super_class) %>%
      dplyr::mutate(banc_nt = tolower(banc_nt)) %>%
      dplyr::distinct(root_id, .keep_all = TRUE)
    results <- results %>%
      dplyr::left_join(banc_nt, by = c("root_888" = "root_id")) %>%
      # Override: neurons assigned to R7/R8 types get histamine as their NT
      dplyr::mutate(banc_nt = ifelse(assigned_cell_type %in% c("R7", "R8"),
                                      "histamine", banc_nt))
  } else {
    results$banc_nt <- NA_character_
    results$super_class <- NA_character_
  }

  results %>%
    dplyr::left_join(target_type_nt, by = c("assigned_cell_type" = "target_cell_type")) %>%
    dplyr::filter(!is.na(banc_nt), banc_nt != "",
                  !is.na(target_consensus_nt), target_consensus_nt != "") %>%
    dplyr::mutate(
      nt_match = banc_nt == target_consensus_nt,
      side = side
    )
}))

if (nrow(nt_data) > 0) {
  # Overall NT agreement
  overall_nt <- nt_data %>%
    dplyr::summarise(n = dplyr::n(), pct = 100 * mean(nt_match))
  message(sprintf("  Overall NT agreement: %.1f%% (%d neurons)", overall_nt$pct, overall_nt$n))

  # NT agreement by side + super_class
  nt_by_sc_side <- nt_data %>%
    dplyr::filter(!is.na(super_class), super_class != "") %>%
    dplyr::group_by(side, super_class) %>%
    dplyr::summarise(n = dplyr::n(),
                     pct_match = 100 * mean(nt_match),
                     .groups = "drop") %>%
    dplyr::filter(n >= 10)

  p_nt_sc <- ggplot2::ggplot(nt_by_sc_side,
      ggplot2::aes(x = reorder(super_class, pct_match), y = pct_match, fill = super_class)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%\n(n=%d)", pct_match, n)),
                       hjust = -0.1, size = 3) +
    ggplot2::scale_fill_manual(values = sc_colors, na.value = "grey60") +
    ggplot2::coord_flip(ylim = c(0, 105)) +
    ggplot2::facet_wrap(~side) +
    ggplot2::labs(x = NULL, y = "NT agreement (%)",
                  title = "Neurotransmitter agreement: BANC prediction vs FAFB type consensus",
                  subtitle = sprintf("Overall: %.1f%% (%s neurons)", overall_nt$pct,
                                     format(overall_nt$n, big.mark = ","))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(file.path(output_dir, "nt_agreement_by_super_class.png"), p_nt_sc,
                  width = 10, height = 5, dpi = 300, bg = "transparent")
  message("  Saved: nt_agreement_by_super_class.png")

  # NT agreement by BANC NT (which neurotransmitters have the most mismatches?)
  nt_by_nt <- nt_data %>%
    dplyr::group_by(side, banc_nt) %>%
    dplyr::summarise(n = dplyr::n(),
                     pct_match = 100 * mean(nt_match),
                     n_mismatch = sum(!nt_match),
                     .groups = "drop") %>%
    dplyr::filter(n >= 5)

  p_nt_nt <- ggplot2::ggplot(nt_by_nt,
      ggplot2::aes(x = reorder(banc_nt, pct_match), y = pct_match)) +
    ggplot2::geom_col(fill = banc_col, alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%\n(%d/%d)", pct_match, n - n_mismatch, n)),
                       hjust = -0.1, size = 3) +
    ggplot2::coord_flip(ylim = c(0, 105)) +
    ggplot2::facet_wrap(~side) +
    ggplot2::labs(x = "BANC predicted NT", y = "NT agreement (%)",
                  title = "NT agreement by BANC neurotransmitter prediction",
                  subtitle = "How often does assigned type's consensus NT match the neuron's predicted NT?") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(file.path(output_dir, "nt_agreement_by_banc_nt.png"), p_nt_nt,
                  width = 10, height = 5, dpi = 300, bg = "transparent")
  message("  Saved: nt_agreement_by_banc_nt.png")
  ######################################
  ### PLOT D3: HISTAMINE POOL        ###
  ######################################

  message("\n=== Plot D3: Histamine neuron assignments ===")

  # Neurons assigned to types with histamine as FAFB consensus NT (R7, R8, etc.)
  histamine_types <- target_type_nt %>%
    dplyr::filter(target_consensus_nt == "histamine") %>%
    dplyr::pull(cell_type)

  # Non-GT neurons assigned to histamine types
  histamine_neurons <- nt_data %>%
    dplyr::filter(target_consensus_nt == "histamine" | assigned_cell_type %in% histamine_types)

  # Also get non-GT neurons whose ASSIGNED type has histamine consensus
  # (including those without banc_nt data)
  histamine_all <- dplyr::bind_rows(lapply(sides, function(side) {
    d <- all_data[[side]]
    d$results %>%
      dplyr::filter(assigned_cell_type %in% histamine_types, !is_anchor) %>%
      dplyr::mutate(side = side,
                    conf_bin = cut(confidence,
                                   breaks = c(-Inf, 0.02, 0.05, 0.1, 0.2, Inf),
                                   labels = c("<0.02", "0.02-0.05", "0.05-0.1", "0.1-0.2", ">0.2")))
  }))

  if (nrow(histamine_all) > 0) {
    # Top assigned types for histamine-pool neurons by confidence
    hist_by_type <- histamine_all %>%
      dplyr::count(assigned_cell_type, conf_bin, side, name = "n") %>%
      dplyr::group_by(side) %>%
      dplyr::mutate(total_side = sum(n)) %>%
      dplyr::ungroup()

    # Keep top types per side
    top_hist_types <- histamine_all %>%
      dplyr::count(assigned_cell_type, sort = TRUE) %>%
      head(15) %>%
      dplyr::pull(assigned_cell_type)

    hist_plot <- hist_by_type %>%
      dplyr::filter(assigned_cell_type %in% top_hist_types) %>%
      dplyr::mutate(assigned_cell_type = factor(assigned_cell_type, levels = rev(top_hist_types)))

    p_hist <- ggplot2::ggplot(hist_plot,
        ggplot2::aes(x = assigned_cell_type, y = n, fill = conf_bin)) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = 1) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~side, scales = "free_x") +
      ggplot2::labs(x = "Assigned cell type", y = "Count (non-GT neurons)",
                    title = sprintf("Histamine-type assignments (non-GT, %d neurons)",
                                    nrow(histamine_all)),
                    subtitle = sprintf("Types with histamine FAFB consensus NT: %s",
                                       paste(histamine_types, collapse = ", ")),
                    fill = "Confidence") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 14, face = "bold"),
                     panel.grid.minor = ggplot2::element_blank())

    ggplot2::ggsave(file.path(output_dir, "histamine_assignments.png"), p_hist,
                    width = 10, height = 6, dpi = 300, bg = "transparent")
    message(sprintf("  Histamine pool: %d non-GT neurons assigned to %s",
                    nrow(histamine_all), paste(histamine_types, collapse = ", ")))
    message("  Saved: histamine_assignments.png")
  } else {
    message("  No histamine-type assignments found")
  }

} else {
  message("  Insufficient NT data for plots")
}

message("\n=== Validation plots complete ===")
message(sprintf("Output: %s/", output_dir))

})
