#' banc-nblast-plot — NBLAST score distribution + per-region matching inventory plots.
#'
#' @section Reads:
#'   - `banc_meta.csv`, per-dataset compiled NBLAST feathers
#'
#' @section Writes:
#'   - `inst/images/nblast/*.png`, `inst/images/cell_types/*.png`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`.

###########################################################
### NBLAST score distribution and matching inventory plots
###
### Generates:
###   1. Per-dataset NBLAST score density plots
###      (MANC, maleCNS, FAFB, hemibrain, mirror)
###   2. Matching inventory pie charts by region
###
### Reads: banc_meta.csv, compiled feather files
### Saves: inst/images/nblast/*.png,
###        inst/images/cell_types/*.png
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: generating NBLAST plots ###")
t_start <- Sys.time()

###########################
### Read data           ###
###########################

banc.meta <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_meta.csv"),
                             col_types = banc.col.types,
                             show_col_types = FALSE)

bc <- banctable_query() %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Read compiled feather files for inventory plots
banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather"))
banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
banc.meta.hemibrain.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))

# Ensure output directories exist
dir.create("inst/images/nblast", showWarnings = FALSE, recursive = TRUE)
dir.create("inst/images/cell_types", showWarnings = FALSE, recursive = TRUE)

###########################
### Helper function     ###
###########################

# Generic NBLAST score distribution plot
plot_nblast_distribution <- function(data, score_col, region_filter = NULL,
                                     dataset_name, filename,
                                     threshold_borderline = 0.2,
                                     threshold_good = 0.3,
                                     threshold_confident = 0.5) {

  if (!is.null(region_filter)) {
    data <- data %>% dplyr::filter(grepl(region_filter, region))
  }

  data <- data %>%
    dplyr::mutate(
      score = as.numeric(.data[[score_col]]),
      assignments_match = dplyr::case_when(
        is.na(score) ~ "no_nblast",
        score < threshold_borderline ~ "bad",
        score >= threshold_confident ~ "confident",
        score >= threshold_good ~ "good",
        score >= threshold_borderline ~ "borderline",
        TRUE ~ "bad"
      )
    ) %>%
    dplyr::mutate(assignments_match = factor(assignments_match,
                                             levels = c("no_nblast", "bad", "borderline", "good", "confident")))

  col.values <- c("no_nblast" = "lightgrey",
                   "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
                   "borderline" = hemibrainr:::hemibrain_bright_colors[["orange"]],
                   "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
                   "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

  match_summary <- data %>%
    dplyr::group_by(assignments_match) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(percentage = count / sum(count) * 100)

  plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
    geom_bar(stat = "identity", width = 0.5) +
    coord_flip() +
    scale_fill_manual(values = col.values) +
    theme_minimal() +
    labs(
      title = sprintf("BANC-%s top NBLAST scores per BANC neuron", dataset_name),
      x = "", y = "percentage", fill = "assignments match"
    ) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))),
              position = position_stack(vjust = 0.5), size = 3.5)

  plot2 <- ggplot(data, aes(x = score, y = (..count..) / nrow(data))) +
    geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
    theme_minimal() +
    labs(
      title = sprintf("normalised density of BANC-%s top NBLAST scores", dataset_name),
      x = sprintf("BANC-%s NBLAST score", dataset_name),
      y = "density", fill = "assignments match"
    ) +
    scale_x_continuous(breaks = seq(0, max(data$score, na.rm = TRUE), by = 0.1)) +
    ggplot2::geom_vline(xintercept = 0,
                        color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
    ggplot2::geom_vline(xintercept = threshold_borderline,
                        color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
    ggplot2::geom_vline(xintercept = threshold_good,
                        color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
    ggplot2::geom_vline(xintercept = threshold_confident,
                        color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
    ggplot2::annotate("text", x = threshold_borderline, y = Inf,
                      label = sprintf("borderline: >%s", threshold_borderline),
                      vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
    ggplot2::annotate("text", x = threshold_good, y = Inf,
                      label = sprintf("good: >%s", threshold_good),
                      vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
    ggplot2::annotate("text", x = threshold_confident, y = Inf,
                      label = sprintf("confident: >%s", threshold_confident),
                      vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

  combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, top = "")
  print(combined_plot)
  ggsave(plot = combined_plot, filename = filename, width = 18, height = 8, dpi = 300)
  message(sprintf("  Saved: %s", filename))
}

###########################
### MANC NBLAST         ###
###########################

plot_nblast_distribution(
  data = banc.meta, score_col = "manc_nblast",
  region_filter = "neck|vnc",
  dataset_name = "MANC",
  filename = "inst/images/nblast/manc_nblast_match_banc.png"
)

###########################
### maleCNS NBLAST      ###
###########################

plot_nblast_distribution(
  data = banc.meta, score_col = "malecns_nblast",
  region_filter = NULL,
  dataset_name = "maleCNS",
  filename = "inst/images/nblast/malecns_nblast_match_banc.png"
)

###########################
### FAFB NBLAST         ###
###########################

plot_nblast_distribution(
  data = banc.meta, score_col = "fafb_nblast",
  region_filter = "brain|neck|optic",
  dataset_name = "FAFB",
  filename = "inst/images/nblast/fafb_nblast_match_banc.png"
)

###########################
### Hemibrain NBLAST    ###
###########################

plot_nblast_distribution(
  data = banc.meta, score_col = "hemibrain_nblast",
  region_filter = "brain|neck ",
  dataset_name = "Hemibrain",
  filename = "inst/images/nblast/hemibrain_nblast_match_banc.png"
)

###########################
### Mirror NBLAST       ###
###########################

plot_nblast_distribution(
  data = banc.meta, score_col = "banc_nblast",
  region_filter = NULL,
  dataset_name = "Mirror",
  filename = "inst/images/nblast/banc_nblast_match_banc.png"
)

##########################
### Matching inventory ###
##########################

# Verified matches pie chart
bc.cts <- bc %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON", status)) %>%
  dplyr::filter(!is.na(root_id), !is.na(supervoxel_id), !is.na(position),
                root_id != "", root_id != "0") %>%
  dplyr::arrange(fafb_cell_type, manc_cell_type, hemibrain_cell_type, cell_type) %>%
  dplyr::filter(!duplicated(root_id)) %>%
  dplyr::mutate(
    fafb_cell_type = dplyr::case_when(
      grepl("auto", fafb_cell_type) ~ NA, is.na(fafb_match) ~ NA, TRUE ~ fafb_cell_type
    ),
    manc_cell_type = dplyr::case_when(
      grepl("auto", manc_cell_type) ~ NA, is.na(manc_match) ~ NA, TRUE ~ manc_cell_type
    ),
    hemibrain_cell_type = dplyr::case_when(
      grepl("auto", hemibrain_cell_type) ~ NA, is.na(hemibrain_match) ~ NA, TRUE ~ hemibrain_cell_type
    )
  ) %>%
  dplyr::filter(!is.na(fafb_match) | !is.na(manc_match) | !is.na(hemibrain_match)) %>%
  dplyr::select(root_id, supervoxel_id, position, fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::arrange(fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::distinct(root_id, fafb_cell_type, manc_cell_type, hemibrain_cell_type, .keep_all = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, fafb_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, manc_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, hemibrain_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!duplicate_flag) %>%
  dplyr::rename(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id, pt_position = position) %>%
  dplyr::select(-duplicate_flag) %>%
  reshape2::melt(id = c("pt_root_id", "pt_supervoxel_id", "pt_position"),
                 value.name = "tag", variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA", "unknown", "Uknown", "fragment", "Fragment", "None", "none", "no_match"),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity", user_id = 355, valid = "t") %>%
  dplyr::mutate(tag2 = gsub("_", " ", tag2)) %>%
  dplyr::filter(!grepl("auto", tag), !grepl("no_match", tag),
                !grepl("auto", tag2), !grepl("no_match", tag2),
                !is.na(tag), !is.na(tag2)) %>%
  dplyr::distinct(pt_root_id, pt_supervoxel_id, pt_position, tag2, tag, user_id, valid)

processed_data <- bc.cts %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(both = any(grepl("fafb", tag2)) & any(grepl("manc", tag2)),
                both2 = any(grepl("fafb", tag2)) & any(grepl("hemibrain", tag2))) %>%
  dplyr::mutate(tag2 = dplyr::case_when(
    both ~ "manual fafb+manc cell type",
    both2 ~ "manual fafb+hemibrain cell type",
    TRUE ~ tag2
  )) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(bc[, c("root_id", "region")], by = c("pt_root_id" = 'root_id')) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(tag2 = ifelse(is.na(tag2), "unmatched", tag2)) %>%
  dplyr::group_by(tag2, region) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(tag2, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

tag2_colors <- scales::hue_pal()(length(unique(processed_data$tag2)))
names(tag2_colors) <- unique(processed_data$tag2)

processed_data <- processed_data %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = tag2)) +
  facet_wrap(vars(region)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), nudge_x = 0.5, size = 2) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tag2_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "match status", title = "distribution of entries by match status")
print(pie_chart)
ggsave(plot = pie_chart,
       filename = "inst/images/cell_types/bc_verified_match_status_distribution.png",
       width = 12, height = 12, dpi = 300)

# High NBLAST score pie chart
nb.thresh <- 0.5
fw.nblast.scores <- banc.meta.fafb.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score >= nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(fafb_cell_type = match_cell_type) %>%
  dplyr::mutate(valid = "f")
mc.nblast.scores <- banc.meta.manc.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score >= nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(manc_cell_type = match_cell_type) %>%
  dplyr::mutate(valid = "f")
hb.nblast.scores <- banc.meta.hemibrain.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score >= nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(hemibrain_cell_type = match_cell_type) %>%
  dplyr::mutate(valid = "f")

nblast.scores <- fw.nblast.scores %>%
  dplyr::full_join(mc.nblast.scores, by = c("pt_root_id", "pt_supervoxel_id", "pt_position")) %>%
  dplyr::full_join(hb.nblast.scores, by = c("pt_root_id", "pt_supervoxel_id", "pt_position")) %>%
  dplyr::filter(!is.na(pt_root_id), !is.na(pt_supervoxel_id), !is.na(pt_position),
                pt_root_id != "", pt_root_id != "0") %>%
  dplyr::select(pt_root_id, pt_supervoxel_id, pt_position,
                fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  reshape2::melt(id = c("pt_root_id", "pt_supervoxel_id", "pt_position"),
                 value.name = "tag", variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA", "unknown", "Uknown", "fragment", "Fragment", "None", "none", "no_match"),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity", tag = paste0(tag, "?"), user_id = 355, valid = "t") %>%
  dplyr::filter(!is.na(tag), !is.na(tag2))

nblast.scores <- tryCatch(
  bancr::banc_updateids(nblast.scores, root.column = "pt_root_id",
                        supervoxel.column = "pt_supervoxel_id",
                        position.column = "pt_position"),
  error = function(e) {
    # CAVE chunkedgraph Redis OOM under bulk updateids (job 38901939, 2026-05-09).
    # Feathers were just refreshed by banc-nblast-compile.R, so pt_root_id is
    # already current; fall through to the unmigrated frame.
    message("  banc_updateids failed (continuing with feather IDs): ",
            conditionMessage(e))
    nblast.scores
  }
)
nblast.scores <- nblast.scores %>%
  dplyr::group_by(pt_root_id, tag2) %>%
  dplyr::filter(!duplicated(pt_position), !duplicated(pt_root_id), pt_root_id != "0") %>%
  dplyr::ungroup()

processed_data2 <- nblast.scores %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(both = any(grepl("fafb", tag2)) & any(grepl("manc", tag2)),
                both2 = any(grepl("fafb", tag2)) & any(grepl("hemibrain", tag2))) %>%
  dplyr::mutate(tag2 = dplyr::case_when(
    both ~ "manual fafb+manc cell type",
    both2 ~ "manual fafb+hemibrain cell type",
    TRUE ~ tag2
  )) %>%
  dplyr::full_join(bc[, c("root_id", "region")], by = c("pt_root_id" = 'root_id')) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(tag2 = ifelse(is.na(tag2), "unmatched", tag2)) %>%
  dplyr::group_by(tag2, region) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(tag2, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

tag2_colors2 <- scales::hue_pal()(length(unique(processed_data2$tag2)))
names(tag2_colors2) <- unique(processed_data2$tag2)

processed_data2 <- processed_data2 %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

pie_chart2 <- ggplot(processed_data2, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = tag2)) +
  facet_wrap(vars(region)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), nudge_x = 0.5, size = 2) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tag2_colors2) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "match status", title = "distribution of entries by match status")
print(pie_chart2)
ggsave(plot = pie_chart2,
       filename = "inst/images/cell_types/bc_high_nblast_match_status_distribution.png",
       width = 12, height = 12, dpi = 300)

message(sprintf("### banc: NBLAST plotting complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
