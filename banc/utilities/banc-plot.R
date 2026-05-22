#' banc-plot — Diagnostic plots for cross-dataset matching + NT predictions.
#'
#' Density / beeswarm / heatmap diagnostics for cross-dataset matching and
#' neurotransmitter prediction sanity checks. Manual exploration script.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/{flywire,manc,fanc}_meta.csv`
#'   - SeaTable `banc_meta`
#'   - `franken_meta()`
#'
#' @section Writes:
#'   - per-block diagnostic plots (consult the script)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`

##############################
### ORGANISE BANC MATCHING ###
##############################
library(dplyr)
library(tidyr)
library(pheatmap)
library(grid)
library(ggbeeswarm)
library(ggpubr)
library(stats)
source("banc/banc-startup.R")

# get meta
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.cts <- unique(fw.meta$cell_type)

# get data
bc.all <- banctable_query() %>%
  dplyr::filter(!grepl("GLIA|NOT_A_NEURON|NOT_NEURON|DELETE|TRACHEA",status))
bc <- subset(bc.all, grepl("ascending|descending", super_class))

# Get franken.meta
franken.meta <- franken_meta()

######################################
#### NEUROTRANSMITTER PREDICTIONS ####
######################################
poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin")

# Create density plot faceted by neurotransmitter_predicted
g.nt <- ggplot(bc.all %>%
                 dplyr::filter(!is.na(neurotransmitter_predicted),
                               neurotransmitter_predicted %in% poss.nts) , 
               aes(x = neurotransmitter_score, fill = neurotransmitter_predicted)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ neurotransmitter_predicted, ncol = 1) + 
  scale_fill_manual(values = paper.cols) +
  labs(
    title = "neurotransmitter score distribution by neuron-level predicted transmitter",
    x = "neurotransmitter score",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )

# Save
ggsave(plot = g.nt,
       filename = "inst/images/nt/banc_neurotransmitter_score_by_neurotransmitter_predicted.png", 
       width = 12, height = 12, dpi = 300)

# Set up the color palette
heatmap.cols <- colorRampPalette(c("grey95", "#FC6882"))
breaks <- seq(from = 0, to = 1, by = 0.01)
h.cols <- heatmap.cols(length(breaks))
regions <- na.omit(unique(bc.all$region))

### Confusion matrix versus GT
poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin")
gt.data <- read_csv("/home/ab714/drosophila_neurotransmitters/gt_sources/bates_2024/202502-gt_data.csv") %>%
  dplyr::arrange(known_nt_confidence,known_nt) %>%
  dplyr::filter(known_nt %in% poss.nts) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE)
for(reg in regions){
  
  # get
  bc.nt.comp <- bc.all %>%
    dplyr::filter(region == reg,
                  !is.na(neurotransmitter_predicted),
                  neurotransmitter_score > 0.5,
                  neurotransmitter_predicted %in% poss.nts) %>%
    dplyr::distinct(root_id, cell_type, top_nt, neurotransmitter_predicted) %>%
    dplyr::left_join(gt.data %>%
                       dplyr::distinct(cell_type, .keep_all = TRUE),
                     by = "cell_type") %>%
    dplyr::filter(!is.na(known_nt),
                  !is.na(neurotransmitter_predicted)) %>%
    dplyr::mutate(in_gt = !is.na(known_nt))
  
  # First, create confusion matrices for all data and GT-only data
  # For all predictions
  all_confusion <- bc.nt.comp %>%
    dplyr::filter(!is.na(neurotransmitter_predicted)) %>%  # Ensure predicted NT is not NA
    dplyr::group_by(neurotransmitter_predicted, known_nt) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      id_cols = neurotransmitter_predicted,
      names_from = known_nt,
      values_from = count,
      values_fill = 0
    )
  
  # Extract row names and convert to matrix format
  row_names_all <- all_confusion$neurotransmitter_predicted
  all_confusion_matrix <- all_confusion %>%
    dplyr::select(-neurotransmitter_predicted) %>%
    as.matrix()
  rownames(all_confusion_matrix) <- row_names_all
  
  # For GT-only predictions
  gt_confusion <- bc.nt.comp %>%
    dplyr::filter(!is.na(neurotransmitter_predicted), in_gt == TRUE) %>%
    dplyr::group_by(neurotransmitter_predicted, known_nt) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      id_cols = neurotransmitter_predicted,
      names_from = known_nt,
      values_from = count,
      values_fill = 0
    )
  
  # Match dimensions to all_confusion_matrix
  gt_confusion_matrix <- matrix(0, nrow = nrow(all_confusion_matrix), ncol = ncol(all_confusion_matrix))
  rownames(gt_confusion_matrix) <- rownames(all_confusion_matrix)
  colnames(gt_confusion_matrix) <- colnames(all_confusion_matrix)
  
  # Fill in the GT data where available
  for (i in 1:nrow(gt_confusion)) {
    pred_nt <- gt_confusion$neurotransmitter_predicted[i]
    row_idx <- which(rownames(gt_confusion_matrix) == pred_nt)
    if (length(row_idx) > 0) {  # If this prediction exists in the all_confusion matrix
      for (col in colnames(gt_confusion)[2:ncol(gt_confusion)]) {
        col_idx <- which(colnames(gt_confusion_matrix) == col)
        if (length(col_idx) > 0) {  # If this known_nt exists in the all_confusion matrix
          gt_confusion_matrix[row_idx, col_idx] <- gt_confusion[[col]][i]
        }
      }
    }
  }
  
  # Calculate row and column totals
  all_row_totals <- rowSums(all_confusion_matrix)
  all_col_totals <- colSums(all_confusion_matrix)
  gt_row_totals <- rowSums(gt_confusion_matrix)
  gt_col_totals <- colSums(gt_confusion_matrix)
  
  # Add totals to the matrices
  all_confusion_with_totals <- cbind(all_confusion_matrix, "total" = all_row_totals)
  all_confusion_with_totals <- rbind(all_confusion_with_totals, "total" = c(all_col_totals, sum(all_confusion_matrix)))
  
  gt_confusion_with_totals <- cbind(gt_confusion_matrix, "total" = gt_row_totals)
  gt_confusion_with_totals <- rbind(gt_confusion_with_totals, "total" = c(gt_col_totals, sum(gt_confusion_matrix)))
  
  # Calculate normalized matrices
  all_norm <- sweep(all_confusion_matrix, 1, all_row_totals, "/")
  all_norm[is.nan(all_norm)] <- 0
  
  gt_norm <- sweep(gt_confusion_matrix, 1, gt_row_totals, "/")
  gt_norm[is.nan(gt_norm)] <- 0
  
  # Create a normalized matrix with totals (setting totals rows/cols to NA for coloring)
  all_norm_with_totals <- matrix(NA, 
                                 nrow = nrow(all_confusion_with_totals), 
                                 ncol = ncol(all_confusion_with_totals))
  all_norm_with_totals[1:(nrow(all_norm_with_totals)-1), 
                       1:(ncol(all_norm_with_totals)-1)] <- all_norm
  rownames(all_norm_with_totals) <- rownames(all_confusion_with_totals)
  colnames(all_norm_with_totals) <- colnames(all_confusion_with_totals)
  
  # Create combined labels with counts and percentages
  number_percent_labels <- matrix("", 
                                  nrow = nrow(all_confusion_with_totals),
                                  ncol = ncol(all_confusion_with_totals))
  
  # Fill in the main part of the matrix (excluding totals)
  for (i in 1:(nrow(all_confusion_matrix))) {
    for (j in 1:(ncol(all_confusion_matrix))) {
      all_count <- all_confusion_matrix[i, j]
      all_percent <- round(all_norm[i, j] * 100, 1)
      
      gt_count <- gt_confusion_matrix[i, j]
      gt_percent <- if(gt_row_totals[i] > 0) round(gt_norm[i, j] * 100, 1) else 0
      
      # Format: "number, percentage (GT: number, percentage)"
      if (all_count > 0) {
        if (gt_count > 0) {
          number_percent_labels[i, j] <- paste0(all_count, ", ", all_percent, "%\n",
                                                "(GT: ", gt_count, ", ", gt_percent, "%)")
        } else {
          number_percent_labels[i, j] <- paste0(all_count, ", ", all_percent, "%")
        }
      }
    }
  }
  
  # Fill in the totals row and column
  for (i in 1:(nrow(all_confusion_matrix))) {
    total_text <- paste0(all_row_totals[i])
    if (gt_row_totals[i] > 0) {
      total_text <- paste0(total_text, "\n(GT: ", gt_row_totals[i], ")")
    }
    number_percent_labels[i, ncol(number_percent_labels)] <- total_text
  }
  
  for (j in 1:(ncol(all_confusion_matrix))) {
    total_text <- paste0(all_col_totals[j])
    if (gt_col_totals[j] > 0) {
      total_text <- paste0(total_text, "\n(GT: ", gt_col_totals[j], ")")
    }
    number_percent_labels[nrow(number_percent_labels), j] <- total_text
  }
  
  # Set the overall total
  total_text <- paste0(sum(all_confusion_matrix))
  if (sum(gt_confusion_matrix) > 0) {
    total_text <- paste0(total_text, "\n(GT: ", sum(gt_confusion_matrix), ")")
  }
  number_percent_labels[nrow(number_percent_labels), ncol(number_percent_labels)] <- total_text
  
  rownames(number_percent_labels) <- rownames(all_confusion_with_totals)
  colnames(number_percent_labels) <- colnames(all_confusion_with_totals)
 
  # Create custom annotation colors for the totals row and column
  annotation_rows <- data.frame(
    is_total = factor(c(rep("no", nrow(all_confusion_matrix)), "yes"))
  )
  rownames(annotation_rows) <- rownames(all_confusion_with_totals)
  
  annotation_cols <- data.frame(
    is_total = factor(c(rep("no", ncol(all_confusion_matrix)), "yes"))
  )
  rownames(annotation_cols) <- colnames(all_confusion_with_totals)
  
  annotation_colors <- list(
    is_total = c(yes = "grey80", no = "white")
  )
  
  # Create the heatmap with custom background colors for totals
  pheatmap(
    mat = all_norm_with_totals,  # Use normalized matrix for coloring
    color = h.cols,
    breaks = breaks,
    show_colnames = TRUE,
    show_rownames = TRUE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    display_numbers = number_percent_labels,  # Show both counts and percentages
    fontsize_number = 8,  # Smaller font to accommodate more text
    number_color = "black",
    main = paste0(reg, " predicted vs known neurotransmitter comparison"),
    border_color = NA,
    annotation_row = annotation_rows,
    annotation_col = annotation_cols,
    annotation_colors = annotation_colors,
    annotation_legend = FALSE,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    na_col = "grey80",
    #filename = sprintf("inst/images/nt/banc_%s_neurotransmitter_prediction_versus_known_comparison.png",reg)
  )
    
}

### Confusion matrix versus FAFB/MANC predictions
for(reg in regions){
  
  # wrangle
  bc.nt.comp <- bc.all %>%
    dplyr::filter(region == reg,
                  !is.na(top_nt), 
                  !is.na(neurotransmitter_predicted),
                  neurotransmitter_predicted %in% poss.nts) %>%
    dplyr::distinct(root_id, cell_type, top_nt, neurotransmitter_predicted) %>%
    dplyr::left_join(gt.data,
                     by = "cell_type") %>%
    dplyr::mutate(in_gt = !is.na(known_nt))
  
  # Create a confusion matrix from bc.nt.comp
  confusion_matrix <- bc.nt.comp %>%
    # Count occurrences of each combination
    dplyr::group_by(top_nt, neurotransmitter_predicted) %>%
    dplyr::summarize(count = n(), .groups = "drop") %>%
    # Convert to wide format for the matrix
    tidyr::pivot_wider(
      id_cols = neurotransmitter_predicted,
      names_from = top_nt,
      values_from = count,
      values_fill = 0
    )
  
  # Extract row names and convert to matrix format
  row_names <- confusion_matrix$neurotransmitter_predicted
  confusion_matrix <- confusion_matrix %>%
    dplyr::select(-neurotransmitter_predicted) %>%
    as.matrix()
  rownames(confusion_matrix) <- row_names
  
  # Calculate row and column totals
  row_totals <- rowSums(confusion_matrix)
  col_totals <- colSums(confusion_matrix)
  
  # Add totals to the matrix
  confusion_matrix_with_totals <- cbind(confusion_matrix, "total" = row_totals)
  confusion_matrix_with_totals <- rbind(confusion_matrix_with_totals, "total" = c(col_totals, sum(confusion_matrix)))
  
  # Calculate proportions for the main confusion matrix (excluding totals)
  confusion_matrix_norm <- sweep(confusion_matrix, 1, row_totals, "/")
  # Replace NaN with 0
  confusion_matrix_norm[is.nan(confusion_matrix_norm)] <- 0
  
  # Create a normalized matrix with totals (setting totals rows/cols to NA for coloring)
  confusion_matrix_norm_with_totals <- matrix(NA, 
                                              nrow = nrow(confusion_matrix_with_totals), 
                                              ncol = ncol(confusion_matrix_with_totals))
  confusion_matrix_norm_with_totals[1:(nrow(confusion_matrix_norm_with_totals)-1), 
                                    1:(ncol(confusion_matrix_norm_with_totals)-1)] <- confusion_matrix_norm
  rownames(confusion_matrix_norm_with_totals) <- rownames(confusion_matrix_with_totals)
  colnames(confusion_matrix_norm_with_totals) <- colnames(confusion_matrix_with_totals)
  
  # Create combined labels with counts and percentages
  number_percent_labels <- matrix("", 
                                  nrow = nrow(confusion_matrix_with_totals),
                                  ncol = ncol(confusion_matrix_with_totals))
  
  # Fill in the main part of the matrix (excluding totals)
  for (i in 1:(nrow(confusion_matrix))) {
    for (j in 1:(ncol(confusion_matrix))) {
      count <- confusion_matrix[i, j]
      if (count > 0) {
        percent <- round(confusion_matrix_norm[i, j] * 100, 1)
        number_percent_labels[i, j] <- paste0(count, "\n(", percent, "%)")
      }
    }
  }
  
  # Fill in the totals row and column
  for (i in 1:(nrow(confusion_matrix))) {
    number_percent_labels[i, ncol(number_percent_labels)] <- row_totals[i]
  }
  for (j in 1:(ncol(confusion_matrix))) {
    number_percent_labels[nrow(number_percent_labels), j] <- col_totals[j]
  }
  number_percent_labels[nrow(number_percent_labels), ncol(number_percent_labels)] <- sum(confusion_matrix)
  
  rownames(number_percent_labels) <- rownames(confusion_matrix_with_totals)
  colnames(number_percent_labels) <- colnames(confusion_matrix_with_totals)
  
  # Create custom annotation colors for the totals row and column
  annotation_rows <- data.frame(
    is_total = factor(c(rep("no", nrow(confusion_matrix)), "yes"))
  )
  rownames(annotation_rows) <- rownames(confusion_matrix_with_totals)
  
  annotation_cols <- data.frame(
    is_total = factor(c(rep("no", ncol(confusion_matrix)), "yes"))
  )
  rownames(annotation_cols) <- colnames(confusion_matrix_with_totals)
  
  annotation_colors <- list(
    is_total = c(yes = "grey80", no = "white")
  )
  
  # Create the heatmap with custom background colors for totals
  pheatmap(
    mat = confusion_matrix_norm_with_totals,  # Use normalized matrix with NA for totals
    color = h.cols,
    breaks = breaks,
    show_colnames = TRUE,
    show_rownames = TRUE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    display_numbers = number_percent_labels,  # Show both raw counts and percentages
    fontsize_number = 10,
    number_color = "black",
    main = paste0(reg, " neurotransmitter prediction comparison"),
    border_color = NA,
    annotation_row = annotation_rows,
    annotation_col = annotation_cols,
    annotation_colors = annotation_colors,
    annotation_legend = FALSE,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    na_col = "grey70"
    #filename = sprintf("inst/images/nt/banc_%s_neurotransmitter_prediction_versus_franken_comparison.png",reg)
  )
}

### Compare against DCV density

# Wrangle
bc.nt.comp <- bc.all %>%
  dplyr::left_join(franken.meta %>%
                     dplyr::distinct(cell_type, .keep_all = TRUE) %>%
                     dplyr::select(cell_type, soma_dcv_density),
                   by = "cell_type") %>%
  dplyr::filter(!is.na(top_nt), 
                !is.na(soma_dcv_density),
                soma_dcv_density>-1,
                !is.na(neurotransmitter_predicted),
                neurotransmitter_predicted %in% poss.nts) %>%
  dplyr::distinct(root_id, cell_type, soma_dcv_density, neurotransmitter_predicted)

# Add log10 transformed and percentile versions of soma_dcv_density
bc.nt.comp <- bc.nt.comp %>%
  dplyr::mutate(
    # Log10 transformation (add small constant to handle zeros)
    soma_dcv_log10 = log10(soma_dcv_density + 1e-10),
    
    # Percentile ranking (0-100) with 4 significant figures
    soma_dcv_percentile = signif(
      stats::ecdf(soma_dcv_density)(soma_dcv_density) * 100, 
      digits = 4
    )
  )

# Calculate summary statistics for each neurotransmitter type using log10 values
nt_stats_log10 <- bc.nt.comp %>%
  dplyr::group_by(neurotransmitter_predicted) %>%
  dplyr::summarize(
    mean_log10 = mean(soma_dcv_log10, na.rm = TRUE),
    median_log10 = stats::median(soma_dcv_log10, na.rm = TRUE),
    count = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_log10))

# Reorder factor levels based on median log10 density
bc.nt.comp <- bc.nt.comp %>%
  dplyr::mutate(
    neurotransmitter_predicted = factor(neurotransmitter_predicted, 
                                        levels = nt_stats_log10$neurotransmitter_predicted)
  )

# Create a visualization with log10 transformed values
p_log10 <- ggplot2::ggplot(bc.nt.comp, 
                           ggplot2::aes(x = neurotransmitter_predicted, y = soma_dcv_log10)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = neurotransmitter_predicted), 
                        alpha = 0.7, 
                        outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.25, color = "black") +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  ggplot2::geom_text(data = nt_stats_log10, 
                     ggplot2::aes(x = neurotransmitter_predicted, 
                                  y = min(bc.nt.comp$soma_dcv_log10) - 0.5,
                                  label = paste0("n=", count)),
                     vjust = 0) +
  ggplot2::scale_fill_manual(values = paper.cols) +  # Use paper.cols for neurotransmitter coloring
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Distribution of log10 Soma DCV Density by Predicted Neurotransmitter",
    x = "Predicted Neurotransmitter",
    y = "log10(Soma DCV Density)",
    fill = "Neurotransmitter"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# Create a visualization with percentiles
p_percentile <- ggplot2::ggplot(bc.nt.comp, 
                                ggplot2::aes(x = neurotransmitter_predicted, y = soma_dcv_percentile)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = neurotransmitter_predicted), 
                        alpha = 0.7, 
                        outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.25, color = "black") +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  ggplot2::scale_fill_manual(values = paper.cols) +  # Use paper.cols for neurotransmitter coloring
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Percentile Distribution of Soma DCV Density by Predicted Neurotransmitter",
    x = "Predicted Neurotransmitter",
    y = "Percentile Rank (0-100)",
    fill = "Neurotransmitter"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# Create density plots for log10 transformed values
p_log10_density <- ggplot2::ggplot(bc.nt.comp, 
                                   ggplot2::aes(x = soma_dcv_log10, fill = neurotransmitter_predicted)) +
  ggplot2::geom_density(alpha = 0.7) +
  ggplot2::facet_wrap(~ neurotransmitter_predicted, scales = "free_y") +
  ggplot2::scale_fill_manual(values = paper.cols) +  # Use paper.cols for neurotransmitter coloring
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Density Distribution of log10 Soma DCV by Neurotransmitter",
    x = "log10(Soma DCV Density)",
    y = "Density",
    fill = "Neurotransmitter"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# Create density plots for percentile values
p_percentile_density <- ggplot2::ggplot(bc.nt.comp, 
                                        ggplot2::aes(x = soma_dcv_percentile, fill = neurotransmitter_predicted)) +
  ggplot2::geom_density(alpha = 0.7) +
  ggplot2::facet_wrap(~ neurotransmitter_predicted, scales = "free_y") +
  ggplot2::scale_fill_manual(values = paper.cols) +  # Use paper.cols for neurotransmitter coloring
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Density Distribution of Soma DCV Percentiles by Neurotransmitter",
    x = "Percentile Rank (0-100)",
    y = "Density",
    fill = "Neurotransmitter"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# Compare original, log10, and percentile distributions
p_compare <- ggplot2::ggplot(bc.nt.comp) +
  # Original values
  ggplot2::geom_histogram(ggplot2::aes(x = soma_dcv_density), 
                          bins = 50, fill = "steelblue", alpha = 0.5) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Distribution of Original Soma DCV Density Values",
    x = "Soma DCV Density",
    y = "Count"
  )

p_compare_log10 <- ggplot2::ggplot(bc.nt.comp) +
  # Log10 values
  ggplot2::geom_histogram(ggplot2::aes(x = soma_dcv_log10), 
                          bins = 50, fill = "darkred", alpha = 0.5) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Distribution of log10 Transformed Soma DCV Density Values",
    x = "log10(Soma DCV Density)",
    y = "Count"
  )

p_compare_percentile <- ggplot2::ggplot(bc.nt.comp) +
  # Percentile values
  ggplot2::geom_histogram(ggplot2::aes(x = soma_dcv_percentile), 
                          bins = 50, fill = "darkgreen", alpha = 0.5) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Distribution of Percentile Ranks of Soma DCV Density",
    x = "Percentile Rank (0-100)",
    y = "Count"
  )

# Perform ANOVA on log10 values
anova_log10 <- stats::aov(soma_dcv_log10 ~ neurotransmitter_predicted, data = bc.nt.comp)
anova_log10_pval <- summary(anova_log10)[[1]]["Pr(>F)"][1,1]

# Perform ANOVA on percentile values
anova_percentile <- stats::aov(soma_dcv_percentile ~ neurotransmitter_predicted, data = bc.nt.comp)
anova_percentile_pval <- summary(anova_percentile)[[1]]["Pr(>F)"][1,1]

# Format p-values using standard format() function
formatted_log10_pval <- ifelse(anova_log10_pval < 0.001, 
                               "p < 0.001", 
                               paste0("p = ", format(anova_log10_pval, digits = 3)))

formatted_percentile_pval <- ifelse(anova_percentile_pval < 0.001, 
                                    "p < 0.001", 
                                    paste0("p = ", format(anova_percentile_pval, digits = 3)))

# Add p-values to plots
p_log10 <- p_log10 + 
  ggplot2::labs(subtitle = paste0("ANOVA ", formatted_log10_pval))

p_percentile <- p_percentile + 
  ggplot2::labs(subtitle = paste0("ANOVA ", formatted_percentile_pval))

# Arrange comparison plots without labels
comparison_distrib <- ggpubr::ggarrange(p_compare, p_compare_log10, p_compare_percentile, 
                                        ncol = 1, nrow = 3)

# Arrange main plots without labels
main_plots <- ggpubr::ggarrange(p_log10, p_percentile, 
                                ncol = 1, nrow = 2)

density_plots <- ggpubr::ggarrange(p_log10_density, p_percentile_density, 
                                   ncol = 1, nrow = 2)

# Display plots
print(comparison_distrib)
print(main_plots)
print(density_plots)

# Save
ggsave(plot = comparison_distrib,
       filename = "inst/images/nt/banc_dcv_transforms.png", 
       width = 24, height = 8, dpi = 300)
ggsave(plot = main_plots,
       filename = "inst/images/nt/banc_nt_by_dcv_boxplots.png", 
       width = 24, height = 8, dpi = 300)
ggsave(plot = density_plots,
       filename = "inst/images/nt/banc_nt_by_dcv_densityplots.png", 
       width = 24, height = 8, dpi = 300)

# Statistical results using log10 values
tukey_log10 <- stats::TukeyHSD(anova_log10)
print("Tukey HSD results for log10 values:")
print(tukey_log10)

# Create a summary table of neurotransmitter comparisons with log10 values
summary_table <- bc.nt.comp %>%
  dplyr::group_by(neurotransmitter_predicted) %>%
  dplyr::summarize(
    n = dplyr::n(),
    mean_dcv = mean(soma_dcv_density, na.rm = TRUE),
    median_dcv = stats::median(soma_dcv_density, na.rm = TRUE),
    mean_log10 = mean(soma_dcv_log10, na.rm = TRUE),
    median_log10 = stats::median(soma_dcv_log10, na.rm = TRUE),
    mean_percentile = mean(soma_dcv_percentile, na.rm = TRUE),
    median_percentile = stats::median(soma_dcv_percentile, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_log10))

print(summary_table)

### Summary of most common confusions

# Get disagreements
bc.nt.comp <- bc.all %>%
  dplyr::filter(!is.na(neurotransmitter_predicted)) %>%
  dplyr::distinct(root_id, cell_type, top_nt, neurotransmitter_predicted) %>%
  dplyr::left_join(gt.data,
                   by = "cell_type") %>%
  dplyr::filter(!is.na(known_nt),
                !is.na(top_nt)) %>%
  dplyr::mutate(top_nt!=neurotransmitter_predicted|neurotransmitter_predicted!=known_nt) %>%
  dplyr::filter(! neurotransmitter_predicted %in% c("histamine","tyramine"))

###############################
#### SENSORY NEURON STATUS ####
###############################

# Process the data
plot_data <- bc.all %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::filter(grepl("sensory",cell_class)|grepl("afferent|sensory",super_class)) %>%
  dplyr::group_by(nerve, status, region) %>% 
  dplyr::mutate(l2_cable_length_um = as.numeric(l2_cable_length_um),
                status = subtract_status(status,c("TRACING_ISSUE","TOO_SMALL"), invert = TRUE)) %>%
  dplyr::mutate(status =  dplyr::case_when(
    proofread == "TRUE" ~ append_status(status,c("PROOFREAD")),
    TRUE ~ append_status(status,c("NOT_PROOFREAD"))
  )) %>%
  dplyr::mutate(status =  dplyr::case_when(
    !is.na(fafb_match) ~ append_status(status,c("MATCHED")),
    !is.na(manc_match) ~ append_status(status,c("MATCHED")),
    !is.na(hemibrain_match) ~ append_status(status,c("MATCHED")),
    is.na(l2_cable_length_um)|l2_cable_length_um<100 ~ "TOO_SMALL",
    TRUE ~ status
  )) %>%
  dplyr::mutate(nerve = gsub("auto:|_exit:|entry:|None|none|_R|_L|,.*","",nerve)) %>%
  dplyr::mutate(nerve = dplyr::case_when(
    is.na(nerve) ~ "unknown",
    nerve=="" ~ "unknown",
    TRUE ~ nerve
  )) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(nerve, region) %>%  
  dplyr::mutate(total = sum(count),
                proportion = count/total,
                label = ifelse(count >= 25, 
                               scales::percent(proportion, accuracy = 0.1),
                               ""))

# Define custom colors
status_colors <- c(
  "MATCHED,PROOFREAD" = "#54BCD1",
  "MATCHED,PROOFREAD,TRACING_ISSUE" = "#1BB6AF",
  "MATCHED,NOT_PROOFREAD" = "#FC6882",
  "MATCHED,NOT_PROOFREAD,TRACING_ISSUE" = "#C70E7B", 
  "PROOFREAD" = "#8FDA04",
  "PROOFREAD,TRACING_ISSUE" =  "#FBBB48",
  "TRACING_ISSUE" = "#EF7C12" ,
  "NOT_PROOFREAD" = "#EE4244",
  "NOT_PROOFREAD,TRACING_ISSUE" = "#C23A4B",
  "TOO_SMALL" = "lightgrey"
)

# Create the plot
g.sens <- ggplot(plot_data, aes(x = nerve, y = count, fill = status)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            size = 3) +
  facet_grid( ~ region, scales = "free_x", space = "free") + 
  scale_fill_manual(values = status_colors) +
  theme_minimal() +
  labs(title = "sensory neuron distribution by nerve and status",
       x = "nerve",
       y = "count",
       fill = "status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(1, "lines"))

# Display the plot
print(g.sens)

# Save
ggsave(plot = g.sens,
       filename = "inst/images/tracing/banc_sensory_neuron_status.png", 
       width = 24, height = 8, dpi = 300)

# Filter the data for sensory neurons
sensory_neurons <- bc.all %>%
  dplyr::filter(grepl("sensory", cell_class) | grepl("afferent|sensory", super_class)) %>%
  dplyr::mutate(l2_cable_length_um = as.numeric(l2_cable_length_um)) %>%
  dplyr::mutate(l2_cable_length_um = dplyr::case_when(
    is.na(l2_cable_length_um) ~ 0,
    l2_cable_length_um > 5000 ~ 5000,
    TRUE ~ l2_cable_length_um
  )) %>%
  dplyr::mutate(
    side = dplyr::case_when(
      is.na(side) ~ "unknown",
      side == "" ~ "unknown",
      !side %in% c("midline", "unknown", "right", "left") ~ "midline",
      TRUE ~ side
    ),
    region = ifelse(is.na(region) | region == "", "unknown", region),
    proofread = dplyr::case_when(
      is.na(proofread) ~ "FALSE",
      proofread == "" ~ "FALSE",
      TRUE ~ proofread
    )
  ) %>%
  dplyr::select(root_id, supervoxel_id, position, l2_cable_length_um, proofread, 
                super_class, nerve, side, region,
                manc_match, manc_nblast_match, manc_nblast, fafb_match, fafb_nblast_match, fafb_nblast)#%>%
  #dplyr::filter(side %in% c("left", "right"))

# Neurons that should be marked as proofread
proofread <- banc_backbone_proofread()
proofread.ids <- unique(as.character(proofread$pt_root_id))
sensory_new_proofread <- sensory_neurons %>%
  dplyr::filter(!root_id %in% proofread.ids,
                (!is.na(fafb_match)|!is.na(manc_match))) %>%
  dplyr::distinct(root_id, supervoxel_id, position) %>%
  dplyr::mutate(proofread = 't',
                user_id = 355,
                valid = "t")
readr::write_csv(sensory_new_proofread, file.path(banc.meta.save.path,"sensory_new_proofread.csv"))
#system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast", banc.meta.save.path, "sensory_new_proofread.csv"))

# Neurons that need to be traced
peripheral.nerves <- banc_peripheral_nerves()
sensory_needs_tracing <- sensory_neurons %>%
  dplyr::filter(!root_id %in% proofread.ids, 
                !is.na(supervoxel_id)) %>%
  dplyr::arrange(l2_cable_length_um) %>%
  dplyr::distinct(root_id, supervoxel_id, position)
readr::write_csv(sensory_needs_tracing, file.path(banc.meta.save.path,"sensory_needs_tracing.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast", banc.meta.save.path, "sensory_needs_tracing.csv"))

# Create the histogram with facets and color by proofread status
g.hist <- ggplot(sensory_neurons, aes(x = l2_cable_length_um, fill = proofread)) +
  geom_histogram(binwidth = 50, position = "stack") +
  geom_vline(xintercept = 100, color = "black", linetype = 2, size = 0.5) +
  facet_grid(side ~ region, scales = "free_y", space = "free") +
  scale_fill_manual(values = c("TRUE" = "#8FDA04", "FALSE" = "#EE4244")) +
  theme_minimal() +
  labs(title = "distribution of cable length for sensory neurons by side and region",
       x = "cable length (µm)",
       y = "count",
       fill = "proofread status") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom")

# Display the plot
print(g.hist)

# Save the plot
ggsave(plot = g.hist,
       filename = "inst/images/tracing/banc_sensory_neuron_cable_length.png", 
       width = 24, height = 8, dpi = 300)

# Process the data for pie chart
processed_data <- plot_data %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

# Calculate the positions for the labels
processed_data <- processed_data %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of sensory neurons by status")

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/tracing/banc_sensory_neuron_status_pie.png", 
       width = 12, height = 12, dpi = 300)

# Process the data for pie chart
processed_data_midbrain <- plot_data %>%
  dplyr::filter(region=="midbrain") %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  ) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart_midbrain <- ggplot(processed_data_midbrain, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of sensory central brain neurons by status")

# Display the plot
print(pie_chart_midbrain)

# Save
ggsave(plot = pie_chart_midbrain,
       filename = "inst/images/tracing/banc_sensory_neuron_status_midbrain_pie.png", 
       width = 12, height = 12, dpi = 300)

# Process the data for pie chart
processed_data_vnc <- plot_data %>%
  dplyr::filter(region%in%c("vnc","ventral_nerve_cord")) %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  ) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart_vnc <- ggplot(processed_data_vnc, 
                        aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, 
                            fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of sensory VNC neurons by status")

# Display the plot
print(pie_chart_vnc)

# Save
ggsave(plot = pie_chart_vnc,
       filename = "inst/images/tracing/banc_sensory_neuron_status_vnc_pie.png", 
       width = 12, height = 12, dpi = 300)

################################
#### EFFERENT NEURON STATUS ####
################################

# Process the data
plot_data <- bc.all %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::filter(grepl("efferent|motor",cell_class)|grepl("motor|efferent",super_class)) %>%
  dplyr::group_by(nerve, status, region) %>% 
  dplyr::mutate(l2_cable_length_um = as.numeric(l2_cable_length_um),
                status = subtract_status(status,c("TRACING_ISSUE","TOO_SMALL"), invert = TRUE)) %>%
  dplyr::mutate(status =  dplyr::case_when(
    proofread == "TRUE" ~ append_status(status,c("PROOFREAD")),
    TRUE ~ append_status(status,c("NOT_PROOFREAD"))
  )) %>%
  dplyr::mutate(status =  dplyr::case_when(
    !is.na(fafb_match) ~ append_status(status,c("MATCHED")),
    !is.na(manc_match) ~ append_status(status,c("MATCHED")),
    !is.na(hemibrain_match) ~ append_status(status,c("MATCHED")),
    is.na(l2_cable_length_um)|l2_cable_length_um<100 ~ "TOO_SMALL",
    TRUE ~ status
  )) %>%
  dplyr::mutate(nerve = gsub("auto:|_exit:|entry:|None|none|_R|_L|,.*","",nerve)) %>%
  dplyr::mutate(nerve = dplyr::case_when(
    is.na(nerve) ~ "unknown",
    nerve=="" ~ "unknown",
    TRUE ~ nerve
  )) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(nerve, region) %>%  
  dplyr::mutate(total = sum(count),
                proportion = count/total,
                label = ifelse(count >= 25, 
                               scales::percent(proportion, accuracy = 0.1),
                               ""))

# Define custom colors
status_colors <- c(
  "MATCHED,PROOFREAD" = "#54BCD1",
  "MATCHED,PROOFREAD,TRACING_ISSUE" = "#1BB6AF",
  "MATCHED,NOT_PROOFREAD" = "#FC6882",
  "MATCHED,NOT_PROOFREAD,TRACING_ISSUE" = "#C70E7B", 
  "PROOFREAD" = "#8FDA04",
  "PROOFREAD,TRACING_ISSUE" =  "#FBBB48",
  "TRACING_ISSUE" = "#EF7C12" ,
  "NOT_PROOFREAD" = "#EE4244",
  "NOT_PROOFREAD,TRACING_ISSUE" = "#C23A4B",
  "TOO_SMALL" = "lightgrey"
)

# Create the plot
g.motors <- ggplot(plot_data, aes(x = nerve, y = count, fill = status)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            size = 3) +
  facet_grid( ~ region, scales = "free_x", space = "free") + 
  scale_fill_manual(values = status_colors) +
  theme_minimal() +
  labs(title = "efferent neuron distribution by nerve and status",
       x = "nerve",
       y = "count",
       fill = "status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(1, "lines"))

# Display the plot
print(g.motors)

# Save
ggsave(plot = g.motors,
       filename = "inst/images/tracing/banc_efferent_neuron_status.png", 
       width = 24, height = 8, dpi = 300)

# Filter the data for efferent neurons
efferent_neurons <- bc.all %>%
  dplyr::filter(grepl("motor|efferent", cell_class) | grepl("motor|efferent", super_class)) %>%
  dplyr::mutate(l2_cable_length_um = as.numeric(l2_cable_length_um)) %>%
  dplyr::mutate(l2_cable_length_um = dplyr::case_when(
    is.na(l2_cable_length_um) ~ 0,
    l2_cable_length_um > 5000 ~ 5000,
    TRUE ~ l2_cable_length_um
  )) %>%
  dplyr::mutate(
    side = dplyr::case_when(
      is.na(side) ~ "unknown",
      side == "" ~ "unknown",
      !side %in% c("midline", "unknown", "right", "left") ~ "midline",
      TRUE ~ side
    ),
    region = ifelse(is.na(region) | region == "", "unknown", region),
    proofread = dplyr::case_when(
      is.na(proofread) ~ "FALSE",
      proofread == "" ~ "FALSE",
      TRUE ~ proofread
    )
  ) %>%
  dplyr::select(root_id, supervoxel_id, position, l2_cable_length_um, proofread, 
                super_class, nerve, side, region,
                manc_match, manc_nblast_match, manc_nblast, fafb_match, fafb_nblast_match, fafb_nblast)#%>%
#dplyr::filter(side %in% c("left", "right"))

# Neurons that should be marked as proofread
proofread.ids <- unique(as.character(proofread$pt_root_id))
efferent_new_proofread <- efferent_neurons %>%
  dplyr::filter(!root_id %in% proofread.ids,
                (!is.na(fafb_match)|!is.na(manc_match))) %>%
  dplyr::distinct(root_id, supervoxel_id, position) %>%
  dplyr::mutate(proofread = 't',
                user_id = 355,
                valid = "t")
readr::write_csv(efferent_new_proofread, file.path(banc.meta.save.path,"efferent_new_proofread.csv"))

# Neurons that need to be traced
peripheral.nerves <- banc_peripheral_nerves()
efferent_needs_tracing <- efferent_neurons %>%
  dplyr::filter(!root_id %in% proofread.ids, 
                !is.na(supervoxel_id)) %>%
  dplyr::arrange(l2_cable_length_um) %>%
  dplyr::distinct(root_id, supervoxel_id, position)
readr::write_csv(efferent_needs_tracing, file.path(banc.meta.save.path,"efferent_needs_tracing.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast", banc.meta.save.path, "efferent_needs_tracing.csv"))

# Create the histogram with facets and color by proofread status
g.hist <- ggplot(efferent_neurons, aes(x = l2_cable_length_um, fill = proofread)) +
  geom_histogram(binwidth = 50, position = "stack") +
  geom_vline(xintercept = 100, color = "black", linetype = 2, size = 0.5) +
  facet_grid(side ~ region, scales = "free_y", space = "free") +
  scale_fill_manual(values = c("TRUE" = "#8FDA04", "FALSE" = "#EE4244")) +
  theme_minimal() +
  labs(title = "distribution of cable length for efferent neurons by side and region",
       x = "cable length (µm)",
       y = "count",
       fill = "proofread status") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom")

# Display the plot
print(g.hist)

# Save the plot
ggsave(plot = g.hist,
       filename = "inst/images/tracing/banc_efferent_neuron_cable_length.png", 
       width = 24, height = 8, dpi = 300)

# Process the data for pie chart
processed_data <- plot_data %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

# Calculate the positions for the labels
processed_data <- processed_data %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of efferent neurons by status")

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/tracing/banc_efferent_neuron_status_pie.png", 
       width = 12, height = 12, dpi = 300)

# Process the data for pie chart
processed_data_midbrain <- plot_data %>%
  dplyr::filter(region=="midbrain") %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  ) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart_midbrain <- ggplot(processed_data_midbrain, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of efferent brain neurons by status")

# Display the plot
print(pie_chart_midbrain)

# Save
ggsave(plot = pie_chart_midbrain,
       filename = "inst/images/tracing/banc_efferent_neuron_status_midbrain_pie.png", 
       width = 12, height = 12, dpi = 300)

# Process the data for pie chart
processed_data_vnc <- plot_data %>%
  dplyr::filter(region%in%c("vnc","ventral_nerve_cord")) %>%
  dplyr::group_by(status) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(status, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  ) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Create the plot
pie_chart_vnc <- ggplot(processed_data_vnc, 
                        aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, 
                            fill = status)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 3) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = status_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "status", title = "distribution of efferent vnc neurons by status")

# Display the plot
print(pie_chart_vnc)

# Save
ggsave(plot = pie_chart_vnc,
       filename = "inst/images/tracing/banc_efferent_neuron_status_vnc_pie.png", 
       width = 12, height = 12, dpi = 300)

################
#### LENGTH ####
################

# Get example fragments
fragment <- "720575941550979134"
small.fragment <- "720575941602628746"
fragment.skel <- banc_read_l2skel(fragment)
small.fragment.skel <- banc_read_l2skel(small.fragment)

# Add log(l2_cable_length_um) to the dataframe
smallest.size <- min(bc.all$l2_cable_length_um[!is.na(bc.all$l2_cable_length_um)&bc.all$l2_cable_length_um>0])
bc.plot <- bc.all %>%
  dplyr::filter(proofread=="TRUE") %>%
  dplyr::mutate(
    l2_cable_length_um = ifelse(is.na(l2_cable_length_um)|l2_cable_length_um==0,smallest.size,l2_cable_length_um),
    log_l2_cable_length = log(l2_cable_length_um),
    region = ifelse(is.na(region),"undetermined",region),
    cell_type = gsub("auto:","",cell_type))

# Define color vectors (replace with your desired hex colors)
region_colors <- c(
  "optic" = hemibrainr:::hemibrain_bright_colors[["marine"]],
  "midbrain" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
  "vnc" = hemibrainr:::hemibrain_bright_colors[["orange"]],
  "neck_connective" = hemibrainr:::hemibrain_bright_colors[["pink"]],
  "undetermined" = "grey70"
)
cell_type_colors <- c(DNa01 = "#AF6125",
                      MBON01 = "pink",
                      pC1a = "green",
                      SA1 = "darkgreen",
                      Mi4 = "#800080", 
                      ORN_DP1m = "#FBBB48",
                      SApp01 = "#008080",
                      KCab = "#D72000",
                      fragment = "black",
                      small_fragment = "grey30")
cell_types_of_interest <- names(cell_type_colors)

# Calculate mean and SD for specified cell types
mean_sd_data <- bc.plot %>%
  dplyr::filter(cell_type %in% cell_types_of_interest) %>%
  dplyr::group_by(cell_type, region) %>%
  dplyr::summarise(mean_length = mean(log_l2_cable_length, na.rm = TRUE),
                   sd_length = sd(log_l2_cable_length, na.rm = TRUE),
                   .groups = "drop")

# Add fragment example
frag_data <- data.frame(
  region = "midbrain",
  cell_type = c("fragment","small_fragment"),
  mean_length = log(c(summary(fragment.skel)$cable.length/1000,summary(small.fragment.skel)$cable.length/1000)),
  sd_length = c(0,0)
)
mean_sd_data <- rbind(mean_sd_data,frag_data)

# Create the plot
p <- ggplot2::ggplot(bc.plot, ggplot2::aes(x = log_l2_cable_length)) +
  ggplot2::geom_histogram(ggplot2::aes(fill = region), 
                          binwidth = 0.5, 
                          color = "white") +
  ggplot2::facet_wrap(~ region, scales = "free_y") +
  ggplot2::scale_fill_manual(values = region_colors, na.value = "gray") +
  ggplot2::labs(title = "distribution of log l2 cable length by region",
                x = "log l2 cable length",
                y = "count") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  # Center the title

# Add mean and SD for specified cell types
g.l2 <- p + ggplot2::geom_vline(data = mean_sd_data, 
                                ggplot2::aes(xintercept = mean_length, color = cell_type),
                                linetype = "solid") +
  ggplot2::geom_vline(xintercept = log(100), 
                      color = "red",
                      linetype = "dashed") +
  ggplot2::geom_errorbarh(data = mean_sd_data,
                          ggplot2::aes(x = mean_length, 
                                       xmin = mean_length - sd_length,
                                       xmax = mean_length + sd_length,
                                       y = 0, color = cell_type),
                          height = 5) +
  ggplot2::scale_color_manual(values = cell_type_colors)

# Display the plot
print(g.l2)

# Save
ggsave(plot = g.l2,
       filename = "inst/images/metrics/banc_log_l2_cable_length_um.png", width = 24, height = 8, dpi = 300)

#################
#### MATCHED ####
#################

# Step 1: Process the data
processed_data <- bc.all %>%
  dplyr::mutate(region = ifelse(is.na(region),"undetermined",region),
                side = ifelse(is.na(side),"undetermined",side),
                proofread = ifelse(is.na(proofread),"FALSE",proofread),
                matched = dplyr::case_when(
                  !is.na(fafb_match) ~"MATCHED",
                  !is.na(manc_match) ~"MATCHED",
                  !is.na(hemibrain_match) ~"MATCHED",
                  grepl("auto",cell_type) ~ "NBLAST",
                  !is.na(cell_type) ~"MATCHED",
                  TRUE ~ "UNMATCHED"
                )) %>%
  dplyr::filter(proofread == "TRUE", region != "undetermined") %>%
  dplyr::mutate(matched = ifelse(is.na(matched),"undetermined",matched),
                region = ifelse(is.na(region),"undetermined",region)) %>%
  dplyr::group_by(matched, region) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  # Calculate positions and labels
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count/sum(count),
    label = paste0(count,"\n(",scales::percent(proportion, accuracy = 0.1), ")"),
    ypos = cumsum(proportion) - proportion/2
  )

# Step 2: Define custom colors
match_colors <- c(
  "MATCHED" = hemibrainr:::hemibrain_bright_colors[["marine"]],
  "UNMATCHED" = 'grey70',
  "NBLAST" = hemibrainr:::hemibrain_bright_colors[["darkyellow"]]
)

# Create the plot
pie_chart <- ggplot(processed_data, aes(x=1, y=proportion, fill=matched)) +
  geom_col(width=1,position="fill", color = "white") +
  ggrepel::geom_text_repel(aes(x = 1.7,label = label), 
            position = position_stack(vjust = 0.5),
            point.size = NA,
            size = 3) +
  coord_polar("y", start=0) + 
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = match_colors) +
  labs(fill = "match status") +
  ggtitle("proportion of made matches for proofread neurons by region") +
  facet_wrap(~region) 

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/metrics/banc_proofread_match_overview.png", 
       width = 8, height = 8, dpi = 300)

###################
#### PROOFREAD ####
###################

# Step 1: Process the data
processed_data <- bc.all %>%
  dplyr::mutate(region = ifelse(is.na(region),"undetermined",region),
                side = ifelse(is.na(side),"undetermined",side),
                proofread = ifelse(is.na(proofread),"FALSE",proofread)
                ) %>%
  dplyr::mutate(region = ifelse(is.na(region),"undetermined",region)) %>%
  dplyr::filter(region != "undetermined") %>%
  dplyr::group_by(proofread, region) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count/sum(count),
    label = paste0(count,"\n(",scales::percent(proportion, accuracy = 0.1), ")"),
    ypos = cumsum(proportion) - proportion/2
  )

# Step 2: Define custom colors
proofread_colors <- c(
  "TRUE" = hemibrainr:::hemibrain_bright_colors[["pink"]],
  "FALSE" = 'lightgrey'
)

# Create the plot
processed_data$region <- factor(processed_data$region, levels = unique(sort(processed_data$region)))
pie_chart <- ggplot(processed_data, aes(x=1, y=proportion, fill=proofread)) +
  geom_col(width=1,position="fill", color = "white") +
  geom_text(aes(x = 1.7,label = label), 
            position = position_stack(vjust = 0.5),
            size = 3) +
  coord_polar("y", start=0) + 
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = proofread_colors) +
  labs(fill = "proofread status") +
  ggtitle("proportion of proofread neurons by region") +
  facet_wrap(~region) 

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/metrics/banc_proofread_overview.png", 
       width = 8, height = 8, dpi = 300)

#################
#### REGIONS ####
#################

# Step 1: Process the data
missing.optic.right <- nrow(subset(fw.meta, side=="right"&super_class=="optic"))
missing.optic.left <- nrow(subset(fw.meta, side=="right"&super_class=="optic"))
missing.optic.right <- missing.optic.right - nrow(subset(bc.all, side=="right"&region=="optic"))
missing.optic.left <- missing.optic.left - nrow(subset(bc.all, side=="left"&region=="optic"))
missing.optic <- data.frame(
  region = "optic",
  side = c(rep("left",missing.optic.left),rep("right",missing.optic.right)),
  proofread = "FALSE"
)
processed_data <- bc.all %>%
  plyr::rbind.fill(missing.optic) %>%
  dplyr::filter(proofread == "TRUE") %>%
  dplyr::mutate(region = ifelse(is.na(region),"undetermined",region),
                side = ifelse(is.na(side),"undetermined",side)) %>%
  dplyr::filter(side != "undetermined", region != "undetermined") %>%
  dplyr::group_by(region, side) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    group = interaction(region,side, sep="_")) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(group,"\n",count," (",scales::percent(proportion, accuracy = 0.1), ")")
  )

# Step 2: Define custom colors (replace with your desired hex codes)
region_colors <- c(
  "optic" = hemibrainr:::hemibrain_bright_colors[["marine"]],
  "midbrain" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
  "vnc" = hemibrainr:::hemibrain_bright_colors[["orange"]],
  "neck_connective" = hemibrainr:::hemibrain_bright_colors[["pink"]],
  "undetermined" = "lightgrey"
)

# Step 3: Create a function to lighten colors
lighten_color <- function(color, factor = 1.4) {
  col_rgb <- col2rgb(color)
  col_rgb <- pmin(col_rgb * factor, 255)
  rgb(t(col_rgb), maxColorValue = 255)
}

# Step 4: Create the color palette
color_palette <- c(
  sapply(region_colors, function(color) color),
  sapply(region_colors, lighten_color, factor = 1.2),
  sapply(region_colors, lighten_color, factor = 1.4)
)
names(color_palette) <- c(
  paste0(names(region_colors), "_left"),
  paste0(names(region_colors), "_right"),
  paste0(names(region_colors), "_undetermined")
)

# Calculate the positions for the labels
processed_data <- processed_data %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05  # This will place the labels just outside the pie
  )

# Step 5: Create the plot
pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = group)) +
  geom_rect() +
  geom_text(aes(label = label, x=4, y=position), 
            nudge_x = .5,
            size = 2) +
  geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette) +  # Assuming you've defined color_palette as before
  xlim(c(0, 4.5)) +  # Adjust this to leave more space for labels
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "region and side", title = "distribution of proofread neurons by region and side")

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/metrics/banc_proofread_overview.png", width = 8, height = 8, dpi = 300)

#################
#### OVERALL ####
#################

# Step 1: Prepare the data for the Sankey diagram
# nuclei <- banc_nuclei()
# nuclei.extra <- subset(nuclei, !nucleus_id%in% bc.all$nucleus_id)
# nuclei.points <- nat::xyzmatrix(nuclei.extra$nucleus_position_nm)
# lrdiffs <- bancr:::banc_lr_position(nuclei.points,units = "nm")
# sides <- ifelse(lrdiffs>0,"right","left")
# nuclei.extra$side <- sides
missing.optic.right <- nrow(subset(fw.meta, side=="right"&super_class=="optic"))
missing.optic.left <- nrow(subset(fw.meta, side=="right"&super_class=="optic"))
missing.optic.right <- missing.optic.right - nrow(subset(bc.all, side=="right"&region=="optic"))
missing.optic.left <- missing.optic.left - nrow(subset(bc.all, side=="left"&region=="optic"))
missing.optic <- data.frame(
  region = "optic",
  side = c(rep("left",missing.optic.left),rep("right",missing.optic.right)),
  proofread = "FALSE"
)
sankey_data <- bc.all %>%
  plyr::rbind.fill(missing.optic) %>%
  dplyr::mutate(region = ifelse(is.na(region),"undetermined",region),
                side = ifelse(is.na(side),"undetermined",side),
                proofread = ifelse(is.na(proofread),"FALSE",proofread),
                matched = dplyr::case_when(
                  !is.na(fafb_match) ~"MATCHED",
                  !is.na(manc_match) ~"MATCHED",
                  !is.na(hemibrain_match) ~"MATCHED",
                  grepl("auto",cell_type) ~ "NBLAST",
                  !is.na(cell_type) ~"MATCHED",
                  TRUE ~ "UNMATCHED"
                )) %>%
  dplyr::filter(side != "undetermined",
                region != "undetermined") %>%
  dplyr::mutate(
    side = paste0(region,"_",side),
    proofread = proofread) %>%
  dplyr::select(region, side, proofread, matched) %>%  # Changed order here
  dplyr::mutate(id = dplyr::row_number()) %>%
  tidyr::pivot_longer(cols = c(region, side, proofread, matched), 
                      names_to = "step", 
                      values_to = "node") %>%
  dplyr::mutate(step = factor(step, levels = c("region", "side", "proofread", "matched"))) %>%  # Changed order here
  dplyr::arrange(id, step) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(next_step = dplyr::lead(step),
                next_node = dplyr::lead(node)) %>%
  dplyr::ungroup() %>%
  dplyr::count(step, next_step, node, next_node, name = "value") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(colour = gsub(".*_","",node),
                label = paste0(node, " (", value, ")"))

# Step 2: Create the Sankey plot
p <- ggplot2::ggplot(sankey_data, ggplot2::aes(x = step, 
                                               next_x = next_step, 
                                               node = node, 
                                               next_node = next_node,
                                               fill = colour,
                                               color = colour,
                                               label = node,
                                               value = value,
                                               label = label)) +
  ggsankey::geom_sankey(flow.alpha = 0.5, 
                        node.color = "black") +
  ggsankey::geom_sankey(ggplot2::aes(value = value), flow.alpha = 0.5, node.color = "black") +
  ggsankey::geom_sankey_label(size = 3, color = "white", fill = "black") +
  scale_fill_manual(values = paper.cols) +
  scale_color_manual(values = paper.cols) +
  ggsankey::theme_sankey(base_size = 18) +
  labs(title = "BANC status sankey diagram",
       x = NULL) +
  theme(legend.position = "none")

# Step 3: Display the plot
print(p)

# Save
ggsave(plot = p,
       filename = "inst/images/metrics/banc_overview_sankey.png", width = 11, height = 6, dpi = 300)

###############
#### FLOW  ####
###############

# Step 1: Process the data
missing.optic.right.intrinsic <- nrow(subset(fw.meta, side=="right"&super_class=="optic")) - nrow(subset(bc.all, side=="right"&region=="optic"&!grepl("sensory",super_class)))
missing.optic.left.intrinsic  <- nrow(subset(fw.meta, side=="right"&super_class=="optic")) - nrow(subset(bc.all, side=="left"&region=="optic"&!grepl("sensory",super_class)))
missing.optic.right.sensory <- nrow(subset(fw.meta, side=="right"&super_class=="optic")) - nrow(subset(bc.all, side=="right"&region=="optic"&grepl("sensory",super_class)))
missing.optic.left.sensory  <- nrow(subset(fw.meta, side=="right"&super_class=="optic")) - nrow(subset(bc.all, side=="left"&region=="optic"&grepl("sensory",super_class)))
missing.optic <- data.frame(
  region = "optic",
  side = c(rep("right",missing.optic.right.intrinsic),rep("right",missing.optic.right.sensory),
           rep("left",missing.optic.left.intrinsic),rep("left",missing.optic.left.sensory)),
  proofread = "FALSE",
  flow = c(rep("optic_intrinsic",missing.optic.right.intrinsic),rep("optic_afferent",missing.optic.right.sensory),
           rep("optic_intrinsic",missing.optic.left.intrinsic),rep("optic_afferent",missing.optic.left.sensory))
)
processed_data <- bc.all %>%
  plyr::rbind.fill(missing.optic) %>%
  dplyr::mutate(flow = ifelse(is.na(flow),"undetermined",flow),
                proofread = ifelse(is.na(proofread),"FALSE",proofread),
                proofread = dplyr::case_when(
                  !is.na(fafb_match) ~"MATCHED",
                  !is.na(manc_match) ~"MATCHED",
                  !is.na(hemibrain_match) ~"MATCHED",
                  grepl("auto",cell_type) ~ "NBLAST",
                  !is.na(cell_type) ~"MATCHED",
                  TRUE ~ proofread
                )) %>%
  dplyr::mutate(super_class=gsub("auto:","",super_class),
                cell_class=gsub("auto:","",cell_class),
                cell_type=gsub("auto:","",cell_type)) %>%
  dplyr::mutate(flow = dplyr::case_when(
    cell_class == "visual_centrifugal" ~ "midbrain_to_optic",
    cell_class == "visual_projection" ~ "optic_to_midbrain",
    super_class == "sensory_ascending" ~ "midbrain_afferent",
    region == "midbrain" & cell_class=="AN" & grepl("sensory", super_class) ~ "midbrain_afferent",
    region == "midbrain" & grepl("sensory|afferent", super_class) ~ "midbrain_afferent",
    region == "optic" & grepl("sensory|afferent", super_class) ~ "optic_afferent",
    region == "vnc" & grepl("sensory|afferent", super_class) ~ "vnc_afferent",
    region == "midbrain" & grepl("motor|efferent", cell_class) ~ "brain_efferent",
    region == "midbrain" & grepl("motor|efferent", super_class) ~ "brain_efferent",
    region == "vnc" & grepl("motor|efferent", cell_class) ~ "vnc_efferent",
    region == "vnc" & grepl("motor|efferent", super_class) ~ "vnc_efferent",
    region == "vnc" & grepl("glia", super_class) ~ "vnc_glia",
    region == "midbrain" & grepl("glia", super_class) ~ "midbrain_glia",
    region == "optic" & grepl("glia", super_class) ~ "optic_glia",
    grepl("ascending|descending", super_class) & cell_class == "AN" ~ "ascending",
    grepl("ascending|descending", super_class) & cell_class == "DN" ~ "descending",
    grepl("ascending|descending", super_class) & cell_class %in% c("descending") ~ "descending",
    grepl("ascending|descending", super_class) & cell_class %in% c("ascending") ~ "ascending",
    region == "vnc" & grepl("intrinsic|central", super_class) ~ "vnc_intrinsic",
    region == "optic" & grepl("intrinsic|central", super_class) ~ "optic_intrinsic",
    region == "midbrain" & grepl("intrinsic|central", super_class) ~ "midbrain_intrinsic",
    region == "optic" ~ "optic_intrinsic",
    region == "vnc" ~ "vnc_intrinsic",
    region == "midbrain" ~ "midbrain_intrinsic",
    TRUE ~ NA
  )) %>%
  dplyr::group_by(proofread, flow) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(flow) %>%
  dplyr::mutate(
    proportion = count/sum(count),
    label = paste0(count,"\n(",scales::percent(proportion, accuracy = 0.1), ")"),
    ypos = cumsum(proportion) - proportion/2
  )

# Create the plot
processed_data$flow <- factor(processed_data$flow, levels = unique(sort(processed_data$flow)))
pie_chart <- ggplot(processed_data, aes(x=1, y=proportion, fill=proofread)) +
  geom_col(width=1,position="fill", color = "white") +
  ggrepel::geom_text_repel(aes(x = 1.7,label = label), 
                           position = position_stack(vjust = 0.5),
                           point.size = NA,
                           segment.colour = NA,
                           size = 3) +
  coord_polar("y", start=0) + 
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = paper.cols) +
  labs(fill = "proofread status") +
  ggtitle("proportion of proofread neurons by flow") +
  facet_wrap(~flow, nrow = 3) 

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/metrics/banc_flow_match_overview.png", 
       width = 12, height = 12, dpi = 300)

#################
#### CREDIT  ####
#################

# Step 1: Process the data
processed_data <- bc.all %>%
  dplyr::mutate(super_class=gsub("auto:","",super_class),
                cell_class=gsub("auto:","",cell_class),
                cell_type=gsub("auto:","",cell_type)) %>%
  dplyr::mutate(super_class = ifelse(is.na(super_class),"undetermined",super_class),
                proofread = ifelse(is.na(proofread),"FALSE",proofread)) %>%
  dplyr::filter(proofread == "TRUE") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_type_source = paste(sort(unique(tolower(unlist(strsplit(cell_type_source,split=","))))),collapse=",")) %>%
  dplyr::ungroup() %>%
  tidyr::separate_rows(cell_type_source, sep = ", ") %>%
  dplyr::group_by(cell_type_source, super_class) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(super_class) %>%
  dplyr::mutate(
    proportion = count/sum(count),
    label = paste0(count,"\n(",scales::percent(proportion, accuracy = 0.1), ")"),
    ypos = cumsum(proportion) - proportion/2
  )

# Create the plot
processed_data$super_class <- factor(processed_data$super_class, levels = unique(sort(processed_data$super_class),decreasing=TRUE))
processed_data$cell_type_source <- factor(processed_data$cell_type_source, levels = unique(sort(processed_data$cell_type_source),decreasing=TRUE))
pie_chart <- ggplot(processed_data, aes(x=1, y=proportion, fill=cell_type_source)) +
  geom_col(width=1,position="fill", color = "white") +
  ggrepel::geom_text_repel(aes(x = 1.7,label = label),
                           position = position_stack(vjust = 0.5),
                           point.size = NA,
                           segment.colour = NA,
                           size = 3) +
  coord_polar("y", start=0) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(fill = "group") +
  ggtitle("proportion of cell type annotations by group") +
  facet_wrap(~super_class, nrow = 3)

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/metrics/banc_match_credit_overview.png", 
       width = 26, height = 12, dpi = 300)

###############
#### CLASS ####
###############
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(grDevices)

# Function to generate LaCroix-like color palette
generate_lacroix_palette <- function(n) {
  lacroix_colors <- c("#F2A5A5", "#F279B2", "#FBCEB1", "#85CDCA", "#B0E0E6", "#D8BFD8", "#FFD700", "#98FB98", "#DDA0DD", "#20B2AA")
  color_ramp <- grDevices::colorRampPalette(lacroix_colors)
  return(color_ramp(n))
}

# Step 1: Process the data
processed_data <- bc.all %>%
  dplyr::mutate(
    cell_class = ifelse(is.na(cell_class), "undetermined", cell_class),
    side = ifelse(is.na(side) | side == "undetermined", "center", side),
    super_class = ifelse(is.na(super_class), "undetermined", super_class)
  ) %>%
  dplyr::mutate(
    cell_class = gsub("auto:", "", cell_class),
    super_class = gsub("auto:", "", super_class),
    side = gsub("auto:", "", side)
  ) %>%
  dplyr::mutate(side = ifelse(side%in%"center","midline",side)) %>%
  dplyr::group_by(super_class, cell_class, side) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(super_class) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = scales::percent(proportion, accuracy = 0.1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cell_class = forcats::fct_reorder(cell_class, count),
    super_class = forcats::fct_reorder(super_class, count, .fun = sum, .desc = TRUE)
  )

# Generate LaCroix-like color palette
n_cell_classes <- length(unique(processed_data$cell_class))
lacroix_palette <- generate_lacroix_palette(n_cell_classes)

# Step 2: Create the grouped and stacked bar chart
bar_chart <- ggplot2::ggplot(processed_data, ggplot2::aes(x = super_class, y = count)) + #fill = cell_class
  ggplot2::geom_bar(stat = "identity", position = "stack", fill = "#85CDCA") +
  # ggplot2::geom_text(ggplot2::aes(label = label), 
  #                    position = ggplot2::position_stack(vjust = 0.5),
  #                    size = 3, color = "white") +
  ggplot2::facet_wrap(~ side, scales = "free_x") +
  ggplot2::coord_flip() +
  #ggplot2::scale_fill_manual(values = lacroix_palette) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 8),
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 16)
  ) +
  ggplot2::labs(
    x = "Super Class",
    y = "Count",
    fill = "Cell Class",
    title = "Distribution of BANC Neurons by Cell Class and Super Class",
    subtitle = "Faceted by Side"
  ) +
  theme(legend.position = "none")

# Display the plot
print(bar_chart)

# Save
ggsave(plot = bar_chart,
       filename = "inst/images/metrics/banc_super_class_overview.png", 
       width = 12, height = 8, dpi = 300)

#####################
#### FAFB NBLAST ####
#####################

# Remove 'auto' prefix and empty entries
bc_cleaned <- bc %>%
  dplyr::filter(!is.na(fafb_nblast_match)) %>%
  dplyr::mutate(
    fafb_nblast = as.numeric(fafb_nblast),
    assignments_match = dplyr::case_when(
      grepl("FAFB_ALT|FAFB_TYPE_CONFLICT",status) ~ "conflict",
      !is.na(fafb_nblast_match) & !is.na(fafb_match) & (fafb_match == fafb_nblast_match) ~ "TRUE",
      !is.na(fafb_nblast_match) & !is.na(fafb_match) & (fafb_match != fafb_nblast_match) ~ "FALSE",
      grepl("TOO_SMALL|TADPOLE",status) ~ "too_small",
      grepl("ISSUE|TRACING",status) ~ "issue",
      grepl("no_match",fafb_png_match) ~"no_match",
      is.na(fafb_match)&!is.na(fafb_nblast_match) ~ "unmatched",
      is.na(fafb_match)&is.na(fafb_nblast_match) ~ "unmatched",
      TRUE ~ "unmatched"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = rev(c("TRUE","FALSE","conflict","too_small","issue","no_match","unmatched"))))

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["pink"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "too_small" = "fragment too small to assess",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "BANC-FAFB matching progress",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of banc_nblast scores
plot2 <- ggplot(bc_cleaned, aes(x = fafb_nblast, 
                                y= (..count..)/nrow(bc_cleaned), 
                                fill = assignments_match)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["pink"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "too_small" = "fragment too small to assess",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-FAFB NBLAST scores",
    x = "BANC-FAFB NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$fafb_nblast, na.rm = TRUE), by = 0.1)) +
  geom_vline(xintercept = 0.3, color = hemibrainr:::hemibrain_bright_colors[["marine"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = 0.3, y = Inf, label = "threshold: 0.3", 
                    vjust = 2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["marine"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/fafb_nblast_match_banc_neck_connective_progress.png", width = 18, height = 8, dpi = 300)

#####################
#### MANC NBLAST ####
#####################

# Remove 'auto' prefix and empty entries
bc_cleaned <- bc %>%
  dplyr::filter(!is.na(manc_nblast_match)) %>%
  dplyr::mutate(
    manc_nblast = as.numeric(manc_nblast),
    assignments_match = dplyr::case_when(
      grepl("MANC_ALT|MANC_TYPE_CONFLICT",status) ~ "conflict",
      !is.na(manc_nblast_match) & !is.na(manc_match) & (manc_match == manc_nblast_match) ~ "TRUE",
      !is.na(manc_nblast_match) & !is.na(manc_match) & (manc_match != manc_nblast_match) ~ "FALSE",
      grepl("TOO_SMALL|TADPOLE",status) ~ "too_small",
      grepl("ISSUE|TRACING",status) ~ "issue",
      grepl("no_match",manc_png_match) ~"no_match",
      is.na(manc_match)&!is.na(manc_nblast_match) ~ "unmatched",
      is.na(manc_match)&is.na(manc_nblast_match) ~ "unmatched",
      TRUE ~ "unmatched"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = rev(c("TRUE","FALSE","conflict","too_small","issue","no_match","unmatched"))))

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["pink"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "too_small" = "fragment too small to assess",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "BANC-MANC matching progress",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of manc_nblast scores
plot2 <- ggplot(bc_cleaned, aes(x = manc_nblast, 
                                y= (..count..)/nrow(bc_cleaned), 
                                fill = assignments_match)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["pink"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "too_small" = "fragment too small to assess",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-MANC NBLAST Scores",
    x = "BANC-MANC NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$manc_nblast, na.rm = TRUE), by = 0.1)) +
  geom_vline(xintercept = 0.3, color = hemibrainr:::hemibrain_bright_colors[["marine"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = 0.3, y = Inf, label = "threshold: 0.3", 
                    vjust = 2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["marine"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/manc_nblast_match_banc_neck_connective.png", width = 18, height = 8, dpi = 300)

#######################
#### mirror NBLAST ####
#######################

# Remove 'auto' prefix and empty entries
bc_cleaned <- bc %>%
  dplyr::filter(!is.na(banc_nblast_match)|!is.na(status)) %>%
  dplyr::mutate(
    banc_nblast = as.numeric(banc_nblast),
    assignments_match = dplyr::case_when(     
      grepl("LR_TYPE_CONFLICT|SIDE_CONFLICT",status) ~ "conflict",
      !is.na(banc_nblast_match) & !is.na(banc_match) & (banc_match == banc_nblast_match) ~ "TRUE",
      !is.na(banc_nblast_match)  & !is.na(banc_match) & (banc_match != banc_nblast_match) ~ "FALSE",
      !is.na(banc_match) ~ "TRUE",
      grepl("TOO_SMALL|TADPOLE",status) ~ "too_small",
      grepl("ISSUE|TRACING|MERGE|INVESTIGATE",status) ~ "issue",
      grepl("no_match",banc_png_match) ~"no_match",
      TRUE ~ "unmatched"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = rev(c("TRUE","FALSE","conflict","too_small","issue","no_match","unmatched"))))

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  as.data.frame()

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["darkyellow"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "too_small" = "fragment too small to assess",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "BANC mirror match  progress",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = ifelse(percentage>1,sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count)),"")), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of banc_nblast scores
plot2 <- ggplot(bc_cleaned, aes(x = banc_nblast, 
                                y= (..count..)/nrow(bc_cleaned), 
                                fill = assignments_match)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["marine"]],
                               "conflict" = hemibrainr:::hemibrain_bright_colors[["pink"]],
                               "no_match" = "grey30",
                               "unmatched" = "lightgrey",
                               "too_small" = hemibrainr:::hemibrain_bright_colors[["orange"]],
                               "issue" = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    labels = c("TRUE" = "top NBLAST match correct", 
                               "FALSE" = "top NBLAST match not correct",
                               "unmatched" = "not yet matched",
                               "no_match" = "match could not be found",
                               "issue" = "tracing issue",
                               "conflict" = "label from another source disagrees")) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC left-right NBLAST Scores",
    x = "BANC left-right NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$banc_nblast, na.rm = TRUE), by = 0.1)) +
  geom_vline(xintercept = 0.3, color = hemibrainr:::hemibrain_bright_colors[["marine"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = 0.3, y = Inf, label = "threshold: 0.3", 
                    vjust = 2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["marine"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/mirror_nblast_match_banc_neck_connective.png", width = 18, height = 8, dpi = 300)

#####################
#### FAFB-match ####
#####################

# Remove 'auto' prefix and empty entries
bc_cleaned <- bc %>%
  dplyr::mutate(
    cell_type = gsub("^auto\\:|\\,.*|\\_.*|\\\n.*|\\*.*", "", cell_type),
    fafb_cell_type = gsub("^auto\\:|\\,.*|\\_.*|\\\n.*|\\*.*", "", fafb_cell_type)
  ) %>%
  dplyr::mutate(fafb_cell_type = dplyr::case_when(
    is.na(fafb_cell_type) ~ fw.meta$cell_type[match(fafb_nblast_match,fw.meta$root_783)],
    fafb_cell_type=="" ~ fw.meta$cell_type[match(fafb_nblast_match,fw.meta$root_783)],
    TRUE ~ fafb_cell_type
  )) %>%
  dplyr::filter(
    proofread=="TRUE",
    !grepl("AN",cell_type),
    cell_type%in%fw.cts,
    !is.na(cell_type) & !is.na(fafb_cell_type),
    cell_type != "" & fafb_cell_type != ""
  ) %>%
  dplyr::mutate(
    correct_assignment = cell_type == fafb_cell_type,
    fafb_nblast = as.numeric(fafb_nblast),
  ) %>%
  dplyr::select(correct_assignment,cell_type,fafb_cell_type,fafb_nblast,proofread)

# Calculate the percentage of correct and incorrect matches
match_summary <- bc_cleaned %>%
  dplyr::group_by(correct_assignment) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)
match_summary_cut <- bc_cleaned %>%
  dplyr::filter(fafb_nblast>=0.3) %>%
  dplyr::group_by(correct_assignment) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of correct and incorrect matches
plot1 <- ggplot(data = match_summary, aes(x = "", y = percentage, fill = correct_assignment)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    labels = c("TRUE" = "correct", "FALSE" = "incorrect")) +
  theme_minimal() +
  labs(
    title = "percentage of correct and incorrect top NBLAST matches",
    x = "",
    y = "percentage",
    fill = "match Result"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 1: Stacked bar plot of correct and incorrect matches
plot2 <- ggplot(data = match_summary_cut, aes(x = "", y = percentage, fill = correct_assignment)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    labels = c("TRUE" = "correct", "FALSE" = "incorrect")) +
  theme_minimal() +
  labs(
    title = "percentage correct among top NBLAST matches >= 0.3",
    x = "",
    y = "percentage",
    fill = "match Result"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of fafb_nblast scores
plot3 <- ggplot(bc_cleaned, aes(x = fafb_nblast, y= (..count..)/nrow(bc_cleaned), 
                                fill = correct_assignment, 
                                group = correct_assignment)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                               "FALSE" = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    labels = c("TRUE" = "correct", "FALSE" = "incorrect")) +
  theme_minimal() +
  labs(
    title = "normalised density of FAFB NBLAST scores",
    x = "FAFB NBLAST Score",
    y = "density",
    fill = "match result"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$fafb_nblast, na.rm = TRUE), by = 0.1)) +
  geom_vline(xintercept = 0.3, color = hemibrainr:::hemibrain_bright_colors[["marine"]], linetype = "dashed", size = 0.5) +
  ggplot2::annotate("text", x = 0.3, y = Inf, label = "threshold: 0.3", 
                    vjust = 2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["marine"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot3, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/fafb_match_banc_neck_connective.png", width = 18, height = 8, dpi = 300)

#####################
#### mirror-match ###
#####################

mirror.correct <- '/n/data1/hms/neurobio/wilson/banc/matching/mirror/correct/'
match.folders <- list.files(mirror.correct, full.names = TRUE)
match.folders <- c(match.folders,'/n/data1/hms/neurobio/wilson/banc/matching/mirror/images/todo/')
mirror.correct.df <- data.frame()
for(mf in match.folders){
  mfn <- basename(mf)
  flist <- sort(list.files(mf, pattern = "png$",recursive = TRUE))
  ids <- gsub(".*_root_id_|_nucleus_id_.*","",flist)
  flist <- flist[!duplicated(ids)]
  count <- length(flist)
  scores <- gsub(".*nblast_score_|_root_id_.*","",flist)
  scores <- as.numeric(scores)
  score.mean <- mean(scores, na.rm = TRUE)
  score.sd <- sd(scores, na.rm = TRUE)
  mirror.correct.df <- rbind(mirror.correct.df, data.frame(match = mfn,
                                                           count = count, 
                                                           score.mean = score.mean,
                                                           score.sd = score.sd))
}

# Reorder factors based on the order in the data
df <- mirror.correct.df

# New order
new_order <- c("1_perfect", "2_confident", "3_good", "4_likely", "5_possible", 
               "6_no_match", "7_does_not_exist", "blue_side_wrong",
               "blue_wrong", "red_wrong", "investigate", "todo")

# Reorder factors based on the new order
df$match <- factor(df$match, levels = new_order)

# Create a new column to differentiate the two divisions of data
df$data_division <- ifelse(df$match %in% c("blue_side_wrong","blue_wrong", "red_wrong",  "investigate", "todo"), "issue", "match")

# Calculate the maximum count for scaling
max_count <- max(df$count, na.rm = TRUE)
total_primary_count <- sum(df$count[df$data_division == "match"])

# Create the plot
g <- ggplot(df, aes(x = match)) +
  geom_hline(yintercept = seq(0, 1, by = 0.1) * max_count, 
             color = "grey80", linetype = "dashed") +
  geom_col(aes(y = count, fill = data_division), width = 0.7, alpha = 0.8) +
  geom_point(aes(y = score.mean * max_count), size = 3, color = hemibrainr:::hemibrain_bright_colors["marine"]) +
  geom_errorbar(aes(y = score.mean * max_count, 
                    ymin = (score.mean - score.sd) * max_count, 
                    ymax = (score.mean + score.sd) * max_count),
                width = 0.2, color = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  scale_fill_manual(values = c("match" = hemibrainr:::hemibrain_bright_colors[["green"]], "issue"  = hemibrainr:::hemibrain_bright_colors[["cerise"]]),
                    name = "") +
  scale_y_continuous(name = "neuron count",
                     sec.axis = sec_axis(~./max_count, name = "blue: normalised NBLAST score", 
                                         breaks = seq(0, 1, by = 0.2),
                                         labels = scales::number_format(accuracy = 0.1))) +
  labs(title = "left-right BANC matches for the neck connective:",
       subtitle = paste("bar height represents no. neurons, points represent mean NBLAST +/- s.d.",
                        "\ntotal reviewed matches:", total_primary_count),
       x = "mirror-match category") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_flip()

# display
print(g)

# save
ggsave(plot = g,
       filename = "inst/images/mirror_match_banc_neck_connective.png", width = 18, height = 8, dpi = 300)





