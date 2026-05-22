### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

############################
### COMPARE CONNECTIVITY ###
############################
source("banc/banc-startup.R")

# Run locally
banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
banc.connectivity.save.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity"
banc.dropbox.connectivity.save.path <- "/Users/abates/HMS Dropbox/Alexander Bates/neuroanat/connectomes"
banc.path <- "/Users/papers/BANC-project/"

# Read meta
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))

############
### BANC ###
############

# Get our best meta data and matching
bc.orig <- banctable_query() 
bc <- bc %>%
  dplyr::mutate(cell_type = ifelse(grepl("auto\\:", cell_type),NA,cell_type),
                fafb_cell_type = ifelse(grepl("auto\\:",fafb_cell_type),NA,cell_type),
                manc_cell_type = ifelse(grepl("auto\\:",manc_cell_type),NA,cell_type))
bc.fafb <- subset(bc, !is.na(fafb_match))
bc.manc <- subset(bc, !is.na(manc_match))
bc.chosen <- bc %>%
  dplyr::filter(region=="neck_connective",
                !is.na(fafb_match),
                !is.na(manc_match))
bc.chosen.ids <- bc.chosen$root_id
bc.fafb.cts <- na.omit(unique(bc.fafb$fafb_cell_type))
bc.fafb.cts <- intersect(bc.fafb.cts,fw.meta$cell_type)
bc.manc.cts <- na.omit(unique(bc.manc$manc_cell_type))
bc.manc.cts <- intersect(bc.manc.cts,mc.meta$cell_type)
bc.cts <- unique(c(bc.fafb.cts, bc.manc.cts))

# Get connective
banc.el.orig <- banc_edgelist(fetch_all_rows = TRUE)
banc.el <- banc.el.orig %>%
  dplyr::rename(count = n) %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(norm = round(count/post_count,6))
banc.el.chosen <- banc.el %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id)) %>%
  dplyr::filter(pre %in% bc.chosen.ids | post %in% bc.chosen.ids) %>%
  dplyr::left_join(bc[,c("root_id","fafb_cell_type","manc_cell_type")], 
                   by = c("pre"="root_id")) %>%
  dplyr::rename(pre_fafb_cell_type = fafb_cell_type,
                pre_manc_cell_type = manc_cell_type) %>%
  dplyr::left_join(bc[,c("root_id","fafb_cell_type","manc_cell_type")], 
                   by = c("post"="root_id")) %>%
  dplyr::rename(post_fafb_cell_type = fafb_cell_type,
                post_manc_cell_type = manc_cell_type) %>%
  dplyr::filter(pre_fafb_cell_type %in% bc.cts | post_fafb_cell_type %in% bc.cts,
                pre_manc_cell_type %in% bc.cts | post_manc_cell_type %in% bc.cts) %>%
  dplyr::select(-post_pt_root_id, -pre_pt_root_id) %>%
  dplyr::arrange(dplyr::desc(count)) %>%
  dplyr::filter(pre!=post, 
                !is.na(pre_fafb_cell_type)|!is.na(pre_manc_cell_type),
                !is.na(post_fafb_cell_type)|!is.na(post_manc_cell_type)) %>%
  dplyr::mutate(dataset = "BANC") %>%
  dplyr::left_join(bc.chosen[,c("root_id","cell_type")], 
                   by = c("pre"="root_id")) %>%
  dplyr::rename(pre_cell_type = cell_type) %>%
  dplyr::left_join(bc.chosen[,c("root_id","cell_type")], 
                   by = c("post"="root_id")) %>%
  dplyr::rename(post_cell_type = cell_type) %>%
  dplyr::mutate(pre_cell_type = coalesce(pre_cell_type, pre_fafb_cell_type, pre_manc_cell_type),
                post_cell_type = coalesce(post_cell_type, post_fafb_cell_type, post_manc_cell_type))
  
############
### FAFB ###
############

# Read from feather files
franken.el <- arrow::read_feather(file.path(banc.connectivity.save.path,"frankenbrain_v.1.5_edgelist_simple.feather"))
franken.meta <- arrow::read_feather(file.path(banc.connectivity.save.path,"frankenbrain_v.1.5_meta.feather"))

# Filter
franken.chosen.cts <- unique(subset(franken.meta, region=="neck_connective" & dataset == "BANC-FAFB-MANC")$cell_type)
franken.el.chosen <- franken.el %>%
  dplyr::left_join(franken.meta[,c("id","cell_type")], by = c("pre"="id")) %>%
  dplyr::rename(pre_cell_type = cell_type) %>%
  dplyr::left_join(franken.meta[,c("id","cell_type")], by = c("post"="id")) %>%
  dplyr::rename(post_cell_type = cell_type) %>%
  dplyr::filter(pre_cell_type %in% bc.cts,
                post_cell_type %in% bc.cts)  %>%
  dplyr::arrange(dplyr::desc(count)) %>%
  dplyr::filter(pre!=post) %>%
  dplyr::mutate(dataset = "franken") 

####################
### Scatter plot ###
####################

franken.el.chosen.ct <- franken.el.chosen %>%
  dplyr::group_by(pre_cell_type, post_cell_type) %>%
  dplyr::group_by(pre_cell_type, post_cell_type) %>%
  dplyr::mutate(franken_count_mean = mean(count),
                franken_norm_mean = mean(norm)) %>%
  dplyr::distinct(pre_cell_type, post_cell_type,franken_count_mean,franken_norm_mean) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(connection = paste0(pre_cell_type,"_",post_cell_type)) %>%
  dplyr::ungroup()

banc.el.chosen.ct <- banc.el.chosen %>%
  dplyr::group_by(pre_cell_type, post_cell_type) %>%
  dplyr::group_by(pre_cell_type, post_cell_type) %>%
  dplyr::mutate(banc_count_mean = mean(count),
                banc_norm_mean = mean(norm)) %>%
  dplyr::distinct(pre_cell_type, post_cell_type,banc_count_mean,banc_norm_mean) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(connection = paste0(pre_cell_type,"_",post_cell_type)) %>%
  dplyr::ungroup()

el.chosen.ct <- dplyr::left_join(
    franken.el.chosen.ct,
    banc.el.chosen.ct %>%
      dplyr::select(banc_count_mean,banc_norm_mean,connection),
    by = "connection"
  ) %>%
  dplyr::filter(!is.na(banc_count_mean)) %>%
  dplyr::distinct(connection, .keep_all = TRUE)

el.chosen.ct.downstream <- el.chosen.ct %>%
  dplyr::filter(pre_cell_type %in% franken.chosen.cts) %>%
  dplyr::mutate(dataset = dplyr::case_when(
    post_cell_type %in% franken.chosen.cts ~ "neck_connective",
    post_cell_type %in% bc.fafb.cts ~ "FAFB",
    post_cell_type %in% bc.manc.cts ~ "MANC",
    TRUE ~ NA
  )) %>%
  dplyr::filter(!is.na(dataset)) %>%
  dplyr::mutate(direction = "downstream")

el.chosen.ct.upstream <- el.chosen.ct %>%
  dplyr::filter(post_cell_type %in% franken.chosen.cts) %>%
  dplyr::mutate(dataset = dplyr::case_when(
    pre_cell_type %in% franken.chosen.cts ~ "neck_connective",
    pre_cell_type %in% bc.fafb.cts ~ "FAFB",
    pre_cell_type %in% bc.manc.cts ~ "MANC",
    TRUE ~ NA
  )) %>%
  dplyr::filter(!is.na(dataset)) %>%
  dplyr::mutate(direction = "upstream")

el.chosen.ct.updown <- rbind(el.chosen.ct.downstream,
                             el.chosen.ct.upstream) %>%
  dplyr::filter(banc_count_mean >= 10, 
                franken_count_mean >= 10) 

# Create the scatter plot
g1 <- ggplot(el.chosen.ct.updown, aes(x = banc_count_mean, y = franken_count_mean)) +
  geom_point(alpha = 0.6) +  # Add points with some transparency
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  #geom_text_repel(aes(label = connection), size = 2, max.overlaps = 5) +  # Add labels for points
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    formula = y ~ x, 
    parse = TRUE,
    label.x = "left",
    label.y = "top",
    size = 3
  ) +
  facet_wrap(direction ~ dataset, scales = "free") +  # Create facets based on dataset
  labs(
    title = "Comparison of Connection Counts: Franken vs BANC",
    subtitle = "Faceted by Dataset",
    x = "BANC Count Mean (log scale)",
    y = "Franken Count Mean (log scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "lightgray", color = NA)
  )

# Save the plot (optional)
ggsave(plot = g1, filename = file.path(banc.path,"images","connectivity","franken_vs_banc_count_andn.png"), 
       width = 16, height = 8, dpi = 300)

# Create the scatter plot
g2 <- ggplot(el.chosen.ct.updown, aes(x = banc_norm_mean, y = franken_norm_mean)) +
  geom_point(alpha = 0.6) +  # Add points with some transparency
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  #geom_text_repel(aes(label = connection), size = 2, max.overlaps = 5) +  # Add labels for points
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    formula = y ~ x, 
    parse = TRUE,
    label.x = "left",
    label.y = "top",
    size = 3
  ) +
  facet_wrap(direction ~ dataset, scales = "free") +  # Create facets based on dataset
  labs(
    title = "Comparison of Connection Norms: Franken vs BANC",
    subtitle = "Faceted by Dataset",
    x = "BANC Norm Mean (log scale)",
    y = "Franken Norm Mean (log scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "lightgray", color = NA)
  )

# Save the plot (optional)
ggsave(plot = g2, filename = file.path(banc.path,"images","connectivity","franken_vs_banc_norm_andn.png"), 
       width = 16, height = 8, dpi = 300)

# Analysis of outliers
outliers <- el.chosen.ct.updown %>%
  dplyr::mutate(norm_diff = banc_norm_mean-franken_norm_mean,
                count_diff = banc_count_mean-franken_count_mean) %>%
  dplyr::arrange(dplyr::desc(norm_diff))

# Create histogram for norm_diff
p1 <- ggplot(outliers, aes(x = norm_diff)) +
  geom_histogram(binwidth = 0.005, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Norm Difference",
       x = "Norm Difference",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Create histogram for count_diff
p2 <- ggplot(outliers, aes(x = count_diff)) +
  geom_histogram(binwidth = 50, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Count Difference",
       x = "Count Difference",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Arrange both plots side by side
g3 <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(filename = file.path(banc.path,"images","connectivity","franken_vs_banc_diff_andn.png"), 
       plot = gridExtra::arrangeGrob(p1, p2, ncol = 2), 
       width = 12, height = 6, dpi = 300)

# What pairs are in the largest norm diff?
largest.outliers <- outliers %>%
  dplyr::filter(norm_diff >= quantile(outliers$norm_diff,0.95)) %>%
  dplyr::left_join(franken.meta[,c("super_class","cell_type")], by = c("post_cell_type"="cell_type")) %>%
  dplyr::rename(post_super_class=super_class) %>%
  dplyr::left_join(franken.meta[,c("super_class","cell_type")], by = c("pre_cell_type"="cell_type")) %>%
  dplyr::rename(pre_super_class=super_class)

# Function to create binned data
create_binned_data <- function(data, cell_type_col) {
  data %>%
    mutate(norm_diff_bin = cut(norm_diff, 
                               breaks = seq(min(norm_diff), max(norm_diff) + 0.05, by = 0.05),
                               include.lowest = TRUE)) %>%
    group_by(norm_diff_bin, direction, dataset, !!sym(cell_type_col)) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    group_by(norm_diff_bin, direction, dataset) %>%
    mutate(total = sum(count),
           proportion = count / total)
}

# Create binned data for pre and post cell types
pre_binned_data <- create_binned_data(largest.outliers, "pre_super_class") %>%
  dplyr::rename(super_class = pre_super_class) %>%
  dplyr::filter(direction=="upstream")
post_binned_data <- create_binned_data(largest.outliers, "post_super_class") %>%
  dplyr::rename(super_class = post_super_class) %>%
  dplyr::filter(direction=="downstream")

# Function to create the plot
create_plot <- function(data, title) {
  ggplot(data, aes(x = norm_diff_bin, y = proportion, fill = super_class)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = title,
         x = "norm_diff (binned)",
         y = "Proportion",
         fill = "Super Class") +
    theme_minimal() +
    facet_grid(direction ~ dataset) +
    coord_flip()
}

# Create plots
pre_plot <- create_plot(pre_binned_data, "Distribution of pre_super_class across norm_diff bins")
post_plot <- create_plot(post_binned_data, "Distribution of post_super_class across norm_diff bins")

# Save the plot (optional)
ggsave(plot = post_plot, filename = file.path(banc.path,"images","connectivity","andn_post_norm_diff_pre_cell_type_distribution.png"), 
       width = 16, height = 8, dpi = 300)
ggsave(plot = pre_plot, filename = file.path(banc.path,"images","connectivity","andn_pre_norm_diff_pre_cell_type_distribution.png"), 
       width = 16, height = 8, dpi = 300)


################################
### Cosine similarity matrix ###
################################
  
# Step 2: Create a connectivity matrix for each cell type
conn_matrix <- franken.el.chosen %>%
  # Create a row for each pre and post connection
  bind_rows(
    select(., cell_type = pre_cell_type, partner = post_cell_type, count),
    select(., cell_type = post_cell_type, partner = pre_cell_type, count)
  ) %>%
  dplyr::filter(!is.na(partner)) %>%
  # Sum the counts for each cell type-partner pair
  group_by(cell_type, partner) %>%
  dplyr::summarise(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  # Create the matrix
  pivot_wider(names_from = partner, values_from = total_count, values_fill = 0) %>%
  dplyr::filter(cell_type %in% bc.cts, !is.na(cell_type)) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Step 3: Calculate cosine similarity
cosine_sim <- lsa::cosine(conn_matrix)

# Step 4: Visualize the cosine similarity matrix
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Create the heatmap
ph <- pheatmap(cosine_sim,
         color = color_palette,
         main = "Cosine Similarity of Cell Type Connectivity",
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename = "cosine_similarity_heatmap.png",
         width = 10,
         height = 8)

# If you want to save the plot
png(filename = file.path(banc.path,"images","cosine_similarity_plot.png"), 
    width = 12, height = 10)
plot(ph)
dev.off()