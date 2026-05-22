### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

######################
### load libraries ###
######################

# Load packages
source("banc/banc-startup.R")
library(uwot)
library(plotly)

# O2
banc.path <- "~/BANC-project/"
banc.connectivity.save.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity"

# Local
banc.path <- "/Users/papers/BANC-project/"
banc.connectivity.save.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity"

####################
### compile data ###
####################

# Read meta data from feather
franken.meta <- arrow::read_feather(file.path(banc.connectivity.save.path,"frankenbrain_v.1.5_meta.feather")) %>%
  dplyr::mutate(region = dplyr::case_when(
    grepl("neck",region) ~ region,
    grepl("^LB|^MX|^MD",hemilineage) ~ "GNG",
    grepl("^aPhN|^MxLbN|^MD|^ON$",nerve) ~ "GNG",
    grepl("^GNG|^DNg",cell_type) ~ "GNG",
    super_class == "visual_projection" ~ "optic",
    super_class == "visual_centrifugal" ~ "optic",
    TRUE ~ region
  ))
seed.neurons <- subset(franken.meta, seed!=""&(grepl("sensory",super_class)|grepl("sensory",cell_class)))$id

# Read influence data
csvs <- list.files(file.path(banc.path,"data/influence/frankenbrain_v1.5/cell_function/"),
                   pattern="csv",
                   full.names = TRUE)
influence.df <- data.frame()
for(csv in csvs){
  if(grepl("motor|unknown|mixed|chemosensory_tactile|proprioceptive_tactile|chemosensory_proprioceptive",csv)){
    next
  }
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = col_character()))
  data <- data[,c(2,4)]
  colnames(data) <- c("id","influence") # influence_norm_unsigned_forward_steady_state
  data$seed <- gsub(".*from_|\\.csv*|_influence.*","",basename(csv))
  influence.df <- rbind(influence.df,data)
}
influence.df$influence <- as.numeric(influence.df$influence)

# Read influence data
pkls <- list.files(file.path(banc.path,"data/cascade/frankenbrain_v1.5/signed/"),
                   pattern=".pkl",
                   full.names = TRUE)
cascade.df <- data.frame()
for(pkl in pkls){
  if(grepl("motor|unknown|mixed|chemosensory_tactile|proprioceptive_tactile|chemosensory_proprioceptive",pkl)){
    next
  }
  data <- py_load_object(pkl)
  data_r <- as.data.frame(data)
  colnames(data_r) <- c("id","distance")
  data_r$seed <- gsub("_to_.*|\\.pkl*","",basename(pkl))
  cascade.df <- rbind(cascade.df,data_r)
}
cascade.df$distance <- as.numeric(cascade.df$distance)
influence.df <- dplyr::left_join(influence.df,cascade.df,by=c("id","seed"))

# Format
influence.meta <- influence.df %>%
  dplyr::left_join(franken.meta %>%
                     dplyr::select(id, region, super_class, hemilineage, 
                                   cell_function, nerve, cell_class, sez_class, cell_sub_class, cell_type, top_nt, input_connections),
                   by = "id") %>%
  dplyr::filter(super_class!="glia", !is.na(super_class)) %>%
  dplyr::mutate(super_class = case_when(
    super_class == "brain_central_other" & !is.na(sez_class) ~ "sez",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(cell_type = ifelse(grepl("KCg-s",cell_type),"KCg-s",cell_type)) %>%
  dplyr::mutate(cell_sub_class = ifelse(is.na(cell_sub_class),gsub("_.*","",cell_type),cell_sub_class)) %>%
  dplyr::mutate(influence_log = 1-log(influence)) 

#############################
### Influence vs. cascade ### 
#############################

# First, let's summarize the data to get mean influence for each distance, super_class, and seed
summarized_data <- influence.meta %>%
  dplyr::filter(influence > 0.0001) %>%
  dplyr::filter(!is.na(distance), !is.na(influence_log)) %>%
  group_by(distance, seed) %>%
  summarise(mean_influence = median(influence_log, na.rm = TRUE),
            .groups = 'drop')

# Make violin plot
g.v <- ggplot(influence.meta %>% dplyr::filter(influence > 0.0001,
                                        !is.na(distance), 
                                        !is.na(influence_log)), 
              aes(x = factor(distance), y = influence_log, fill = as.character(distance))) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_point(data = summarized_data, 
             aes(y = mean_influence),
             position = position_dodge(width = 0.9),
             color = "black", size = 1) +
  facet_wrap(~ seed, scales = "free_y") +
  scale_fill_manual(values = paper.cols) +
  labs(title = "Influence by Distance and Super Class, Faceted by Seed",
       x = "Distance",
       y = "Influence (log10 scale)",
       fill = "Super Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal") +
  # scale_y_log10(
  #   breaks = scales::trans_breaks("log10", function(x) 10^x),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) +
  guides(fill = guide_legend(nrow = 1))

# Save
ggsave(plot = g.v, 
       filename = file.path(banc.path, "images", "influence", "influence_cascade_comparison.png"), 
       width = 12, height = 6, dpi = 300)

# Get unique super_classes
super_classes <- unique(influence.meta$super_class)

# Function to create density plot for a single super_class
create_density_plot <- function(data, super_class) {
  
  # Prepare data for the overall distribution
  overall_data <- data %>% 
    dplyr::filter(influence > 0.0001,
                  !is.na(distance), 
                  !is.na(influence)) %>%
    dplyr::mutate(super_class = "other", distance = -1) 
  
  # Combine with the specific super_class data
  plot_data <- bind_rows(
    data %>% 
      dplyr::filter(influence > 0.0001) %>%
      dplyr::filter(!is.na(distance), !is.na(influence), super_class == !!super_class),
    overall_data
  ) %>%
    dplyr::mutate(distance = as.character(distance))
  values <- c(rep("solid",length(unique(plot_data$distance))))
  names(values) <- c(unique(plot_data$distance))
  values[["-1"]] <- "dashed"
  
  # plot
  ggplot(plot_data, 
         aes(x = influence, 
             color = distance, 
             group = interaction(distance, super_class))) +
    geom_density(aes(linetype = super_class), linewidth = 1) +
    facet_wrap(~ seed, scales = "free_y") +
    scale_color_manual(values = paper.cols) +
    #scale_linetype_manual(values = values) +
    labs(title = paste("Influence Distribution by Distance for", super_class),
         subtitle = "Faceted by Seed (Dashed line: Overall distribution)",
         x = "Influence (log10 scale)",
         y = "Density",
         color = "Distance",
         linetype = "Distribution") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "horizontal") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1))
}

# Create and save a density plot for each super_class
for (sc in super_classes) {
  g.d <- create_density_plot(influence.meta, sc)
  
  # Save the plot
  ggsave(plot = g.d, 
         filename = file.path(banc.path, "images", "cascade", 
                              paste0("influence_distribution_by_distance_", 
                                     gsub(" ", "_", sc), ".png")),
         width = 16, height = 10, dpi = 300)
}

#############################
### Influence progression ### 
#############################

# Prepare the data
binned_data <- influence.meta %>%
  dplyr::filter(influence_log <= 15) %>%
  dplyr::mutate(influence_bin = cut(influence_log, 
                        breaks = seq(min(influence_log), max(influence_log) + 0.01, by = 0.01),
                        include.lowest = TRUE,
                        labels = seq(min(influence_log), max(influence_log), by = 0.01))) %>%
  dplyr::group_by(seed, influence_bin, super_class) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  complete(seed, influence_bin, super_class, fill = list(count = 0))

# Create the plot
g.progress <- ggplot(binned_data, aes(x = influence_bin, y = count, fill = super_class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ seed) +
  scale_fill_viridis_d() +
  labs(title = "distribution of super classes across influence values",
       x = "1-log(influence), binned",
       y = "no. neurons",
       fill = "super class") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text.x = element_blank()) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
  scale_fill_manual(values=paper.cols)

# Save
ggsave(plot = g.progress, 
       filename = file.path(banc.path,"images","influence","influence_super_class_distribution.png"), 
       width = 16, height = 8, dpi = 300)

# Create the plot
g.dens.brain <- ggplot(influence.meta %>%
                   dplyr::filter(region %in% c("midbrain","sez")), 
                 aes(x = influence_log, color = super_class)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ seed, scales = "free_y") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(title = "distribution of super classes across influence values",
       x = "1-log(influence), binned",
       y = "density",
       fill = "super class") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) +
  xlim(-2.5,12) +
  scale_color_manual(values=paper.cols)

# Save
ggsave(plot = g.dens.brain, 
       filename = file.path(banc.path,"images","influence","influence_midbrain_super_class_density.png"), 
       width = 16, height = 8, dpi = 300)

# Create the plot
g.dens.optic <- ggplot(influence.meta %>%
                         dplyr::filter(region %in% c("optic","optic_lobes")), 
                       aes(x = influence_log, color = super_class)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ seed, scales = "free_y") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(title = "distribution of super classes across influence values",
       x = "1-log(influence), binned",
       y = "density",
       fill = "super class") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) +
  xlim(-2.5,12) +
  scale_color_manual(values=paper.cols)

# Save
ggsave(plot = g.dens.optic, 
       filename = file.path(banc.path,"images","influence","influence_optic_super_class_density.png"), 
       width = 16, height = 8, dpi = 300)

# Create the plot
g.dens.vnc <- ggplot(influence.meta %>%
                         dplyr::filter(region %in% c("neck_connective","vnc")), 
                       aes(x = influence_log, color = super_class)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ seed, scales = "free_y") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(title = "distribution of super classes across influence values",
       x = "1-log(influence), binned",
       y = "density",
       fill = "super class") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) +
  xlim(-2.5,12) +
  scale_color_manual(values=paper.cols)

# Save
ggsave(plot = g.dens.vnc, 
       filename = file.path(banc.path,"images","influence","influence_vnc_super_class_density.png"), 
       width = 16, height = 8, dpi = 300)

#########################
### Influence heatmap ### 
#########################

# Turn into matrix
heatmap_matrix <- reshape2::acast(data = influence.meta %>% 
                                    dplyr::filter(influence_log<15), 
                                  formula = cell_type ~ seed, 
                                  value.var = "influence_log",
                                  fun.aggregate = function(x) mean(x,na.rm=TRUE))
heatmap_matrix <- scale(heatmap_matrix)

# Create annotation for cell types (rows)
cell_type_annotation <- influence.meta %>%
  dplyr::distinct(cell_type, super_class) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_type)) %>%
  column_to_rownames("cell_type")

# Ensure row names in annotation match matrix
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix), , drop = FALSE]

# Create annotation for seeds (columns)
seed_annotation <- franken.meta %>%
  dplyr::filter(!is.na(cell_function), !is.na(region)) %>%
  dplyr::distinct(cell_function, region) %>%
  dplyr::distinct(seed = cell_function, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_function)) %>%
  dplyr::mutate(cell_function = paste0(cell_function,"_influence")) %>%
  column_to_rownames("cell_function") 

# Ensure column names in annotation match matrix
seed_annotation <- seed_annotation[colnames(heatmap_matrix), , drop = FALSE]

# Create a color palette from blue to white to red
color_palette <- colorRampPalette(c("red","coral", "white", "skyblue","blue"))(100)

# Create breaks centered at 0
min_val <- min(heatmap_matrix, na.rm = T)
max_val <- max(heatmap_matrix, na.rm =T)
range <- max(abs(min_val), abs(max_val))
breaks <- seq(-range, range, length.out = 101)

# Help merge
.merge_hclust <- function(hclist) {
  #-- Merge
  d <- as.dendrogram(hclist[[1]])
  for (i in 2:length(hclist)) {
    d <- merge(d, as.dendrogram(hclist[[i]]))
  }
  as.hclust(d)
}

# Define the hclust_semisupervised and .merge_hclust functions
hclust_semisupervised <- function(data, groups, 
                                  dist_method = "euclidean",
                                  dist_p = 2, 
                                  hclust_method = "complete") {
  hclist <- lapply(groups, function (group) {
    hclust(dist(data[group,], 
                method = dist_method, 
                p = dist_p), 
           method = hclust_method)
  })
  hc <- .merge_hclust(hclist)
  data_reordered <- data[unlist(groups),]
  return(list(data = data_reordered, hclust = hc))
}

# Group cell types by super_class
groups <- split(rownames(cell_type_annotation), cell_type_annotation$super_class)

# Apply semi-supervised clustering
clustering_result <- hclust_semisupervised(data = heatmap_matrix, 
                                           groups = groups, 
                                           dist_method = "euclidean", 
                                           hclust_method = "ward.D2")

# Reorder the matrix and annotation based on the clustering result
heatmap_matrix_normalized <- clustering_result$data
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix_normalized), , drop = FALSE]

# First, create the annotation_colors list using your paper.cols vector
annotation_colors <- list(
  super_class = paper.cols[names(paper.cols) %in% unique(cell_type_annotation$super_class)],
  region = paper.cols[names(paper.cols) %in% unique(seed_annotation$region)],
  flow = paper.cols[names(paper.cols) %in% unique(seed_annotation$flow)],
  seed = paper.cols[names(paper.cols) %in% unique(seed_annotation$seed)]
)

# Create the heatmap
pheatmap(heatmap_matrix_normalized,
         color = color_palette,
         breaks = breaks,
         annotation_row = cell_type_annotation,
         annotation_col = seed_annotation,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         cluster_rows = clustering_result$hclust,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_row = 6,
         treeheight_row = 0,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         cutree_cols = length(unique(influence.meta$region)),
         filename = file.path(banc.path, "images", "influence", "influence_heatmap.png"),
         width = 15,
         height = 12)

#######################
### Cascade heatmap ### 
#######################

# Turn into matrix
heatmap_matrix <- reshape2::acast(data = influence.meta %>% 
                                    dplyr::filter(!is.na(distance)), 
                                  formula = cell_type ~ seed, 
                                  value.var = "distance",
                                  fun.aggregate = function(x) mean(x,na.rm=TRUE))

# Create annotation for cell types (rows)
cell_type_annotation <- influence.meta %>%
  dplyr::filter(!is.na(cell_type),!is.na(distance)) %>%
  dplyr::distinct(cell_type, super_class) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  column_to_rownames("cell_type")

# Ensure row names in annotation match matrix
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix), , drop = FALSE]

# Create annotation for seeds (columns)
seed_annotation <- franken.meta %>%
  dplyr::filter(!is.na(cell_function), !is.na(region)) %>%
  dplyr::distinct(cell_function, region) %>%
  dplyr::distinct(seed = cell_function, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_function)) %>%
  column_to_rownames("cell_function") 

# Ensure column names in annotation match matrix
seed_annotation <- seed_annotation[colnames(heatmap_matrix), , drop = FALSE]

# Create a color palette from blue to white to red
color_palette <- colorRampPalette(c("red","coral", "white", "skyblue","blue"))(100)

# Create breaks centered at 0
min_val <- min(heatmap_matrix, na.rm = T)
max_val <- max(heatmap_matrix, na.rm =T)
breaks <- seq(min_val, max_val, length.out = 101)

# Group cell types by super_class
groups <- split(rownames(cell_type_annotation), cell_type_annotation$super_class)

# Apply semi-supervised clustering
clustering_result <- hclust_semisupervised(data = heatmap_matrix, 
                                           groups = groups, 
                                           dist_method = "euclidean", 
                                           hclust_method = "ward.D2")

# Reorder the matrix and annotation based on the clustering result
heatmap_matrix_normalized <- clustering_result$data
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix_normalized), , drop = FALSE]

# First, create the annotation_colors list using your paper.cols vector
annotation_colors <- list(
  super_class = paper.cols[names(paper.cols) %in% unique(cell_type_annotation$super_class)],
  region = paper.cols[names(paper.cols) %in% unique(seed_annotation$region)],
  flow = paper.cols[names(paper.cols) %in% unique(seed_annotation$flow)],
  seed = paper.cols[names(paper.cols) %in% unique(seed_annotation$seed)]
)

# Create the heatmap
pheatmap(heatmap_matrix_normalized,
         color = color_palette,
         breaks = breaks,
         annotation_row = cell_type_annotation,
         annotation_col = seed_annotation,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         cluster_rows = clustering_result$hclust,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_row = 6,
         treeheight_row = 0,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         cutree_cols = length(unique(influence.meta$region)),
         filename = file.path(banc.path, "images", "cascade", "cascade_heatmap.png"),
         width = 15,
         height = 12)


##################
### KC heatmap ### 
##################

# Turn into matrix
heatmap_matrix <- reshape2::acast(data = influence.meta %>% 
                                    dplyr::filter(influence_log<12,
                                                  grepl("^KC",cell_type)), 
                                  formula = id ~ seed, 
                                  value.var = "influence",
                                  fun.aggregate = function(x) mean(x,na.rm=TRUE))
heatmap_matrix <- t(apply(heatmap_matrix, 1, function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE)))

# Create annotation for cell types (rows)
cell_type_annotation <- influence.meta %>%
  dplyr::filter(!is.na(cell_type),grepl("^KC",cell_type)) %>%
  dplyr::mutate(cell_type = ifelse(grepl("KCg-s",cell_type),"KCg-s",cell_type)) %>%
  dplyr::distinct(id, cell_type) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  column_to_rownames("id")

# Ensure row names in annotation match matrix
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix), , drop = FALSE]

# Create annotation for seeds (columns)
seed_annotation <- franken.meta %>%
  dplyr::filter(!is.na(cell_function), !is.na(region)) %>%
  dplyr::distinct(cell_function, region) %>%
  dplyr::distinct(seed = cell_function, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_function)) %>%
  dplyr::mutate(cell_function = paste0(cell_function,"_influence")) %>%
  column_to_rownames("cell_function") 

# Ensure column names in annotation match matrix
seed_annotation <- seed_annotation[colnames(heatmap_matrix), , drop = FALSE]

# Group cell types by cell_type
groups <- split(rownames(cell_type_annotation), cell_type_annotation$cell_type)

# Apply semi-supervised clustering
clustering_result <- hclust_semisupervised(data = heatmap_matrix, 
                                           groups = groups, 
                                           dist_method = "euclidean", 
                                           hclust_method = "ward.D2")

# Reorder the matrix and annotation based on the clustering result
heatmap_matrix_normalized <- clustering_result$data
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix_normalized), , drop = FALSE]

# First, create the annotation_colors list using your paper.cols vector
annotation_colors <- list(
  cell_type = paper.cols[names(paper.cols) %in% unique(cell_type_annotation$cell_type)],
  region = paper.cols[names(paper.cols) %in% unique(seed_annotation$region)],
  flow = paper.cols[names(paper.cols) %in% unique(seed_annotation$flow)],
  seed = paper.cols[names(paper.cols) %in% unique(seed_annotation$seed)]
)

# Logarithmic color scale
epsilon <- min(heatmap_matrix_normalized[heatmap_matrix_normalized > 0],na.rm=TRUE) / 10
heatmap_matrix_log <- log10(heatmap_matrix_normalized + epsilon)
min_val <- min(heatmap_matrix_log,na.rm=TRUE)
max_val <- max(heatmap_matrix_log,na.rm=TRUE)
breaks <- 10^seq(min_val, max_val, length.out = 101)
colors <- viridis::viridis(100)
colors <- colorRampPalette(rev(c("red","coral", "white", "skyblue","blue")))(100)

# Create the heatmap
pheatmap(heatmap_matrix_normalized,
         color = colors,
         breaks = breaks,
         annotation_row = cell_type_annotation,
         annotation_col = seed_annotation,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         cluster_rows = clustering_result$hclust,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_row = 6,
         treeheight_row = 0,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         cutree_cols = length(unique(influence.meta$region)),
         filename = file.path(banc.path, "images", "influence", "influence_kc_heatmap.png"),
         width = 15,
         height = 12)

##################
### CX heatmap ### 
##################

# Turn into matrix
heatmap_matrix <- reshape2::acast(data = influence.meta %>% 
                                    dplyr::filter(influence_log<12,
                                                  !is.na(distance),
                                                  grepl("^CX",cell_class)), 
                                  formula = id ~ seed, 
                                  value.var = "distance",
                                  fun.aggregate = function(x) mean(x,na.rm=TRUE))
#heatmap_matrix <- t(apply(heatmap_matrix, 1, function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE)))

# Create annotation for cell types (rows)
cell_type_annotation <- influence.meta %>%
  dplyr::filter(!is.na(cell_type),
                !grepl("^DM",cell_sub_class),
                grepl("^CX_input",cell_class)) %>%
  dplyr::distinct(id, cell_sub_class) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  column_to_rownames("id")

# Ensure row names in annotation match matrix
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix), , drop = FALSE]

# Create annotation for seeds (columns)
seed_annotation <- franken.meta %>%
  dplyr::filter(!is.na(cell_function), !is.na(region)) %>%
  dplyr::distinct(cell_function, region) %>%
  dplyr::distinct(seed = cell_function, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_function)) %>%
  dplyr::mutate(cell_function = paste0(cell_function,"_influence")) %>%
  column_to_rownames("cell_function") 

# Ensure column names in annotation match matrix
seed_annotation <- seed_annotation[colnames(heatmap_matrix), , drop = FALSE]

# Group cell types by cell_type
groups <- split(rownames(cell_type_annotation), cell_type_annotation$cell_sub_class)

# Apply semi-supervised clustering
clustering_result <- hclust_semisupervised(data = heatmap_matrix, 
                                           groups = groups, 
                                           dist_method = "euclidean", 
                                           hclust_method = "ward.D2")

# Reorder the matrix and annotation based on the clustering result
heatmap_matrix_normalized <- clustering_result$data
cell_type_annotation <- cell_type_annotation[rownames(heatmap_matrix_normalized), , drop = FALSE]

# First, create the annotation_colors list using your paper.cols vector
annotation_colors <- list(
  cell_sub_class = paper.cols[names(paper.cols) %in% unique(cell_type_annotation$cell_sub_class)],
  region = paper.cols[names(paper.cols) %in% unique(seed_annotation$region)],
  flow = paper.cols[names(paper.cols) %in% unique(seed_annotation$flow)],
  seed = paper.cols[names(paper.cols) %in% unique(seed_annotation$seed)]
)

# Logarithmic color scale
epsilon <- min(heatmap_matrix_normalized[heatmap_matrix_normalized > 0],na.rm=TRUE) / 10
heatmap_matrix_log <- log10(heatmap_matrix_normalized + epsilon)
min_val <- min(heatmap_matrix_log,na.rm=TRUE)
max_val <- max(heatmap_matrix_log,na.rm=TRUE)
breaks <- 10^seq(min_val, max_val, length.out = 101)
colors <- viridis::viridis(100)
colors <- colorRampPalette(c("red","coral", "white", "skyblue","blue"))(100)

# Create the heatmap
pheatmap(heatmap_matrix_normalized,
         color = colors,
         breaks = breaks,
         annotation_row = cell_type_annotation,
         annotation_col = seed_annotation,
         #annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         cluster_rows = clustering_result$hclust,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_row = 6,
         treeheight_row = 0,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         cutree_cols = length(unique(influence.meta$region)),
         filename = file.path(banc.path, "images", "influence", "influence_cx_heatmap.png"),
         width = 15,
         height = 12)

############
### UMAP ### 
############

# Turn into matrix
influence.m <- reshape2::acast(data = influence.df %>% dplyr::filter(!id %in% seed.neurons), 
                               formula = id ~ seed, 
                               value.var = "influence",
                               fun.aggregate = mean)

# Represent as UMAP
set.seed(42)  # for reproducibility
umap_result <- uwot::umap(influence.m, n_neighbors = 15, min_dist = 0.1, n_components = 2)

# Create a data frame with UMAP coordinates
umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  id = rownames(umap_result)
) %>% 
  dplyr::left_join(franken.meta %>%
                     dplyr::select(id, region,super_class, hemilineage, cell_function, nerve, cell_class, cell_sub_class, cell_type, input_connections),
                   by = "id") %>%
  dplyr::filter(input_connections>=120)
umap_df$cell_class = umap_df$super_class

# Consolidate cell class
classes <- sort(table(umap_df$cell_class), decreasing = TRUE)
chosen.classes <- names(classes[classes>50])
umap_df$cell_class[!umap_df$cell_class%in%chosen.classes] <- "other"

# Shave off outliers
umap_df <- umap_df %>%
  dplyr::filter(UMAP1 < quantile(UMAP1,.999),
                UMAP1 > quantile(UMAP1,.001)) %>%
  dplyr::filter(UMAP2 < quantile(UMAP2,.999),
                UMAP2 > quantile(UMAP2,.001))

# Create a color map using paper.cols
color_map <- paper.cols[names(paper.cols) %in% unique(umap_df$cell_class)]

# If there are any cell classes without a corresponding color in paper.cols,
# assign them a default color (e.g., grey)
missing_classes <- setdiff(unique(umap_df$cell_class), names(color_map))
if(length(missing_classes) > 0) {
  color_map[missing_classes] <- "#CCCCCC"  # light grey for missing colors
}

# Create the base layer with all points in grey
p <- plotly::plot_ly(umap_df, 
                     x = ~UMAP1, y = ~UMAP2, type = 'scatter', mode = 'markers',
                     marker = list(size = 5, opacity = 0.3, color = 'lightgrey'),
                     hoverinfo = 'none',
                     showlegend = FALSE,
                     name = 'all') %>%
  plotly::layout(title = "UMAP representation of neurons",
                 xaxis = list(title = "UMAP1"),
                 yaxis = list(title = "UMAP2"))

# Add colored layer for each cell type
for (entry in unique_cell_types) {
  subset_df <- umap_df[umap_df$cell_class == entry, ]
  p <- plotly::add_trace(p,
                         data = subset_df,
                         x = ~UMAP1, y = ~UMAP2, type = 'scatter', mode = 'markers',
                         marker = list(size = 5, opacity = 0.5, color = color_map[entry]),
                         hoverinfo = 'text',
                         name = entry,
                         showlegend = TRUE,
                         text = ~paste("cell_type:", cell_type, 
                                       "<br>cell_sub_class:", cell_sub_class,
                                       "<br>cell_class:", cell_class))
}

# Display the plot
p

# Save the plot as an HTML file
htmlwidgets::saveWidget(p, file.path(banc.path, "images/plotly/influence_umap_frankenbrain.html"))
webshot::webshot(file.path(banc.path, "images/plotly/influence_umap_frankenbrain.html"), 
        file.path(banc.path, "images/plotly/influence_umap_frankenbrain.png"))

######################
### UMAP cell type ### 
######################

# Turn into matrix
influence.df.ct <- influence.df %>%
  dplyr::left_join(franken.meta %>%
                     dplyr::select(id, region,super_class, hemilineage, cell_function, nerve, cell_class, cell_sub_class, cell_type, input_connections),
                   by = "id") %>%
  dplyr::filter(!id %in% seed.neurons, super_class != "sensory") %>%
  dplyr::group_by(cell_type, seed) %>%
  dplyr::mutate(infuence = sum(influence, na.rm = TRUE),
                input_connections = sum(input_connections, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cell_type, seed, .keep_all = TRUE) %>%
  dplyr::filter(input_connections>=120)

# Make matrix
influence.m <- reshape2::acast(data = influence.df.ct, 
                               formula = cell_type ~ seed, 
                               value.var = "influence")

# Represent as UMAP
set.seed(42)  # for reproducibility
umap_result <- uwot::umap(influence.m, 
                          n_neighbors = 15, 
                          min_dist = 0.1, 
                          n_components = 2)

# Create a data frame with UMAP coordinates
umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  cell_type = rownames(umap_result)
) %>% 
  dplyr::left_join(franken.meta %>%
                     dplyr::select(cell_type, region, super_class, hemilineage, cell_function, nerve, cell_class, cell_sub_class),
                   by = "cell_type") %>%
  dplyr::distinct(cell_type, .keep_all = TRUE)
umap_df$cell_class = umap_df$super_class

# Consolidate cell class
classes <- sort(table(umap_df$cell_class), decreasing = TRUE)
chosen.classes <- names(classes[classes>10])
umap_df$cell_class[!umap_df$cell_class%in%chosen.classes] <- "other"

# Shave off outliers
umap_df <- umap_df %>%
  dplyr::filter(UMAP1 < quantile(UMAP1,.999),
                UMAP1 > quantile(UMAP1,.001)) %>%
  dplyr::filter(UMAP2 < quantile(UMAP2,.999),
                UMAP2 > quantile(UMAP2,.001))

# Create a color map using paper.cols
color_map <- paper.cols[names(paper.cols) %in% unique(umap_df$cell_class)]

# If there are any cell classes without a corresponding color in paper.cols,
# assign them a default color (e.g., grey)
missing_classes <- setdiff(unique(umap_df$cell_class), names(color_map))
if(length(missing_classes) > 0) {
  color_map[missing_classes] <- "#CCCCCC"  # light grey for missing colors
}

# Create the base layer with all points in grey
p <- plotly::plot_ly(umap_df, 
                     x = ~UMAP1, y = ~UMAP2, type = 'scatter', mode = 'markers',
                     marker = list(size = 5, opacity = 0.3, color = 'lightgrey'),
                     hoverinfo = 'none',
                     showlegend = FALSE,
                     name = 'all') %>%
  plotly::layout(title = "UMAP representation of neurons",
                 xaxis = list(title = "UMAP1"),
                 yaxis = list(title = "UMAP2"))

# Add colored layer for each cell type
for (entry in unique_cell_types) {
  subset_df <- umap_df[umap_df$cell_class == entry, ]
  p <- plotly::add_trace(p,
                         data = subset_df,
                         x = ~UMAP1, y = ~UMAP2, type = 'scatter', mode = 'markers',
                         marker = list(size = 5, opacity = 0.5, color = color_map[entry]),
                         hoverinfo = 'text',
                         name = entry,
                         showlegend = TRUE,
                         text = ~paste("cell_type:", cell_type, 
                                       "<br>cell_sub_class:", cell_sub_class,
                                       "<br>cell_class:", cell_class))
}

# Display the plot
p

# Save the plot as an HTML file
htmlwidgets::saveWidget(p, file.path(banc.path, "images/plotly/influence_cell_type_umap_frankenbrain.html"))
webshot::webshot(file.path(banc.path, "images/plotly/influence_cell_type_umap_frankenbrain.html"), 
                 file.path(banc.path, "images/plotly/influence_cell_type_umap_frankenbrain.png"))


# ############
# ### PHATE ### 
# ############
# 
# # Load required libraries
# library(phateR)
# library(dplyr)
# library(plotly)
# library(RColorBrewer)
# library(caret) 
# 
# # Represent as PHATE
# set.seed(42)  # for reproducibility
# phate_result <- phate(influence.m, ndim = 2, k = 15)
# 
# # Create a data frame with PHATE coordinates
# phate_df <- data.frame(
#   PHATE1 = phate_result$embedding[,1],
#   PHATE2 = phate_result$embedding[,2],
#   id = rownames(influence.m)
# ) %>% 
#   dplyr::left_join(franken.meta %>%
#                      dplyr::select(id, region, hemilineage, cell_function, nerve, cell_class, cell_sub_class, cell_type),
#                    by = "id")
# 
# # Consolidate cell class
# classes <- sort(table(phate_df$cell_class), decreasing = TRUE)
# chosen.classes <- names(classes[classes>50])
# phate_df$cell_class[!phate_df$cell_class%in%chosen.classes] <- "other"
# 
# # Shave off outliers
# phate_df <- phate_df %>%
#   dplyr::filter(PHATE1 < quantile(PHATE1,.999),
#                 PHATE1 > quantile(PHATE1,.001)) %>%
#   dplyr::filter(PHATE2 < quantile(PHATE2,.999),
#                 PHATE2 > quantile(PHATE2,.001))
# 
# # Create a color palette that can handle more categories
# unique_cell_types <- unique(phate_df$cell_class)
# n_types <- length(unique_cell_types)
# if (n_types <= 12) {
#   # If 12 or fewer categories, use Set3 which has 12 colors
#   colors <- RColorBrewer::brewer.pal(n_types, "Set3")
# } else {
#   # If more than 12 categories, use a continuous color palette
#   colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_types)
# }
# 
# # Create a named vector for color mapping
# color_map <- setNames(colors, unique_cell_types)
# 
# # Create the base plot
# p <- plotly::plot_ly(phate_df, 
#                      x = ~PHATE1, y = ~PHATE2, type = 'scatter', mode = 'markers',
#                      marker = list(size = 5, opacity = 0.3, color = 'lightgrey'),
#                      hoverinfo = 'none',
#                      showlegend = FALSE,
#                      name = 'all') %>%
#   plotly::layout(title = "PHATE representation of neurons",
#                  xaxis = list(title = "PHATE1"),
#                  yaxis = list(title = "PHATE2"))
# 
# # Add colored layer for each cell type
# for (entry in unique_cell_types) {
#   subset_df <- phate_df[phate_df$cell_class == entry, ]
#   p <- plotly::add_trace(p,
#                          data = subset_df,
#                          x = ~PHATE1, y = ~PHATE2, type = 'scatter', mode = 'markers',
#                          marker = list(size = 5, opacity = 0.5, color = color_map[entry]),
#                          hoverinfo = 'text',
#                          name = entry,
#                          showlegend = TRUE,
#                          text = ~paste("cell_type:", cell_type, 
#                                        "<br>cell_sub_class:", cell_sub_class,
#                                        "<br>cell_class:", cell_class))
# }
# 
# # Calculate feature importance and project onto PHATE space
# scaled_data <- scale(influence.m)
# feature_importance <- apply(scaled_data, 2, function(x) {
#   cor(x, phate_result$embedding[,1], method = "spearman")^2 +
#     cor(x, phate_result$embedding[,2], method = "spearman")^2
# })
# 
# # Select top features (e.g., top 10)
# top_features <- names(sort(feature_importance, decreasing = TRUE)[1:10])
# 
# for (feature in top_features) {
#   feature_values <- scaled_data[, feature]
#   
#   # Calculate direction based on correlation
#   cor_x <- cor(feature_values, phate_result$embedding[,1])
#   cor_y <- cor(feature_values, phate_result$embedding[,2])
#   
#   # Normalize to get unit vector
#   magnitude <- sqrt(cor_x^2 + cor_y^2)
#   unit_x <- cor_x / magnitude
#   unit_y <- cor_y / magnitude
#   
#   # Scale the length of the arrow (adjust 0.1 as needed)
#   scale_factor <- 0.1 * max(abs(phate_result$embedding))
#   end_x <- unit_x * scale_factor
#   end_y <- unit_y * scale_factor
#   
#   p <- add_annotations(p,
#                        x = end_x, y = end_y,
#                        ax = end_x, ay = end_y,
#                        text = feature,
#                        showarrow = TRUE,
#                        arrowhead = 2,
#                        arrowsize = 1,
#                        arrowwidth = 2,
#                        arrowcolor = "red")
# }
# 
# # Display the plot
# p
# 
# # Save the plot as an HTML file
# htmlwidgets::saveWidget(p, "~/BANC-project/images/influence_phate_frankenbrain_cell_class.html")

###########
### PCA ###
###########

# # Perform PCA
# pca_result <- prcomp(influence.m, scale. = FALSE)
# 
# # Extract the first two principal components
# pca_df <- data.frame(
#   PC1 = pca_result$x[,1],
#   PC2 = pca_result$x[,2],
#   id = rownames(influence.m)
# ) %>% 
#   dplyr::left_join(franken.meta %>%
#                      dplyr::select(id, region,super_class, hemilineage, cell_function, nerve, cell_class, cell_sub_class, cell_type, input_connections),
#                    by = "id") %>%
#   dplyr::filter(input_connections>=120)
# 
# # Calculate the proportion of variance explained by each PC
# var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
# 
# # Create the scatter plot
# p <- plotly::plot_ly(pca_df, x = ~PC1, y = ~PC2, type = 'scatter', mode = 'markers',
#                      marker = list(size = 5, opacity = 0.5),
#                      hoverinfo = 'text',
#                      text = ~paste("Neuron ID:", id)) %>%
#   plotly::layout(
#     title = "PCA of neurons",
#     xaxis = list(title = paste0("PC1 (", round(var_explained[1]*100, 1), "%)")),
#     yaxis = list(title = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"))
#   )
# 
# # Add arrows for original dimensions
# for(i in 1:ncol(influence.m)) {
#   if(grepl("unknown",colnames(influence.m)[i])){
#     next
#   }
#   p <- plotly::add_trace(
#     p,
#     x = c(0, pca_result$rotation[i,1]),
#     y = c(0, pca_result$rotation[i,2]),
#     mode = 'lines+text',
#     line = list(color = 'red'),
#     text = colnames(influence.m)[i],
#     textposition = 'top center',
#     showlegend = FALSE
#   )
# }

# Display the plot
#p