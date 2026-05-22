### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

######################
### load libraries ###
######################

# Load packages
source("banc/banc-startup.R")
library(plotly)

######################
### influence plot ###
######################

# Meta
banc.connectivity.save.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity"
franken.meta <- arrow::read_feather(file.path(banc.connectivity.save.path,"frankenbrain_v.1.2_meta.feather")) %>%
  dplyr::mutate(region = dplyr::case_when(
    grepl("neck",region) ~ region,
    grepl("^LB|^MX|^MD",hemilineage) ~ "GNG",
    grepl("^aPhN|^MxLbN|^MD|^ON$",nerve) ~ "GNG",
    super_class == "visual_projection" ~ "optic",
    super_class == "visual_centrifugal" ~ "optic",
    super_class %in% c("ascending_neuron","descending_neuron") ~ "neck_connective",
    TRUE ~ region
  ))

# Read data
csvs <- list.files("~/BANC-project/data/influence/frankenbrain_v1.1/unsigned/seed/",
                   pattern="csv",
                   full.names = TRUE)
csvs <- csvs[grepl("chemosensory",csvs)]
csvs <- csvs[!grepl("CvC",csvs)]
influence.df <- data.frame()
for(csv in csvs){
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = readr::col_character()))
  data <- data[,3:5]
  colnames(data) <- c("id","distance","class")
  data$seed <- gsub(".*from_|\\.csv*","",basename(csv))
  influence.df <- rbind(influence.df,data)
}
influence.df$distance <- as.numeric(influence.df$distance)
influence.df$influence <- 1/influence.df$distance
seed.neurons <- subset(franken.meta, seed%in%influence.df$seed)$id
gng.neurons <- subset(franken.meta, !seed%in%influence.df$seed & region=="GNG" & !grepl("sensory",super_class))$id
meta.gng <- subset(franken.meta, id %in% gng.neurons)
meta.gng$cell_class[meta.gng$cell_class==""] = meta.gng$super_class[meta.gng$cell_class==""]

# Read data
csvs <- list.files("~/BANC-project/data/influence/frankenbrain_v1.1/unsigned/seed/",
                   pattern="csv",
                   full.names = TRUE)
csvs <- csvs[grepl("chemosensory",csvs)]
csvs <- csvs[!grepl("CvC",csvs)]
influence.df <- data.frame()
for(csv in csvs){
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = readr::col_character()))
  data <- data[,3:5]
  colnames(data) <- c("id","distance","class")
  data$seed <- gsub(".*from_|\\.csv*","",basename(csv))
  influence.df <- rbind(influence.df,data)
}
influence.df$distance <- as.numeric(influence.df$distance)
influence.df$influence <- 1/influence.df$distance
seed.neurons <- subset(franken.meta, seed%in%influence.df$seed)$id
gng.neurons <- subset(franken.meta, !seed%in%influence.df$seed & region=="GNG" & !grepl("sensory",super_class))$id
meta.gng <- subset(franken.meta, id %in% gng.neurons)
meta.gng$cell_class[meta.gng$cell_class==""] = meta.gng$super_class[meta.gng$cell_class==""]

# Read data
csvs <- list.files("~/BANC-project/data/influence/frankenbrain_v1.1/unsigned/seed/",
                   pattern="csv",
                   full.names = TRUE)
csvs <- csvs[grepl("chemosensory",csvs)]
csvs <- csvs[!grepl("CvC",csvs)]
influence.df <- data.frame()
for(csv in csvs){
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = readr::col_character()))
  data <- data[,3:5]
  colnames(data) <- c("id","distance","class")
  data$seed <- gsub(".*from_|\\.csv*","",basename(csv))
  influence.df <- rbind(influence.df,data)
}
influence.df$distance <- as.numeric(influence.df$distance)
influence.df$influence <- 1/influence.df$distance
seed.neurons <- subset(franken.meta, seed%in%influence.df$seed)$id
gng.neurons <- subset(franken.meta, !seed%in%influence.df$seed & region=="GNG" & !grepl("sensory",super_class))$id
meta.gng <- subset(franken.meta, id %in% gng.neurons)
meta.gng$cell_class[meta.gng$cell_class==""] = meta.gng$super_class[meta.gng$cell_class==""]

# Show strongly influenced neurons, broken down by nerve
influence.df.strong <- influence.df %>%
  dplyr::group_by(id) %>%
  dplyr::filter(any(influence>0.5)) %>%
  dplyr::ungroup()
influence.strong.m <- reshape2::acast(data = influence.df.strong %>% dplyr::filter(!id %in% seed.neurons), 
                               formula = id ~ seed, 
                               value.var = "influence")
influence.strong.gng <- influence.strong.m[rownames(influence.strong.m)%in%gng.neurons,]

# Turn into matrix
influence.m <- reshape2::acast(data = influence.df %>% dplyr::filter(!id %in% seed.neurons), 
                               formula = id ~ seed, 
                               value.var = "influence")
influence.m <- influence.m[rownames(influence.m)%in%gng.neurons,]
influence.m.df <- as.data.frame(influence.m)
influence.m.df$id <- rownames(influence.m)

# Load required libraries
library(phateR)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(caret) 

# Represent as PHATE
set.seed(42)  # for reproducibility
phate_result <- phate(influence.m, ndim = 2, k = 15)

# Create a data frame with PHATE coordinates
phate_df <- data.frame(
  PHATE1 = phate_result$embedding[,1],
  PHATE2 = phate_result$embedding[,2],
  id = rownames(influence.m)
) %>% 
  dplyr::left_join(franken.meta %>%
                     dplyr::select(id, region, hemilineage, cell_function, nerve, cell_class, cell_sub_class, cell_type),
                   by = "id") %>%
  dplyr::left_join(influence.m.df,
                   by = "id")

# Consolidate cell class
classes <- sort(table(phate_df$cell_class), decreasing = TRUE)
chosen.classes <- names(classes[classes>50])
phate_df$cell_class[!phate_df$cell_class%in%chosen.classes] <- "other"

# Shave off outliers
phate_df <- phate_df %>%
  dplyr::filter(PHATE1 < quantile(PHATE1,.999),
                PHATE1 > quantile(PHATE1,.001)) %>%
  dplyr::filter(PHATE2 < quantile(PHATE2,.999),
                PHATE2 > quantile(PHATE2,.001))

# Create a color palette that can handle more categories
unique_cell_types <- unique(phate_df$cell_class)
n_types <- length(unique_cell_types)
if (n_types <= 12) {
  # If 12 or fewer categories, use Set3 which has 12 colors
  colors <- RColorBrewer::brewer.pal(n_types, "Set3")
} else {
  # If more than 12 categories, use a continuous color palette
  colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_types)
}

# Create a named vector for color mapping
color_map <- setNames(colors, unique_cell_types)

# Create the base plot
p <- plotly::plot_ly(phate_df, 
                     x = ~PHATE1, y = ~PHATE2, type = 'scatter', mode = 'markers',
                     marker = list(size = 5, opacity = 0.3),
                     color = ~chemosensory__ProLNleft,
                     hoverinfo = 'none',
                     showlegend = FALSE,
                     name = 'all') %>%
  plotly::layout(title = "PHATE representation of neurons",
                 xaxis = list(title = "PHATE1"),
                 yaxis = list(title = "PHATE2"))

# Add colored layer for each cell type
for (entry in unique_cell_types) {
  subset_df <- phate_df[phate_df$cell_class == entry, ]
  p <- plotly::add_trace(p,
                         data = subset_df,
                         x = ~PHATE1, y = ~PHATE2, type = 'scatter', mode = 'markers',
                         marker = list(size = 5, opacity = 0.5, color = color_map[entry]),
                         hoverinfo = 'text',
                         name = entry,
                         showlegend = TRUE,
                         text = ~paste("cell_type:", cell_type, 
                                       "<br>cell_sub_class:", cell_sub_class,
                                       "<br>cell_class:", cell_class))
}

# Calculate feature importance and project onto PHATE space
scaled_data <- scale(influence.m)
feature_importance <- apply(scaled_data, 2, function(x) {
  cor(x, phate_result$embedding[,1], method = "spearman")^2 +
    cor(x, phate_result$embedding[,2], method = "spearman")^2
})

# Select top features (e.g., top 10)
top_features <- names(sort(feature_importance, decreasing = TRUE)[1:10])

for (feature in top_features) {
  feature_values <- scaled_data[, feature]
  
  # Calculate direction based on correlation
  cor_x <- cor(feature_values, phate_result$embedding[,1])
  cor_y <- cor(feature_values, phate_result$embedding[,2])
  
  # Normalize to get unit vector
  magnitude <- sqrt(cor_x^2 + cor_y^2)
  unit_x <- cor_x / magnitude
  unit_y <- cor_y / magnitude
  
  # Scale the length of the arrow (adjust 0.1 as needed)
  scale_factor <- 0.1 * max(abs(phate_result$embedding))
  end_x <- unit_x * scale_factor
  end_y <- unit_y * scale_factor
  
  p <- add_annotations(p,
                       x = end_x, y = end_y,
                       ax = end_x, ay = end_y,
                       text = feature,
                       showarrow = TRUE,
                       arrowhead = 2,
                       arrowsize = 1,
                       arrowwidth = 2,
                       arrowcolor = "red")
}

# Display the plot
p

# Show heatmap with GNG neurons
influence.gng <- influence.m[rownames(influence.m)%in%gng.neurons,]
influence.gng <- t(apply(influence.gng, 1, function(x) (x - min(x)) / (max(x) - min(x))))

# Plot heatmap
p.brain <- pheatmap(influence.gng[sort(rownames(influence.gng)),sort(colnames(influence.gng))], 
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    #color = colorRampPalette(c("grey95","#FC6882"))(100),
                    breaks = seq(0, max(influence.gng), length.out = 101),
                    #annotation_row = row_ann,
                    #annotation_col = col_ann,
                    #annotation_colors = list(cell_sub_class = class_colors),
                    main = paste("heatmap for brain GNG neurons"),
                    fontsize = 10             # Adjust font size if names are crowded
)

##################
### get labels ###
##################

# Get meta data
bc <- banctable_query() %>%
  dplyr::select(-starts_with("_")) %>%
  dplyr::arrange(cell_type)

# franken meta
meta <- arrow::read_feather(file.path(banc.connectivity.save.path,"frankenbrain_v.1.2_meta.feather"))

# fafb/manc meta
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types)) %>%
  dplyr::arrange(cell_type)
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types)) %>%
  dplyr::arrange(cell_type)

# Different gustatory partitions
meta.gust <- meta %>%
  dplyr::filter(flow %in% c("brain_afferent","vnc_afferent"),
                grepl("gustatory|chemosensory",cell_function)) %>%
  dplyr::mutate(cell_sub_class = case_when(
    grepl("AbN",nerve) ~ "AbN",
    dataset=="MANC" ~ gsub("_R.*|_L.*","",nerve),
    dataset=="cross-matched" ~ gsub("_R.*|_L.*","",nerve),
    cell_class!="gustatory" ~ cell_class,
    cell_sub_class=="" ~ "unknown",
    !is.na(cell_sub_class) ~ cell_sub_class,
    is.na(cell_sub_class) ~ "unknown",
    TRUE ~ "unknown"
  )) %>%
  dplyr::mutate(cell_sub_class = ifelse(cell_sub_class=="","unknown",cell_sub_class))
fw.meta.gust <- fw.meta %>%
  dplyr::filter(flow=="afferent",
                grepl("gustatory|chemosensory",cell_function)) %>%
  dplyr::mutate(cell_sub_class = ifelse(cell_class!="gustatory",cell_class,cell_sub_class),
                cell_sub_class = ifelse(is.na(cell_sub_class),"unknown",cell_sub_class)) %>%
  dplyr::mutate(cell_sub_class = ifelse(is.na(cell_sub_class),"unknown",cell_sub_class)) %>%
  dplyr::mutate(cell_type = ifelse(is.na(cell_type),"unknown",cell_type))
mc.meta.gust <- mc.meta %>%
  dplyr::filter(flow == "afferent",
                grepl("gustatory|chemosensory",cell_function)) %>%
  dplyr::mutate(nerve = ifelse(is.na(nerve),"unknown",nerve),
                cell_sub_class = gsub("_R.*|_L.*","",nerve)) %>%
  dplyr::mutate(cell_sub_class = ifelse(is.na(cell_sub_class),"unknown",cell_sub_class),
                cell_sub_class = ifelse(grepl("AbN",cell_sub_class),"AbN",cell_sub_class))
bc.gusts <- bc %>%
  dplyr::filter(fafb_match %in% fw.meta.gust$root_783|manc_match %in% mc.meta.gust$bodyid) %>%
  dplyr::select(-cell_sub_class) %>%
  dplyr::left_join(mc.meta.gust[,c("bodyid","cell_sub_class")], by = c("manc_match"="bodyid")) %>%
  dplyr::left_join(fw.meta.gust[,c("root_783","cell_sub_class")], by = c("fafb_match"="root_783")) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    is.na(cell_sub_class.x) ~ cell_sub_class.y,
    is.na(cell_sub_class.y) ~ cell_sub_class.x,
    !is.na(cell_sub_class.x) ~ cell_sub_class.x,
    !is.na(cell_sub_class.x) ~ cell_sub_class.y,
    TRUE ~ "unknown"
  )) %>%
  dplyr::mutate(cell_sub_class = ifelse(is.na(cell_sub_class),"unknown",cell_sub_class)) %>%
  dplyr::select(-cell_sub_class.x,-cell_sub_class.y)
bc.ids.brain <- bc.gusts$root_id[bc.gusts$region=="midbrain"]
bc.ids.vnc <- bc.gusts$root_id[bc.gusts$region=="vnc"]
bc.ids.neck <- bc.gusts$root_id[bc.gusts$region=="neck_connective"]

# Get brain sub classes
bc.ids.bitter <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="bitter"]
bc.ids.accessory_pharyngeal_nerve_sensory_group1 <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="accessory_pharyngeal_nerve_sensory_group1"]
bc.ids.accessory_pharyngeal_nerve_sensory_group2 <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="accessory_pharyngeal_nerve_sensory_group2"]
bc.ids.pharyngeal_nerve_sensory_group2 <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="pharyngeal_nerve_sensory_group2"]
bc.ids.enteric_gustatory <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="pharyngeal_nerve_sensory_group2"]
bc.ids.low_salt <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="low-salt"]
bc.ids.sugar_water <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="sugar/water"]
bc.ids.taste_peg <- bc.gusts$root_id[bc.gusts$region=="midbrain"&bc.gusts$cell_sub_class=="taste peg"]
bc.ids.brain.MesoLN <- bc.gusts$root_id[grepl("SA",bc.gusts$manc_cell_type)&bc.gusts$cell_sub_class=="MesoLN"]
bc.ids.brain.MetaLN <- bc.gusts$root_id[grepl("SA",bc.gusts$manc_cell_type)&bc.gusts$cell_sub_class=="MetaLN"]
bc.ids.brain.ProLN <- bc.gusts$root_id[grepl("SA",bc.gusts$manc_cell_type)&bc.gusts$cell_sub_class=="ProLN"]

# Get vnc sub classes
bc.ids.AbN <- bc.gusts$root_id[bc.gusts$cell_sub_class%in%c("AbN","AbN2","AbN3","AbN4","AbNT")]
bc.ids.ADMN <- bc.gusts$root_id[bc.gusts$cell_sub_class=="ADMN"]
bc.ids.MesoLN <- bc.gusts$root_id[bc.gusts$cell_sub_class=="MesoLN"]
bc.ids.MetaLN <- bc.gusts$root_id[bc.gusts$cell_sub_class=="MetaLN"]
bc.ids.ProLN <- bc.gusts$root_id[bc.gusts$cell_sub_class=="ProLN"]

# Gustatory focused ascending neurons

# Gustatory focused descending neurons

# Gustatory focused motor neurons

#################
### BANC plot ###
#################

# Get BANC meshes
mesh.obj <- file.path(banc.obj.save.path,paste0(unique(bc.gusts$root_id),".obj"))
mesh.obj <- mesh.obj[file.exists(mesh.obj)]
bc.meshes <- nat::as.neuronlist(pbapply::pblapply(mesh.obj,function(x) {
  n = readobj::read.obj(x, convert.rgl = TRUE)[[1]]
  n$id = gsub("\\.obj$","",basename(x))
  n
}))
names(bc.meshes) <- unlist(sapply(bc.meshes, function(x) x$id))  
bc.meshes <- nat::nlapply(bc.meshes,
                          Rvcg::vcgQEdecim,
                          percent = 0.1)
bc.meshes.brain <- bc.meshes[names(bc.meshes)%in%bc.ids.brain]
bc.meshes.neck <- bc.meshes[names(bc.meshes)%in%bc.ids.neck]
bc.meshes.vnc <- bc.meshes[names(bc.meshes)%in%bc.ids.vnc]

# Plot flange region
gng <- subset(banc_brain_neuropils.surf,"GNG")

# Plot
g.anat <- ggplot2::ggplot() +
  geom_neuron(x = banc_neuropil.surf,
              cols = c("grey50","grey75"),
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.1) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.bitter],
              cols = c("#8B0000", "#FF4500"),  # Dark red to orange-red
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.low_salt],
              cols = c("#4B0082", "#9400D3"),  # Indigo to violet
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.sugar_water],
              cols = c("#006400", "#32CD32"),  # Dark green to lime green
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.taste_peg],
              cols = c("#FF1493", "#FF69B4"),  # Deep pink to hot pink
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.enteric_gustatory],
              cols = c("#FF4500", "#FFA500"),  # Orange-red to orange
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.accessory_pharyngeal_nerve_sensory_group1],
              cols = c("#1E90FF", "#87CEEB"),  # Dodger blue to sky blue
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.accessory_pharyngeal_nerve_sensory_group2],
              cols = c("#DAA520", "#FFD700"),  # Goldenrod to gold
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.AbN],
              cols = c("#008B8B", "#20B2AA"),  # Dark cyan to light sea green
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.ADMN],
              cols = c("#800080", "#9932CC"),  # Purple to dark orchid
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.MesoLN],
              cols = c("#2E8B57", "#3CB371"),  # Sea green to medium sea green
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.MetaLN],
              cols = c("#A52A2A", "#DC143C"),  # Brown to crimson
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  geom_neuron(bc.meshes[names(bc.meshes)%in%bc.ids.ProLN],
              cols = c("#000080", "#4169E1"),  # Navy blue to royal blue
              rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
              alpha = 0.15) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::guides(fill = "none", color = "none") +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(hjust = 0, size = 8,
                                                    face = "bold",
                                                    colour = "black"),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(0, 0, 0, 0),
                 panel.spacing = ggplot2::unit(0, "cm"),
                 panel.border = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank()) +
  ggplot2::labs(title = '')

# Save whle BANC plot
ggsave(g.anat,
       filename = "inst/images/gustatory_network/banc_gustatory_neurons.png",
       height = 8, width = 8)

# Plot for just the brain
g.front <- ggplot2::ggplot() +
  geom_neuron(x = gng,
              cols = c("grey50","grey75"),
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.1) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.bitter]),
              cols = c("#8B0000", "#FF4500"),  # Dark red to orange-red
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.low_salt]),
              cols = c("#4B0082", "#9400D3"),  # Indigo to violet
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.sugar_water]),
              cols = c("#006400", "#32CD32"),  # Dark green to lime green
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.taste_peg]),
              cols = c("#FF1493", "#FF69B4"),  # Deep pink to hot pink
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.enteric_gustatory]),
              cols = c("#FF4500", "#FFA500"),  # Orange-red to orange
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.accessory_pharyngeal_nerve_sensory_group1]),
              cols = c("#1E90FF", "#87CEEB"),  # Dodger blue to sky blue
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.accessory_pharyngeal_nerve_sensory_group2]),
              cols = c("#DAA520", "#FFD700"),  # Goldenrod to gold
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.brain.MesoLN]),
              cols = c("#2E8B57", "#3CB371"),  # Sea green to medium sea green
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.brain.MetaLN]),
              cols = c("#A52A2A", "#DC143C"),  # Brown to crimson
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(invert = TRUE,bc.meshes[names(bc.meshes)%in%bc.ids.brain.ProLN]),
              cols = c("#000080", "#4169E1"),  # Navy blue to royal blue
              rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
              alpha = 0.15) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::guides(fill = "none", color = "none") +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(hjust = 0, size = 8,
                                                    face = "bold",
                                                    colour = "black"),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(0, 0, 0, 0),
                 panel.spacing = ggplot2::unit(0, "cm"),
                 panel.border = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank()) +
  ggplot2::labs(title = '')

# Save brain BANC plot
ggsave(g.front,
       filename = "inst/images/gustatory_network/banc_front_gustatory_neurons.png",
       height = 8, width = 8)

# Plot for VNC
g.vnc <- ggplot2::ggplot() +
  geom_neuron(x = banc_vnc_neuropil.surf,
              cols = c("grey50","grey75"),
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.1) +
  geom_neuron(banc_decapitate(bc.meshes[names(bc.meshes)%in%bc.ids.AbN]),
              cols = c("#008B8B", "#20B2AA"),  # Dark cyan to light sea green
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(bc.meshes[names(bc.meshes)%in%bc.ids.ADMN]),
              cols = c("#800080", "#9932CC"),  # Purple to dark orchid
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(bc.meshes[names(bc.meshes)%in%bc.ids.MesoLN]),
              cols = c("#2E8B57", "#3CB371"),  # Sea green to medium sea green
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(bc.meshes[names(bc.meshes)%in%bc.ids.MetaLN]),
              cols = c("#A52A2A", "#DC143C"),  # Brown to crimson
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.15) +
  geom_neuron(banc_decapitate(bc.meshes[names(bc.meshes)%in%bc.ids.ProLN]),
              cols = c("#000080", "#4169E1"),  # Navy blue to royal blue
              rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
              alpha = 0.15) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::guides(fill = "none", color = "none") +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(hjust = 0, size = 8,
                                                    face = "bold",
                                                    colour = "black"),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(0, 0, 0, 0),
                 panel.spacing = ggplot2::unit(0, "cm"),
                 panel.border = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank()) +
  ggplot2::labs(title = '')

# Save VNC BANC plot
ggsave(g.vnc,
       filename = "inst/images/gustatory_network/banc_vnc_gustatory_neurons.png",
       height = 8, width = 8)

#####################
### NBLAST groups ###
#####################


######################################
### Determine second order neurons ###
######################################
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)

# Cosine similarity matrix
df <- readr::read_csv("/home/ab714/BANC-project/data/connectivity_embedding/gustatory/gustatory_output_cosine_distances.csv", 
                      col_types = cols(
                        .default = col_number(),
                        pre_id = col_character())) %>%
  as.data.frame()
rownames(df) <- df$pre_id
distance_matrix_full <- as.matrix(df[,-1])

# Create a collapsed distance matrix
distance_matrix <- distance_matrix_full
meta.gust$anatomy_group <- paste0(meta.gust$cell_type,"_",meta.gust$cell_sub_class)
rownames(distance_matrix) <- meta.gust$anatomy_group[match(rownames(distance_matrix),meta.gust$id)]
colnames(distance_matrix) <- meta.gust$anatomy_group[match(colnames(distance_matrix),meta.gust$id)]
distance_matrix = apply(distance_matrix, 2, function(i) tapply(i, rownames(distance_matrix), median, na.rm = TRUE))
distance_matrix = t(apply(t(distance_matrix), 2, function(i) tapply(i, colnames(distance_matrix), median, na.rm = TRUE)))

# Prepare row and column annotations
row_annotation <- meta.gust %>% 
  dplyr::distinct(pre_id = anatomy_group, dataset, cell_sub_class) %>%
  dplyr::filter(pre_id %in% rownames(distance_matrix)) %>%
  dplyr::arrange(match(pre_id, rownames(distance_matrix)))

# Prepare column annotation similarly
col_annotation <- meta.gust %>% 
  dplyr::distinct(pre_id = id, dataset, cell_sub_class) %>%
  dplyr::filter(pre_id %in% colnames(distance_matrix)) %>%
  dplyr::arrange(match(pre_id, colnames(distance_matrix)))

# Color palettes
class_colors <- paper.cols[unique(row_annotation$cell_sub_class)]

# Filter distance matrix for the specific dataset
dataset_row_annotation <- row_annotation %>% 
  dplyr::filter(dataset %in% c("FAFB","cross-matched"))

# Filter distance matrix rows and columns
subset_matrix <- distance_matrix[dataset_row_annotation$pre_id, 
                                 dataset_row_annotation$pre_id]

# Create row annotation for this dataset
row_ann <- dataset_row_annotation %>%
  dplyr::select(cell_sub_class) %>%
  as.data.frame()
rownames(row_ann) <- dataset_row_annotation$pre_id

# Create column annotation for this dataset
col_ann <- dataset_row_annotation %>%
  select(cell_sub_class) %>%
  as.data.frame()
rownames(col_ann) <- dataset_row_annotation$pre_id

# Plot heatmap
p.brain <- pheatmap(subset_matrix[sort(rownames(subset_matrix)),sort(colnames(subset_matrix))], 
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    color = colorRampPalette(c("#FC6882","grey95"))(100),
                    breaks = seq(0, max(subset_matrix), length.out = 101),
                    annotation_row = row_ann,
                    annotation_col = col_ann,
                    annotation_colors = list(cell_sub_class = class_colors),
                    main = paste("heatmap for brain gustatory neurons"),
                    fontsize = 10             # Adjust font size if names are crowded
)

# Create a collapsed distance matrix
distance_matrix <- distance_matrix_full
meta.gust$anatomy_group <- paste0(meta.gust$cell_type,"_",meta.gust$cell_sub_class)
rownames(distance_matrix) <- meta.gust$anatomy_group[match(rownames(distance_matrix),meta.gust$id)]
colnames(distance_matrix) <- meta.gust$anatomy_group[match(colnames(distance_matrix),meta.gust$id)]
distance_matrix = apply(distance_matrix, 2, function(i) tapply(i, rownames(distance_matrix), median, na.rm = TRUE))
distance_matrix = t(apply(t(distance_matrix), 2, function(i) tapply(i, colnames(distance_matrix), median, na.rm = TRUE)))

# Prepare row and column annotations
row_annotation <- meta.gust %>% 
  dplyr::select(pre_id = anatomy_group, dataset, cell_sub_class) %>%
  dplyr::distinct(pre_id, cell_sub_class, .keep_all = TRUE) %>%
  dplyr::filter(pre_id %in% rownames(distance_matrix)) %>%
  dplyr::arrange(match(pre_id, rownames(distance_matrix)))

# Prepare column annotation similarly
col_annotation <- meta.gust %>% 
  dplyr::select(pre_id = anatomy_group, dataset, cell_sub_class) %>%
  dplyr::distinct(pre_id, cell_sub_class, .keep_all = TRUE) %>%
  dplyr::filter(pre_id %in% colnames(distance_matrix)) %>%
  dplyr::arrange(match(pre_id, colnames(distance_matrix)))

# Color palettes
class_colors <- scales::hue_pal()(length(unique(row_annotation$cell_sub_class)))
names(class_colors) <- unique(row_annotation$cell_sub_class)

# Filter distance matrix for the specific dataset
dataset_row_annotation <- row_annotation %>% 
  dplyr::filter(dataset %in% c("MANC","cross-matched"))

# Filter distance matrix rows and columns
subset_matrix <- distance_matrix[dataset_row_annotation$pre_id, 
                                 dataset_row_annotation$pre_id]

# Create row annotation for this dataset
row_ann <- dataset_row_annotation %>%
  dplyr::select(cell_sub_class) %>%
  as.data.frame()
rownames(row_ann) <- dataset_row_annotation$pre_id

# Create column annotation for this dataset
col_ann <- dataset_row_annotation %>%
  select(cell_sub_class) %>%
  as.data.frame()
rownames(col_ann) <- dataset_row_annotation$pre_id

# Plot heatmap
p.vnc <- pheatmap(subset_matrix[sort(rownames(subset_matrix)),sort(colnames(subset_matrix))], 
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  show_rownames = TRUE,
                  show_colnames = TRUE,
                  color = colorRampPalette(c("#FC6882","grey95"))(100),
                  breaks = seq(0, max(subset_matrix), length.out = 101),
                  annotation_row = row_ann,
                  annotation_col = col_ann,
                  annotation_colors = list(cell_sub_class = class_colors),
                  main = "heatmap vnc chemosensors",
                  fontsize = 10
)