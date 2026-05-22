### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

######################
### load libraries ###
######################

# Load packages
library(bancr)
library(tidyverse)
library(ggraph)
library(igraph)
library(dplyr)
library(tidygraph)
library(plotly)
library(arrow)
library(pheatmap)
library(tidyr)
library(furrr)
library(future)
library(parallel)
library(progressr)

# Data sources on files1
fafb.meta.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_meta.feather"
fafb.edgelist.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_edgelist.feather"
hemibrain.meta.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/hemibrain_v.1.2.1_meta.feather"
hemibrain.edgelist.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/hemibrain_v.1.2.1_edgelist.feather"
manc.meta.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/manc_1.2.1_meta.feather"
manc.edgelist.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/manc_1.2.1_edgelist.feather"
banc.meta.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/banc_meta.feather"

# Data sources on data1
fafb.meta.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_meta.feather"
fafb.edgelist.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_edgelist.feather"
hemibrain.meta.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/hemibrain_v.1.2.1_meta.feather"
hemibrain.edgelist.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/hemibrain_v.1.2.1_edgelist.feather"
manc.meta.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/manc_1.2.1_meta.feather"
manc.edgelist.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/manc_1.2.1_edgelist.feather"
banc.meta.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity/banc_meta.feather"

# output
manc.images <- "inst/images/manc/"

################
### get data ###
################

# Get the meta data, cell types, etc.
### NOTE: if this is slow, it could be because you are reading from remote, copy and store locally for faster performance ###
mc.meta <- arrow::read_feather(manc.meta.path) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::arrange(dplyr::desc(top_nt_p)) %>%
  dplyr::mutate(top_nt=top_nt[1]) %>%
  dplyr::ungroup()

# Read the edgelist
mc.elist <- arrow::read_feather(manc.edgelist.path)

################################
### choose  analysis problem ###
################################

# Let's make a plot of direct sensory to motor connections
mc.sensory <- mc.meta %>%
  dplyr::filter(grepl("CB.FB",cell_type),
                cell_class %in% c("afferent","sensory"))
mc.motor <- mc.meta %>%
  dplyr::filter(grepl("motor",cell_type),
                cell_class  %in% c("efferent","motor"))

# First, collect the IDs we're interested in
sensory_ids <- mc.sensory$bodyid
motor_ids <- mc.motor$bodyid

###################################
### filter edgelist for analysis ###
###################################

### NOTE: if this is slow, it could be because you are reading from remote, copy and store locally for faster performance ###

# Filter the edgelist for sensory -> motor connections
mc.elist.th <- mc.elist %>%
  dplyr::filter(count >= 5,
                pre %in% sensory_ids,
                post %in% motor_ids)

############
### plot ###
############

# First, let's merge the edgelist with meta data
mc.elist.merged <- mc.elist.th %>%
  left_join(mc.meta %>% select(bodyid, cell_type, top_nt), by = c("pre" = "bodyid")) %>%
  rename(pre_cell_type = cell_type, pre_top_nt = top_nt) %>%
  left_join(mc.meta %>% select(bodyid, cell_type, top_nt), by = c("post" = "bodyid")) %>%
  rename(post_cell_type = cell_type, post_top_nt = top_nt)

# Reshape and collapse the data by cell_type
mc.elist.collapsed <- mc.elist.merged %>%
  mutate(connection_type = paste(pre_label, "to", post_label)) %>%
  group_by(pre_cell_type, post_cell_type, connection_type) %>%
  summarise(count = sum(count), 
            pre_top_nt = first(pre_top_nt),
            post_top_nt = first(post_top_nt),
            .groups = 'drop')

# Heatmap plotting function
create_heatmap <- function(data, conn_type, filename) {
  mat_data <- data %>%
    filter(connection_type == conn_type) %>%
    select(pre_cell_type, post_cell_type, count)
  
  # Create matrix and collapse by row and column names
  mat <- mat_data %>%
    pivot_wider(names_from = post_cell_type, values_from = count, values_fill = 0) %>%
    column_to_rownames("pre_cell_type") %>%
    as.matrix()
  
  # Collapse matrix by taking mean of duplicate rows and columns
  mat <- aggregate(mat, by = list(rownames(mat)), FUN = sum)
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mat <- t(aggregate(t(mat), by = list(colnames(mat)), FUN = sum))
  colnames(mat) <- mat[1,]
  mat <- mat[-1,]
  
  # Convert matrix to numeric
  mat.names <-  dimnames(mat)
  mat <- apply(mat, 2, as.numeric)
  dimnames(mat) <- mat.names
  
  # Format numbers for display
  mat_display <- matrix(sprintf("%.2f", mat), nrow = nrow(mat))
  rownames(mat_display) <- rownames(mat)
  colnames(mat_display) <- colnames(mat)
  
  # Create annotation data frames
  row_ann <- data %>%
    filter(connection_type == conn_type) %>%
    select(pre_cell_type, pre_top_nt) %>%
    distinct() %>%
    group_by(pre_cell_type) %>%
    summarise(pre_top_nt = paste(unique(pre_top_nt), collapse = "/")) %>%
    column_to_rownames("pre_cell_type")
  
  col_ann <- data %>%
    filter(connection_type == conn_type) %>%
    select(post_cell_type, post_top_nt) %>%
    distinct() %>%
    group_by(post_cell_type) %>%
    summarise(post_top_nt = paste(unique(post_top_nt), collapse = "/")) %>%
    column_to_rownames("post_cell_type")
  
  # Ensure annotations match matrix dimensions
  row_ann <- row_ann[rownames(mat), , drop = FALSE]
  col_ann <- col_ann[colnames(mat), , drop = FALSE]
  
  # Create color palette for top_nt
  nt_colors <- c(
    acetylcholine = "#EF7C12",
    glutamate = "#8FDA04",
    gaba = "#1BB6AF",
    serotonin = "#FBBB48",
    dopamine = "#F4E3C7",
    octopamine = "#C70E7B",
    unknown = "grey"
  )
  ann_colors <- list(pre_top_nt = nt_colors, post_top_nt = nt_colors)
  
  # Create the heatmap
  heatmap_plot <- pheatmap(mat,
                           main = conn_type,
                           annotation_row = row_ann,
                           annotation_col = col_ann,
                           annotation_colors = ann_colors,
                           display_numbers = mat_display,
                           number_color = "black",
                           clustering_method = "ward.D2",
                           fontsize_number = 8,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = TRUE,  # This ensures row names are displayed
                           fontsize_row = 8,      # Adjust this value to change row name font size
                           angle_col = 45)        # Angle column labels for better readability
  
  # Save the plot
  ggsave(filename, heatmap_plot, width = 12, height = 10, units = "in")
}

# Create heatmaps for each connection type
for (conn_type in unique(mc.elist.collapsed$connection_type)) {
  filename <- file.path(fafb.images, paste0("heatmap_collapsed_", gsub(" ", "_", conn_type), ".png"))
  create_heatmap(mc.elist.collapsed, conn_type, filename)
}