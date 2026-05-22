### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Wrangling Mo's connectivity data with dn taxonomy

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)
library(R.matlab)

# Make sure all functions query BANC and not FAFB
choose_banc()

###############################
### get data from BANCTABLE ###
###############################

bc <- banctable_query()
meta <- bc %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
    is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
    TRUE ~ cell_type
  ),
  top_nt = gsub("auto\\:","",top_nt),
  cell_sub_class = gsub("auto\\:","",cell_sub_class),
  cell_class = gsub("auto\\:","",cell_class),
  super_class = gsub("auto\\:","",super_class)) %>%
  dplyr::filter(region %in% c("neck_connective"),
                !is.na(side), 
                super_class == "descending")

###############################
### load connectivity data ####
###############################

os_info <- Sys.info()
if (os_info['sysname'] == "Darwin") {
  dn_dir <- "/Users/papers/BANC-project/data/connectivity_embedding/descending/"
} else if (os_info['sysname'] == "Windows") {
  dn_dir <- "C:/Users/papers/BANC-project/data/connectivity_embedding/descending/"
} else {
  print(paste("You are using:", os_info['sysname']))
}
rm(os_info)

# pull metadata associated to DNs, in particular some features of the inputs/outputs
dn_metadata <- read_csv(str_c(dn_dir, "descending_annotated_metadata.csv"),
                        col_types = cols(.default = col_guess(), `id...1` = col_character(), `id...2` = col_character()))
# For each DN a pairwise cosine distance between them is generated using:
#  1) all inputs (collapsed by celltype)
dn_input_cos_dist <- read_csv(str_c(dn_dir, "descending_input_cosine_distances.csv"),
                                      col_types = cols(.default = col_guess(), `post_id` = col_character()))
dn_input_umap <- read_csv(str_c(dn_dir, "descending_input_umap_embeddings.csv"),
                          col_types = cols(.default = col_guess(), `post_id` = col_character()))
#  2) all inputs and outputs (collapsed by celltype)
dn_joint_cos_dist <- read_csv(str_c(dn_dir, "descending_joint_cosine_distances.csv"),
                            col_types = cols(.default = col_guess(), `post_id` = col_character()))
dn_joint_umap <- read_csv(str_c(dn_dir, "descending_joint_umap_embeddings.csv"),
                          col_types = cols(.default = col_guess(), `post_id` = col_character()))
#  3) all outputs (collapsed by celltype)
dn_output_cos_dist <- read_csv(str_c(dn_dir, "descending_output_cosine_distances.csv"),
                                       col_types = cols(.default = col_guess(), `pre_id` = col_character()))
dn_output_umap <- read_csv(str_c(dn_dir, "descending_output_umap_embeddings.csv"),
                           col_types = cols(.default = col_guess(), `pre_id` = col_character()))

# The profiles of the input/output synapses are:
dn_input_synapse_profiles <- read_csv(str_c(dn_dir, "descending_input_synapse_count_profiles.csv"),
                                      col_types = cols(.default = col_guess(), `post_id` = col_character()))
dn_output_synapse_profiles <- read_csv(str_c(dn_dir, "descending_output_synapse_count_profiles.csv"),
                                       col_types = cols(.default = col_guess(), `pre_id` = col_character()))
rm(dn_dir)

#########################
### define functions ####
#########################

# define function to update ids per matrix
update_rows_column_names <- function(latest_id, synapse_type, input_dist) {
  
  # update ids (it assumes the matrix is symmetric and that it uses the same order for all matrices)
  column_names <- c(synapse_type, latest_id)
  
  # pass updated ids to distance matrices
  input_dist[[synapse_type]] <- latest_id
  colnames(input_dist) <- column_names
  
  # reduce matrix to unique values
  input_dist <- input_dist[match(unique(latest_id), latest_id), ]
  input_dist <- input_dist[, c(1, match(unique(latest_id), colnames(input_dist)))]
  
  return(list(input_dist = input_dist))
  
}
update_rows_names <- function(latest_id, synapse_type, input_dist) {
  
  # pass updated ids to distance matrices
  input_dist[[synapse_type]] <- latest_id

  # reduce matrix to unique values
  input_dist <- input_dist[match(unique(latest_id), latest_id), ]

  return(list(input_dist = input_dist))
  
}

# define function to select particular ids
get_subset_row_columns <- function(selected_id, synapse_type, input_dist) {
  
  # update ids (it assumes the matrix is symmetric and that it uses the same order for all matrices)
  latest_id <- input_dist[[synapse_type]]

  # reduce matrix to values matching selected_id
  input_dist <- input_dist[match(selected_id, latest_id), ]
  input_dist <- input_dist[, c(1, match(selected_id, colnames(input_dist)))]
  
  return(list(input_dist = input_dist))
  
}
get_subset_row <- function(selected_id, synapse_type, input_dist) {
  
  # update ids (it assumes the matrix is symmetric and that it uses the same order for all matrices)
  latest_id <- input_dist[[synapse_type]]
  
  # reduce matrix to values matching selected_id
  input_dist <- input_dist[match(selected_id, latest_id), ]

  return(list(input_dist = input_dist))
  
}

# define function to go from table to data format to plot heatmaps
get_linear_matrix_with_dn_type <- function(sorted_matrix, synapse_type, dn_subset) {
  
  # Convert matrix to long format while preserving full numeric precision
  # Remove pre_id column
  numeric_matrix <- as.matrix(sorted_matrix[, -c(1)])
  rownames(numeric_matrix) <- sorted_matrix[[synapse_type]]
  
  matrix_long <- as.data.table(numeric_matrix, keep.rownames = "Neuron1")
  matrix_long <- melt(matrix_long, id.vars = "Neuron1",
                      variable.name = "Neuron2",
                      value.name = "CosineDistance")
  
  # Add dn_type information
  matrix_long <- matrix_long %>%
    mutate(Class1 = dn_subset$dn_type[match(Neuron1, dn_subset$root_id)],
           Class2 = dn_subset$dn_type[match(Neuron2, dn_subset$root_id)])
  
  # Create a sorting order that groups by dn_type and maintains neuron order within each type
  neuron_order <- matrix_long %>%
    mutate(dn_type = dn_subset$dn_type[match(Neuron1, dn_subset$root_id)]) %>%
    arrange(dn_type, Neuron1) %>%
    pull(Neuron1) %>%
    unique()
  
  # Convert Neuron1 and Neuron2 to factors with the custom order
  matrix_long <- matrix_long %>%
    mutate(
      Neuron1 = factor(Neuron1, levels = neuron_order),
      Neuron2 = factor(Neuron2, levels = neuron_order)
    )  
  
  return(list(matrix_long = matrix_long))
  
}

# function to save R variable as MatLab variable
save_matlab_vars <- function(distance_vector, umap_dims, datadir, matfilename) {

  writeMat(con = str_c(datadir, matfilename, ".mat"), 
         cos_dist = distance_vector$CosineDistance, 
         dn_types_rows = as.numeric(distance_vector$Class1),
         dn_types_columns = as.numeric(distance_vector$Class2),
         umap_dim_one = as.numeric(umap_dims[['UMAP dim 1']]),
         umap_dim_two = as.numeric(umap_dims[['UMAP dim 2']]),
         dn_type = as.numeric(umap_dims$dn_type))
  
}

###################################################
### extract, sort and format data for plotting ####
###################################################

# sort dns by dn_type, remove na, update root_id, and keep unique
dn_subset <- meta %>%
  arrange(as.numeric(dn_type)) %>%
  filter(!is.na(dn_type))
dn_subset$root_id <- banc_latestid(dn_subset$root_id)
dn_subset <- dn_subset %>%
  distinct(root_id, .keep_all = TRUE)

# sort neuron ids by dn_type
sorted_dns <- dn_subset %>%
  arrange(as.numeric(dn_type)) %>%
  pull(root_id)
# get latest ids
sorted_dns <- banc_latestid(sorted_dns)
# use only unique values
sorted_dns <- sorted_dns[match(unique(sorted_dns), sorted_dns)]

# get latest ids from distance matrices and update matrix column and rows
#   given that tables/matrices have the same order of post_id, you only 
#   need to do this once and then replace them all
latest_id <- banc_latestid(dn_input_cos_dist$post_id)

dn_input_cos_dist <- update_rows_column_names(latest_id, "post_id", dn_input_cos_dist)$input_dist
dn_joint_cos_dist <- update_rows_column_names(latest_id, "post_id", dn_joint_cos_dist)$input_dist
dn_output_cos_dist <- update_rows_column_names(latest_id, "pre_id", dn_output_cos_dist)$input_dist

dn_input_umap <- update_rows_names(latest_id, "post_id",
                                   dn_input_umap)$input_dist
dn_joint_umap <- update_rows_names(latest_id, "post_id",
                                   dn_joint_umap)$input_dist
dn_output_umap <- update_rows_names(latest_id, "pre_id",
                                    dn_output_umap)$input_dist
dn_input_synapse_profiles <- update_rows_names(latest_id, "post_id",
                                               dn_input_synapse_profiles)$input_dist
dn_output_synapse_profiles <- update_rows_names(latest_id, "pre_id",
                                               dn_output_synapse_profiles)$input_dist

# get intersection of neurons
sel_sorted_dns <- intersect(sorted_dns, colnames(dn_input_cos_dist))

dn_input_cos_dist <- get_subset_row_columns(sel_sorted_dns, "post_id", 
                                               dn_input_cos_dist)$input_dist
dn_joint_cos_dist <- get_subset_row_columns(sel_sorted_dns, "post_id", 
                                               dn_joint_cos_dist)$input_dist
dn_output_cos_dist <- get_subset_row_columns(sel_sorted_dns, "pre_id", 
                                                dn_output_cos_dist)$input_dist

dn_input_umap <- get_subset_row(sel_sorted_dns, "post_id",
                                   dn_input_umap)$input_dist
dn_joint_umap <- get_subset_row(sel_sorted_dns, "post_id",
                                dn_joint_umap)$input_dist
dn_output_umap <- get_subset_row(sel_sorted_dns, "pre_id",
                                    dn_output_umap)$input_dist

dn_input_synapse_profiles <- get_subset_row(sel_sorted_dns, "post_id",
                                            dn_input_synapse_profiles)$input_dist
dn_output_synapse_profiles <- get_subset_row(sel_sorted_dns, "pre_id",
                                             dn_output_synapse_profiles)$input_dist

# update length of dn_subset data
dn_subset <- dn_subset %>%
  filter(root_id %in% sel_sorted_dns) %>%
  distinct(root_id, .keep_all = TRUE)

# update dn_metadata data
# Note: dn_metadata$id...1 and dn_metadata$id...2 are identical
dn_metadata$id...1 <- latest_id
dn_metadata$id...2 <- latest_id
dn_metadata <- get_subset_row(sel_sorted_dns, "id...1",
                              dn_metadata)$input_dist
dn_metadata$dn_type <- as.numeric(dn_subset$dn_type)

# collapse over dn_type
dn_metadata_meanpertype <- dn_metadata %>%
  group_by(dn_type) %>%
  summarize(
    id...1 = list(id...1),  # Keep `root_id` as a list
    across(where(is.numeric) & !c(id...1), mean, na.rm = TRUE),  # Calculate mean for numeric columns except `root_id`
    .groups = "drop"  # Ungroup the result
  )

# generate plot compatible variables (it also adds DN type to matrix)
sorted_dn_input_cos_dist_mat <- get_linear_matrix_with_dn_type(dn_input_cos_dist, 
                                                            "post_id", dn_subset)$matrix_long
sorted_dn_joint_cos_dist_mat <- get_linear_matrix_with_dn_type(dn_joint_cos_dist, 
                                                            "post_id", dn_subset)$matrix_long
sorted_dn_output_cos_dist_mat <- get_linear_matrix_with_dn_type(dn_output_cos_dist, 
                                                            "pre_id", dn_subset)$matrix_long

# add dn type to UMAP
dn_input_umap$dn_type <- dn_subset$dn_type
dn_output_umap$dn_type <- dn_subset$dn_type
dn_joint_umap$dn_type <- dn_subset$dn_type

###################################################
### save data in R and MatLab compatible files ####
###################################################

# Define directory where to save data
#   It assumes that the current directory is "/bancpipeline/"
datadir <- str_c(getwd(), "/data/dn/dn_connectivity")
if (!dir.exists(datadir)) {
  dir.create(datadir, recursive = TRUE)
}

# save a table of cosine distance matrices
write.csv(dn_input_cos_dist, file = str_c(datadir, "dn_input_cos_dist.csv"), 
          row.names = FALSE)
write.csv(dn_joint_cos_dist, file = str_c(datadir, "dn_joint_cos_dist.csv"), 
          row.names = FALSE)
write.csv(dn_output_cos_dist, file = str_c(datadir, "dn_output_cos_dist.csv"), 
          row.names = FALSE)

# save a table of metadata
write.csv(dn_metadata_meanpertype[, !names(dn_metadata_meanpertype) %in% "id...1"], 
          file = str_c(datadir, "dn_metadata_meanpertype.csv"), 
          row.names = FALSE)

# Save relevant variables in rds format
saveRDS(sorted_dn_input_cos_dist_mat, file = str_c(datadir, "dn_input_cos_dist.rds"))
saveRDS(sorted_dn_joint_cos_dist_mat, file = str_c(datadir, "dn_joint_cos_dist.rds"))
saveRDS(sorted_dn_output_cos_dist_mat, file = str_c(datadir, "dn_output_cos_dist.rds"))

saveRDS(dn_input_umap, file = str_c(datadir, "dn_input_umap.rds"))
saveRDS(dn_joint_umap, file = str_c(datadir, "dn_joint_umap.rds"))
saveRDS(dn_output_umap, file = str_c(datadir, "dn_output_umap.rds"))

saveRDS(dn_metadata_meanpertype, file = str_c(datadir, "dn_metadata_meanpertype.rds"))

# save matrix in matlab format
save_matlab_vars(sorted_dn_input_cos_dist_mat, dn_input_umap, datadir, "dn_input_cos_dist")
save_matlab_vars(sorted_dn_joint_cos_dist_mat, dn_joint_umap, datadir, "dn_joint_cos_dist")
save_matlab_vars(sorted_dn_output_cos_dist_mat, dn_output_umap, datadir, "dn_output_cos_dist")