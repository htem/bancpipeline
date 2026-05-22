### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Plot connectivity data with dn taxonomy

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(R.matlab)
library(stringr)

###############################
### load connectivity data ####
###############################

# Define directory where to save data
#   It assumes that the current directory is "/bancpipeline/"
datadir <- str_c(getwd(), "/data/dn_connectivity/")

sorted_dn_input_cos_dist_mat <- readRDS(str_c(datadir, "dn_input_cos_dist.rds"))
sorted_dn_joint_cos_dist_mat <- readRDS(str_c(datadir, "dn_joint_cos_dist.rds"))
sorted_dn_output_cos_dist_mat <- readRDS(str_c(datadir, "dn_output_cos_dist.rds"))

dn_input_umap <- readRDS(str_c(datadir, "dn_input_umap.rds"))
dn_joint_umap <- readRDS(str_c(datadir, "dn_joint_umap.rds"))
dn_output_umap <- readRDS(str_c(datadir, "dn_output_umap.rds"))

dn_metadata_meanpertype <- readRDS(str_c(datadir, "dn_metadata_meanpertype.rds"))

plot_dist_heatmap <- function(input_matrix, figtitle, figname, figure.dir) {
  
  # set figure settings
  width_in_inches <- 10
  height_in_inches <- 10
  dpi <- 50
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, figname, ".png")
  png(filename = output_file, width = width, height = height, res = dpi)

  # Plot with the neuron labels
  p <- ggplot(input_matrix, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
    geom_tile() +
    scale_fill_gradient(
      low = "blue",
      high = "white",
      limits = c(0.8, 1)) +
    labs(title = figtitle,
         x = "DN Type",
         y = "DN Type",
         fill = "Cosine Distance") +
    theme(axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
    # Replace x and y-axis labels with neuron dn_type
    #scale_x_discrete(labels = input_matrix$Class1) +
    #scale_y_discrete(labels = input_matrix$Class1)
  
  # Print the plot object to the graphics device
  print(p)
  
  dev.off()

}

plot_umap_xy <- function(input_matrix, figtitle, figname, figure.dir){

  # set figure settings
  width_in_inches <- 10
  height_in_inches <- 10
  dpi <- 50
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, figname, ".png")
  png(filename = output_file, width = width, height = height, res = dpi)
    
  p <- ggplot(input_matrix, aes(x = `UMAP dim 1`, y = `UMAP dim 2`, color = dn_type)) +
    geom_point(size = 2) +  # Adjust point size as needed
    labs(title = figtitle,
         x = "UMAP Dimension 1",
         y = "UMAP Dimension 2",
         color = "Neuron Type") +
    theme_minimal()

  # Print the plot object to the graphics device
  print(p)
  
  dev.off()
    
}

#############################################################################
### Plot heatmap of dn cosine distances based on input, output, or joint ####
#############################################################################

# heat map of input cosine distances
figure.dir <- str_c(getwd(), "/figures/dn_connectivity_dist/")
figname <- "input_cos_dist_heatmap_linear_r"
figtitle <- "DN Input Distance Matrix Sorted by Class"

plot_dist_heatmap(sorted_dn_input_cos_dist_mat, figtitle, figname, figure.dir)

# heat map of output cosine distances
figname <- "output_cos_dist_heatmap_linear_r"
figtitle <- "DN Output Distance Matrix Sorted by Class"

plot_dist_heatmap(sorted_dn_output_cos_dist_mat, figtitle, figname, figure.dir)

# heat map of joint cosine distances
figname <- "joint_cos_dist_heatmap_linear_r"
figtitle <- "DN Joint Distance Matrix Sorted by Class"

plot_dist_heatmap(sorted_dn_joint_cos_dist_mat, figtitle, figname, figure.dir)

#####################################################################################
### Plot UMAP projection of dn cosine distances based on input, output, or joint ####
#####################################################################################

# first and second dimension of input cosine distances UMAP
figure.dir <- str_c(getwd(), "/figures/dn_connectivity_dist/")
figname <- "input_cos_dist_umap_r"
figtitle <- "DN Input Distance by DN morphological type"

plot_umap_xy(dn_input_umap, figtitle, figname, figure.dir)

# first and second dimension of output cosine distances UMAP
figure.dir <- str_c(getwd(), "/figures/dn_connectivity_dist/")
figname <- "output_cos_dist_umap_r"
figtitle <- "DN Output Distance by DN morphological type"

plot_umap_xy(dn_output_umap, figtitle, figname, figure.dir)

# first and second dimension of joint cosine distances UMAP
figure.dir <- str_c(getwd(), "/figures/dn_connectivity_dist/")
figname <- "input_cos_dist_umap_r"
figtitle <- "DN Joint Distance by DN morphological type"

plot_umap_xy(dn_joint_umap, figtitle, figname, figure.dir)

#####################################################################################
### Plot UMAP projection of dn cosine distances based on input, output, or joint ####
#####################################################################################

# heatmaps - input

# Reshape the data into long format for columns with "input" in their names
heatmap_data <- summary_df %>%
  pivot_longer(
    cols = contains("input_super") & !contains("dataset"),
    #cols = !contains("dataset"),# Select columns with "input" in their names
    names_to = "Metric",       # New column for the column names
    values_to = "Value"        # New column for the values
  )

# Create the heatmap
ggplot(heatmap_data, aes(x = Metric, y = dn_type, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white") +  # Customize color scale
  labs(
    title = "Heatmap of Input Metrics by Neuron Type",
    x = "Metric",
    y = "Neuron Type",
    fill = "Value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# heatmaps - output

# Reshape the data into long format for columns with "input" in their names
heatmap_data <- summary_df %>%
  pivot_longer(
    cols = contains("output_super") & !contains("dataset"),
    #cols = !contains("dataset"),# Select columns with "input" in their names
    names_to = "Metric",       # New column for the column names
    values_to = "Value"        # New column for the values
  )

# Create the heatmap
ggplot(heatmap_data, aes(x = Metric, y = dn_type, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white") +  # Customize color scale
  labs(
    title = "Heatmap of Output Metrics by Neuron Type",
    x = "Metric",
    y = "Neuron Type",
    fill = "Value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))