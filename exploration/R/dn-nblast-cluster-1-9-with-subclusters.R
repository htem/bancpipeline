### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Cluster DNs, use various distances and explore the pattern of splitting

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)

#########################
### define functions ####
#########################

plot_dendrogram_v1 <- function(hc, meta, skt, clus.cut.height, subclus.cut.height, figure.dir) {
  
  # Get groups/clusters
  groups <- cutree(hc, h = clus.cut.height)
  
  # generate a dendrogram
  hc.col2 <- as.dendrogram(hc)
  
  # Reorder clusters to have the leftmost as cluster 1
  # Get the order of the labels in the dendrogram
  order_of_labels <- order.dendrogram(hc.col2)
  
  # Get a vector of group labels ordered by the dendrogram
  ordered_group_labels <- groups[order_of_labels]
  
  # Find unique groups in their left-to-right order
  unique_groups <- unique(ordered_group_labels)
  
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
  
  for(cluster in unique(groups)){
    
    print(paste("running cluster #", cluster))
    
    # get ordered cluster number (from left to right)
    cluster.ordered <- unique_groups[cluster]
    
    # Get members 
    members <- names(groups[groups == cluster.ordered])
    neurons <- skt[members]
    
    # Plot new sub-dendrogram
    nb.all.sub <- nb.all.full[members, members]
    hc.sub <- nhclust(scoremat = nb.all.sub, method = "ward.D2")
    hc.sub$height <- sqrt(hc.sub$height)
    hc.sub.col <- dendroextras::color_clusters(hc.sub, h = subclus.cut.height)
    plot(hc.sub.col, labels = T)
    abline(h = subclus.cut.height)
    
    temp_labels <- paste0(meta[match(labels(hc.sub.col), meta$root_id), "cell_type"], 
                          "_", meta[match(labels(hc.sub.col), meta$root_id), "side"])
    hc.sub.col <- dendextend::set_labels(hc.sub.col, temp_labels)
    
    width_in_inches <- 20
    height_in_inches <- 10
    dpi <- 300
    width <- width_in_inches * dpi
    height <- height_in_inches * dpi
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "-0_subdend.png")
    png(filename = output_file, width = width, height = height, res = dpi)
    
    par(mar = c(3, 1, 1, 10))
    plot(hc.sub.col, labels = T, horiz = TRUE)
    par(mar = c(5, 4, 4, 2) + 0.1)
    dev.off()
    
    # Plot neurons
    clear3d()
    open3d(windowRect = c(20, 30, 1000, 1000))
    banc_view()
    plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
    plot3d(hc.sub, db = neurons, lwd = 0.5, soma = 2500, h = subclus.cut.height)
    
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "-0_skt.png")
    snapshot3d(output_file)
    close3d()
    
    # Get groups/clusters
    subgroups <- cutree(hc.sub, h = subclus.cut.height)
    
    # generate a dendrogram
    hc.col3 <- as.dendrogram(hc.sub)
    
    # Reorder clusters to have the leftmost as cluster 1
    # Get the order of the labels in the dendrogram
    sub_order_of_labels <- order.dendrogram(hc.col3)
    
    # Get a vector of group labels ordered by the dendrogram
    sub_ordered_group_labels <- subgroups[sub_order_of_labels]
    
    # Find unique groups in their left-to-right order
    unique_sub_groups <- unique(sub_ordered_group_labels)
    
    for(subcluster in unique(subgroups)){
      
      print(paste("running subcluster #", subcluster))
      
      # get ordered cluster number (from left to right)
      subcluster.ordered <- unique_sub_groups[subcluster]
      
      # Get members 
      submembers <- names(subgroups[subgroups == subcluster.ordered])
      subneurons <- skt[submembers]
      
      # Plot dendrograms
      figure.suffix <- str_c(cluster, "-", subcluster)
      plot_selected_cluster_dend(hc.col3, meta, subgroups, subcluster.ordered,
                                 figure.dir, figure.suffix)
      # Plot neurons
      plot_selected_cluster_neurons(subneurons, figure.dir, figure.suffix)     
      
    }
    
  }
}

plot_selected_cluster_neurons <- function(subneurons, figure.dir, figure.suffix) {
  
  # Plot neurons per selected cluster
  clear3d()
  open3d(windowRect = c(20, 30, 1000, 1000))
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(subneurons, lwd = 0.5, soma = 2500)
  
  output_file <- str_c(figure.dir, 'DN_nblast_clus_', figure.suffix,  "_skt.png")
  snapshot3d(output_file)
  close3d()
}

plot_selected_cluster_dend <- function(hc, meta, groups, cluster.ordered,
                                       figure.dir, figure.suffix) {
  
  # define colors to use for cluster to be plotted
  other_color <- "#000000"
  target_color <- "#FF0000"
  
  # plot selected cluster in full dendrogram
  temp_hc <- as.dendrogram(hc)
  temp_col <- ifelse(groups == cluster.ordered, target_color, other_color)
  temp_col <- temp_col[order.dendrogram(temp_hc)]
  temp_col <- factor(temp_col, unique(temp_col))
  temp_labels <- paste0(meta[match(labels(temp_hc), meta$root_id), "cell_type"], "_", 
                        meta[match(labels(temp_hc), meta$root_id), "side"])
  temp_hc <- temp_hc %>% 
    color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
    dendextend::set("labels_colors", as.character(temp_col))
  temp_hc <- dendextend::set_labels(temp_hc,temp_labels)
  
  width_in_inches <- 20
  height_in_inches <- 10
  dpi <- 100
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, 'DN_nblast_clus_', figure.suffix, "_subdend.png")
  png(filename = output_file, width = width, height = height, res = dpi)
  
  par(mar = c(3, 1, 1, 10))
  plot(temp_hc, labels = T, horiz = TRUE)
  dev.off()
}

#################
### load data ###
#################

# Define data directory
datadir <- str_c(getwd(), "/data/dn/")

# Load relevant variables
meta <- readRDS(str_c(datadir, "meta_dn.rds"))
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
nb.all.full <- readRDS(str_c(datadir, "nblast_score_all_dn.rds"))

####################################
### Cluster using NBLAST scores  ###
####################################

# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height <- sqrt(hc$height)

# cluster
clus.cut.height <- 1.9
subclus.cut.height <- 1.2
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_with_subclus_dist_1-9_1-2/")
plot_dendrogram_v1(hc, meta, l2, clus.cut.height, subclus.cut.height, figure.dir)

clus.cut.height <- 1.9
subclus.cut.height <- 1.4
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_with_subclus_dist_1-9_1-4/")
plot_dendrogram_v1(hc, meta, l2, clus.cut.height, subclus.cut.height, figure.dir)