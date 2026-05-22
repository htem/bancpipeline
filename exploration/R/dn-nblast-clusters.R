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

#################
### load data ###
#################

# define data directory
datadir <- str_c(getwd(), "/data/dn/")

# load relevant variables
meta <- readRDS(str_c(datadir, "meta_dn.rds"))
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
nb.all.full <- readRDS(str_c(datadir, "nblast_score_all_dn.rds"))

#########################
### define functions ####
#########################

plot_dendrogram_v1 <- function(hc, meta, skt, cut.height, figure.dir) {
  
  # Get groups/clusters
  groups <- cutree(hc, h = cut.height)
  
  # generate a dendrogram
  hc.col2 <- as.dendrogram(hc)
  
  # Reorder clusters to have the leftmost as cluster 1
  # Get the order of the labels in the dendrogram
  order_of_labels <- order.dendrogram(hc.col2)
  
  # Get a vector of group labels ordered by the dendrogram
  ordered_group_labels <- groups[order_of_labels]
  
  # Find unique groups in their left-to-right order
  unique_groups <- unique(ordered_group_labels)
  
  # define colors to use for cluster to be plotted
  other_color <- "#000000"
  target_color <- "#FF0000"
  
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
  
  for(cluster in unique(groups)){
    
    # get ordered cluster number (from left to right)
    cluster.ordered <- unique_groups[cluster]
    
    # Get members 
    members <- names(groups[groups==cluster])
    members <- names(groups[groups == cluster.ordered])
    neurons <- skt[members]
    
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
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "_fulldend.png")
    png(filename = output_file, width = width, height = height, res = dpi)
    
    plot(temp_hc, ylab = "Height", leaflab = "none")
    abline(h=cut.height)
    dev.off()
    
    # Plot new sub-dendrogram
    nb.all.sub <- nb.all.full[members,members]
    nb.all.sub <- nb.all.full[members, members]
    hc.sub <- nhclust(scoremat = nb.all.sub, method = "ward.D2")
    hc.sub$height <- sqrt(hc.sub$height)
    cut.height.sub <- 1 
    hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
    plot(hc.sub.col, labels = T)
    abline(h=cut.height.sub)
    
    temp_labels <- paste0(meta[match(labels(hc.sub.col), meta$root_id), "cell_type"], 
                          "_", meta[match(labels(hc.sub.col), meta$root_id), "side"])
    hc.sub.col <- dendextend::set_labels(hc.sub.col, temp_labels)
    
    dpi <- 300
    width <- width_in_inches * dpi
    height <- height_in_inches * dpi
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "_subdend.png")
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
    plot3d(hc.sub, db=neurons, lwd=0.5, soma = 2500, h=cut.height.sub)
    
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "_skt.png")
    snapshot3d(output_file)
    close3d()
    
  }
}

####################################
### cluster using NBLAST scores  ###
####################################

# cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# flatten
hc$height <- sqrt(hc$height)

# try coarse clustering
cut.height <- 3
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_dist_3/")
plot_dendrogram_v1(hc, meta, l2, cut.height, figure.dir)

cut.height <- 2.4
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_dist_2-4/")
plot_dendrogram_v1(hc, meta, l2, cut.height, figure.dir)

cut.height <- 1.9
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_dist_1-9/")
plot_dendrogram_v1(hc, meta, l2, cut.height, figure.dir)

cut.height <- 1.5
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_dist_1-5/")
plot_dendrogram_v1(hc, meta, l2, cut.height, figure.dir)