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
library(nat.nblast)
library(dendextend)

# Make sure all functions query BANC and not FAFB
# choose_banc()

################
### get data ###
################

# Get meta data
# make sure you run: banctable_set_token(user = "", pwd = "") to generate banc token
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
  dplyr::filter(region%in%c("neck_connective"),
                !is.na(side), 
                super_class=="descending")
# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)

# Get all ids
ids <- unique(meta$root_id)
ids.left <- unique(subset(meta, side=="left")$root_id)
l2dps <- banc_read_l2dp(ids)
# l2dps[ids.left] <- banc_mirror(l2dps[ids.left])

# Get skeletons for visualisation
# Get skeletons for visualization (units in nanometers)
l2 <- banc_read_l2skel(ids)
# l2[ids.left] <- banc_mirror(l2[ids.left])

# mirror neurons to right hemisphere (mirror require units in nanometers)
l2[ids.left] <- banc_mirror(l2[ids.left])

# Make dotprops for NBLASTing (inaccurate but fast to acquire) (units in microns)
l2dps <- dotprops(l2/1000)

# Re-root to soma where this is known
# banc.roots <- banc_roots()
banc.roots <- bancr:::banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots, estimate = TRUE)
rm(banc.roots)

# Save the state of this R session, to avoid doing the step above next time you load this file
# save.image()

########################
### visualizing data ###
########################

# plot all neurons (mirrored to right side)
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(l2, lwd=0.5, soma = 5000)

figure.dir <- str_c(here::here(), "/figures/")
snapshot3d(str_c(figure.dir, "allDNs_skt.png"))
close3d()

# plot example neuron (DNg34)
DNg34.ids <- subset(meta, cell_type == "DNg34")
DNg34.ids <- DNg34.ids$root_id

open3d(windowRect = c(20, 30, 1000, 1000))
wire3d(banc_neuropil.surf/1e3, alpha = 0.1, col='lightgrey')
plot3d(l2dps[DNg34.ids[1]], lwd=3, soma = 5000, col = 'red')
plot3d(l2dps[DNg34.ids[2]], lwd=100, soma = 5000, col = 'blue')

# Check NBLAST distance for example neuron
# DNg34 <- l2dps[DNg34.ids]
# nb.all.sub <- nat.nblast::nblast_allbyall(DNg34, normalisation = "mean")

###############
### NBLAST  ###
###############

# Run NBLAST
nb.all.full <- nat.nblast::nblast_allbyall(l2dps, normalisation = "mean")

# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height <- sqrt(hc$height)

# Color clusters

# Split data into clusters and color clusters (each with a different color)
cut.height <- 1.5
hc.col <- dendroextras::color_clusters(hc, h = cut.height)
# Generate DN names
DN_labels <- paste0(meta[match(labels(hc.col), meta$root_id), "cell_type"], "_", 
                    meta[match(labels(hc.col), meta$root_id), "side"])
hc.col <- dendextend::set_labels(hc.col, DN_labels)

# Plot
plot(hc.col, labels = T)
abline(h = cut.height)
# Set custom labels
hc.col <- dendextend::set_labels(hc.col, DN_labels)

# See neurons
clear3d()
banc_view()
plot3d(hc, db = l2, h = cut.height, soma = 5000)

##############################
### visualizing clustering ###
##############################

# figure settings
width_in_inches <- 20
height_in_inches <- 10
dpi <- 100
width <- width_in_inches * dpi
height <- height_in_inches * dpi

# 1) plot and save dendrogram without labels
figure.dir <- str_c(here::here(), "/figures/")
png(filename = str_c(figure.dir, "AllDNs_clus_fulldend.png"), width = width, height = height, res = dpi)
plot(hc.col, ylab = "Height", leaflab = "none")
# to add labels run:
# plot(hc.col, , labels = T)
abline(h = cut.height)
dev.off()

# 2) custom plot clusters (select one in red and the rest in black)
hc.col2 <- as.dendrogram(hc)

# Get groups/clusters
groups <- cutree(hc, h = cut.height)

other_color <- "#000000"
target_color <- "#FF0000"
target_group <- 1
temp_col <- ifelse(groups == target_group, target_color, other_color)
temp_col <- temp_col[order.dendrogram(hc.col2)]
temp_col <- factor(temp_col, unique(temp_col))

hc.col2 <- hc.col2 %>% 
  color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
  dendextend::set("labels_colors", as.character(temp_col))

# Get groups
# Set custom labels
DN_labels <- paste0(meta[match(labels(hc.col2), meta$root_id), "cell_type"], "_", 
                    meta[match(labels(hc.col2), meta$root_id), "side"])
hc.col2 <- dendextend::set_labels(hc.col2, DN_labels)

# plot and save dendrogram without labels
figure.dir <- str_c(here::here(), "/figures/")
png(filename = str_c(figure.dir, "AllDNs_clus_", as.character(target_group), "_fulldend.png"), 
    width = width, height = height, res = dpi)
plot(hc.col2, ylab = "Height", leaflab = "none")
# to add labels run:
# plot(hc.col2, , labels = T)
abline(h = cut.height)
dev.off()

##############################################################################
### Split into clusters and plot each cluster (ordered from left to right) ###
##############################################################################

# Get groups/clusters
groups <- cutree(hc, h = cut.height)
table(groups)

###########################
### Manual Cell typing  ###
###########################
# Reorder clusters to have the leftmost as cluster 1
# Get the order of the labels in the dendrogram
order_of_labels <- order.dendrogram(hc.col2)

# Get a vector of group labels ordered by the dendrogram
ordered_group_labels <- groups[order_of_labels]

# Find unique groups in their left-to-right order
unique_groups <- unique(ordered_group_labels)

choices <- list()
other_color <- "#000000"
target_color <- "#FF0000"
figure.dir <- str_c(here::here(), "/figures/DN_nblast/")

if (!dir.exists(figure.dir)) {
  dir.create(figure.dir, recursive = TRUE)
}

#choices <- list()
for(cluster in unique(groups)){
  
  # get ordered cluster number (from left to right)
  cluster.ordered <- unique_groups[cluster]
  
  # Get members 
  members <- names(groups[groups==cluster])
  members <- names(groups[groups == cluster.ordered])
  neurons <- l2[members]
  
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
  
  # Get choices
  # db <- nlscan(neurons, col = "black", lwd = 3, soma = 5000)
  # choices[[cluster]] <- db
  
  # Pause for user
  #message(paste(choices,collapse=", "))
  #readline("prompt: press any key to continue to the next cluster ")
  #message(paste(choices,collapse=", "))
  #readline("prompt: press any key to continue to the next cluster ")
  
}

# Record manually entries here