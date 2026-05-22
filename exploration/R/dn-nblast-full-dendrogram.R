### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Plot full NBLAST-based dendrogram and neurons colorcoded by cluster (using a arbitrary initial distance of 1.5)

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

# Define data directory
datadir <- str_c(getwd(), "/data/dn/")

# Save relevant variables
meta <- readRDS(str_c(datadir, "meta_dn.rds"))
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
l2dps <-readRDS(str_c(datadir, "l2dps_dn.rds"))
nb.all.full <- readRDS(str_c(datadir, "nblast_score_all_dn.rds"))

####################################
### Cluster using NBLAST scores  ###
####################################

# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height <- sqrt(hc$height)

# Split data into clusters and color clusters (each with a different color)
cut.height <- 1.9
hc.col <- dendroextras::color_clusters(hc, h = cut.height)

# Generate DN names
DN_labels <- paste0(meta[match(labels(hc.col), meta$root_id), "cell_type"], "_", 
                    meta[match(labels(hc.col), meta$root_id), "side"])
hc.col <- dendextend::set_labels(hc.col, DN_labels)

##############################
### visualizing clustering ###
##############################

# Define figure directory
figure.dir <- str_c(getwd(), "/figures/")

# figure settings
width_in_inches <- 20
height_in_inches <- 10
dpi <- 100
width <- width_in_inches * dpi
height <- height_in_inches * dpi

# 1) plot and save dendrogram without labels
png(filename = str_c(figure.dir, "AllDNs_clus_fulldend.png"), width = width, height = height, res = dpi)
plot(hc.col, ylab = "Height", leaflab = "none")
# to add labels run:
# plot(hc.col, , labels = T)
abline(h = cut.height)
dev.off()

# 2) See neurons
# Plot neurons
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(hc, db = l2, lwd = 0.5, soma = 2500, h = cut.height)
snapshot3d(str_c(figure.dir, "AllDNs_skt.png"))
close3d()