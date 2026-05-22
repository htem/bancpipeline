### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Script to explore oviDNs cluster (cluster 15)
#   it plots the anatomy of cell types withing this cluster

######################
### load libraries ###
######################

library(writexl)
library(natverse)
library(tidyverse)
library(bancr)
library(nat.nblast)
library(tidyverse)
library(dendextend)

#########################
### define functions ####
#########################

# Plot neurons
plot_sel_skt <- function(ids2plot, skt, figure.dir, figure.name) {
  
  # get neurons skeletons
  ids2plot <- banc_latestid(ids2plot)
  neurons <- skt[ids2plot]
  
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
  
  # Plot neurons
  clear3d()
  open3d(windowRect = c(20, 30, 1000, 1000))
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(neurons, lwd = 0.5, soma = 2500)
  
  if (length(figure.name) > 0) {
    output_file <- str_c(figure.dir, figure.name, ".png")
    snapshot3d(output_file)
  }


}

###################
### get abdDNs ####
###################

# pull abdDN cluster of interest
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
                super_class == "descending") %>%
  dplyr::filter(dn_type == "15") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
  
# load skeletons
datadir <- str_c(getwd(), "/data/dn/")
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
names(l2) <- banc_latestid(names(l2))

# get abdDNs ids
abdDNs.ids <- meta$root_id

# define directory where to save data
if (!dir.exists(datadir)) {
  dir.create(datadir, recursive = TRUE)
}

# save relevant variables for clustering and plotting
saveRDS(meta, file = str_c(datadir, "meta_abd_dn.rds"))

############################
### plot abdDNs anatomy ####
############################

# set figure directory name
figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/")

# plot neurons to confirm morphology
plot_sel_skt(meta$root_id, l2, figure.dir, c())

# plot each cell type
cell_type <- unique(meta$cell_type)
for (i in cell_type) {
  plot_sel_skt(meta$root_id[grepl(i, meta$cell_type)], l2, figure.dir, i)
}

# check neurons in BANC
bancsee(banc_ids = meta$root_id[grepl("DNpe047", meta$cell_type)], 
        open = TRUE)
bancsee(banc_ids = meta$root_id[grepl("DNpe034", meta$cell_type)], 
        open = TRUE)
bancsee(banc_ids = meta$root_id[grepl("DNpe046", meta$cell_type)], 
        open = TRUE)
bancsee(banc_ids = meta$root_id[grepl("oviDNb", meta$cell_type)], 
        open = TRUE)
bancsee(banc_ids = meta$root_id[grepl("oviDNa_a", meta$cell_type)], 
        open = TRUE)

# check celltype of interest across datasets
#cell_type2use <- "DNpe047"
cell_type2use <- "oviDNa_b"

# Find neuron in manc
fafb.id2use <- fafb.meta$root_783[grepl(cell_type2use, fafb.meta$cell_type)]
# Find neuron in fafb
manc.id2use <- manc.meta$bodyid[grepl(cell_type2use, manc.meta$cell_type)]
# Find neuron in banc
banc.id2use <- meta$root_id[grepl(cell_type2use, meta$cell_type)]

# plot all neurons
bancsee(banc_ids = banc.id2use,
        fafb_ids = fafb.id2use, open = TRUE)

bancsee(banc_ids = meta$root_id[grepl("DNpe046", meta$cell_type)],
        fafb_ids = fafb.id2use, open = TRUE)