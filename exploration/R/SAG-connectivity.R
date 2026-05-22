### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Script to explore pC1 inputs

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

#####################
### get metadata ####
#####################

# pull BANC metadata
bc <- banctable_query()
bc.meta <- bc %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
    is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
    TRUE ~ cell_type
  ),
  top_nt = gsub("auto\\:","",top_nt),
  cell_sub_class = gsub("auto\\:","",cell_sub_class),
  cell_class = gsub("auto\\:","",cell_class),
  super_class = gsub("auto\\:","",super_class)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# pull fafb metadata

# Get the meta data, cell types, etc.
fafb.meta <- arrow::read_feather(file.path(banc.connectivity.save.path, "flywire_783_meta.feather"))

##############################
### Explore putative SAGs ####
##############################

# sag ids
sag.ids <- bc.meta$root_id[grepl("AN_SMP_2", bc.meta$cell_type)]

# Find neuron in manc
cell_type2use <- "SAG"
fafb.id2use <- fafb.meta$root_783[grepl(cell_type2use, fafb.meta$hemibrain_type)]

# plot all neurons
bancsee(banc_ids = sag.ids,
        fafb_ids = fafb.id2use, open = TRUE)

###########################
### explore SAGs inputs ####
###########################

# get banc edgelist from local
datadir <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/data/dn_connectivity/"
banc.el <- readRDS(str_c(datadir, "banc_edgelist.rds"))

# generate metadata for pre/post synaptic partners
banc.pre <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.pre) <- paste0("pre_", colnames(banc.pre))
banc.post <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.post) <- paste0("post_", colnames(banc.post))

# add metadata to synaptic partners
banc.el <- banc.el %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id)) %>%
  # generate normalized input weight
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(post_count = sum(n),
                pre_norm = n/post_count) %>%
  dplyr::ungroup() %>%
  # generate normalized output weight
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(pre_count = sum(n),
                post_norm = n/post_count) %>%
  dplyr::ungroup() %>%  
  dplyr::left_join(banc.pre,
                   by = c("pre" = "pre_root_id")) %>%
  dplyr::left_join(banc.post,
                   by = c("post" = "post_root_id")) %>%
  dplyr::rename(count = n)

# extract inputs to SAGs
norm_syn_thres <- 0.0025
abs_syn_thres <- 5

sag.inputs <- banc.el %>% 
  dplyr::filter(post_pt_root_id %in% sag.ids) %>% 
  dplyr::filter(count > abs_syn_thres) %>%
  arrange(desc(pre_norm))

unique(sag.inputs$pre_cell_type)

sag.inputs.ans.skt <- banc_read_l2skel(unique(sag.inputs$pre))

# plot top cell types
celltype2plot <- "IN"
celltype2plot <- "EN"
celltype2plot <- "MN"

ids2plot <- unique(sag.inputs$pre[grepl(celltype2plot, sag.inputs$pre_cell_type)])

clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(sag.inputs.ans.skt[ids2plot], lwd = 0.5, soma = 2500)

bancsee(banc_ids = ids2plot,
        fafb_ids = fafb.id2use, open = TRUE)

#######################################
### check outputs of AN_multi_107  ####
#######################################

# extract outputs to AN_multi_107
abs_syn_thres <- 5

AN_multi_107.outpus <- banc.el %>%
  dplyr::filter(pre_pt_root_id %in% AN_multi_107.id) %>%
  dplyr::filter(count > abs_syn_thres) %>%
  arrange(desc(post_norm))

AN_multi_107.outpus.skt <- banc_read_l2skel(unique(AN_multi_107.outpus$post))

# plot some cell types
# celltype2plot <- "CL313"
# celltype2plot <- "CB0666"
# celltype2plot <- "auto:DNp13"
celltype2plot <- "AVLP569"
celltype2plot <- "CB1161"
celltype2plot <- "auto:CL123,CRE061"
celltype2plot <- "auto:vpoEN"
celltype2plot <- "AVLP"
celltype2plot <- "CB"
celltype2plot <- "CL"

ids2plot <- unique(AN_multi_107.outpus$post[grepl(celltype2plot, AN_multi_107.outpus$post_cell_type)])
#ids2plot <- AN_multi_107.outpus$post[19]

clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(AN_multi_107.outpus.skt[ids2plot], lwd = 0.5, soma = 2500)