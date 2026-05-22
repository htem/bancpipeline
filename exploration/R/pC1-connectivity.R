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
### Explore putative pC1s ####
##############################

cell_type2use <- "pC1"

# Find neuron in manc
fafb.id2use <- fafb.meta$root_783[grepl(cell_type2use, fafb.meta$hemibrain_type)]

# Find neuron in banc
banc.id2use <- bc.meta$root_id[grepl(cell_type2use, bc.meta$cell_type)]

# plot all neurons
bancsee(banc_ids = banc.id2use,
        fafb_ids = fafb.id2use, open = TRUE)

################################################
### get edgelist and populate with metadata ####
################################################

# get banc edgelist
# get from local
setwd("C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/data/dn_connectivity/")
datadir <- getwd()
banc.el <- readRDS(str_c(datadir, "banc_edgelist.rds"))

# get from remote
# banc.el <- banc_edgelist()
# saveRDS(banc.el, file = str_c(datadir, "banc_edgelist.rds"))

banc.pre <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.pre) <- paste0("pre_", colnames(banc.pre))
banc.post <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.post) <- paste0("post_", colnames(banc.post))

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

###########################
### explore pC1 inputs ####
###########################

# pC1 ids
pC1.ids <- bc.meta$root_id[grepl("pC1", bc.meta$cell_type)]

# extract inputs to pC1s
norm_syn_thres <- 0.0025
abs_syn_thres <- 5

pC1.inputs <- banc.el %>% 
  dplyr::filter(post_pt_root_id %in% pC1.ids) %>% 
  dplyr::filter(pre_norm > norm_syn_thres, count > abs_syn_thres)

# display top AN inputs

pC1.inputs.ans <- pC1.inputs %>% 
  dplyr::filter(grepl("AN", pC1.inputs$pre_cell_type)) %>%
  arrange(desc(pre_norm))

unique(pC1.inputs.ans$pre_cell_type)
# "AN00A006"
# "AN_FLA_SMP_2"
# "AN_multi_107"
# "AN_multi_82"
# "AN05B103"
# "AN05B101"
# "AN_4_None"
# "AN01A055"

pC1.inputs.ans.skt <- banc_read_l2skel(unique(pC1.inputs.ans$pre))

# plot top cell types
celltype2plot <- "AN_"
ids2plot <- unique(pC1.inputs.ans$pre[grepl(celltype2plot, pC1.inputs.ans$pre_cell_type)])
  
clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(pC1.inputs.ans.skt[ids2plot], lwd = 0.5, soma = 2500)

celltype2plot <- "AN_FLA_SMP_2"
ids2plot <- unique(pC1.inputs.ans$pre[grepl(celltype2plot, pC1.inputs.ans$pre_cell_type)])

clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(pC1.inputs.ans.skt[ids2plot], lwd = 0.5, soma = 2500)

#####################################
### Compare pC1 input ANs to SAG ####
#####################################

AN.celltype <- "AN00A006"
AN.celltype <- "AN_FLA_SMP_2"
AN.celltype <- "AN_multi_82"
AN.celltype <- "AN05B103"
AN.celltype <- "AN05B101"
AN.celltype <- "AN_4_None"
AN.celltype <- "AN01A055"

# Find neuron in manc
fafb.id2use <- fafb.meta$root_783[grepl("SAG", fafb.meta$hemibrain_type)]

# Find neuron in banc
banc.id2use <- unique(pC1.inputs.ans$pre[grepl(AN.celltype, pC1.inputs.ans$pre_cell_type)])

# plot all neurons
bancsee(banc_ids = banc.id2use,
        fafb_ids = fafb.id2use, open = TRUE)

###########################
### explore pC1 output ####
###########################

# extract inputs to pC1s
norm_syn_thres <- 0.0025
abs_syn_thres <- 5

pC1.outputs <- banc.el %>% 
  dplyr::filter(pre_pt_root_id %in% pC1.ids) %>% 
  dplyr::filter(count > abs_syn_thres)

# display top DN outputs

pC1.outputs.dns <- pC1.outputs %>% 
  dplyr::filter(grepl("DN", pC1.outputs$post_cell_type)) %>%
  arrange(desc(post_norm))

unique(pC1.outputs.dns$post_cell_type)

# major outputs:
# "DNp37" very strong! (3-14%)
# "DNp46", "DNp60" AND "DNp62" mid (1-2%)
# "DNp68", "DNae001", "DNa11", "DNg13", "DNpe053"

########################
### find other pC1s ####
########################

# get pC1 metadata
pC1.meta <- bc.meta %>%
  dplyr::filter(root_id %in% pC1.ids) %>%
  dplyr::select(root_id, supervoxel_id, cell_type) %>%
  dplyr::mutate(supervoxel_id = as.character(supervoxel_id))
pC1.meta$root_id <- banc_latestid(pC1.meta$root_id)
pC1.meta <- as_tibble(pC1.meta)

# get csv with nblast matches
csv.name <- "C:/Users/Diego/Downloads/brain_and_nerve_cord_nblast_banc_mirror_nblast.csv"
data <- readr::read_csv(csv.name, col_types = readr:::cols(.default = col_character()))
data <- data %>%
  dplyr::filter(pt_supervoxel_id %in% unique(pC1.meta$supervoxel_id)) %>%
  dplyr::mutate(pt_supervoxel_id = as.character(pt_supervoxel_id))
  
# update ids
data$pt_root_id <- banc_latestid(data$pt_root_id)
data$match_root_id  <- banc_latestid(data$match_root_id)

pC1.matches <- data %>%
  dplyr::left_join(pC1.meta, by = c("pt_supervoxel_id" = "supervoxel_id"))

print(str_c(pC1.matches$match_root_id, "   ", pC1.matches$cell_type))