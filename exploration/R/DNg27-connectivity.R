### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Find neurons connected to DNg27 - to find corazonin

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

#####################
### get metadata ####
#####################

bancpipeline.dir <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/"
setwd(bancpipeline.dir)

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

# load from local
datadir <- paste0(bancpipeline.dir, "data/abd_dn_an/")
bc.meta <- readRDS(str_c(datadir, "meta_banc.rds"))
# saveRDS(bc.meta, file = str_c(datadir, "meta_banc.rds"))

################################################
### get edgelist and populate with metadata ####
################################################

# load banc edgelist from local
datadir <- paste0(bancpipeline.dir, "data/dn_connectivity/")
banc.el <- readRDS(str_c(datadir, "banc_edgelist.rds"))

# collect metadata for pre and postsynaptic partners
banc.pre <- bc.meta %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.pre) <- paste0("pre_", colnames(banc.pre))
banc.post <- bc.meta %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(banc.post) <- paste0("post_", colnames(banc.post))

# append metadata to edgelist
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

###########################################
### define DNg27 and get inputs/outputs ####
###########################################

# DNg27 ids
DNg27 <- bc.meta %>%
  dplyr::filter(grepl("DNg27", cell_type))
DNg27.ids <- DNg27$root_id

# extract inputs to DNg27
norm_syn_thres <- 0.0025
abs_syn_thres <- 5

DNg27.inputs <- banc.el %>% 
  dplyr::filter(post_pt_root_id %in% DNg27.ids) %>% 
  dplyr::filter(pre_norm > norm_syn_thres, count > abs_syn_thres) %>% 
  arrange(desc(pre_norm))
  
# extract outputs of DNg27
abs_syn_thres <- 5

DNg27.outputs <- banc.el %>% 
  dplyr::filter(pre_pt_root_id %in% DNg27.ids) %>% 
  dplyr::filter(count > abs_syn_thres) %>% 
  arrange(desc(post_norm))
  
# print outputs
DNg27.outputs.brain <- DNg27.outputs %>%
  dplyr::filter(!grepl("vnc", post_region)) %>% 
  arrange(desc(post_norm))

print(str_c(DNg27.outputs$post, "__",
            round(DNg27.outputs$post_norm*100, digits = 3), "%__",
            DNg27.outputs$post_cell_type, "%__",
            DNg27.outputs$post_region))

print(str_c(DNg27.outputs.brain$post, "__",
            round(DNg27.outputs.brain$post_norm*100, digits = 3), "%__",
            DNg27.outputs.brain$post_cell_type, "%__",
            DNg27.outputs.brain$post_region))

# save results
datadir <- paste0(bancpipeline.dir, "data/vigette_idea_1/")
saveRDS(pC1.inputs, file = str_c(datadir, "pC1_inputs.rds"))
saveRDS(pC1.outputs, file = str_c(datadir, "pC1_outputs.rds"))