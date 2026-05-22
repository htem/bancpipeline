### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

#############################################
### Influence scores for selected neurons ###
#############################################

# load required libraries
library(bancr)
library(tidyverse)
library(progress)
library(dplyr)
library(tidyr)
library(pheatmap)

# define main directory paths
# Run locally
if (.Platform$OS.type == "windows") {
  bancpipeline.dir <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/"
} else if (.Platform$OS.type == "unix") {
  bancpipeline.dir <- "/Users/diegopinedo/Dropbox (Personal)/LabScripts/ImportedTools/bancpipeline/"
}
setwd(bancpipeline.dir)

# setup some basic variables
# usethis::edit_file(paste0(bancpipeline.dir, "analysis/R/abd-startup.R"))
source(paste0(bancpipeline.dir, "analysis/R/abd-startup.R"))

setwd(bancpipeline.dir)

# load metadata from frankenbrain
# it creates: 
#   franken.orig
#   franken.meta
#   franken.meta.pre
#   franken.meta.post
# usethis::edit_file(paste0(banc.path, "R/startup/franken-meta.R"))
source(paste0(banc.path, "R/startup/franken-meta.R"))
rm(franken.meta_duplicates, franken.meta_unique, duplicated_ids)

# load metadata from banc
# it creates: 
#   bc.orig
#   banc.meta
#   banc.meta.pre
#   banc.meta.post
# usethis::edit_file(paste0(banc.path, "R/startup/banc-meta.R"))
source(paste0(banc.path, "R/startup/banc-meta.R"))

# update ids for ANs and DNs in frankenbrain
franken.ids <- unique(franken.meta$banc_id[grepl("descending_neuron", franken.meta$cell_class) |
                                       grepl("ascending_neuron", franken.meta$cell_class)])
franken.ids <- franken.ids[!is.na(franken.ids)]

franken.lids <- banc_latestid(franken.ids)
# find indexes of ids to be updated
id.indices <- match(franken.meta$banc_id, franken.ids)
valid.indices <- which(!is.na(id.indices))
# replace ids with latest ids
franken.meta$banc_id[valid.indices] <- franken.ids[id.indices[valid.indices]]
rm(franken.ids, franken.lids, id.indices, valid.indices)

# update cell type name:
franken.meta$cell_type[grepl("DNxn177", franken.meta$cell_type)] <- "DNpe053"

# fix some names in frankenbrain
# fix region name
franken.meta$region[grepl("o,p,t,i,c,_,l,o,b,e,s", franken.meta$region)] <- "optic_lobes"
# fix cell_function name
franken.meta$cell_function[grepl("o,c,e,l,l,a,r", franken.meta$cell_function)] <- "ocellar"
franken.meta$cell_function[grepl("v,i,s,u,a,l,_,c,h,r,o,m,a,t,i,c", franken.meta$cell_function)] <- "visual_chromatic"
franken.meta$cell_function[grepl("v,i,s,u,a,l,_,a,c,h,r,o,m,a,t,i,c", franken.meta$cell_function)] <- "visual_achromatic"
# fix super_class name
franken.meta$super_class[grepl("o,p,t,i,c", franken.meta$super_class)] <- "optic"
franken.meta$super_class[grepl("s,e,n,s,o,r,y", franken.meta$super_class)] <- "sensory"
# fix cell_class name
franken.meta$cell_class[grepl("v,i,s,u,a,l", franken.meta$cell_class)] <- "visual"
franken.meta$cell_class[grepl("o,p,t,i,c,_,l,o,b,e,s", franken.meta$cell_class)] <- "optic_lobes"
franken.meta$cell_class[grepl("M,E,>,L,A", franken.meta$cell_class)] <- "ME>LA"
franken.meta$cell_class[grepl("L,A,>,M,E", franken.meta$cell_class)] <- "LA>ME"
franken.meta$cell_class[grepl("L,A", franken.meta$cell_class)] <- "LA"

# load all influence scores
# it creates: efferent.ids and motor.meta
#             influence.meta that contains influence scores with metadata
#                 temporary variables: 
#                     influence.df, influence.df.log (normalized after log)
#                     rev.influence.df, rev.influence.coarse.df
#             main types of influence scores: 
#                 forward_coarse (includes body parts, and modalitites:
#                                 olfactory|hygrosensory|thermosensory|visual|
#                                 chemosensory|proprioceptive|tactile|visual_projection)
#                 forward_fine (all from seed)
#                 reverse_fine (reverse from AN and DN types)
#                 reverse_coarse ()
# usethis::edit_file(paste0(banc.path, "R/startup/franken-seeds.R"))
source(paste0(banc.path, "R/startup/franken-seeds.R"))

# usethis::edit_file(paste0(banc.path, "R/startup/banc_build_influence.R"))
source(paste0(banc.path, "R/startup/banc_build_influence.R"))
rm(influence.master, seed.groups, seed.sets, con, table_name)

# loads influence scores
# it creates:
#   andn.meta
#   influence.neck.meta
#     post_ is metadata for each id
#     pre_ is metadata for each seed (by pre_cell_type)
#     excluding post_super_class: glia, afferent|sensory, or NA
#     define post_super_class: cases of sez
#     what do you select with this: 
#       !(grepl("forward", resolution) & grepl("neck_connecive", post_region)) | 
#       !(grepl("reverse", resolution) & !grepl("neck_connecive", post_region))
# usethis::edit_file(paste0(banc.path, "R/startup/banc-andn-influence.R"))
source(paste0(banc.path, "R/startup/banc-andn-influence.R"))

# usethis::edit_file(paste0(banc.path, "R/startup/banc-influence.R"))
source(paste0(banc.path, "R/startup/banc-influence.R"))

# Aggregate influences across sides and across JO
influence.meta <- influence.meta %>%
  dplyr::mutate(seed = dplyr::case_when(
    grepl("JO_",seed) ~ "JON",
    grepl("_right$|_r$|_left$|_l$",seed) ~ gsub("_right$|_r$|_left$|_l$","",seed),
    TRUE ~ seed
  )) %>%
  dplyr::mutate(seed = dplyr::case_when(
    grepl("_right$|_left$",seed) ~ gsub("_right$|_left$","",seed),
    TRUE ~ seed
  )) %>%
  dplyr::group_by(seed, resolution, id) %>%
  dplyr::mutate(influence = sum(influence, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(seed, resolution) %>%
  dplyr::mutate(influence_scaled = as.vector((influence - median(influence, na.rm = TRUE))/stats::mad(influence, na.rm = TRUE))) %>%
  dplyr::ungroup()

# load pC1 inputs and outputs
temp.directory <- paste0(bancpipeline.dir, "data/vigette_idea_1/")
pC1.inputs <- readRDS(str_c(temp.directory, "pC1_inputs.rds"))
pC1.outputs <- readRDS(str_c(temp.directory, "pC1_outputs.rds"))
rm(temp.directory)

#################
### plot pC1s ###
#################

pC1types <- unique(pC1.inputs$post_cell_type)
pC1types <- sort(pC1types)
  
# select type of influence scores
influence2use = "influence_scaled"

# define figure specs
width.siz <- 5
height.siz <- c(10, 20)
dpi.png <- 100

# define figure name and directory
figure.dir <- "C:/Users/Diego/Desktop/tempfiles/"
figure.name <- "pC1s"

plot_sensory_influence_from_ct(pC1types, influence2use, influence.meta, 
                               scaled_heatmap_palette3, scaled_heatmap_breaks3,
                               figure.dir, figure.name, width.siz, height.siz, dpi.png)

####################################################
### plot top 10 pC1s inputs/outputs per pC1 type ###
####################################################

# print top inputs to pC1s per type
pC1.inputs.coarse <- get_summary_inputs(pC1.inputs)

# define figure name and directory
figure.dir <- "C:/Users/Diego/Desktop/tempfiles/"

# Note:
#   oviEN is SMP550
#   "PAL02", and "SLP212a", "SLP212b", "SLP212c" are other oviDN inputs

for (ct2use in unique(pC1.inputs.coarse$post_cell_type)) {
  
  figure.name <- paste0(ct2use, "_and_inputs")
  ct2plot <- unique(pC1.inputs.coarse$pre_cell_type[grepl(ct2use, pC1.inputs.coarse$post_cell_type)])
  ct2plot <- ct2plot[ct2plot != ct2use & !is.na(ct2plot)]
  ct2plot <- c(ct2use, ct2plot)

  plot_sensory_influence_from_ct(ct2plot[1:10], influence2use, influence.meta, 
                                 scaled_heatmap_palette3, scaled_heatmap_breaks3,
                                 figure.dir, figure.name, width.siz, height.siz, dpi.png)
  
}

# print top outputs to pC1s per type
pC1.outputs.coarse <- get_summary_outputs(pC1.outputs)

# Note:
#   aIPgs are SMP555 (aIPg_a), CB1127 (aIPg_b), CB1877(aIPg_c)

# define figure name and directory
figure.dir <- "C:/Users/Diego/Desktop/tempfiles/"

for (ct2use in unique(pC1.outputs.coarse$pre_cell_type)) {
  
  figure.name <- paste0(ct2use, "_and_outputs")
  ct2plot <- unique(pC1.outputs.coarse$post_cell_type[grepl(ct2use, pC1.outputs.coarse$pre_cell_type)])
  ct2plot <- ct2plot[ct2plot != ct2use & !is.na(ct2plot)]
  ct2plot <- c(ct2use, ct2plot)
  
  plot_sensory_influence_from_ct(ct2plot[1:10], influence2use, influence.meta, 
                                 scaled_heatmap_palette3, scaled_heatmap_breaks3,
                                 figure.dir, figure.name, width.siz, height.siz, dpi.png)
  
  plot_reverse_influence_from_ct(ct2plot[1:10], influence2use, influence.meta, 
                                 scaled_heatmap_palette3, scaled_heatmap_breaks3,
                                 figure.dir, figure.name, width.siz, height.siz, dpi.png)
  
  }

#############################################
### plot pC1 input ANs and pC1 output DNs ###
#############################################

# evaluate input ANs to pC1s
# display input ANs found vs available for sensory influence analysis
pC1.inputs.ans <- pC1.inputs %>% 
  dplyr::filter(grepl("AN", pC1.inputs$pre_cell_type)) %>% 
  arrange(desc(pre_norm))

# cell types to search for
query.ct <- unique(pC1.inputs.ans$pre_cell_type[
  pC1.inputs.ans$pre_pt_root_id %in% unique(pC1.inputs.ans$pre_pt_root_id)])
# cell types found
existing.ct <- unique(pC1.inputs.ans$pre_cell_type[
  pC1.inputs.ans$pre_pt_root_id %in% unique(influence.meta$id)])

# overwrite cell_type for some banc_id
ids2edit <- pC1.inputs.ans %>%
  dplyr::distinct(pre_pt_root_id, .keep_all = TRUE) %>%
  dplyr::select(pre_pt_root_id, pre_cell_type)
id_indices <- match(influence.meta$id, ids2edit$pre_pt_root_id)
valid_indices <- which(!is.na(id_indices))
# Replace IDs with a match
influence.meta$cell_type[valid_indices] <- ids2edit$pre_cell_type[id_indices[valid_indices]]
rm(ids2edit, id_indices, valid_indices)

# display
cat("Cell types to look for:\n", query.ct, "\n",
    "Found:\n", existing.ct, "\n",
    "Not found:\n", setdiff(query.ct, existing.ct), "\n")

# evaluate output DNs from pC1s
# display output DNs found vs available for sensory influence analysis
pC1.outputs.dns <- pC1.outputs %>% 
  dplyr::filter(grepl("DN", pC1.outputs$post_cell_type)) %>% 
  arrange(desc(post_norm))

# cell types to search for
query.ct <- unique(pC1.outputs.dns$post_cell_type[
  pC1.outputs.dns$post_pt_root_id %in% unique(pC1.outputs.dns$post_pt_root_id)])
# cell types found
existing.ct <- unique(pC1.outputs.dns$post_cell_type[
  pC1.outputs.dns$post_pt_root_id %in% unique(influence.meta$id)])

# overwrite cell_type for some banc_id
ids2edit <- pC1.outputs.dns %>%
  dplyr::distinct(post_pt_root_id, .keep_all = TRUE) %>%
  dplyr::select(post_pt_root_id, post_cell_type)
id_indices <- match(influence.meta$id, ids2edit$post_pt_root_id)
valid_indices <- which(!is.na(id_indices))
# Replace IDs with a match
influence.meta$cell_type[valid_indices] <- ids2edit$post_cell_type[id_indices[valid_indices]]
rm(ids2edit, id_indices, valid_indices)

# display
cat("Cell types to look for:\n", query.ct, "\n",
    "Found:\n", existing.ct, "\n",
    "Not found:\n", setdiff(query.ct, existing.ct), "\n")

# select type of influence scores
influence2use = "influence_scaled"

# define figure specs
width.siz <- 5
height.siz <- c(10, 20)
dpi.png <- 100

# define figure name and directory
figure.dir <- "C:/Users/Diego/Desktop/tempfiles/"
figure.name <- "ANs_pC1in"

plot_sensory_influence_from_id(unique(pC1.inputs.ans$pre_pt_root_id), 
                              influence2use, influence.meta, 
                              scaled_heatmap_palette3, scaled_heatmap_breaks3,
                              figure.dir, figure.name, width.siz, height.siz, dpi.png)

# define figure name and directory
figure.dir <- "C:/Users/Diego/Desktop/tempfiles/"
figure.name <- "DNs_pC1out"

plot_sensory_influence_from_id(unique(pC1.outputs.dns$post_pt_root_id), 
                               influence2use, influence.meta, 
                               scaled_heatmap_palette3, scaled_heatmap_breaks3,
                               figure.dir, figure.name, width.siz, height.siz, dpi.png)

plot_reverse_influence_from_id(unique(pC1.outputs.dns$post_pt_root_id), 
                               influence2use, influence.meta, 
                               scaled_heatmap_palette3, scaled_heatmap_breaks3,
                               figure.dir, figure.name, width.siz, height.siz, dpi.png)

#######################################################
### Select influence scores for neurons of interest ###
#######################################################

# define cell types to explore
celltype2use <- c(# oviDNs
          "oviDNa_a", "oviDNa_b", "oviDNb",
          # pC1 main outputs
         "DNp37", "DNp60", "DNp62", "DNp46", "DNp68",
         "DNa11", "DNae001", "DNg13", "DNpe053")

##############################
### plot cell type anatomy ###
##############################

# identify and plot cells using celltype
celltype2use.meta <- banc.meta %>%
  dplyr::select(root_id, cell_type) %>%
  dplyr::filter(cell_type %in% celltype2use) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
celltype2use.skt <- banc_read_l2skel(unique(celltype2use.meta$root_id))

figure.dir <- "C:/Users/Diego/Desktop/tempfiles/anatomy/"

for(ct2p in unique(celltype2use.meta$cell_type)){
  
  ids2plot <- celltype2use.meta$root_id[grepl(ct2p, celltype2use.meta$cell_type)]
  figure.name <- paste0("skt_", ct2p)
  
  plot_sel_skt(ids2plot, celltype2use.skt, figure.dir, figure.name)
  
}

# identify and plot cells using id
ids2use <- banc_latestid(unique(pC1.inputs.ans$pre))
celltype2use.meta <- banc.meta %>%
  dplyr::select(root_id, cell_type) %>%
  dplyr::filter(root_id %in% ids2use) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
celltype2use.skt <- banc_read_l2skel(unique(celltype2use.meta$root_id))

figure.dir <- "C:/Users/Diego/Desktop/tempfiles/anatomy/"

for(ct2p in unique(celltype2use.meta$cell_type)){
  
  ids2plot <- celltype2use.meta$root_id[grepl(ct2p, celltype2use.meta$cell_type)]
  figure.name <- paste0("skt_", ct2p)
  
  plot_sel_skt(ids2plot, celltype2use.skt, figure.dir, figure.name)
  
}

#################
### functions ###
#################

# Extract matrix with influence scores of selected cell types
get_influence_matrix_from_celltype <- function(celltype2use, resolution2use, influence2use, influence.meta) {

  # update some custom cell types
  celltype2use <- manual_celltype_renaming(celltype2use)
  
  # Get matrix of scores
  influence_matrix <- influence.meta %>%
    # choose particular resolution type
    dplyr::filter(resolution == resolution2use,
                  !is.na(influence)) %>%
    # get just selected cell types   
    dplyr::mutate(cell_type = manual_celltype_renaming(cell_type)) %>%
    dplyr::filter(cell_type %in% celltype2use) %>%
    # reshape the data from a long format to a wide format
    reshape2::dcast(seed ~ cell_type,
                    # it computes the median across cells of the same cell type
                    fun.aggregate = median,
                    # choose influence score to use
                    value.var = influence2use,
                    # missing data is filled with 10
                    fill = 10)
  
  # Set row names and remove the seed column
  rownames(influence_matrix) <- influence_matrix$seed
  influence_matrix$seed <- NULL
  
  # Convert to matrix
  influence_matrix <- as.matrix(influence_matrix)

  # Reorder columns to match the order in celltype2use
  order_indices <- match(celltype2use, colnames(influence_matrix))
  valid_indices <- order_indices[!is.na(order_indices)]
  cat("Cell types not found:\n", paste(setdiff(celltype2use, colnames(influence_matrix)), collapse = "\n"), "\n")
  influence_matrix <- influence_matrix[, valid_indices, drop = FALSE]
  
  return(influence_matrix = influence_matrix)
  
}

get_influence_matrix_from_id <- function(id2use, resolution2use, influence2use, influence.meta) {
  
  # Get matrix of scores
  influence.meta <- influence.meta %>%
    # choose particular resolution type
    dplyr::filter(resolution == resolution2use,
                  !is.na(influence)) %>%
    # get just selected ids   
    dplyr::filter(id %in% id2use) %>%
    dplyr::mutate(cell_type = manual_celltype_renaming(cell_type))
  
  # make mini influence.meta
  id_and_celltype <- influence.meta %>%
    dplyr::select(id, cell_type) %>%
    dplyr::distinct(id,.keep_all = TRUE)
  
  influence_matrix <- influence.meta %>%
    # reshape the data from a long format to a wide format
    reshape2::dcast(seed ~ cell_type,
                    # it computes the median across cells of the same cell type
                    fun.aggregate = median,
                    # choose influence score to use
                    value.var = influence2use,
                    # missing data is filled with 10
                    fill = 10)
  
  # Set row names and remove the seed column
  rownames(influence_matrix) <- influence_matrix$seed
  influence_matrix$seed <- NULL
  
  # Convert to matrix
  influence_matrix <- as.matrix(influence_matrix)
  
  # Reorder columns to match the order in id2use
  cat("IDs types not found:\n", paste(id_and_celltype$cell_type[!id_and_celltype$id %in% id2use]), collapse = "\n")

  return(influence_matrix = influence_matrix)
  
}

# Generate annotations for brain regions to populate rows
generate_region_annotation <- function(franken.meta, influence_matrix2use) {

  # Region annotations
  region_annotation <- franken.meta %>%
    dplyr::filter(!is.na(cell_function)) %>%
    dplyr::select(cell_function, region) %>%
    dplyr::distinct(cell_function, .keep_all = TRUE) %>%
    as.data.frame()
  rownames(region_annotation) <- region_annotation$cell_function
  region_annotation <- region_annotation[rownames(region_annotation) %in% rownames(influence_matrix2use),]
  region_annotation$cell_function <- NULL

  return(region_annotation = region_annotation)
  
}

# Generate cell type annotations for selected neurons to populate columns
generate_celltype_annotation <- function(influence.meta, celltype2use, influence_matrix2use) {
  
  # Get cell types for each pre_root_id
  col_annotation <- influence.meta %>%
    dplyr::filter(!is.na(cell_type)) %>%
    dplyr::select(cell_type, super_class) %>%
    dplyr::filter(cell_type %in% celltype2use) %>%
    dplyr::distinct(cell_type, .keep_all = TRUE) %>%
    as.data.frame()
  rownames(col_annotation) <- as.character(col_annotation$cell_type)
  col_annotation <- col_annotation[rownames(col_annotation) %in% colnames(influence_matrix2use),]
  col_annotation <- data.frame(cell_type = as.character(col_annotation$cell_type))
  rownames(col_annotation) <- col_annotation$cell_type
 
  return(col_annotation = col_annotation)
   
}

# Generate heatmap
simple_plot <- function(influence_matrix2use, annotation2use_row, 
                        annotation2use_column, colors2use, breaks2use, 
                        figure.dir, figure.name, width.siz, height.siz, dpi.png){
  
  # Calculate optimal cell dimensions based on desired total dimensions
  num_rows <- nrow(influence_matrix2use)
  num_cols <- ncol(influence_matrix2use)
  
  cellwidth_optimal <- (0.4 * width.siz) / num_cols
  cellheight_optimal <- (0.4 * height.siz) / num_rows
  
  # deal with na
  influence_matrix2use[is.na(influence_matrix2use)] <- 0
  influence_matrix2use[is.nan(influence_matrix2use)] <- 0
  #influence_matrix2use[is.infinite(influence_matrix2use)] <- 0
  
  # generate heatmap
  # if (length(annotation2use_row) > 0 && nchar(annotation2use_row[1]) > 0) {
  #  
  #   ph4 <- pheatmap(
  #     influence_matrix2use,
  #     cluster_rows = TRUE,
  #     cluster_cols = TRUE,
  #     color = colors2use,
  #     breaks = breaks2use,
  #     show_rownames = TRUE,
  #     show_colnames = TRUE,
  #     main = "",
  #     fontsize = 8,
  #     cellwidth = cellwidth_optimal * dpi.png,
  #     cellheight = cellheight_optimal * dpi.png,
  #     cellwidth = cellwidth_optimal * dpi.png,
  #     cellheight = cellheight_optimal * dpi.png,
  #     annotation_row = annotation2use_row,
  #     annotation_col = annotation2use_column,
  #     annotation_legend = FALSE,
  #     legend = TRUE,
  #     width = unit(width.siz, "inches"),
  #     height = unit(height.siz, "inches"),
  #     treeheight_row = 5, 
  #     treeheight_col = 5
  #   )
  #    
  # } else {
    
    ph4 <- pheatmap(
      influence_matrix2use,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colors2use,
      breaks = breaks2use,
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = "",
      fontsize = 8,
      legend = TRUE,
      width = unit(width.siz, "inches"),
      height = unit(height.siz, "inches"),
      treeheight_row = 5, 
      treeheight_col = 5
    )
    
  #}
  
  if (length(figure.name) > 0 && nchar(figure.name[1]) > 0) {

    if (!dir.exists(figure.dir)) {
      dir.create(figure.dir, recursive = TRUE)
    } 
    
    output_file <- str_c(figure.dir, figure.name, ".png")
    png(filename = output_file, width = width.siz, height = height.siz, units = "in", res = dpi.png)
    
    print(ph4)
    
    dev.off()
    
  } else {

    print(ph4)
        
  }

}

# Plot skeletons of selected cell types
plot_sel_skt <- function(ids2plot, skt, figure.dir, figure.name) {
  
  # get neurons skeletons
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
  
  close3d()
  
}

# Plot sensory influences
plot_sensory_influence_from_ct <- function(celltype2use, influence2use, influence2use.meta, 
                        scaled_heatmap_palette3, scaled_heatmap_breaks3, 
                        figure.dir, figure.name, width.siz, height.siz, dpi.png){
  
  # pull coarse scores
  resolution2use = "forward_coarse"
  influence_matrix_f_coarse <- get_influence_matrix_from_celltype(celltype2use, resolution2use, influence2use, influence2use.meta)
  
  # edit influence names and reorder
  rownames(influence_matrix_f_coarse) <- gsub("\\bendocrine\\b", "vnc_endocrine", rownames(influence_matrix_f_coarse))
  influence_matrix_f_coarse <- reorder_forward_coarse_scores(influence_matrix_f_coarse)
  
  # pull fine scores
  resolution2use = "forward_fine"
  influence_matrix_f_fine <- get_influence_matrix_from_celltype(celltype2use, resolution2use, influence2use, influence2use.meta)
  
  # reorder coarse sensory influences
  influence_matrix_f_fine <- reorder_forward_fine_scores(influence_matrix_f_fine)
  
  # plot heatmaps
  figure.name_1 <- paste0(figure.name, "_sensory_coarse")
  simple_plot(influence_matrix_f_coarse, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_1, width.siz, height.siz[1], dpi.png)
  
  figure.name_2 <- paste0(figure.name, "_sensory_fine")
  simple_plot(influence_matrix_f_fine, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_2, width.siz, height.siz[2], dpi.png)
  
}

plot_sensory_influence_from_id <- function(id2use, influence2use, influence2use.meta, 
                                            scaled_heatmap_palette3, scaled_heatmap_breaks3, 
                                            figure.dir, figure.name, width.siz, height.siz, dpi.png){
  
  # pull coarse scores
  resolution2use = "forward_coarse"
  influence_matrix_f_coarse <- get_influence_matrix_from_id(id2use, resolution2use, influence2use, influence2use.meta)
  
  # edit influence names and reorder
  rownames(influence_matrix_f_coarse) <- gsub("\\bendocrine\\b", "vnc_endocrine", rownames(influence_matrix_f_coarse))
  influence_matrix_f_coarse <- reorder_forward_coarse_scores(influence_matrix_f_coarse)
  
  # pull fine scores
  resolution2use = "forward_fine"
  influence_matrix_f_fine <- get_influence_matrix_from_id(id2use, resolution2use, influence2use, influence2use.meta)
 
  # reorder coarse sensory influences
  influence_matrix_f_fine <- reorder_forward_fine_scores(influence_matrix_f_fine)
  
  # plot heatmaps
  figure.name_1 <- paste0(figure.name, "_sensory_coarse")
  simple_plot(influence_matrix_f_coarse, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_1, width.siz, height.siz[1], dpi.png)
  
  figure.name_2 <- paste0(figure.name, "_sensory_fine")
  simple_plot(influence_matrix_f_fine, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_2, width.siz, height.siz[2], dpi.png)
  
}

# Plot motor influences
plot_reverse_influence_from_ct <- function(celltype2use, influence2use, influence2use.meta, 
                                           scaled_heatmap_palette3, scaled_heatmap_breaks3, 
                                           figure.dir, figure.name, width.siz, height.siz, dpi.png){
  
  # pull coarse scores
  resolution2use = "reverse_coarse"
  influence_matrix_f_coarse <- get_influence_matrix_from_celltype(celltype2use, resolution2use, influence2use, influence2use.meta)
  
  # edit influence names and reorder
  influence_matrix_f_coarse <- reorder_reverse_coarse_scores(influence_matrix_f_coarse)
  
  # pull fine scores
  resolution2use = "reverse_fine"
  influence_matrix_f_fine <- get_influence_matrix_from_celltype(celltype2use, resolution2use, influence2use, influence2use.meta)
  
  # reorder coarse sensory influences
  influence_matrix_f_fine_1 <- get_reverse_fine_scores_group_1(influence_matrix_f_fine)
  influence_matrix_f_fine_2 <- get_reverse_fine_scores_group_2(influence_matrix_f_fine)
  
  # plot heatmaps
  figure.name_1 <- paste0(figure.name, "_reverse_fine_1")
  simple_plot(influence_matrix_f_fine_1, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_1, width.siz, height.siz[2], dpi.png)
  
  figure.name_2 <- paste0(figure.name, "_reverse_fine_2")
  simple_plot(influence_matrix_f_fine_2, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_2, width.siz, height.siz[2], dpi.png)
  
}

plot_reverse_influence_from_id <- function(id2use, influence2use, influence2use.meta, 
                                           scaled_heatmap_palette3, scaled_heatmap_breaks3, 
                                           figure.dir, figure.name, width.siz, height.siz, dpi.png){
  
  # pull coarse scores
  resolution2use = "reverse_coarse"
  influence_matrix_f_coarse <- get_influence_matrix_from_id(id2use, resolution2use, influence2use, influence2use.meta)
  
  # edit influence names and reorder
  influence_matrix_f_coarse <- reorder_reverse_coarse_scores(influence_matrix_f_coarse)
  
  # pull fine scores
  resolution2use = "reverse_fine"
  influence_matrix_f_fine <- get_influence_matrix_from_id(id2use, resolution2use, influence2use, influence2use.meta)
  
  # reorder coarse sensory influences
  influence_matrix_f_fine_1 <- get_reverse_fine_scores_group_1(influence_matrix_f_fine)
  influence_matrix_f_fine_2 <- get_reverse_fine_scores_group_2(influence_matrix_f_fine)
  
  # plot heatmaps
  figure.name_1 <- paste0(figure.name, "_reverse_coarse")
  simple_plot(influence_matrix_f_coarse, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_1, width.siz, height.siz[1], dpi.png)
  
  figure.name_2 <- paste0(figure.name, "_reverse_fine_1")
  simple_plot(influence_matrix_f_fine_1, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_2, width.siz, height.siz[2], dpi.png)

  figure.name_3 <- paste0(figure.name, "_reverse_fine_2")
  simple_plot(influence_matrix_f_fine_2, c(""), c(""),
              scaled_heatmap_palette3, scaled_heatmap_breaks3, 
              figure.dir, figure.name_3, width.siz, height.siz[2], dpi.png)
  
}
# Manual editing of seed names
manual_celltype_renaming <- function(cell_type) {
  
  old_celltype_name <- c("SMP550", "SMP555", "CB1127", "CB1877", "CL278", "VES061")
  new_celltype_name <- c("oviEN", "aIPg_a", "aIPg_b", "aIPg_c", "pIP5", "pIP5")
    
  for (i in 1:length(old_celltype_name)) {
    cell_type <- gsub(old_celltype_name[i], new_celltype_name[i], cell_type)
  }

  return(cell_type = cell_type)
  
}

reorder_forward_coarse_scores <- function(influence_matrix){
  
  # custom row order for coarse sensory influences
  row_order <- c("visual_achromatic", "visual_chromatic", "olfactory", "thermosensory", 
                 "hygrosensory", "JO", "eye_bristle", "head_bristle", "ocellar", "neck", 
                 "brain_endocrine", "vnc_endocrine", "gustatory_brain", "ciberial_mechanosensory", "pharynx", 
                 "enteric_gustatory_brain", "sensory_ascending", "ascending", "chemosensory", "tactile", "proprioceptive", 
                 "prothoracic_chordotonal_organ", "front_leg", "middle_leg", "hind_leg", 
                 "wing", "haltere", "notum", "abdomen")
  
  # Reorder the rows of the matrix
  order_indices <- match(row_order, rownames(influence_matrix))
  influence_matrix <- influence_matrix[na.omit(order_indices), , drop = FALSE]
  
  return(influence_matrix = influence_matrix)
  
}

reorder_forward_fine_scores <- function(influence_matrix){
  
  # custom row order for coarse sensory influences
  row_order <- c("visual_achromatic", "visual_chromatic", "visual_centrifugal", "visual_projection",
                 "olfactory", "olfactory_NA", "heating", "cold", "cooling", "evaporative_cooling", "humid", "moist", 
                 "JON", "mechanosensory", "eye_bristle", "head_bristle", "proboscis_bristle", "ocellar", "ocellar_center",
                 "endocrine", # which one is this brain or vnc?
                 "taste_peg", "low_salt", "bitter", "sugar_water", "dry", 
                 "ciberial_mechanosensory_neuron", "accessory_pharyngeal_nerve_sensory_group1", "accessory_pharyngeal_nerve_sensory_group2", 
                 "enteric_gustatory", # change to "enteric_gustatory_brain"
                 "sensory_ascending", 
                 "sensory_ascending_pr_n",
                 # Prosternal nerve (PrN) A slender nerve that projects anteriorly from the ventral nerve cord (VNC), 
                 # medial to the base of the dorsal prothoracic nerve to the prosternal sense organ (Power, 1948)
                 "sensory_ascending_pro_ln", "sensory_ascending_meso_ln", # pro/meso/metathoracic leg nerve
                 "sensory_ascending_d_meta_n",
                 # dorsal mesothoracic nerve (DMetaN), houses sensory axons from halteres
                 "sensory_ascending_admn", "sensory_ascending_pdmn", 
                 # anterior/posterior dorsal mesothoracic nerve, the wing nerve joins ADMN (which houses wing sensory fibers)
                 # pdmn houses motor neurons to jump muscle and wing power muscles
                 "sensory_ascending_ab_n_4", # 4th abdominal nerve (AbN1-4)
                 "chemosensory_pro_ln", "chemosensory_meso_ln", "chemosensory_meta_ln", # pro/meso/metathoracic leg nerve
                 "chemosensory_admn", 
                 "chemosensory_ab_n_2", "chemosensory_ab_n_3", "chemosensory_ab_n_4", # 2nd-4th abdominal nerve (AbN2-4)
                 "chemosensory_ab_nt", # abdominal nerve trunk
                 "SA_VTV_pro_meso_meta", "SA_VTV_PDMN", # SA? Ventral median tract of ventral cervical fasciculus (VTV)
                 "tactile_admn", "tactile_pdmn", 
                 "tactile_v_pro_n", "tactile_d_pro_n", "tactile_d_meta_n", 
                 # ventral prothoracic nerve, and dorsal prothoracic/mesothoracic nerve
                 "tactile_pro_ln", "tactile_meso_ln", "tactile_meta_ln",
                 "proprioceptive_pro_cn", # prothoracic chordotonal nerve
                 "proprioceptive_pro_an",  # prothoracic accessory nerve
                 "proprioceptive_pro_ntbd", # nerve to be determine?
                 "proprioceptive_pr_n", "proprioceptive_v_pro_n", "proprioceptive_d_pro_n", "proprioceptive_d_meta_n", 
                 "proprioceptive_pro_ln", "proprioceptive_meso_ln", "proprioceptive_meta_ln", # pro/meso/metathoracic leg nerve
                 "proprioceptive_admn", "proprioceptive_pdmn", "proprioceptive_ab_n_3", "proprioceptive_ab_nt", 
                 "haltere")
  
  # Reorder the rows of the matrix
  order_indices <- match(row_order, rownames(influence_matrix))
  influence_matrix <- influence_matrix[na.omit(order_indices), , drop = FALSE]
  
  # custom editing of names
  rownames(influence_matrix) <- gsub("\\bendocrine\\b", "brain_endocrine", rownames(influence_matrix))
  
  return(influence_matrix = influence_matrix)
  
}

reorder_reverse_coarse_scores <- function(influence_matrix){
  
  # custom row order for coarse sensory influences
  row_order <- c("eye_motor_neuron", "antennal_motor_neuron", "neck_motor_neuron", # neck motor neuron head
                 "nm_motor_neuron", # neck motor neuron vnc
                 "proboscis_motor_neuron", 
                 "haustellum_motor_neuron", 
                 # haustellum: a part of the proboscis of a Drosophila fly that extends and flexes to help the fly feed
                 "salivary_motor_neuron", "ingestion_motor_neuron", "crop_motor_neuron", 
                 "lateral_NSC", "medial_NSC", "SEZ-NSC",
                 # neurosecretory neurons:
                 #  lateral: CRZ, DH31, ITP, unknown
                 #  medial: DILP, DH44, DMS, unknown
                 #  SEZ: hugin, CAPA
                 "fl_motor_neuron", "ml_motor_neuron", "hl_motor_neuron", # fore/mid/hind leg motor neuron
                 "wm_motor_neuron", # wing motor neuron
                 "hm_motor_neuron", # haltere motor neuron
                 "ad_motor_neuron", # abdomen motor neuron
                 "ad_encodrine_neuron", # abdomen endocrine neuron
                 "xm_motor_neuron", # unknown motor neurons?
                 "EN", # endocrine_vnc
                 "EA" # efferent_ascending
                 )
  
  # Reorder the rows of the matrix
  order_indices <- match(row_order, rownames(influence_matrix))
  influence_matrix <- influence_matrix[na.omit(order_indices), , drop = FALSE]
  
  # custom editing of names
  newnames <- rownames(influence_matrix)
  newnames <- gsub("nm_motor_neuron", "brain_neck_motor_neuron", newnames)
  newnames <- gsub("neck_motor_neuron", "vnc_neck_motor_neuron", newnames)
  newnames <- gsub("wm_motor_neuron", "wing_motor_neuron", newnames)
  newnames <- gsub("hm_motor_neuron", "haltere_motor_neuron", newnames)
  newnames <- gsub("fl_motor_neuron", "foreleg_motor_neuron", newnames)
  newnames <- gsub("ml_motor_neuron", "midleg_motor_neuron", newnames)
  newnames <- gsub("hl_motor_neuron", "hindleg_motor_neuron", newnames)
  newnames <- gsub("ad_motor_neuron", "abd_motor_neuron", newnames)
  newnames <- gsub("ad_encodrine_neuron", "abd_encodrine", newnames)
  newnames <- gsub("\\bEN\\b", "vnc_endocrine", newnames)
  newnames <- gsub("\\bEA\\b", "efferent_ascending", newnames)
  newnames <- gsub("_motor_neuron", "_mn", newnames)
  
  rownames(influence_matrix) <- newnames
  
  return(influence_matrix = influence_matrix)
  
}

get_reverse_fine_scores_group_1 <- function(influence_matrix){
  
  # group endocrine and efferent ascending
  row_order <- c(# eye motor neuron
               "CB0804", "CB0901", 
               # antennal motor neuron
               "CB0723", "CB0750", "CB0810", "CB0873", "CB0886",
               # neck motor neurons
               "neck_motor_neuron", "FNM2",
               "CB0918", "CB0706", "CB0705", "CB0835", "CB0899", "CB0913", "CB0831", "CB0916", "CB0004", "CB0838", 
               # other neck motor neuron
               "MNnm03", "MNnm07", "MNnm08", "MNnm09", "MNnm10", "MNnm11", 
               "MNnm12", "MNnm13", "MNnm14", 
               # ingestion related motor neurons
               "ingestion_motor_neuron", 
               "CB0915", "CB0700", "CB0701", "CB0703", "CB0708", "CB0715", "CB0728", "CB0769", "CB0904", "CB0914", "PRW.143", 
               # salivary MN
               "CB0836", 
               # crop MN
               "CB0764", "CEM",
               # proboscis MN
               "CB0720", "CB0762", "CB0783", "CB0789", "CB0845", "CB0858", "CB0861", "CB0871", "CB0911", 
               # haustellum MN
               "CB0875", "CB0882", 
               # endocrine brain
               "PI", "Hugin-RG", "DNES1", "DNES2", "DNES3", "CAPA",
               # corazonin brain
               "CB0905", "CB0726",
               # endocrine at VNC
               "mesVUM-MJ_midline", "PSI", 
               # VNC endocrine
               "EN21X001", "EN00B003_midline", "EN00B002_midline", "EN00B004_midline", 
               "EN00B027_midline", "EN00B021_midline", "EN00B020_midline", "ENXXX012", 
               "EN00B024_midline", "EN00B012_midline", "EN00B010_midline", "EN27X010", 
               "EN00B009_midline", "EN00B014_midline", "EN00B007_midline", "EN00B008_midline",
               "EN00B026_midline", "EN00B015_midline", "EN00B001_midline", "EN00B023_midline", 
               "EN00B018_midline", "EN00B011_midline", "EN00B017_midline", "EN00B013_midline", 
               "EN00B016_midline", "ENXXX226", "EN00B025_midline", "ENXXX286", "EN00B019_midline", 
               "ENXXX128", 
               # VNC efferent ascending (endocrine ascending)
               "EA00B006_midline", "EA27X006", "EAXXX079", "EA00B022_midline", "EA06B010")

  # Reorder the rows of the matrix
  order_indices <- match(row_order, rownames(influence_matrix))
  influence_matrix <- influence_matrix[na.omit(order_indices), , drop = FALSE]
  
  return(influence_matrix = influence_matrix)
  
}

get_reverse_fine_scores_group_2 <- function(influence_matrix){
  
  row_order <- c(# indirect flight muscles
    "DLMn a, b", "DLMn c-f", "DVMn 1a-c", "DVMn 2a, b", "DVMn 3a, b", "hDVM MN", 
    # direct flight muscles
    "i1 MN", "i2 MN", "hi1 MN", "hi2 MN", "iii1 MN", "hiii2 MN", "iii3 MN", # iii2 is missing? are i and hi equivalent?
    "hg1 MN", "hg2 MN", "hg3 MN", "hg4 MN", "b1 MN", "b2 MN", "b3 MN", 
    # tension wing motor neurons
    "tp1 MN", "tp2 MN", "tpn MN", "ps1 MN", "ps2 MN", 
    # other wing motor neurons?
    "MNwm35", "MNwm36", "ADNM1", "ADNM2", 
    # haltere motor neurons
    "MNhm03", "MNhm42", "MNhm43", 
    # leg swing muscles
    "Tergopleural/Pleural promotor MN", "Sternal anterior rotator MN", 
    "Sternal adductor MN", "Tergotr. MN", "Sternotrochanter MN", "Tr extensor MN", 
    "Fe reductor MN", "Ti extensor MN", "Ta levator MN", 
    "STTMm", "TTMn", 
    # leg stance muscles
    "Pleural remotor/abductor MN", "Sternal posterior rotator MN", "Tr flexor MN", 
    "Acc. tr flexor MN", "Ti flexor MN", "Acc. ti flexor MN", "Ta depressor MN", 
    # other leg mn
    "MNfl10", "MNml10", "MNml29", "MNml76", "MNml77", "MNml78", "MNml79", "MNml82", "MNml84",
    "MNml80", "MNml81", "MNml83", "MNml85", "MNml86", "MNhl01", "MNhl02", "MNhl29", "MNhl59", 
    "MNhl60", "MNhl61", "MNhl62", "MNhl63", "MNhl64", "MNhl65", "MNhl66", "MNhl68", "MNhl67", 
    "MNhl69", "MNhl70", "MNhl71", "MNhl72", "MNhl73", "MNhl74", "MNhl75", "MNhl87", "MNhl88", 
    "ltm MN", "ltm2-femur MN", "ltm1-tibia MN",
    # abdominal motor neurons
    # Note some of these MN are endocrine too
    "MNad01", "MNad02", "MNad03", "MNad04", "MNad05", "MNad06", "MNad07", "MNad08", "MNad09", 
    "MNad10", "MNad11", "MNad12", "MNad13", "MNad14", "MNad15", "MNad16", "MNad17", "MNad18", 
    "MNad19", "MNad20", "MNad21", "MNad22", "MNad23", "MNad24", "MNad25", "MNad26", "MNad27", 
    "MNad28", "MNad29", "MNad30", "MNad31", "MNad32", "MNad33", "MNad34", "MNad35", "MNad36", 
    "MNad37", "MNad38", "MNad39", "MNad40", "MNad41", "MNad42", "MNad43", "MNad44", "MNad45", 
    "MNad46", "MNad47", "MNad48", "MNad49", "MNad50", "MNad51", "MNad52", "MNad53", "MNad54", 
    "MNad55", "MNad56", "MNad57", "MNad58", "MNad59", "MNad61", "MNad62", "MNad63", "MNad64", 
    "MNad65", "MNad66", "MNad67", "MNad68", "MNad69", 
    # unknown
    "MNxm01", "MNxm02", "MNxm03")
  
  # Reorder the rows of the matrix
  order_indices <- match(row_order, rownames(influence_matrix))
  influence_matrix <- influence_matrix[na.omit(order_indices), , drop = FALSE]
  
  return(influence_matrix = influence_matrix)

}

# display summary connectivity
get_summary_inputs <- function(input.meta){
  
  # simplify inputs
  input.meta.coarse <- input.meta %>%
    dplyr::mutate(pre_cell_type = gsub("auto:", "", pre_cell_type)) %>%
    # sum weights of all ids of the same cell type per post synaptic neuron
    dplyr::filter(post != pre) %>%
    dplyr::group_by(pre_cell_type, post) %>%
    dplyr::mutate(pre_cell_type_count = n()) %>%
    dplyr::mutate(pre_norm = sum(pre_norm)) %>%
    dplyr::ungroup() %>%
    # select only unique values of pre cell type per post ids
    dplyr::distinct(pre_cell_type, post, .keep_all = TRUE) %>%
    # average weights of same pre cell type across neurons of post_cell_type
    dplyr::group_by(pre_cell_type, post_cell_type) %>%
    dplyr::mutate(post_cell_count = n()) %>%
    dplyr::mutate(pre_norm = sum(pre_norm)/post_cell_count) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pre_cell_type, post_cell_type, .keep_all = TRUE) %>%
    arrange(desc(pre_norm))
  
  for (celltype2use in unique(input.meta.coarse$post_cell_type)) {
    
    print("*********************************")
    print(paste0("Running cell type: ", celltype2use))
    print("*********************************")
    
    print(paste0(
      input.meta.coarse$pre_cell_type[grepl(celltype2use, input.meta.coarse$post_cell_type)], 
      "  ", 
      round(input.meta.coarse$pre_norm[grepl(celltype2use, input.meta.coarse$post_cell_type)], digits = 3) * 100, 
      "%  in = ",
      input.meta.coarse$pre_cell_type_count[grepl(celltype2use, input.meta.coarse$post_cell_type)],
      " on = ",
      input.meta.coarse$post_cell_count[grepl(celltype2use, input.meta.coarse$post_cell_type)]
    ))
    
    print("*********************************")
    
  }

  return(input.meta.coarse = input.meta.coarse)
  
}

get_summary_outputs <- function(output.meta){
  
  # simplify inputs
  output.meta.coarse <- output.meta %>%
    dplyr::mutate(post_cell_type = gsub("auto:", "", post_cell_type)) %>%
    # sum weights of all ids of the same cell type per post synaptic neuron
    dplyr::filter(pre != post) %>%
    dplyr::group_by(post_cell_type, pre) %>%
    dplyr::mutate(post_cell_type_count = n()) %>%
    dplyr::mutate(post_norm = sum(post_norm)) %>%
    dplyr::ungroup() %>%
    # select only unique values of post cell type per pre ids
    dplyr::distinct(post_cell_type, pre, .keep_all = TRUE) %>%
    # average weights of same post cell type across neurons of pre_cell_type
    dplyr::group_by(post_cell_type, pre_cell_type) %>%
    dplyr::mutate(pre_cell_count = n()) %>%
    dplyr::mutate(post_norm = sum(post_norm)/pre_cell_count) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(post_cell_type, pre_cell_type, .keep_all = TRUE) %>%
    arrange(desc(post_norm))
  
  for (celltype2use in unique(output.meta.coarse$pre_cell_type)) {
    
    print("*********************************")
    print(paste0("Running cell type: ", celltype2use))
    print("*********************************")
    
    print(paste0(
      output.meta.coarse$post_cell_type[grepl(celltype2use, output.meta.coarse$pre_cell_type)], 
      "  ", 
      round(output.meta.coarse$post_norm[grepl(celltype2use, output.meta.coarse$pre_cell_type)], digits = 3) * 100, 
      "%  in = ",
      output.meta.coarse$post_cell_type_count[grepl(celltype2use, output.meta.coarse$pre_cell_type)],
      " on = ",
      output.meta.coarse$pre_cell_count[grepl(celltype2use, output.meta.coarse$pre_cell_type)]
    ))
    
    print("*********************************")
    
  }
  
  return(output.meta.coarse = output.meta.coarse)
  
}