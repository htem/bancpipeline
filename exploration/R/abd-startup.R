### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# startup code

#################
##### PATHS #####
#################

# Run locally
if (.Platform$OS.type == "windows") {
  
  banc.meta.save.path <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/data/abd_dn_an"
  banc.path <- "C:/Users/papers/BANC-project/"
  banc.meta.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/meta"
  banc.connectivity.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/connectivity"
  banc.influence.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/influence"
  banc.dropbox.connectivity.save.path <- "C:/Users/Diego/Dropbox (Personal)/neuroanat/connectomes"
  banc.dropbox.influence.save.path <- "C:/Users/Diego/Dropbox (Personal)/neuroanat/influence"
  
} else if (.Platform$OS.type == "unix") {

  banc.meta.save.path <- "/Users/diegopinedo/Dropbox (Personal)/LabScripts/ImportedTools/bancpipeline/data/abd_dn_an"
  banc.path <- "/Users/papers/BANC-project/"
  banc.meta.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/meta"
  banc.connectivity.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/connectivity"
  banc.influence.save.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/influence"
  banc.dropbox.connectivity.save.path <- "C:/Users/Diego/Dropbox (Personal)/neuroanat/connectomes"
  banc.dropbox.influence.save.path <- "C:/Users/Diego/Dropbox (Personal)/neuroanat/influence"
  
}

setwd(banc.path)

# for reference check:
# usethis::edit_file(paste0(banc.path, "R/startup/banc-startup.R"))
# usethis::edit_file(paste0(banc.path, "R/startup/banc-functions.R"))

source("R/startup/banc-functions.R")
# it creates functions:
#   cosine_similarity_sparse
#   .merge_hclust
#   hclust_semisupervised
#   adjust_color_brightness

###################
##### OPTIONS #####
###################
options(pillar.sigfig=15)

#####################
##### LIBRARIES #####
#####################

# load required libraires
options(scipen = 999)
library(bancr)

# SQL column types
inf.col.types <- readr::cols(
  id = readr::col_character(),
  is_seed = readr::col_logical(),
  .default = readr::col_number())

# SQL column types
banc.col.types <- readr::cols(
  .default = readr::col_character(),
  cleft_segid  = readr::col_character(),
  centroid_x = readr::col_number(),
  centroid_y = readr::col_number(),
  centroid_z = readr::col_number(),
  bbox_bx = readr::col_number(),
  bbox_by = readr::col_number(),
  bbox_bz = readr::col_number(),
  bbox_ex = readr::col_number(),
  bbox_ey = readr::col_number(),
  bbox_ez = readr::col_number(),
  presyn_segid = readr::col_character(),
  postsyn_segid  = readr::col_character(),
  presyn_x = readr::col_integer(),
  presyn_y = readr::col_integer(),
  presyn_z = readr::col_integer(),
  postsyn_x = readr::col_integer(),
  postsyn_y = readr::col_integer(),
  postsyn_z = readr::col_integer(),
  clefthash = readr::col_number(),
  partnerhash = readr::col_integer(),
  size = readr::col_number(),
  l2_root = readr::col_number(),
  l2_nodes = readr::col_number(),
  l2_segments = readr::col_number(),
  l2_branchpoints = readr::col_number(),
  l2_endpoints = readr::col_number(),
  l2_cable_length = readr::col_number(),
  l2_n_trees = readr::col_number(),
  nb = readr::col_number(),
  score = readr::col_number(),
  hemibrain_nblast = readr::col_number(),
  fafb_nblast = readr::col_number(),
  manc_nblast = readr::col_number(),
  X = readr::col_number(),
  Y = readr::col_number(),
  Z = readr::col_number(),
  x = readr::col_number(),
  y = readr::col_number(),
  z = readr::col_number(),
  dcv_density = readr::col_number(),
  dcv_count = readr::col_integer(),
  acetylcholine= readr::col_number() , 
  glutamate= readr::col_number() , 
  gaba= readr::col_number() , 
  glycine= readr::col_number() , 
  dopamine= readr::col_number() , 
  serotonin= readr::col_number() , 
  octopamine= readr::col_number() , 
  tyramine= readr::col_number() , 
  histamine= readr::col_number() , 
  nitric_oxide= readr::col_number() ,
  `allatostatin-a`= readr::col_number() , 
  `allatostatin-c`= readr::col_number() , 
  amnesiac= readr::col_number() , 
  bursicon= readr::col_number() , 
  capability= readr::col_number() , 
  ccap= readr::col_number() , 
  ccha1= readr::col_number() , 
  cnma= readr::col_number() , 
  corazonin= readr::col_number() , 
  darc1= readr::col_number() , 
  dh31= readr::col_number() , 
  dh331= readr::col_number() , 
  dh44= readr::col_number() , 
  dilp2= readr::col_number() , 
  dilp3= readr::col_number() , 
  dilp5= readr::col_number() , 
  dnpf= readr::col_number() , 
  drosulfakinin= readr::col_number() , 
  eclosion_hormone= readr::col_number() , 
  fmrf= readr::col_number() , 
  fmrfa= readr::col_number() , 
  hugin= readr::col_number() ,
  itp= readr::col_number() , 
  leucokinin= readr::col_number() , 
  mip= readr::col_number() , 
  myosuppressin= readr::col_number() , 
  myosupressin= readr::col_number() , 
  natalisin= readr::col_number() , 
  negative= readr::col_number() , 
  neuropeptide= readr::col_number() , 
  neuropeptides= readr::col_number() , 
  npf= readr::col_number() , 
  nplp1= readr::col_number() , 
  orcokinin= readr::col_number() , 
  pdf= readr::col_number() , 
  proctolin= readr::col_number() , 
  sifamide= readr::col_number() , 
  snpf= readr::col_number() , 
  space_blanket= readr::col_number() , 
  tachykinin= readr::col_number() , 
  trissin= readr::col_number(),
  input_connections = readr::col_number(),
  output_connections = readr::col_number(),
  total_connections = readr::col_number(),
  UMAP1 = readr::col_number(),
  UMAP2 = readr::col_number(),
  mitochondria = readr::col_number() , 
  mitochondria_volume = readr::col_number()
)

###############
### Colours ###
###############

paper.cols.df <- read.csv("settings/paper_colours.csv")
paper.cols <- paper.cols.df$hex
names(paper.cols) <- paper.cols.df$label
class.order <- unique(names(paper.cols))

# Create a function to generate n colors
cerise_limon_base <- c("#EE5B32", "#F6B83C", "#4BA747", "#5BB6E4", "#7C378A")
cerise_limon_palette <- grDevices::colorRampPalette(cerise_limon_base)
scale_color_cerise_limon <- function(n) {
  n <- sort(na.omit(n))
  values <- cerise_limon_palette(length(n))
  names(values) <- n
  scale_color_manual(values = values)
}
scale_fill_cerise_limon <- function(n) {
  n <- sort(na.omit(n))
  values <- cerise_limon_palette(length(n))
  names(values) <- n
  scale_fill_manual(values = values)
}

#####################
### COLOUR SCALES ###
#####################

# Colour breaks
n_breaks <- 200

# Logarithmic color scale
connection_heatmap_breaks <- seq(0, 
                                 500, 
                                 length.out = 101)
connection_heatmap_palette <- colorRampPalette(rev(c("grey20", "grey35", "grey30", "grey40", "grey60", "white")))(n_breaks-1)

# Create cosine color palette
cosine_heatmap_breaks <- seq(0, 1, length.out = n_breaks)
cosine_heatmap_palette <- colorRampPalette(c("#007BC3", "white",  "#D72000"))(n_breaks - 1)

# # Logarithmic color scale
# influence_heatmap_breaks <- expm1(seq(min(influence.meta$influence_log,na.rm = TRUE), 
#                                       max(influence.meta$influence_log,na.rm = TRUE), 
#                                       length.out = 101))
# influence_heatmap_palette <- colorRampPalette(rev(c("grey30", "grey50", "grey60", "white")))(n_breaks-1)

# Colours for cascade
cascade_heatmap_breaks <- seq(0, 
                              10, 
                              length.out = n_breaks)
cascade_heatmap_palette <- colorRampPalette(c("white", "black"))(n_breaks-1)

# Create scaled color palette
scaled_heatmap_breaks <- seq(-5, 5, length.out = n_breaks)
scaled_heatmap_palette <- colorRampPalette(c("#0054c3","#007BC3", "#58b6ed","grey90","#f08c7a","#D72000","#D72000"))(n_breaks - 1)

# Create scaled color palette 2
scaled_heatmap_breaks2 <- seq(-10, 10, length.out = n_breaks)
scaled_heatmap_palette2 <- colorRampPalette(c("#0054c3","#007BC3", "#58b6ed","grey90","#f08c7a","#D72000","#D72000"))(n_breaks - 1)

# Create scaled color palette 3
scaled_heatmap_breaks3 <- seq(-1, 500, length.out = n_breaks)
scaled_heatmap_palette3 <- colorRampPalette(c("grey80","#f08c7a","#D72000","#D72000"))(n_breaks - 1)

# Create scaled color palette 4
scaled_heatmap_breaks4 <- seq(-20, 20, length.out = n_breaks)
scaled_heatmap_palette4 <- colorRampPalette(c("#0054c3","#007BC3", "#58b6ed","grey90","#f08c7a","#D72000","#D72000"))(n_breaks - 1)

# Create scaled color palette 5
scaled_heatmap_breaks5 <- seq(-1, 100, length.out = n_breaks)
scaled_heatmap_palette5 <- colorRampPalette(c("grey80","#f08c7a","#D72000","#D72000"))(n_breaks - 1)

##############################################
##### Functional annotations for neurons #####
##############################################

cns.functions <- bancr::banctable_query(sql = "SELECT * FROM functions") %>%
  dplyr::select(-starts_with("_")) %>%
  dplyr::filter(!is.na(cell_type))

# Read DN data output
umap.dn.df <- read_csv(file = "data/banc_dn_functional_classes.csv", col_types = banc.col.types) %>%
  dplyr::mutate(cluster = gsub("SD_|ED_","DN_",cluster),
                cluster = gsub("SA_|EA_","AN_",cluster)) %>%
  dplyr::mutate(clusterno = gsub(".*_","",cluster))
umap.an.df <- read_csv(file = "data/banc_an_functional_classes.csv", col_types = banc.col.types) %>%
  dplyr::mutate(cluster = gsub("SD_|ED_","DN_",cluster),
                cluster = gsub("SA_|EA_","AN_",cluster)) %>%
  dplyr::mutate(clusterno = gsub(".*_","",cluster))
classes.dn.df <- read_csv(file = "data/banc_dn_functional_classes_by_neuron.csv", col_types = banc.col.types) %>%
  dplyr::mutate(cluster = gsub("SD_|ED_","DN_",cluster),
                cluster = gsub("SA_|EA_","AN_",cluster)) %>%
  dplyr::mutate(clusterno = gsub(".*_","",cluster))
classes.an.df <- read_csv(file = "data/banc_an_functional_classes_by_neuron.csv", col_types = banc.col.types) %>%
  dplyr::mutate(cluster = gsub("SD_|ED_","DN_",cluster),
                cluster = gsub("SA_|EA_","AN_",cluster)) %>%
  dplyr::mutate(clusterno = gsub(".*_","",cluster))

#########################
##### ANNOUNCEMENTS #####
#########################

# Messages for debugging
message("R_MAX_VSIZE: ", Sys.getenv("R_MAX_VSIZE"))
message(".libPaths: ", print(.libPaths()))
message("bancr: ", packageVersion("bancr"))
print("##### SESSION INFO #####")
print(sessionInfo())
print("##### SESSION INFO #####")