### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Script to look at top downstream targets of JO-F to find potential aDNs

######################
### load libraries ###
######################

library(fafbseg)
library(writexl)
library(natverse)
library(tidyverse)
library(bancr)
library(tidyverse)
library(nat.nblast)

# fetch canned connectivity *and* cell type data 
choose_segmentation('flywire31')
download_flywire_release_data('all')

# to check version run:
flywire_connectome_data_version()

# set search parameters
cell_meta = flytable_meta()
JO_F = subset(cell_meta, cell_type == "JO-F")
min_thresh = 5

JO_F.out <- flywire_partner_summary2(JO_F, partners = 'out', threshold = min_thresh, summarise = T)
JO_F.outDN <- JO_F.out %>%
  dplyr::filter(grepl("DN", JO_F.out$cell_type))

# load banc skeletons and metadata
# define data directory
datadir <- str_c(getwd(), "/data/dn/")

# load relevant variables
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
banc_meta <- readRDS(str_c(datadir, "meta_dn.rds"))

ids <- banc_meta$root_id
unique_ids <- unique(ids)
unique_idx <- match(unique_ids, ids)
unique_dn_names <- banc_meta$cell_type[unique_idx]
skt <- l2[unique_ids]

figure.dir <- str_c(getwd(), "/figures/dn_temp/JO_F_outputs/")

plot_output_dns(JO_F.outDN$cell_type, JO_F.outDN$weight, skt, unique_dn_names, figure.dir)

plot_output_dns <- function(input_celltype, input_weight, skt, skt_celltype, figure.dir) {

  unique_ct <- unique(input_celltype)
  unique_idx_ct <- match(unique_ct, input_celltype)
  unique_weight <- input_weight[unique_idx_ct]

  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
    
  k <- 1
  
  for(cell2run in unique_ct){
    
    print(paste("running celltype #", cell2run))
    
    # Get members 
    neurons <- l2[skt_celltype %in% cell2run]
    
    # Plot neurons
    clear3d()
    open3d(windowRect = c(20, 30, 1000, 1000))
    banc_view()
    plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
    plot3d(neurons, lwd = 0.5, soma = 2500)
    
    output_file <- str_c(figure.dir, "output_", as.character(k), 
                         "_weight_", as.character(unique_weight[unique_ct %in% cell2run])
                         , "_", cell2run, ".png")
    snapshot3d(output_file)
    close3d()
    
    k <- k + 1
    
  }
  
}

# Notes: based on morphology DNge011 and DNge012 look like aDN1, not clear how aDN2 looks like
#   (those are two grooming DNs mentioned in hampel et al 2015)
#   see https://ngl.flywire.ai/?local_id=11ff56af1243058c4d6a839eb9a4e987
#   this looks like aBN1 https://codex.flywire.ai/app/cell_details?data_version=783&root_id=720575940630907434