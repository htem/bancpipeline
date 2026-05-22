### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Load banctable and save DN metadata, DN skeletons & dotprops and all-to-all nblast scores

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)

choose_banc()

#########################
### define functions ####
#########################

assign_nblast_cluster <- function(hc, cut.height, bc) {
  
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
  
  for(cluster in unique(groups)){
    
    # get ordered cluster number (from left to right)
    cluster.ordered <- unique_groups[cluster]
    
    # Get members 
    members <- names(groups[groups == cluster.ordered])

    # update IDs names
    members <- banc_latestid(members)
    
    # assign auto clusters
    bc$auto_dn_type[bc$root_id %in% members] <- as.character(cluster)
    bc$dn_type[bc$root_id %in% members] <- as.character(cluster)
    
  }
  
  #subset bc to just DNs before uploading to seatable
  dn_types <- subset(bc, super_class == "descending")
  dn_types <- as.data.frame(dn_types)
  dn_types$auto_dn_type[is.na(dn_types$auto_dn_type)] = ""
  banctable_update_rows(df = dn_types[,c("_id", "root_id", "auto_dn_type", "dn_type")], 
                        table = "banc_meta", base = "banc_meta")
  
}

###############################
### get data from BANCTABLE ###
###############################

# Get meta data
# make sure you run: banctable_set_token(user = "", pwd = "") to generate banc token
bc <- banctable_query()

# update ID names
bc$root_id <- banc_latestid(bc$root_id)

# define data directory
datadir <- str_c(getwd(), "/data/dn/")

# Load relevant variables
nb.all.full <- readRDS(str_c(datadir, "nblast_score_all_dn.rds"))

# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height <- sqrt(hc$height)

cut.height <- 1.9
assign_nblast_cluster(hc, cut.height, bc)