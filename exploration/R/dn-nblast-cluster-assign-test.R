### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Double check that all DNs have assigned type

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)

#########################
### define functions ####
#########################

plot_neurons_percluster <- function(dn_type, unique_ids, l2, figure.dir) {
  
  # get groups/clusters
  groups <- as.numeric(dn_type)
  na_n <- sum(is.na(groups))
  
  # remove na
  unique_ids <- unique_ids[!is.na(groups)]
  groups <- groups[!is.na(groups)]
  print(str_c('removing ', na_n, ' NA'))
  
  unique_groups <- sort(unique(groups))
  
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
  
  for(cluster in unique_groups){
    
    print(paste("running cluster #", cluster))
    
    # Get members 
    selids <- unique_ids[groups == cluster]
    neurons <- l2[selids]
    
    # Plot neurons
    clear3d()
    open3d(windowRect = c(20, 30, 1000, 1000))
    banc_view()
    plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
    plot3d(neurons, lwd = 0.5, soma = 2500)
    
    output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "_test_skt.png")
    snapshot3d(output_file)
    close3d()
    
  }
  
}

test_skeletons <- function(ids2use, l2) {
  
  # get unique ids
  unique_ids <- unique(ids2use)
  
  ids2fix <- character()
  
  for(id in unique_ids){
    
    # Attempt to find and plot skeleton
    successful <- tryCatch({
      if (id %in% names(l2)) {
        plot3d(l2[id])
      } else {
        ids2fix <- c(ids2fix, id)
        stop("ID not found in data list.")
      }
      TRUE  # Indicate success
    }, error = function(e) {
      message("An error occurred while processing ID ", id, ": ", conditionMessage(e))
      FALSE  # Indicate failure
    })
    
  }
  
  return(list(ids2fix = ids2fix))
  
}

#################
### load data ###
#################

# get meta data
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
  dplyr::filter(region %in% c("neck_connective"),
                !is.na(side), 
                super_class == "descending")

# define data directory
datadir <- str_c(getwd(), "/data/dn/")

# load relevant variables
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
names(l2) <- banc_latestid(names(l2))

# set directory name
figure.dir <- str_c(getwd(), "/figures/dn_temp/dn_dist_1-9/")

# get ID names
ids <- meta$root_id
ids <- banc_latestid(ids)
unique_ids <- unique(ids)
# deal with ids that change (are deleted somehow?)
unique_ids <- intersect(unique_ids, names(l2))
unique_idx <- match(unique_ids, ids)

ids2fix <- test_skeletons(unique_ids, l2)$ids2fix

# get dn type for selected unique IDs
dn_type <- meta$dn_type
dn_type <- dn_type[unique_idx]

# get celltype
cell_type <- meta$cell_type[unique_idx]
super_class <- meta$super_class[unique_idx]
  
plot_neurons_percluster(dn_type, unique_ids, l2, figure.dir)

###########################
## plot selected neurons ##
###########################

neurons <- l2[unique_ids[(cell_type %in% "aSP22")]]
neurons <- l2[unique_ids[(cell_type %in% "AN_4_None")]]

# neurons that are not DNs
non_DNs <- !grepl("DN", cell_type)
neurons <- l2[unique_ids[non_DNs]]
cell_type[non_DNs]
super_class[non_DNs]

# Notes
# "720575941625215387" not descending (5)
# "720575941470536459" is mislabeled as AN (7)
# "720575941498142177" is mislabeled as AN (9)
# "720575941515072003" is mislabeled as AN (11)
# "720575941644313761" is mislabeled as AN (12)
# "720575941510361298" is mislabeled as AN (13)

# Plot neurons
clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(neurons[11], lwd = 0.5, soma = 2500)