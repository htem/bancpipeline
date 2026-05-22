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
library(progress)
library(tidyverse)
library(ggplot2)

# Make sure all functions query BANC and not FAFB
# choose_banc()

#########################
### define functions ####
#########################

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

###############################
### get data from BANCTABLE ###
###############################

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
meta$root_id <- banc_latestid(meta$root_id)

# get all ids
ids <- unique(meta$root_id)
ids <- banc_latestid(ids)
ids.left <- unique(subset(meta, side == "left")$root_id)
ids.left <- banc_latestid(ids.left)

# get skeletons for visualization (units in nanometers)
l2 <- banc_read_l2skel(ids)

# optional: test that you can read all skeletons
ids2fix <- test_skeletons(ids, l2)$ids2fix

# mirror neurons to right hemisphere (mirror require units in nanometers)
l2[ids.left] <- banc_mirror(l2[ids.left])

# make dotprops for NBLASTing (inaccurate but fast to acquire) (units in microns)
l2dps <- dotprops(l2/1000)

# re-root to soma where this is known
banc.roots <- bancr:::banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots, estimate = TRUE)
rm(banc.roots)

# define directory where to save data
datadir <- str_c(getwd(), "/data/dn/")
if (!dir.exists(datadir)) {
  dir.create(datadir, recursive = TRUE)
}

# save relevant variables for clustering and plotting
saveRDS(meta, file = str_c(datadir, "meta_dn.rds"))
saveRDS(l2, file = str_c(datadir, "l2_dn.rds"))
saveRDS(l2dps, file = str_c(datadir, "l2dps_dn.rds"))

##############################
### generate NBLAST scores ###
##############################

# run NBLAST
nb.all.full <- nat.nblast::nblast_allbyall(l2dps, normalisation = "mean")

# save relevant variables for clustering and plotting
saveRDS(nb.all.full, file = str_c(datadir, "nblast_score_all_dn.rds"))

rm(datadir)