### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

######################
### load libraries ###
######################

# Load packages
library(bancr)
library(tidyverse)
library(nat.nblast)

################
### get data ###
################

# Get meta data
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
  dplyr::filter(region%in%c("neck_connective"),
                !is.na(side), 
                super_class=="ascending")
choose_banc()
# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)
ids <- unique(meta$root_id)
ids <- banc_latestid(ids)
ids.left <- unique(subset(meta, side=="left")$root_id)
ids <- na.omit(ids)
ids <- ids[ids!="0"]
l2dps <- banc_read_l2dp(ids)
#l2dps[ids.left] <- banc_mirror(l2dps[ids.left])

# Get skeletons for visualisation
l2 <- banc_read_l2skel(ids)
l2[ids.left] <- banc_mirror(l2[ids.left])

# Re-root to soma where this is known
bancr:::banc.roots <- banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots, estimate = TRUE)

# Save the state of this R session, to avoid doing the step above next time you load this file
save.image()

###############
### NBLAST  ###
###############

# Run NBLAST
nb.all.full <- nat.nblast::nblast_allbyall(l2dps, normalisation = "mean")

#separate out right hand
nblast.right <- nb.all.full[!rownames(nb.all.full) %in% ids.left, 
                                   !colnames(nb.all.full) %in% ids.left]

# Print the dimensions of the new matrix to verify
print(dim(nb.all.full))
save.image()
# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height=sqrt(hc$height)

# cluster ---------------------------------------------------------


# Color clusters
cut.height <- 2.2 # needs optimisation
hc.col <- dendroextras::color_clusters(hc, h = cut.height)
hc.col <- dendextend::set_labels(hc.col,  paste0(meta[match(labels(hc.col),meta$root_id),"cell_type"],"_",meta[match(labels(hc.col),meta$root_id),"side"]))

# Plot
plot(hc.col, labels = T)
abline(h=cut.height)

# See neurons
clear3d()
banc_view()
plot3d(hc, db = l2, h = cut.height, soma = 5000)

# Get groups
groups <- cutree(hc, h = cut.height)
table(groups)
i=1
for(cluster in unique(groups)){
  
  # Get members 
  members <- names(groups[groups==cluster])
  neurons <- l2[members]
  
  # add to an_type
  bc$an_type[bc$root_id %in% members] <- as.character(cluster)
  
  # Plot neurons
  clear3d()
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(neurons, lwd=0.5, soma = 2500, col="black")
  
  # Use sprintf to create a filename with the current value of i
  filename <- sprintf('group_%d.png', i)
  rgl.snapshot(filename)
  i=i+1
}


# Upload to Seatable! -----------------------------------------------------
# 
#subset bc to just ANs before uploading to seatable
an_types <- subset(bc, super_class=="ascending")
an_types <- as.data.frame(an_types)
an_types[is.na(an_types)] = ""
ids <- data.frame(`_id`=an_types$`_id`)
banctable_update_rows(df = an_types[,c("_id", "an_type")], table = "banc_meta", base = "banc_meta")

# Update - FROM ALEX
bc.update <- bc %>%
  dplyr::filter(root_id != "0") %>%
  dplyr::mutate(across(where(is.character), ~replace_na(., "")),
                across(where(is.factor), ~as.factor(replace_na(as.character(.), ""))),
                across(where(is.numeric), ~replace_na(., 0)),
                across(where(is.logical), ~.))
bc.update.present <- subset(bc.update, bc.update$`_id`!="")
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update.present[,c("_id", "an_type")]), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

#Alex Stop here

# NOTE - GROUP NUMBERS DO NOT MATCH BC$AN_TYPE DUE TO DUPLICATE NEURONS
###########################
### Manual Cell typing  ###
###########################

choices <- list()
for(cluster in unique(groups)){
  
  # Get members 
  members <- names(groups[groups==cluster])
  neurons <- l2[members]
  
  # add to an_type
  bc$an_type[bc$root_id %in% members] <- cluster
  
  # Plot new sub-dendrogram
  nb.all.sub <- nb.all.full[members,members]
  hc.sub <- nhclust(scoremat = nb.all.sub, method = "ward.D2")
  hc.sub$height <- sqrt(hc.sub$height)
  cut.height.sub <- 1.1 
  hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
  plot(hc.sub.col, labels = T)
  abline(h=cut.height.sub)
  
  # Plot neurons
  clear3d()
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(hc.sub, db=neurons, lwd=0.5, soma = 2500, h=cut.height.sub)
  
  # Get choices
  db <- nlscan(neurons, col = "black", lwd = 3, soma = 5000)
  choices[[cluster]] <- db
  
  # Pause for user
  message(paste(choices,collapse=", "))
  readline("prompt: press any key to continue to the next cluster ")
  
}

 


# Record manually entries here