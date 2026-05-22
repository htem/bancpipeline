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

################
### get data ###
################

# Get meta data
bc <- banctable_query()
meta <- bc %>%
  # dplyr::mutate(cell_type = dplyr::case_when(
  #   is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
  #   is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
  #   TRUE ~ cell_type
  # ),
  # top_nt = gsub("auto\\:","",top_nt),
  # cell_sub_class = gsub("auto\\:","",cell_sub_class),
  # cell_class = gsub("auto\\:","",cell_class),
  # super_class = gsub("auto\\:","",super_class)) %>%
  dplyr::filter(#region%in%c("neck_connective"),
                #grepl("^SA|^Sx",cell_type)|grepl("sensory|afferent",super_class)|grepl("sensory|afferent",cell_class),
                is.na(root_position_nm)&(is.na(nucleus_id)|nucleus_id=="0")
                )

# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)
ids <- unique(meta$root_id)
ids <- ids[!is.na(ids)]
l2 <- banc_read_l2skel(ids)

# Re-root to soma where this is known
# banc.roots <- bancr:::banc_roots()
# l2 <- banc_reroot(l2, roots = banc.roots, estimate = TRUE)

#############################
### NBAST to find nerves  ###
#############################

# Det dotprops
#l2dps <- banc_read_l2dp(ids)
l2dps <- dotprops(l2/1000, OmitFailures = TRUE)
# Run NBLAST
nb.all.full <- nat.nblast::nblast_allbyall(l2dps, normalisation = "mean")
# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2", maxneurons = Inf)
# Flatten
# hc$height=sqrt(hc$height)
# Color clusters
cut.height <- 5 # needs optimisation
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

#############################
### manually choose roots ###
#############################

choices <- list()
nopen3d()
for(cluster in unique(groups)){
  
  # Get members 
  members <- names(groups[groups==cluster])
  neurons <- l2[members]
  
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
  plot3d(hc.sub, db=neurons, lwd=0.1, soma = 2500, h=cut.height.sub)
  
  # Get choices
  neurons <- add_field_seq(neurons,names(neurons),field="id")
  points <- nlapply(neurons, function(x){
    df <- as.data.frame(nat::xyzmatrix(x$d))[nat::endpoints(x),]
    df$root_id <- x$id
    df
  })
  points <- do.call(rbind,points)
  chosen.points <- select_points(points = xyzmatrix(points))
  spheres3d(chosen.points,col="purple",radius=5000)
  
  # Add
  chosen <- as.data.frame(chosen.points)
  chosen <- dplyr::left_join(chosen,points,by=c("X","Y","Z")) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  choices[[cluster]] <- chosen
  print(knitr::kable(chosen))
  spheres3d(xyzmatrix(chosen),col="black",radius=5000)
  
  # Pause for user
  readline("prompt: press any key to continue to the next cluster ")
}
choices.points <- do.call(rbind,choices)

############################
### check that IDs match ###
############################

choices.points$root_position_nm <- apply(xyzmatrix(choices.points),1,hemibrainr:::paste_coords)
choices.points$root_position_nm  <- gsub("\\(|\\)","",choices.points$root_position_nm)
choices.points$root_position <- apply(banc_nm2raw(choices.points$root_position_nm),1,hemibrainr:::paste_coords)
choices.points$root_position  <- gsub("\\(|\\)","",choices.points$root_position)

###############################
### manually calculate side ###
###############################

lrdiffs <- bancr:::banc_lr_position(xyzmatrix(choices.points),units = "nm")
sides <- ifelse(lrdiffs>0,"right","left")
choices.points$side <- sides

#########################
### update banc table ###
#########################
update.df <- left_join(choices.points, meta[,c("_id","root_id")], by = "root_id")
update.df <- subset(update.df, !is.na(side))
update.df[is.na(update.df)] <- ""
update.df <- subset(update.df, update.df$`_id`!="")
banctable_update_rows(base = 'banc_meta', 
                      table = "banc_meta", 
                      df = update.df, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

###################################
### change root position column ###
###################################
bc.roots <- banctable_query("SELECT _id, root_id, root_position_nm, root_position from banc_meta")
missed <- is.na(bc.roots$root_position_nm)&!is.na(bc.roots$root_position)
if(any(missed)){
  bc.roots$root_position_nm[missed] <- apply(banc_nm2raw(bc.roots$root_position[missed]),1,hemibrainr:::paste_coords)
}
bc.roots$root_position <- apply(banc_nm2raw(bc.roots$root_position_nm),1,hemibrainr:::paste_coords)
bc.roots$root_position[bc.roots$root_position=="(NA, NA, NA)"] <- ''
bc.roots$root_position_nm[bc.roots$root_position_nm=="(NA, NA, NA)"] <- ''
bc.roots$root_position_nm  <- gsub("\\(|\\)","",bc.roots$root_position_nm)
bc.roots$root_position  <- gsub("\\(|\\)","",bc.roots$root_position)
banctable_update_rows(base = 'banc_meta', 
                      table = "banc_meta", 
                      df = bc.roots, 
                      append_allowed = FALSE, 
                      chunksize = 1000)