### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Script to explore oviDNs cluster

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

display_ordered_inputs_or_outputs_per_cell_type <- function(neuron.meta, input.meta, synapse.type, dn_type) {
  
  # simplify metadata of input neurons
  if (synapse.type == "pre") {
    
    input.meta <- input.meta %>% 
      dplyr::select(pre, pre_cell_type, pre_top_nt, pre_side, norm, post, post_cell_type) %>%
      dplyr::rename(rootid = pre) %>%
      dplyr::rename(cell_type = pre_cell_type) %>%
      dplyr::rename(top_nt = pre_top_nt) %>%
      dplyr::rename(side = pre_side) %>%
      dplyr::rename(weight = norm) %>%
      dplyr::rename(dnid = post) %>%
      dplyr::rename(dn_cell_type = post_cell_type)

  } else if (synapse.type == "post") {
    
    input.meta <- input.meta %>% 
      dplyr::select(post, post_cell_type, post_top_nt, post_side, pre_norm, pre, pre_cell_type) %>%
      dplyr::rename(rootid = post) %>%
      dplyr::rename(cell_type = post_cell_type) %>%
      dplyr::rename(top_nt = post_top_nt) %>%
      dplyr::rename(side = post_side) %>%
      dplyr::rename(weight = pre_norm) %>%
      dplyr::rename(dnid = pre) %>%
      dplyr::rename(dn_cell_type = pre_cell_type)

  }
  
  # get unique cell types
  cell_type <- neuron.meta$cell_type

  # sort inputs per cell type
  for(current_ct in unique(cell_type)){
    
    # get subset
    weta_temp <- input.meta %>%
      dplyr::filter(dn_cell_type %in% current_ct)
    max_weight <- max(weta_temp$weight)
    
    if (length(dn_type) > 0 && nchar(dn_type[1]) > 0) {
      weta_temp <- weta_temp %>% 
        dplyr::filter(
          Reduce(`|`, lapply(dn_type, function(pattern) grepl(pattern, cell_type)))
        )
    }
    
    # sort in ascending order and generate a count of downstream neurons connected to (of the same cell type)
    #   and a mean norm input value (across downstream targets)
    weta_temp <- weta_temp %>%
      arrange(desc(weight)) %>%
      dplyr::group_by(rootid) %>%
      dplyr::mutate(n_count = length(weight),
                    weight_mean = sum(weight)/n_count) %>%
      dplyr::ungroup()
    
    weta_temp <- weta_temp %>%
      dplyr::distinct(rootid, .keep_all = TRUE) %>%
      arrange(desc(weight_mean))  %>%
      mutate(weight_mean = round(weight_mean, digits = 3)) %>%
      mutate(weight_mean = weight_mean * 100)
    max_weight <- round(max_weight, digits = 3) * 100

    print(paste0("Running cell type : ", current_ct, " max weight = ", max_weight))
    print(weta_temp[, c("weight_mean", "n_count", "cell_type", "side")], n = 40)
    
  }
  
}

plot_dendrogram_v1 <- function(nb_scores, cut.height, skt, input.meta, synapse.type, dn_type, figure.dir, figure.name) {
  
  # simplify metadata of input neurons
  if (synapse.type == "pre") {
    
    input.meta <- input.meta %>% 
      dplyr::select(pre, pre_cell_type, pre_top_nt, pre_side, norm, post_cell_type) %>%
      dplyr::distinct(pre, .keep_all = TRUE) %>%
      dplyr::rename(rootid = pre) %>%
      dplyr::rename(cell_type = pre_cell_type) %>%
      dplyr::rename(top_nt = pre_top_nt) %>%
      dplyr::rename(side = pre_side) %>%
      dplyr::rename(weight = norm) %>%
      dplyr::rename(dn_cell_type = post_cell_type)
    
    figure.dir.sub <- str_c(figure.dir, "inputs/")
    
  } else if (synapse.type == "post") {
    
    input.meta <- input.meta %>% 
      dplyr::select(post, post_cell_type, post_top_nt, post_side, pre_norm, pre_cell_type) %>%
      dplyr::distinct(post, .keep_all = TRUE) %>%
      dplyr::rename(rootid = post) %>%
      dplyr::rename(cell_type = post_cell_type) %>%
      dplyr::rename(top_nt = post_top_nt) %>%
      dplyr::rename(side = post_side) %>%
      dplyr::rename(weight = pre_norm) %>%
      dplyr::rename(dn_cell_type = pre_cell_type)

    figure.dir.sub <- str_c(figure.dir, "outputs/")
    
  }
 
  if (length(dn_type) > 0 && nchar(dn_type[1]) > 0) {
    input.meta <- input.meta %>% 
      dplyr::filter(
        Reduce(`|`, lapply(dn_type, function(pattern) grepl(pattern, cell_type)))
      )
    
    nb_scores <- nb_scores[input.meta$rootid, input.meta$rootid]
    
  }
   
  # create target directory  
  if (!dir.exists(figure.dir.sub)) {
    dir.create(figure.dir.sub, recursive = TRUE)
  }  
  
  # change variable name
  # confirm count variable to use instead of norm
  
  # simplify nt names
  input.meta$top_nt[input.meta$top_nt == "acetylcholine"] <- "ACh"
  input.meta$top_nt[input.meta$top_nt == "glutamate"] <- "Glut"
  input.meta$top_nt[input.meta$top_nt == "serotonin"] <- "5-HT"
  
  # create figure directory
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir, recursive = TRUE)
  }
  
  # cluster
  hc <- nhclust(scoremat = nb_scores, method = "ward.D2")
  
  # flatten
  hc$height <- sqrt(hc$height)
  
  # Get groups/clusters
  groups <- cutree(hc, h = cut.height)
  
  # generate a dendrogram
  hc.col2 <- as.dendrogram(hc)

  # plot  dendrogram
  width_in_inches <- 20
  height_in_inches <- 10
  dpi <- 100
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, figure.name, ".png")
  png(filename = output_file, width = width, height = height, res = dpi)
  
  plot(hc.col2, ylab = "Height", leaflab = "none")
  abline(h = cut.height)
  dev.off()
  
  # Reorder clusters to have the leftmost as cluster 1
  # Get the order of the labels in the dendrogram
  order_of_labels <- order.dendrogram(hc.col2)
  
  # Get a vector of group labels ordered by the dendrogram
  ordered_group_labels <- groups[order_of_labels]
  
  # Find unique groups in their left-to-right order
  unique_groups <- unique(ordered_group_labels)
  
  # define colors to use for cluster to be plotted
  other_color <- "#000000"
  target_color <- "#FF0000"
  
  for(cluster in unique(groups)){
    
    # get ordered cluster number (from left to right)
    cluster.ordered <- unique_groups[cluster]
    
    # Get members 
    members <- names(groups[groups == cluster.ordered])
    neurons <- skt[members]
    
    # plot selected cluster in full dendrogram
    temp_hc <- as.dendrogram(hc)
    temp_col <- ifelse(groups == cluster.ordered, target_color, other_color)
    temp_col <- temp_col[order.dendrogram(temp_hc)]
    temp_col <- factor(temp_col, unique(temp_col))
    temp_hc <- temp_hc %>% 
      color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
      dendextend::set("labels_colors", as.character(temp_col))
    
    width_in_inches <- 20
    height_in_inches <- 10
    dpi <- 100
    width <- width_in_inches * dpi
    height <- height_in_inches * dpi
    output_file <- str_c(figure.dir.sub, 'nblast_clus_', cluster, "_fulldend.png")
    png(filename = output_file, width = width, height = height, res = dpi)
    
    plot(temp_hc, ylab = "Height", leaflab = "none")
    abline(h = cut.height)
    dev.off()
    
    # Plot new sub-dendrogram
    nb_scores.sub <- nb_scores[members, members]
    hc.sub <- nhclust(scoremat = nb_scores.sub, method = "ward.D2")
    hc.sub$height <- sqrt(hc.sub$height)
    cut.height.sub <- 1 
    hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
    
    # replace labels
    idx <- match(labels(hc.sub.col), input.meta$rootid)
    
    new_labels <- paste0(input.meta$cell_type[idx], 
                         "_", input.meta$top_nt[idx],
                         "_", input.meta$side[idx],
                         "_", round(input.meta$weight[idx], digits = 3)*100, "%")
    
    hc.sub.col <- dendextend::set_labels(hc.sub.col, new_labels)
    
    dpi <- 300
    width <- width_in_inches * dpi
    height <- height_in_inches * dpi
    output_file <- str_c(figure.dir.sub, 'nblast_clus_', cluster, "_subdend.png")
    png(filename = output_file, width = width, height = height, res = dpi)
    
    par(mar = c(3, 1, 1, 10))
    plot(hc.sub.col, labels = T, horiz = TRUE)
    par(mar = c(5, 4, 4, 2) + 0.1)
    dev.off()
    
    # Plot neurons
    clear3d()
    open3d(windowRect = c(20, 30, 1000, 1000))
    banc_view()
    plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
    plot3d(hc.sub, db = neurons, lwd = 0.5, soma = 2500, h = cut.height.sub)
    
    output_file <- str_c(figure.dir.sub, 'nblast_clus_', cluster, "_skt.png")
    snapshot3d(output_file)
    close3d()
    
  }
  
}

# delete temporary variables
rm(nb_scores, cut.height, skt, input.meta, synapse.type, dn_type)
rm(figure.dir.sub, hc, groups, hc.col2)

plot_dendrogram_v2 <- function(nb_scores, skt, cell_type2use, input.meta, synapse.type, figure.dir) {
  
  # simplify metadata of input neurons
  if (synapse.type == "pre") {
    
    input.meta <- input.meta %>% 
      dplyr::select(pre, pre_cell_type, pre_top_nt, pre_side, norm, post, pre_type, post_cell_type) %>%
      dplyr::rename(rootid = pre) %>%
      dplyr::rename(cell_type = pre_cell_type) %>%
      dplyr::rename(top_nt = pre_top_nt) %>%
      dplyr::rename(side = pre_side) %>%
      dplyr::rename(weight = norm) %>%
      dplyr::rename(dnid = post) %>%
      dplyr::rename(dn_cell_type = post_cell_type) %>%
      dplyr::rename(conntype = pre_type)
    
  } else if (synapse.type == "post") {
    
    input.meta <- input.meta %>% 
      dplyr::select(post, post_cell_type, post_top_nt, post_side, pre_norm, pre, pre_type, pre_cell_type) %>%
      dplyr::rename(rootid = post) %>%
      dplyr::rename(cell_type = post_cell_type) %>%
      dplyr::rename(top_nt = post_top_nt) %>%
      dplyr::rename(side = post_side) %>%
      dplyr::rename(weight = pre_norm) %>%
      dplyr::rename(dnid = pre) %>%
      dplyr::rename(dn_cell_type = pre_cell_type) %>%
      dplyr::rename(conntype = pre_type)
    
  }
  
  # simplify nt names
  input.meta$top_nt[input.meta$top_nt == "acetylcholine"] <- "ACh"
  input.meta$top_nt[input.meta$top_nt == "glutamate"] <- "Glut"
  input.meta$top_nt[input.meta$top_nt == "serotonin"] <- "5-HT"
  
  # sort inputs per cell type
  for(current_ct in unique(cell_type2use)){
    
    # print cell type
    print(current_ct)
    
    # get inputs to selected cell type
    subset_input_meta <- input.meta %>%
      dplyr::filter(dn_cell_type %in% current_ct)
    
    # get clusters present in this subset
    conntype <- subset_input_meta$conntype
    
    if (synapse.type == "pre") {
      figure.dir.type <- str_c(figure.dir, "inputs_to_", current_ct, "/")
    } else if (synapse.type == "post") {
      figure.dir.type <- str_c(figure.dir, "outputs_to_", current_ct, "/")
    }
    
    # create figure directory
    if (!dir.exists(figure.dir.type)) {
      dir.create(figure.dir.type, recursive = TRUE)
    }
    
    for(current_conntype in unique(conntype)){
      
      # get neurons
      members <- subset_input_meta %>%
        dplyr::filter(conntype == current_conntype)
      
      # if just one neuron duplicate
      if (length(members$rootid) == 1) {
        
        # edit rows
        members <- rbind(members, members[1, ])
        
        # replace labels
        new_labels <- paste0(members$cell_type, 
                             "_", members$top_nt,
                             "_", members$side,
                             "_", round(members$weight, digits = 3)*100, "%",
                             "_", members$cell_type,
                             "_", members$side, "_*single*")
      } else {
        
        # replace labels
        new_labels <- paste0(members$cell_type, 
                             "_", members$top_nt,
                             "_", members$side,
                             "_", round(members$weight, digits = 3)*100, "%",
                             "_", members$cell_type,
                             "_", members$side)
        
      }
      
      # make and plot sub-dendrogram
      nb_scores.sub <- nb_scores[members$rootid, members$rootid]
      hc.sub <- nhclust(scoremat = nb_scores.sub, method = "ward.D2")
      hc.sub$height <- sqrt(hc.sub$height)
      cut.height.sub <- 1 
      hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
      hc.sub.col <- dendextend::set_labels(hc.sub.col, new_labels)
      
      width_in_inches <- 20
      height_in_inches <- 10
      dpi <- 100
      width <- width_in_inches * dpi
      height <- height_in_inches * dpi
      output_file <- str_c(figure.dir.type, 'nblast_clus_', current_conntype, "_subdend_all.png")
      png(filename = output_file, width = width, height = height, res = dpi)
      
      par(mar = c(3, 1, 1, 20))
      plot(hc.sub.col, labels = T, horiz = TRUE)
      par(mar = c(5, 4, 4, 2) + 0.1)
      dev.off()
      
      # make and plot sub-dendrogram of unique
      members <- members %>%
        dplyr::group_by(rootid) %>%
        dplyr::mutate(n_count = length(weight),
                      weight_mean = sum(weight)/n_count) %>%
        dplyr::ungroup() %>%
        mutate(weight_mean = round(weight_mean, digits = 3)) %>%
        dplyr::distinct(rootid, .keep_all = TRUE)
      
      members <- members %>% 
        mutate(weight_mean = weight_mean * 100)

      # get neuron skeletons
      neurons <- skt[members$rootid]
      
      # if just one neuron duplicate
      if (length(members$rootid) == 1) {
        
        # edit rows
        members <- rbind(members, members[1, ])
        
        # replace labels
        new_labels <- paste0(members$cell_type, 
                             "_", members$top_nt,
                             "_", members$side,
                             "_", members$weight_mean, "%",
                             "_", members$n_count, "n", "_*single*")
      } else {
        
        # replace labels
        new_labels <- paste0(members$cell_type, 
                             "_", members$pre_top_nt,
                             "_", members$pre_side,
                             "_", members$weight_mean, "%",
                             "_", members$n_count, "n")
        
      }
      
      nb_scores.sub <- nb_scores[members$rootid, members$rootid]
      hc.sub <- nhclust(scoremat = nb_scores.sub, method = "ward.D2")
      hc.sub$height <- sqrt(hc.sub$height)
      cut.height.sub <- 1 
      hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
      hc.sub.col <- dendextend::set_labels(hc.sub.col, new_labels)
      
      width_in_inches <- 20
      height_in_inches <- 10
      dpi <- 100
      width <- width_in_inches * dpi
      height <- height_in_inches * dpi
      output_file <- str_c(figure.dir.type, 'nblast_clus_', current_conntype, "_subdend.png")
      png(filename = output_file, width = width, height = height, res = dpi)
      
      par(mar = c(3, 1, 1, 20))
      plot(hc.sub.col, labels = T, horiz = TRUE)
      par(mar = c(5, 4, 4, 2) + 0.1)
      dev.off()
      
      # plot neurons
      clear3d()
      open3d(windowRect = c(20, 30, 1000, 1000))
      banc_view()
      plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
      
      if (length(names(neurons)) == 1) {
        
        plot3d(neurons, lwd = 0.5, soma = 2500)
        
      } else {
        
        plot3d(hc.sub, db = neurons, lwd = 0.5, soma = 2500, h = cut.height.sub)
        
      }
      
      output_file <- str_c(figure.dir.type, 'nblast_clus_', current_conntype, "_skt.png")
      snapshot3d(output_file)
      close3d()
      
    }
    
  }
  
}

gen_clusters_v1 <- function(nb_scores, cut.height) {
  
  # cluster
  hc <- nhclust(scoremat = nb_scores, method = "ward.D2")
  
  # flatten
  hc$height <- sqrt(hc$height)
  
  # Get groups/clusters
  groups <- cutree(hc, h = cut.height)
  
  # generate a dendrogram
  hc.col2 <- as.dendrogram(hc)
  
  # Get the order of the labels in the dendrogram
  order_of_labels <- order.dendrogram(hc.col2)
  
  # Get a vector of group labels ordered by the dendrogram
  ordered_group_labels <- groups[order_of_labels]
  
  # Find unique groups in their left-to-right order
  unique_groups <- unique(ordered_group_labels)
  
  clus.per.id <- data.frame(root_id = unique(names(groups)), 
                            input_type = "NA", 
                            stringsAsFactors = FALSE)
  
  for(cluster in unique(groups)){
    
    # get ordered cluster number (from left to right)
    cluster.ordered <- unique_groups[cluster]
    
    # Get members 
    members <- names(groups[groups == cluster.ordered])
    
    # assign auto clusters
    clus.per.id$input_type[clus.per.id$root_id %in% members] <- as.character(cluster)
    
  }
  
  return(list(clus.per.id = clus.per.id))
  
}

###################
### get abdDNs ####
###################

# pull abdDN cluster of interest
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
                super_class == "descending") %>%
  dplyr::filter(dn_type == "15") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# load skeletons
datadir <- str_c(getwd(), "/data/dn/")

# set figure directory name
figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/")

# get abdDNs ids
abdDNs.ids <- meta$root_id

######################################
### get abdDNs outputs and inputs ####
######################################

# vignette exploration
el <- banc_edgelist()

# populate edgelist with some metadata from bc and add norm input synapses
bc.pre <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(bc.pre) <- paste0("pre_", colnames(bc.pre))
bc.post <- bc %>% 
  dplyr::select(root_id, region, cell_class, cell_sub_class, cell_type, top_nt, side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
colnames(bc.post) <- paste0("post_", colnames(bc.post))
el.norm <- el %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id)) %>%
  # generate normalized input weight
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(post_count = sum(n),
                norm = n/post_count) %>%
  dplyr::ungroup() %>%
  # generate normalized output weight
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(pre_count = sum(n),
                pre_norm = n/post_count) %>%
  dplyr::ungroup() %>%  
  dplyr::left_join(bc.pre,
                   by = c("pre" = "pre_root_id")) %>%
  dplyr::left_join(bc.post,
                   by = c("post" = "post_root_id")) %>%
  dplyr::rename(count = n)

# abs synapse threshold (values to try 10, 25, 50)
abs_syn_thres <- 10
# norm synapse threshold (values to try 0.0025 0.005 0.01)
norm_syn_thres <- 0.0025

# extract inputs and outputs to abdDNs
abdDN.inputs <- el.norm %>% 
  dplyr::filter(post_pt_root_id %in% abdDNs.ids) %>% 
  dplyr::filter(norm > norm_syn_thres, count > abs_syn_thres)

# abs synapse threshold (values to try 10, 25, 50)
abs_syn_thres <- 10
abdDN.outputs <- el.norm %>% 
  dplyr::filter(pre_pt_root_id %in% abdDNs.ids) %>%
  dplyr::filter(count > abs_syn_thres)

# display table of top 40 inputs
display_ordered_inputs_or_outputs_per_cell_type(meta, abdDN.inputs, "pre", c(""))
display_ordered_inputs_or_outputs_per_cell_type(meta, abdDN.inputs, "pre", c("DN", "AN"))

# display table of top 40 outputs
display_ordered_inputs_or_outputs_per_cell_type(meta, abdDN.outputs, "post", c(""))
display_ordered_inputs_or_outputs_per_cell_type(meta, abdDN.outputs, "post", c("DN", "AN"))

# check if a particular cell_type is in there
sum(grepl("pal02", unique(abdDN.inputs$pre_cell_type), ignore.case = TRUE))

# save relevant variables for clustering and plotting
saveRDS(abdDN.inputs, file = str_c(datadir, "meta_abddn_inputs.rds"))
saveRDS(abdDN.outputs, file = str_c(datadir, "meta_abddn_outputs.rds"))

########################################################
### get Outputs and inputs skt, dp and nblast score ####
########################################################

# get skeletons and dotprops of:
#   input neurons
skts.inputs <- banc_read_l2skel(unique(abdDN.inputs$pre))
ids.left <- unique(subset(abdDN.inputs, pre_side == "left")$pre)
skts.inputs[ids.left] <- banc_mirror(skts.inputs[ids.left])
dps.inputs <- dotprops(skts.inputs/1000)

#   output neurons
skts.outputs <- banc_read_l2skel(unique(abdDN.outputs$post))
ids.left <- unique(subset(abdDN.outputs, post_side == "left")$post)
skts.outputs[ids.left] <- banc_mirror(skts.outputs[ids.left])
dps.outputs <- dotprops(skts.outputs/1000)

# nblast inputs and outputs
nb.inputs.full <- nat.nblast::nblast_allbyall(dps.inputs, normalisation = "mean")
nb.output.full <- nat.nblast::nblast_allbyall(dps.outputs, normalisation = "mean")

# define directory where to save data
if (!dir.exists(datadir)) {
  dir.create(datadir, recursive = TRUE)
}

# save relevant variables for clustering and plotting
saveRDS(skts.inputs, file = str_c(datadir, "l2_abddn_inputs.rds"))
saveRDS(skts.outputs, file = str_c(datadir, "l2_abddn_outputs.rds"))
saveRDS(dps.inputs, file = str_c(datadir, "l2dps_abddn_inputs.rds"))
saveRDS(dps.outputs, file = str_c(datadir, "l2dps_abddn_outputs.rds"))
saveRDS(nb.inputs.full, file = str_c(datadir, "nblast_score_all_abddn_inputs.rds"))
saveRDS(nb.output.full, file = str_c(datadir, "nblast_score_all_abddn_outputs.rds"))

########################################################
### Classify / visualize output/input neurons ####
########################################################

# manual fix of cell_type oviDNa_b and DNpe046
oviDNa_b.new <- c("720575941563681799", "720575941654228180")
abdDN.inputs$post_cell_type[abdDN.inputs$post %in% oviDNa_b.new] <- "oviDNa_b"
abdDN.outputs$pre_cell_type[abdDN.outputs$pre %in% oviDNa_b.new] <- "oviDNa_b"
meta$cell_type[meta$root_id %in% oviDNa_b.new] <- "oviDNa_b"

# plot dendrogram
cut.height <- 1.7

figure.name <- "abdDNs_inputs_fulldend"
plot_dendrogram_v1(nb.inputs.full, cut.height, skts.inputs, abdDN.inputs, "pre", "", figure.dir, figure.name)

figure.name <- "abdDNs_outputs_fulldend"
plot_dendrogram_v1(nb.output.full, cut.height, skts.outputs, abdDN.outputs, "post", "", figure.dir, figure.name)

# just ANs and DNs
figure.name <- "abdDNs_inputs_ANDNdend"
figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/ANDN/")
plot_dendrogram_v1(nb.inputs.full, cut.height, skts.inputs, abdDN.inputs, "pre", c("AN", "DN"), figure.dir, figure.name)

figure.name <- "abdDNs_outputs_ANDNdend"
figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/ANDN/")
plot_dendrogram_v1(nb.output.full, cut.height, skts.outputs, abdDN.outputs, "post", c("AN", "DN"), figure.dir, figure.name)

# assign types to inputs and add to abdDN.inputs
clus.per.id <- gen_clusters_v1(nb.inputs.full, cut.height)$clus.per.id
abdDN.inputs <- abdDN.inputs %>%
  dplyr::left_join(clus.per.id, by = c("pre"="root_id")) %>%
  dplyr::rename(pre_type = input_type)

clus.per.id <- gen_clusters_v1(nb.output.full, cut.height)$clus.per.id
abdDN.outputs <- abdDN.outputs %>%
  dplyr::left_join(clus.per.id, by = c("post"="root_id")) %>%
  dplyr::rename(pre_type = input_type)

# plot input types of each cell type separately
figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/allinputs/")
plot_dendrogram_v2(nb.inputs.full, skts.inputs, unique(meta$cell_type), abdDN.inputs, "pre", figure.dir)

figure.dir <- str_c(getwd(), "/figures/dn_temp/abdDNs/alloutputs/")
plot_dendrogram_v2(nb.output.full, skts.outputs, unique(meta$cell_type), abdDN.outputs, "post", figure.dir)


# plot ANs
celltype2use <- "DNpe047"
n_indeces <- unique(abdDN.inputs$pre[grepl("AN", abdDN.inputs$pre_cell_type) & grepl(celltype2use, abdDN.inputs$post_cell_type)])

temp.dir <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/figures/dn_temp/abdDNs/"
clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(skts.inputs[n_indeces], lwd = 0.5, soma = 2500)

output_file <- str_c(temp.dir, celltype2use, "_input_AN_skt.png")
snapshot3d(output_file)
close3d()




celltype2use <- "DNpe047"
n_indeces <- unique(abdDN.outputs$post[grepl("AN", abdDN.outputs$post_cell_type) & grepl(celltype2use, abdDN.outputs$pre_cell_type)])

temp.dir <- "C:/Users/Diego/Dropbox (HMS)/LabScripts/ImportedTools/bancpipeline/figures/dn_temp/abdDNs/"
clear3d()
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(skts.outputs[n_indeces], lwd = 0.5, soma = 2500)

output_file <- str_c(temp.dir, celltype2use, "_output_AN_skt.png")
snapshot3d(output_file)
close3d()


# load influence score and check which IDs are matched

is.dir <- "C:/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/Rachel_request/"
is.sel <- read.csv(str_c(is.dir, "JO_influence.csv"), colClasses = c(
  "matrix_index" = "numeric",
  "id" = "character",
  "JO_influence_norm_unsigned_forward_steady_state" = "numeric"))

intersect(banc_latestid(is.sel$id), meta$root_id)

#sort(meta$cell_type) 
#save csvs of cells to use for influence scores