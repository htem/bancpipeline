### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# get neuron to neuropil synapses

######################
### load libraries ###
######################

library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(R.matlab)
library(stringr)
library(progress)

###############################
### load connectivity data ####
###############################

# define directory where to save data
datadir <- str_c(getwd(), "/data/dn/")

# load relevant variables for clustering and plotting
l2 <- readRDS(str_c(datadir, "l2_dn.rds"))
meta <- readRDS(str_c(datadir, "meta_dn.rds"))

#########################
### define functions ####
#########################

sort_append_collapse_by_dn_type <- function(input_id_to_np_count, dn_subset) {
  
  # get ids from input_id_to_np_count
  latest_id <- banc_latestid(rownames(input_id_to_np_count))
  rownames(input_id_to_np_count) <- latest_id
  
  # reduce size of dn_subset to match neurons from input_id_to_np_count
  dn_subset <- dn_subset %>%
    filter(root_id %in% latest_id) %>%
    distinct(root_id, .keep_all = TRUE)
  
  # sort input_id_to_np_count based on dn_subset$root_id
  input_id_to_np_count <- input_id_to_np_count[match(dn_subset$root_id, rownames(input_id_to_np_count)), ]
  
  # test if the two arrays are identical
  if (identical(dn_subset$root_id, rownames(input_id_to_np_count))) {
    print('correct match')
  }
  
  # add dn type field
  input_id_to_np_count$dn_type <- as.numeric(dn_subset$dn_type)
  
  # collapse over dn_type
  combine_by_dn_type <- input_id_to_np_count %>%
    group_by(dn_type) %>%
    summarize(
      # Calculate mean for numeric columns
      across(where(is.numeric), mean, na.rm = TRUE),
      # Ungroup the result
      .groups = "drop"
    )
  
  return(list(combine_by_dn_type = combine_by_dn_type))
  
}

##################################
### neuron to Neuropil overlap ###
##################################

# add synapses to each neuron
names(l2) <- banc_latestid(names(l2))
l2_syn = banc_add_synapses(l2, OmitFailures = TRUE)

# 1) brain df

# make progress bar
total_iterations <- length(banc_brain_neuropils.surf$RegionList) * length(names(l2))
pb <- progress_bar$new(
  format = "[:bar] :percent Complete. Elapsed: :elapsed. ETA: :eta",
  total = total_iterations,
  clear = FALSE,
  width = 60
)

# initialize brain neuropil overlap variable
df_brain = data.frame()

# evaluate overlap of neurons to each neuropil in banc_brain_neuropils.surf$RegionList
for (region in banc_brain_neuropils.surf$RegionList) {
  np = subset(banc_brain_neuropils.surf, region)
  
  for (id in names(l2)) {
    
    tryCatch({
      
      # select input synapses (postsynaptic) to selected neuron (from neuronlist)
      #   0 is presynaptic 1 is postsynaptic
      a = xyzmatrix(subset(l2_syn[[id]]$connectors, prepost == 1))
      # evaluate overlap with selected neuropil
      b = pointsinside(a, np)
      # select non-empty XYZ points synapses contained in this neuropil
      c = a[b, ]
      
      # calculate the number of synapses contained in this neuropil
      if (is.null(c)) {
        count = 0
        norm = 0
      } else if (is.null(nrow(c))) {
        count = 1
        norm = 1 / nrow(a)
      } else if (nrow(c) == 0) {
        count = 0
        norm = 0
      } else {
        count = nrow(c)
        norm = nrow(c) / nrow(a)
      }
      
      # append this synapse count to variable res
      #   (which contains fields: neuropil, neuron, and count)
      res = data.frame(neuropil = region, neuron = id, count = count)
      # append res to df_brain
      df_brain = rbind(df_brain, res)
      
    }, error = function(e) {
      message(sprintf("\nError processing neuron %s in region %s: %s", id, region, e$message))
      # Optionally, you can add a placeholder row to df_brain to indicate the error
      # df_brain = rbind(df_brain, data.frame(neuropil = region, neuron = id, count = NA))
    })
    
    pb$tick()
    
  }
  
}

# 2) vnc df

# make progress bar
total_iterations <- length(banc_vnc_neuropils.surf$RegionList) * length(names(l2))
pb <- progress_bar$new(
  format = "[:bar] :percent Complete. Elapsed: :elapsed. ETA: :eta",
  total = total_iterations,
  clear = FALSE,
  width = 60
)

# initialize brain neuropil overlap variable
df_vnc = data.frame()

# evaluate overlap of neurons to each neuropil in banc_vnc_neuropils.surf$RegionList
for (region in banc_vnc_neuropils.surf$RegionList) {
  
  np = subset(banc_vnc_neuropils.surf, region)
  
  for (id in names(l2)) {
    
    tryCatch({
      
      # select input synapses (presynaptic) to selected neuron (from neuronlist)
      #   0 is presynaptic 1 is postsynaptic
      a = xyzmatrix(subset(l2_syn[[id]]$connectors, prepost == 0))
      # evaluate overlap with selected neuropil
      b = pointsinside(a, np)
      # select non-empty XYZ points synapses contained in this neuropil
      c = a[b, ]
      
      # calculate the number of synapses contained in this neuropil
      if (is.null(c)) {
        count = 0
        norm = 0
      } else if (is.null(nrow(c))) {
        count = 1
        norm = 1 / nrow(a)
      } else if (nrow(c) == 0) {
        count = 0
        norm = 0
      } else {
        count = nrow(c)
        norm = nrow(c) / nrow(a)
      }
      
      # append this synapse count to variable res
      #   (which contains fields: neuropil, neuron, and count)
      res = data.frame(neuropil = region, neuron = id, count = count)
      # append res to df_brain
      df_vnc = rbind(df_vnc, res)
      
    }, error = function(e) {
      message(sprintf("\nError processing neuron %s in region %s: %s", id, region, e$message))
      # Optionally, you can add a placeholder row to df_brain to indicate the error
      # df_brain = rbind(df_brain, data.frame(neuropil = region, neuron = id, count = NA))
    })
    
    pb$tick()
    
  }
  
}

# define directory where to save data
datadir <- str_c(getwd(), "/data/dn_synapses")
if (!dir.exists(datadir)) {
  dir.create(datadir, recursive = TRUE)
}

# save relevant variables
saveRDS(df_vnc, file = str_c(datadir, "df_vnc.rds"))
saveRDS(df_brain, file = str_c(datadir, "df_brain.rds"))
saveRDS(l2_syn, file = str_c(datadir, "l2_syn_dns.rds"))

# Add synapse number per neuropil

# remove rows of neurons without any synapse in brain or vnc neuropils
non_zero_brain = df_brain[df_brain$count > 0,]  
non_zero_brain$neuropil <- sapply(non_zero_brain$neuropil, function(x) paste0("b_", x))

non_zero_vnc = df_vnc[df_vnc$count > 0,]
non_zero_vnc$neuropil <- sapply(non_zero_vnc$neuropil, function(x) paste0("v_", x))

# combine brain and vnc neuropils
combined_df <- bind_rows(non_zero_brain, non_zero_vnc)
combined_df$neuron_num <- 1

# collapse synapse number over neuron and neuropil into a matrix of neuron x neuropil
combine_by_neuron <- combined_df %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(count)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron <- t(combine_by_neuron)
colnames(combine_by_neuron) <- combine_by_neuron[1,]
combine_by_neuron <- combine_by_neuron[-1,]
combine_by_neuron <- data.frame(combine_by_neuron)
combine_by_neuron <- mutate_all(combine_by_neuron, function(x) as.numeric(as.character(x)))

# reshape the data and calculate the sums per neuropil
region_map <- combine_by_neuron %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") %>%
  group_by(b_type, v_type) %>%
  summarise(sum_value = sum(b_value + v_value)) %>%
  pivot_wider(names_from = v_type, values_from = sum_value)

# add normalized synapse number per neuropil

# remove rows of neurons without any synapse in brain or vnc neuropils
# norm should be # of synapses / total number of synapses in that neuron in the brain
# add normalized values
df_brain <- df_brain %>%
  group_by(neuron) %>%
  mutate(total = sum(count)) %>%
  mutate(norm = count / total) %>%
  ungroup()
df_vnc <- df_vnc %>%
  group_by(neuron) %>%
  mutate(total = sum(count)) %>%
  mutate(norm = count / total) %>%
  ungroup()

non_zero_brain_norm = subset(df_brain, norm > 0)  
non_zero_brain_norm$neuropil <- sapply(non_zero_brain_norm$neuropil, function(x) paste0("b_", x))

non_zero_vnc_norm = subset(df_vnc, norm > 0) 
non_zero_vnc_norm$neuropil <- sapply(non_zero_vnc_norm$neuropil, function(x) paste0("v_", x))

# combine brain and vnc
combined_df_norm <- bind_rows(non_zero_brain_norm, non_zero_vnc_norm)
combined_df_norm$neuron_num <- 1

# collapse over neuron and neuropil
combine_by_neuron_norm <- combined_df_norm %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(norm)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron_norm <- t(combine_by_neuron_norm)
colnames(combine_by_neuron_norm) <- combine_by_neuron_norm[1,]
combine_by_neuron_norm <- combine_by_neuron_norm[-1,]
combine_by_neuron_norm <-data.frame(combine_by_neuron_norm)
combine_by_neuron_norm <- mutate_all(combine_by_neuron_norm, function(x) as.numeric(as.character(x)))

# reshape the data and calculate the sums
region_map_norm <- combine_by_neuron_norm %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") %>%
  group_by(b_type, v_type) %>%
  summarise(sum_value = sum(b_value + v_value)) %>%
  pivot_wider(names_from = v_type, values_from = sum_value)

# add neuron number per neuropil

#select cutoff for number of synapses in a neuropil for a neuron to "count" as in that neuropil
brain_cutoff = df_brain[df_brain$count > 50,]  
brain_cutoff$neuropil <- sapply(brain_cutoff$neuropil, function(x) paste0("b_", x))

vnc_cutoff = df_vnc[df_vnc$count > 50,]
vnc_cutoff$neuropil <- sapply(vnc_cutoff$neuropil, function(x) paste0("v_", x))

# combine brain and vnc
combined_df_neuron <- bind_rows(brain_cutoff, vnc_cutoff)
combined_df_neuron$neuron_num <- 1

combine_by_neuron_num <- combined_df_neuron %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(neuron_num)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron_num <- t(combine_by_neuron_num)
colnames(combine_by_neuron_num) <- combine_by_neuron_num[1,]
combine_by_neuron_num <- combine_by_neuron_num[-1,]
combine_by_neuron_num <- data.frame(combine_by_neuron_num)
combine_by_neuron_num <- mutate_all(combine_by_neuron_num, function(x) as.numeric(as.character(x)))

# reshape the data and calculate the sums
region_map_num <- combine_by_neuron_num %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") 
region_map_num$both <- region_map_num$b_value+region_map_num$v_value
region_number <- subset(region_map_num, both == 2)

# save all matrices
saveRDS(combine_by_neuron, file = str_c(datadir, "combine_by_neuron.rds"))
saveRDS(region_map, file = str_c(datadir, "region_map.rds"))
saveRDS(combine_by_neuron_norm, file = str_c(datadir, "combine_by_neuron_norm.rds"))
saveRDS(region_map_norm, file = str_c(datadir, "region_map_norm.rds"))
saveRDS(combine_by_neuron_num, file = str_c(datadir, "combine_by_neuron_num.rds"))
saveRDS(region_map_num, file = str_c(datadir, "region_map_num.rds"))

rm(datadir)

###############################
### get data from BANCTABLE ###
###############################

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

#####################################################
## collapse synapses per neuropil and per dn type ###
#####################################################

# if only running this part seprataly load the requiered variables:
# datadir <- str_c(getwd(), "/data/dn_synapses/")
# combine_by_neuron <- readRDS(str_c(datadir, "combine_by_neuron.rds"))
# combine_by_neuron_norm <- readRDS(str_c(datadir, "combine_by_neuron_norm.rds"))
# combine_by_neuron_num <- readRDS(str_c(datadir, "combine_by_neuron_num.rds"))

# sort dns by dn_type, remove na, update root_id, and keep unique
dn_subset <- meta %>%
  arrange(as.numeric(dn_type)) %>%
  filter(!is.na(dn_type))
dn_subset$root_id <- banc_latestid(dn_subset$root_id)
dn_subset <- dn_subset %>%
  distinct(root_id, .keep_all = TRUE)

# sort ids, append dn_type and collapse by dn type
combine_by_dn_type <- sort_append_collapse_by_dn_type(
  combine_by_neuron, dn_subset)$combine_by_dn_type

combine_by_dn_type_norm <- sort_append_collapse_by_dn_type(
  combine_by_neuron_norm, dn_subset)$combine_by_dn_type

combine_by_dn_type_num <- sort_append_collapse_by_dn_type(
  combine_by_neuron_num, dn_subset)$combine_by_dn_type

# define data directory
datadir <- str_c(getwd(), "/data/dn_synapses/")

# save as rds
saveRDS(combine_by_dn_type, file = str_c(datadir, "combine_by_dn_type.rds"))
saveRDS(combine_by_dn_type_norm, file = str_c(datadir, "combine_by_dn_type_norm.rds"))
saveRDS(combine_by_dn_type_num, file = str_c(datadir, "combine_by_dn_type_num.rds"))

# save a table of cosine distance matrices
write.csv(combine_by_dn_type, file = str_c(datadir, "combine_by_dn_type.csv"), 
          row.names = FALSE)
write.csv(combine_by_dn_type_norm, file = str_c(datadir, "combine_by_dn_type_norm.csv"), 
          row.names = FALSE)
write.csv(combine_by_dn_type_num, file = str_c(datadir, "combine_by_dn_type_num.csv"), 
          row.names = FALSE)

rm(datadir)