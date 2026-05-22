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
library(parallel)
library(progress)
library(tidyverse)
library(ggplot2)

# Make sure all functions query BANC and not FAFB
choose_banc()

# get data ----------------------------------------------------------------



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

# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)
ids <- unique(meta$root_id)
ids <- banc_latestid(ids)
ids<- unique(ids)

# Get skeletons for visualisation
l2 <- banc_read_l2skel(ids)
#l2[ids.left] <- banc_mirror(l2[ids.left])

# Re-root to soma where this is known
banc.roots <- bancr:::banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots)

# add synapses
y = banc_add_synapses(l2[1:500], OmitFailures=TRUE)
y2 = banc_add_synapses(l2[501:1200], OmitFailures=TRUE)
y3 = banc_add_synapses(l2[1201:1777], OmitFailures=TRUE)
# Save the state of this R session, to avoid doing the step above next time you load this file
save.image()
y_total <- c(y, y2, y3)
y<- y_total

# brain df -------------------------------------------------------
total_iterations <- length(banc_brain_neuropils.surf$RegionList) * length(names(l2))
pb <- progress_bar$new(
  format = "[:bar] :percent Complete. Elapsed: :elapsed. ETA: :eta",
  total = total_iterations,
  clear = FALSE,
  width = 60
)

df_brain = data.frame()

for (region in banc_brain_neuropils.surf$RegionList) {
  np = subset(banc_brain_neuropils.surf, region)
  for (id in names(l2)) {
    tryCatch({
      l = l2[[id]]
      #0 is presynaptic 1 is postsynaptic
      a = xyzmatrix(subset(y[[id]]$connectors, prepost == 0))
      b = pointsinside(a, np)
      c = a[b, ]
      
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
      
      res = data.frame(neuropil = region, neuron = id, count = count)
      df_brain = rbind(df_brain, res)
      
    }, error = function(e) {
      message(sprintf("\nError processing neuron %s in region %s: %s", id, region, e$message))
      # Optionally, you can add a placeholder row to df_brain to indicate the error
      # df_brain = rbind(df_brain, data.frame(neuropil = region, neuron = id, count = NA))
    })
    
    pb$tick()
  }
}


# vnc df --------------------------------------------------------------
total_iterations <- length(banc_vnc_neuropils.surf$RegionList) * length(names(l2))
pb <- progress_bar$new(
  format = "[:bar] :percent Complete. Elapsed: :elapsed. ETA: :eta",
  total = total_iterations,
  clear = FALSE,
  width = 60
)

df_vnc = data.frame()

for (region in banc_vnc_neuropils.surf$RegionList) {
  np = subset(banc_vnc_neuropils.surf, region)
  for (id in names(l2)) {
    tryCatch({
      l = l2[[id]]
      #0 is presynaptic 1 is postsynaptic
      a = xyzmatrix(subset(y[[id]]$connectors, prepost == 1))
      b = pointsinside(a, np)
      c = a[b, ]
      
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
      
      res = data.frame(neuropil = region, neuron = id, count = count)
      df_vnc = rbind(df_vnc, res)
      
    }, error = function(e) {
      message(sprintf("\nError processing neuron %s in region %s: %s", id, region, e$message))
      # Optionally, you can add a placeholder row to df_brain to indicate the error
      # df_brain = rbind(df_brain, data.frame(neuropil = region, neuron = id, count = NA))
    })
    
    pb$tick()
  }
}


# add synapses
# add normalized synapses
# multiply normalized 

# Add synapses ------------------------------------------------------------------

#select cutoff
non_zero_brain = df_brain[df_brain$count>0,]  
non_zero_brain$neuropil <- sapply(non_zero_brain$neuropil, function(x) paste0("b_", x))

non_zero_vnc = df_vnc[df_vnc$count>0,]
non_zero_vnc$neuropil <- sapply(non_zero_vnc$neuropil, function(x) paste0("v_", x))

# combine brain and vnc
combined_df <- bind_rows(non_zero_brain, non_zero_vnc)
combined_df$neuron_num <- 1

# collapse over neuron and neuropil
combine_by_neuron <- combined_df %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(count)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron <- t(combine_by_neuron)
colnames(combine_by_neuron) <- combine_by_neuron[1,]
combine_by_neuron<-combine_by_neuron[-1,]
combine_by_neuron <-data.frame(combine_by_neuron)
combine_by_neuron <- mutate_all(combine_by_neuron, function(x) as.numeric(as.character(x)))


# Reshape the data and calculate the sums
region_map <- combine_by_neuron %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") %>%
  group_by(b_type, v_type) %>%
  summarise(sum_value = sum(b_value + v_value)) %>%
  pivot_wider(names_from = v_type, values_from = sum_value)

##### PLOT UNSCALED #####
# Reshape the data
df_long <- region_map %>%
  pivot_longer(cols = -b_type, names_to = "variable", values_to = "value")

# Create the heatmap
ggplot(df_long, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropils by adding synapses AN",
       x = "vnc regions",
       y = "brain regions",
       fill = "number of synapses (added)")

##### PLOT VNC SCALED #####
# Normalize columns to sum to 1
normalized_region_number <- region_map %>%
  ungroup() %>%  # Remove any grouping
  mutate(across(where(is.numeric), ~ ./sum(., na.rm = TRUE))) %>%
  group_by(!!!groups(region_map))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil adding synapses AN VNC SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of vnc neuropil output ")
##### PLOT BRAIN SCALED #####
normalized_region_number <- region_map %>%
  ungroup() %>%  # Remove any grouping
  rowwise() %>%  # Operate row by row
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%  # Remove rowwise grouping
  mutate(across(where(is.numeric), ~ . / row_sum)) %>%
  select(-row_sum) %>%  # Remove the row_sum column
  group_by(!!!groups(region_map))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil with summed synapses AN brain SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of brain AN  input ")


# Add normalized synapses -------------------------------------------------
#select cutoff
# norm should be # of synapses / total number of synapses in that neuron in the brain
#add normalized values
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
combine_by_neuron <- combined_df_norm %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(norm)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron <- t(combine_by_neuron)
colnames(combine_by_neuron) <- combine_by_neuron[1,]
combine_by_neuron<-combine_by_neuron[-1,]
combine_by_neuron <-data.frame(combine_by_neuron)
combine_by_neuron <- mutate_all(combine_by_neuron, function(x) as.numeric(as.character(x)))


# Reshape the data and calculate the sums
region_map <- combine_by_neuron %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") %>%
  group_by(b_type, v_type) %>%
  summarise(sum_value = sum(b_value + v_value)) %>%
  pivot_wider(names_from = v_type, values_from = sum_value)

##### PLOT UNSCALED #####
# Reshape the data
df_long <- region_map %>%
  pivot_longer(cols = -b_type, names_to = "variable", values_to = "value")

# Create the heatmap
ggplot(df_long, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropils by adding normalized synapses AN",
       x = "vnc regions",
       y = "brain regions",
       fill = "normalized # synapses")


##### PLOT VNC SCALED #####
# Normalize columns to sum to 1
normalized_region_number <- region_map %>%
  ungroup() %>%  # Remove any grouping
  mutate(across(where(is.numeric), ~ ./sum(., na.rm = TRUE))) %>%
  group_by(!!!groups(region_map))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil adding normalized synapses AN VNC SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of vnc neuropil output ")

##### PLOT BRAIN SCALED #####
# Normalize ROWS to sum to 1
normalized_region_number <- region_map %>%
  ungroup() %>%  # Remove any grouping
  rowwise() %>%  # Operate row by row
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%  # Remove rowwise grouping
  mutate(across(where(is.numeric), ~ . / row_sum)) %>%
  select(-row_sum) %>%  # Remove the row_sum column
  group_by(!!!groups(region_map))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil with summed normalized synapses brain region SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of brain AN  input ")




# Multiply VNC norm x Brain number ----------------------------------------

result <- non_zero_brain_norm %>%
  # Select and rename columns from brain dataframe
  select(neuron, brain_neuropil = neuropil, brain_count = count) %>%
  # Join with VNC dataframe
  full_join(
    non_zero_vnc_norm %>%
      select(neuron, vnc_neuropil = neuropil, vnc_norm = norm),
    by = "neuron"
  ) %>%
  # Calculate the product only for matching rows
  mutate(product = ifelse(!is.na(brain_count) & !is.na(vnc_norm), brain_count * vnc_norm, NA)) %>%
  # Remove rows where either brain_count or vnc_norm is NA
  filter(!is.na(brain_count) & !is.na(vnc_norm)) %>%
  # Select final columns in desired order
  select(neuron, vnc_neuropil, brain_neuropil, product)


reshaped_result <- result %>%
  # Group by vnc_neuropil and brain_neuropil
  group_by(vnc_neuropil, brain_neuropil) %>%
  # Sum the products for each group
  summarize(total_product = sum(product, na.rm = TRUE), .groups = "drop") %>%
  # Reshape the data to have brain regions as columns
  pivot_wider(
    names_from = brain_neuropil,
    values_from = total_product,
    values_fill = 0  # Fill NA values with 0
  ) %>%
  # Move vnc_neuropil to be the first column
  select(vnc_neuropil, everything())

##### PLOT UNSCALED #####

# Reshape the data
df_long <- reshaped_result %>%
  pivot_longer(cols = -vnc_neuropil, names_to = "brain_neuropil", values_to = "value")


# Create the heatmap
ggplot(df_long, aes(x = vnc_neuropil, y = brain_neuropil, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropils by normalized vnc*brain AN",
       x = "vnc regions",
       y = "brain regions",
       fill = "normalized vnc * brain synapses (summed across neurons)")

##### PLOT VNC SCALED #####
# for this dataframe, need to normalize rows to sum to one (vnc is row not column)
# Normalize rows to sum to 1
normalized_region_number <- reshaped_result %>%
  ungroup() %>%  # Remove any grouping
  rowwise() %>%  # Operate row by row
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%  # Remove rowwise grouping
  mutate(across(where(is.numeric), ~ . / row_sum)) %>%
  select(-row_sum) %>%  # Remove the row_sum column
  group_by(!!!groups(reshaped_result))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = vnc_neuropil, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropils by normalized vnc*brain AN VNC SCALED",
       x = "vnc regions",
       y = "brain regions",
       fill = "normalized vnc * brain synapses (summed across neurons)")


##### PLOT BRAIN SCALED #####
# Normalize columns to sum to 1
normalized_region_number <- reshaped_result %>%
  ungroup() %>%  # Remove any grouping
  mutate(across(where(is.numeric), ~ ./sum(., na.rm = TRUE))) %>%
  group_by(!!!groups(reshaped_result))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = vnc_neuropil, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropils by normalized vnc*brain AN BRAIN SCALED",
       x = "vnc regions",
       y = "brain regions",
       fill = "Brain normalized (vnc * brain synapses) (summed across neurons)")



# Neuron number -----------------------------------------------------------
#select cutoff for number of synapses in a neuropil for a neuron to "count" as in that neuropil
brain_cutoff = df_brain[df_brain$count>50,]  
brain_cutoff$neuropil <- sapply(brain_cutoff$neuropil, function(x) paste0("b_", x))

vnc_cutoff = df_vnc[df_vnc$count>50,]
vnc_cutoff$neuropil <- sapply(vnc_cutoff$neuropil, function(x) paste0("v_", x))


# follow-ups --------------------------------------------------------------
# combine brain and vnc
combined_df_neuron <- bind_rows(brain_cutoff, vnc_cutoff)
combined_df_neuron$neuron_num <- 1
# this is what i can search for interesting neurons
combined_df_total <- full_join(brain_cutoff, vnc_cutoff, by = "neuron", suffix = c("_brain", "_vnc"))

### ANENNAL LOBE ###
# look at AL neurons
al_ans <- subset(combined_df_total, grepl("_AL_", neuropil_brain) & grepl("MANC", neuropil_vnc))
al_an_ids <- unique(al_ans$neuron)

# Suppose 'neuron_ids' is your list of neuron IDs, and 'new_df' is the new data frame
banc_al_an <- subset(bc, root_id %in% al_an_ids)

# plot all al an neurons
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
al_r = subset(banc_brain_neuropils.surf, "ITO_midbrain_AL_R")
al_l = subset(banc_brain_neuropils.surf, "ITO_midbrain_AL_L")
plot3d(al_r, alpha=0.1)
plot3d(al_l, alpha=0.1)
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in al_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500) 
}


### MUSHROOM ###
# look at MB calyx neurons
mb_ans <- subset(combined_df_total, grepl("MB", neuropil_brain) & grepl("MANC", neuropil_vnc))
mb_an_ids <- unique(mb_ans$neuron)
mb_r = subset(banc_brain_neuropils.surf, "ITO_midbrain_MB_CA_R")
mb_l = subset(banc_brain_neuropils.surf, "ITO_midbrain_MB_CA_L")
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(mb_r, alpha=0.1)
plot3d(mb_l, alpha=0.1)
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in mb_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500) 
}

### CX ###
# look at FB  neurons
fb_ans <- subset(combined_df_total, grepl("FB", neuropil_brain) & grepl("MANC", neuropil_vnc))
fb_an_ids <- unique(fb_ans$neuron)
fb = subset(banc_brain_neuropils.surf, "ITO_midbrain_FB")
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(fb, alpha=0.1)
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in fb_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500) 
}

# LOOK AT PB
pb_ans <- subset(combined_df_total, grepl("PB", neuropil_brain) & grepl("MANC", neuropil_vnc))
pb_an_ids <- unique(pb_ans$neuron)

### MVAC ###
# look at mvac to avlp and saddel and wedge
mvac_avlp_ans <- subset(combined_df_total, grepl("mVAC", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("AVLP", neuropil_brain))
mvac_avpl_an_ids <- unique(mvac_avlp_ans$neuron)

mvac_sad_ans <- subset(combined_df_total, grepl("mVAC", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("SAD", neuropil_brain))
mvac_sad_an_ids <- unique(mvac_sad_ans$neuron)

mvac_wed_ans <- subset(combined_df_total, grepl("mVAC", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("WED", neuropil_brain))
mvac_wed_an_ids <- unique(mvac_wed_ans$neuron)

mvac_ans <- subset(combined_df_total, grepl("mVAC", neuropil_vnc) & grepl("MANC", neuropil_vnc))
mvac_an_ids <- unique(mvac_ans$neuron)

# plot all mvac to saddel, avlp, and wedge neurons
# plot all al an neurons
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
avlp_r = subset(banc_brain_neuropils.surf, "ITO_midbrain_AVLP_R")
avlp_l = subset(banc_brain_neuropils.surf, "ITO_midbrain_AVLP_L")
plot3d(avlp_l, alpha = 0.1)
plot3d(avlp_r, alpha = 0.1)
wed_r = subset(banc_brain_neuropils.surf, "ITO_midbrain_WED_R")
wed_l = subset(banc_brain_neuropils.surf, "ITO_midbrain_WED_L")
plot3d(wed_l, alpha = 0.1)
plot3d(wed_r, alpha = 0.1)
sad = subset(banc_brain_neuropils.surf, "ITO_midbrain_SAD")
plot3d(sad, alpha = 0.1)
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in mvac_avpl_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="green") 
}
for (neuron_id in mvac_sad_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="blue") 
}
for (neuron_id in mvac_wed_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="red") 
}

#look in banc data for type
banc_mvac_avlp_an <- subset(bc, root_id %in% mvac_avpl_an_ids)
banc_mvac_an <- subset(bc, root_id %in% mvac_an_ids)

### HALTERES ### 
# look at haltere to ips
haltere_an <- subset(combined_df_total, grepl("HTct", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("IPS", neuropil_brain))
haltere_an_ids <- unique(haltere_an$neuron)
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in haltere_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="pink") 
}

### ABD NEUROMERE ###
abd_an <- subset(combined_df_total, grepl("ANm", neuropil_vnc) & grepl("MANC", neuropil_vnc))
abd_an_ids <- unique(abd_an$neuron)
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

# Loop over each neuron ID in neuron_ids
for (neuron_id in abd_an_ids) {
  # Filter data for the current neuron ID
  plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="purple") 
}
### SUBDIVIDE GNG ###
t1_gng <- subset(combined_df_total, grepl("LNp_T1", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("GNG", neuropil_brain))
t1_gng_ids <- unique(t1_gng$neuron)
gng <-subset(combined_df_total, grepl("GNG", neuropil_brain))
gng <- unique(gng$neuron)
wing_gng <- subset(combined_df_total, grepl("WTct", neuropil_vnc) & grepl("MANC", neuropil_vnc) & grepl("GNG", neuropil_brain))
wing_gng_ids <- unique(wing_gng$neuron)
#plot leg and wing
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
#plot3d(banc_brain_neuropils.surf$Regions$ITO_midbrain_CRE_L, alpha=0.2, col="lightblue")
# Load ggplot2 for plotting
library(ggplot2)

##

# # Loop over each neuron ID in neuron_ids
# for (neuron_id in t1_gng_ids) {
#   # Filter data for the current neuron ID
#   plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="orange") 
# }
# for (neuron_id in wing_gng_ids) {
#   # Filter data for the current neuron ID
#   plot3d(y[neuron_id], lwd=0.5, soma = 2500, col="lightblue") 
# }


# get sensory nerves upstream mvac ans ------------------------------------
mvac_an_ids <- banc_latestid(mvac_an_ids)
banc_mvac_ans <- subset(bc, root_id %in% mvac_an_ids)
nmvac_AN <- length(banc_mvac_ans)
inputs_ans <- data.frame()

for (l in 1:nmvac_AN){
  print(paste('Working on', banc_mvac_ans[l,]$root_id, '...'))
  inputs_an = banc_partner_summary(banc_mvac_ans[l,]$root_id, partners = 'in', threshold = 10)
  inputs_ans = rbind(inputs_ans, inputs_an)
}

an_input_vnc <- subset(bc, root_id %in% inputs_ans$pre_id)
ok_celllabels <- c("ascending, innervates leg, sensory neuron", "ascending, sensory neuron", "wing sensory neuron", "sensory neuron", "sensory_neuron")
an_input_vnc <- subset(an_input_vnc, cell_class %in% ok_celllabels)
sn_upstream_ans <- data.frame(root_id = an_input_vnc$root_id)
# reshaping to plot heatmaps ----------------------------------------------



combine_by_neuron_num <- combined_df_neuron %>%
  group_by(neuron, neuropil) %>%
  summarise(total = sum(neuron_num)) %>%
  pivot_wider(names_from = neuron, values_from = total, values_fill = 0)
combine_by_neuron_num <- t(combine_by_neuron_num)
colnames(combine_by_neuron_num) <- combine_by_neuron_num[1,]
combine_by_neuron_num<-combine_by_neuron_num[-1,]
combine_by_neuron_num <-data.frame(combine_by_neuron_num)
combine_by_neuron_num <- mutate_all(combine_by_neuron_num, function(x) as.numeric(as.character(x)))


# Reshape the data and calculate the sums
region_map <- combine_by_neuron_num %>%
  # Exclude the id column since it's not needed
  pivot_longer(cols = starts_with("b_"), names_to = "b_type", values_to = "b_value") %>%
  pivot_longer(cols = starts_with("v_"), names_to = "v_type", values_to = "v_value") 
region_map$both <- region_map$b_value+region_map$v_value
region_number <-subset(region_map, both==2)
# remove GNG to visualize better
# region_number <- subset(region_number, b_type != "b_GNG")

region_number <- region_number%>%
  group_by(b_type, v_type) %>%
  summarise(sum_value = sum(b_value + v_value)) %>%
  pivot_wider(names_from = v_type, values_from = sum_value)


#na to 0
region_number <- region_number %>% mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate_all(~ . / 2)

# Reshape the data
df_long <- region_number %>%
  pivot_longer(cols = -b_type, names_to = "variable", values_to = "value")

##### PLOT UNSCALED #####
ggplot(df_long, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil with innervating neurons AN",
       x = "VNC regions",
       y = "brain regions",
       fill = "number of neurons")

##### PLOT VNC SCALED #####
# Normalize columns to sum to 1
normalized_region_number <- region_number %>%
  ungroup() %>%  # Remove any grouping
  mutate(across(where(is.numeric), ~ ./sum(., na.rm = TRUE))) %>%
  group_by(!!!groups(region_number))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil with innervating neurons AN VNC SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of vnc neuropil output ")

##### PLOT BRAIN SCALED #####

# Normalize ROWS to sum to 1
normalized_region_number <- region_number %>%
  ungroup() %>%  # Remove any grouping
  rowwise() %>%  # Operate row by row
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%  # Remove rowwise grouping
  mutate(across(where(is.numeric), ~ . / row_sum)) %>%
  select(-row_sum) %>%  # Remove the row_sum column
  group_by(!!!groups(region_number))  # Reapply the original grouping

# Reshape the data into long format
heatmap_data <- normalized_region_number %>%
  pivot_longer(cols = where(is.numeric), 
               names_to = "variable", 
               values_to = "value")

# Create the heatmap
ggplot(heatmap_data, aes(x = variable, y = b_type, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # This uses a color-blind friendly palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Heatmap of Neuropil with innervating neurons AN brain SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of brain AN  input ")