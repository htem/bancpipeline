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
                super_class=="descending")

# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)
ids <- unique(meta$root_id)

# Get skeletons for visualisation
l2 <- banc_read_l2skel(ids)

# Re-root to soma where this is known
banc.roots <- bancr:::banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots)

# add synapses
y = banc_add_synapses(l2, OmitFailures=TRUE)

# Save the state of this R session, to avoid doing the step above next time you load this file
save.image()

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
      df_vnc = rbind(df_vnc, res)
      
    }, error = function(e) {
      message(sprintf("\nError processing neuron %s in region %s: %s", id, region, e$message))
      # Optionally, you can add a placeholder row to df_brain to indicate the error
      # df_brain = rbind(df_brain, data.frame(neuropil = region, neuron = id, count = NA))
    })
    
    pb$tick()
  }
}


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
combine_by_neuron <- combine_by_neuron[-1,]
combine_by_neuron <- data.frame(combine_by_neuron)
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
  labs(title = "Heatmap of Neuropils by adding synapses DN",
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
  labs(title = "Heatmap of Neuropil adding synapses DN VNC SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of vnc neuropil INPUT ")

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
  labs(title = "Heatmap of Neuropil with summed synapses DN brain SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of brain DN output ")


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
  labs(title = "Heatmap of Neuropils by adding normalized synapses DN",
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
  labs(title = "Heatmap of Neuropil adding normalized synapses DN VNC SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of vnc neuropil input ")

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
       fill = "proportion of brain DN output ")




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
  labs(title = "Heatmap of Neuropils by normalized vnc*brain DN",
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
  labs(title = "Heatmap of Neuropils by normalized vnc*brain DN VNC SCALED",
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
  labs(title = "Heatmap of Neuropils by normalized vnc*brain DN BRAIN SCALED",
       x = "vnc regions",
       y = "brain regions",
       fill = "Brain normalized (vnc * brain synapses) (summed across neurons)")



# Neuron number -----------------------------------------------------------
#select cutoff for number of synapses in a neuropil for a neuron to "count" as in that neuropil
brain_cutoff = df_brain[df_brain$count>50,]  
brain_cutoff$neuropil <- sapply(brain_cutoff$neuropil, function(x) paste0("b_", x))

vnc_cutoff = df_vnc[df_vnc$count>50,]
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
  labs(title = "Heatmap of Neuropil with innervating neurons DN",
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
  labs(title = "Heatmap of Neuropil with innervating neurons DN VNC SCALED",
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
  labs(title = "Heatmap of Neuropil with innervating neurons DN brain SCALED",
       x = "VNC regions",
       y = "brain regions",
       fill = "proportion of brain DN input")