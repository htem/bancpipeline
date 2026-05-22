### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Wrangling Mo's connectivity data with an taxonomy

# Load packages
library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)


# Make sure all functions query BANC and not FAFB
choose_banc()

# get data
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


# load connectivity data
an_metadata <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_annotated_metadata.csv", col_types = cols(.default = col_guess(), `id...1` = col_character(), `id...2` = col_character()))
an_input_cosine_distances <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_input_cosine_distances.csv", col_types = cols(.default = col_guess(), `post_id` = col_character()))
an_input_synapse_profiles <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_input_synapse_count_profiles.csv", col_types = cols(.default = col_guess(), `post_id` = col_character()))
an_input_umap <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_input_umap_embeddings.csv", col_types = cols(.default = col_guess(), `post_id` = col_character()))
an_joint_cosine <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_joint_cosine_distances.csv", col_types = cols(.default = col_guess(), `post_id` = col_character()))
an_joint_umap <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_joint_umap_embeddings.csv", col_types = cols(.default = col_guess(), `post_id` = col_character()))
an_output_cosine_distances <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_output_cosine_distances.csv", col_types = cols(.default = col_guess(), `pre_id` = col_character()))
an_output_synapse_profiles <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_output_synapse_count_profiles.csv", col_types = cols(.default = col_guess(), `pre_id` = col_character()))
an_output_umap <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/ascending/ascending_output_umap_embeddings.csv", col_types = cols(.default = col_guess(), `pre_id` = col_character()))

# input synapses by an type -----------------------------------------------


# plot input cosine distances ordered by an type

#make sure ids are up to date
an_input_cosine_distances$post_id <- banc_latestid(an_input_cosine_distances$post_id)
colnames(an_input_cosine_distances) <- banc_latestid(colnames(an_input_cosine_distances))
colnames(an_input_cosine_distances)[1] <- "post_id"

# Sort neuron IDs by an_type and remove nas
sorted_neurons <- bc %>%
  arrange(an_type) %>%
  filter(!is.na(an_type)) %>%
  pull(root_id)
sorted_neurons <- banc_latestid(sorted_neurons)

# Filter the distance matrix to include only valid neurons
# Filter rows based on `sorted_neurons`
filtered_input_cosine_distances <- an_input_cosine_distances[an_input_cosine_distances$post_id %in% sorted_neurons, ]

# Filter columns, but retain the first column `post_id`
filtered_input_cosine_distances <- filtered_input_cosine_distances[, 
                                                                   colnames(filtered_input_cosine_distances) %in% c("post_id", sorted_neurons)]

filtered_sorted_neurons <- sorted_neurons[sorted_neurons %in% filtered_input_cosine_distances$post_id]
an_subset <- subset(bc, root_id %in% filtered_sorted_neurons)

sorted_matrix <- filtered_input_cosine_distances[
  match(filtered_sorted_neurons, filtered_input_cosine_distances$post_id), 
  c("post_id", filtered_sorted_neurons)
]


# Convert matrix to long format while preserving full numeric precision
numeric_matrix <- as.matrix(sorted_matrix[, -c(1)])  # Remove pre_id column
rownames(numeric_matrix) <- sorted_matrix$post_id

matrix_long <- as.data.table(numeric_matrix, keep.rownames = "Neuron1")
matrix_long <- melt(matrix_long, id.vars = "Neuron1", variable.name = "Neuron2", value.name = "CosineDistance")

# Add class information
matrix_long <- matrix_long %>%
  mutate(Class1 = bc$an_type[match(Neuron1, bc$root_id)],
         Class2 = bc$an_type[match(Neuron2, bc$root_id)])

# Create a sorting order that groups by an_type and maintains neuron order within each type
neuron_order <- matrix_long %>%
  mutate(an_type = bc$an_type[match(Neuron1, bc$root_id)]) %>%
  arrange(an_type, Neuron1) %>%
  pull(Neuron1) %>%
  unique()

# Convert Neuron1 and Neuron2 to factors with the custom order
matrix_long <- matrix_long %>%
  mutate(
    Neuron1 = factor(Neuron1, levels = neuron_order),
    Neuron2 = factor(Neuron2, levels = neuron_order)
  )

# Plot with the neuron labels
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white") +
  labs(title = "AN Input Distance Matrix Sorted by Class",
       x = "Neuron ID",
       y = "Neuron ID",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Ensure a mapping from neuron IDs to an_type
neuron_to_type <- matrix_long %>%
  distinct(Neuron1, Class1) %>%
  rename(Neuron = Neuron1, an_type = Class1)

# Plot the heatmap with an type labels
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, max(matrix_long$CosineDistance, na.rm = TRUE))) +
  labs(title = "AN Input Distance Matrix Sorted by Class",
       x = "AN Type",
       y = "AN Type ",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank()) +
  # Replace x and y-axis labels with neuron an_type
  scale_x_discrete(labels = neuron_to_type$an_type) +
  scale_y_discrete(labels = neuron_to_type$an_type)


# now plot with dividers for type
# Create a table with unique neurons per class
unique_neurons <- matrix_long %>%
  distinct(Class1, Neuron1)  # Keep unique Class1-Neuron1 combinations
# Calculate cumulative positions for boundaries
x_boundaries <- unique_neurons %>%
  group_by(Class1) %>%
  summarise(position = n()) %>%         # Count unique neurons in each class
  mutate(cumulative_position = cumsum(position))  # Calculate cumulative totals
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, max(matrix_long$CosineDistance, na.rm = TRUE))) +
  labs(title = "AN Input Distance Matrix Sorted by AN Type",
       x = "AN",
       y = "AN",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # Add class separation lines
  geom_vline(data = x_boundaries, aes(xintercept = cumulative_position + 0.5), color = "red", linetype = "dashed") +
  geom_hline(data = x_boundaries, aes(yintercept = cumulative_position + 0.5), color = "red", linetype = "dashed")

# output synapses by an type -----------------------------------------------


# plot output cosine distances ordered by an type

#make sure ids are up to date
an_output_cosine_distances$pre_id <- banc_latestid(an_output_cosine_distances$pre_id)
colnames(an_output_cosine_distances) <- banc_latestid(colnames(an_output_cosine_distances))
colnames(an_output_cosine_distances)[1] <- "pre_id"


# Filter the distance matrix to include only valid neurons
# Filter rows based on `sorted_neurons`
filtered_output_cosine_distances <- an_output_cosine_distances[an_output_cosine_distances$pre_id %in% sorted_neurons, ]

# Filter columns, but retain the first column `post_id`
filtered_output_cosine_distances <- filtered_output_cosine_distances[, 
                                                                   colnames(filtered_output_cosine_distances) %in% c("pre_id", sorted_neurons)]

filtered_sorted_neurons <- sorted_neurons[sorted_neurons %in% filtered_output_cosine_distances$pre_id]
an_subset <- subset(bc, root_id %in% filtered_sorted_neurons)

sorted_matrix <- filtered_output_cosine_distances[
  match(filtered_sorted_neurons, filtered_output_cosine_distances$pre_id), 
  c("pre_id", filtered_sorted_neurons)
]


# Convert matrix to long format while preserving full numeric precision
numeric_matrix <- as.matrix(sorted_matrix[, -c(1)])  # Remove pre_id column
rownames(numeric_matrix) <- sorted_matrix$pre_id

matrix_long <- as.data.table(numeric_matrix, keep.rownames = "Neuron1")
matrix_long <- melt(matrix_long, id.vars = "Neuron1", variable.name = "Neuron2", value.name = "CosineDistance")

# Add class information
matrix_long <- matrix_long %>%
  mutate(Class1 = bc$an_type[match(Neuron1, bc$root_id)],
         Class2 = bc$an_type[match(Neuron2, bc$root_id)])

# Create a sorting order that groups by an_type and maintains neuron order within each type
neuron_order <- matrix_long %>%
  mutate(an_type = bc$an_type[match(Neuron1, bc$root_id)]) %>%
  arrange(an_type, Neuron1) %>%
  pull(Neuron1) %>%
  unique()

# Convert Neuron1 and Neuron2 to factors with the custom order
matrix_long <- matrix_long %>%
  mutate(
    Neuron1 = factor(Neuron1, levels = neuron_order),
    Neuron2 = factor(Neuron2, levels = neuron_order)
  )

# Plot with the neuron labels
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  labs(title = "AN Output Distance Matrix Sorted by Class",
       x = "Neuron ID",
       y = "Neuron ID",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Ensure a mapping from neuron IDs to an_type
neuron_to_type <- matrix_long %>%
  distinct(Neuron1, Class1) %>%
  rename(Neuron = Neuron1, an_type = Class1)

# Plot the heatmap with an type labels
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, max(matrix_long$CosineDistance, na.rm = TRUE))) +
  labs(title = "AN Output Distance Matrix Sorted by Class",
       x = "AN Type",
       y = "AN Type ",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank()) +
  # Replace x and y-axis labels with neuron an_type
  scale_x_discrete(labels = neuron_to_type$an_type) +
  scale_y_discrete(labels = neuron_to_type$an_type)


# now plot with dividers for type
# Create a table with unique neurons per class
unique_neurons <- matrix_long %>%
  distinct(Class1, Neuron1)  # Keep unique Class1-Neuron1 combinations
# Calculate cumulative positions for boundaries
x_boundaries <- unique_neurons %>%
  group_by(Class1) %>%
  summarise(position = n()) %>%         # Count unique neurons in each class
  mutate(cumulative_position = cumsum(position))  # Calculate cumulative totals
ggplot(matrix_long, aes(x = Neuron1, y = Neuron2, fill = CosineDistance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, max(matrix_long$CosineDistance, na.rm = TRUE))) +
  labs(title = "AN Output Distance Matrix Sorted by AN Type",
       x = "AN ID",
       y = "AN ID",
       fill = "Cosine Distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # Add class separation lines
  geom_vline(data = x_boundaries, aes(xintercept = cumulative_position + 0.5), color = "red", linetype = "dashed") +
  geom_hline(data = x_boundaries, aes(yintercept = cumulative_position + 0.5), color = "red", linetype = "dashed")


# now looking at umap -----------------------------------------------------
# input
an_input_umap$post_id <- banc_latestid(an_input_umap$post_id)
an_type <- data.frame('post_id'=bc$root_id, 'an_type' = bc$an_type)
an_input_umap <- an_input_umap %>%
  left_join(an_type, by = "post_id", relationship = "many-to-many")

# Plot UMAP with points colored by an_type
ggplot(an_input_umap, aes(x = `UMAP dim 1`, y = `UMAP dim 2`, color = an_type)) +
  geom_point(size = 2) +  # Adjust point size as needed
  labs(title = "AN Input Plot Colored by AN Type",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Neuron Type") +
  theme_minimal()

# output
an_output_umap$pre_id <- banc_latestid(an_output_umap$pre_id)
an_type <- data.frame('pre_id'=bc$root_id, 'an_type' = bc$an_type)
an_output_umap <- an_output_umap %>%
  left_join(an_type, by = "pre_id", relationship = "many-to-many")

# Plot UMAP with points colored by an_type
ggplot(an_output_umap, aes(x = `UMAP dim 1`, y = `UMAP dim 2`, color = an_type)) +
  geom_point(size = 2) +  # Adjust point size as needed
  labs(title = "AN Output Plot Colored by AN Type",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Neuron Type") +
  theme_minimal()

# both
an_joint_umap$post_id <- banc_latestid(an_joint_umap$post_id)
an_type <- data.frame('post_id'=bc$root_id, 'an_type' = bc$an_type)
an_joint_umap <- an_joint_umap %>%
  left_join(an_type, by = "post_id", relationship = "many-to-many")

# Plot UMAP with points colored by an_type
ggplot(an_joint_umap, aes(x = `UMAP dim 1`, y = `UMAP dim 2`, color = an_type)) +
  geom_point(size = 2) +  # Adjust point size as needed
  labs(title = "AN Joint Plot Colored by AN Type",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Neuron Type") +
  theme_minimal()


# class descriptions from metadata ----------------------------------------
# update IDs
an_metadata$id...1 <- banc_latestid(an_metadata$id...1)

# join with an type
an_type <- data.frame('id...1'=bc$root_id, 'an_type' = bc$an_type)
an_metadata <- an_metadata %>%
  left_join(an_type, by = 'id...1', relationship = "many-to-many")

#collapse over an_type
summary_df <- an_metadata %>%
  group_by(an_type) %>%
  summarize(
    id...1 = list(id...1),  # Keep `root_id` as a list
    across(where(is.numeric) & !c(id...1), mean, na.rm = TRUE),  # Calculate mean for numeric columns except `root_id`
    .groups = "drop"  # Ungroup the result
  )

# View the result
print(summary_df)

library(ggplot2)

# Create the bar chart
library(ggplot2)

# Assume summary_df is your dataframe
# Generate bar plots for all numeric columns except root_ids
for (col_name in colnames(summary_df)[sapply(summary_df, is.numeric)]) {
  # Create the plot
  p <- ggplot(summary_df, aes(x = an_type, y = .data[[col_name]])) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(
      title = paste("Mean", col_name, "by Type"),
      x = "Neuron Type",
      y = paste("Mean", col_name)
    ) +
    theme_minimal()
  
  ggsave(filename = paste0("plot_", col_name, ".png"), plot = p, width = 6, height = 4)
  
  # Print the plot (or save it)
  print(p)  # Use ggsave() here if saving the plots
}


# heatmaps - input

# Reshape the data into long format for columns with "input" in their names
heatmap_data <- summary_df %>%
  pivot_longer(
    cols = contains("input_super") & !contains("dataset"),
    #cols = !contains("dataset"),# Select columns with "input" in their names
    names_to = "Metric",       # New column for the column names
    values_to = "Value"        # New column for the values
  )

# Create the heatmap
ggplot(heatmap_data, aes(x = Metric, y = an_type, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white") +  # Customize color scale
  labs(
    title = "Heatmap of Input Metrics by Neuron Type",
    x = "Metric",
    y = "Neuron Type",
    fill = "Value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# heatmaps - output

# Reshape the data into long format for columns with "input" in their names
heatmap_data <- summary_df %>%
  pivot_longer(
    cols = contains("output_super") & !contains("dataset"),
    #cols = !contains("dataset"),# Select columns with "input" in their names
    names_to = "Metric",       # New column for the column names
    values_to = "Value"        # New column for the values
  )

# Create the heatmap
ggplot(heatmap_data, aes(x = Metric, y = an_type, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white") +  # Customize color scale
  labs(
    title = "Heatmap of Output Metrics by Neuron Type",
    x = "Metric",
    y = "Neuron Type",
    fill = "Value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))