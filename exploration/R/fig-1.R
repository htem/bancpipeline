### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# For testing out figures for Figure 3/4 of the banc paper

library(dplyr)
library(purrr)
library(readr)
library(bancr)
library(tidyr)
library(data.table)
library(tibble)  # for column_to_rownames

# AN OR DN
dn = "YES"

# import nb.all.full
#nb.all.full<- readRDS("/Users/sophiarenauld/Documents/GitHub/bancpipeline/data/nb.all.dns.rds")
#nb.all.full <- nb.all.full.dn
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

if (dn == "YES") {
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
}


# recent nb.all.full names
row_names <- rownames(nb.all.full)
# Apply the 'banc_id_isrecent()' function to each row name
row_names <- banc_latestid(row_names)

# Rename the rows of the dataframe
rownames(nb.all.full) <- row_names
colnames(nb.all.full) <- row_names

# Function to remove duplicate rows or columns
remove_duplicates <- function(mat) {
  # Get unique row/column names
  unique_names <- unique(rownames(mat))
  
  # Keep only the first occurrence of each unique row/column
  mat <- mat[match(unique_names, rownames(mat)), ]
  mat <- mat[, match(unique_names, colnames(mat))]
  
  return(mat)
}

# Remove duplicate rows and columns
similarity_matrix <- remove_duplicates(nb.all.full)

# Print the dimensions of the new matrix
print(dim(similarity_matrix))

an_working <- subset(bc, root_id %in% rownames(similarity_matrix))
an_working <- an_working[!duplicated(an_working$root_id), ]
#an_working$an_type
column_to_keep <- ifelse(dn == "YES", "dn_type", "an_type")
an_working <- an_working[, c("root_id", column_to_keep)]

if (dn == "YES") {
  an_working$an_type <- an_working$dn_type
  an_working <- an_working %>% select(-dn_type)
}
  


nb_df <- as.data.frame(similarity_matrix)
nb_df$root_id <- rownames(similarity_matrix)

merged_df <- nb_df %>% left_join(an_working, by = "root_id")
collapsed_df <- merged_df %>%
  group_by(an_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            root_id = first(root_id)) %>%
  ungroup()


neuron_to_type <- an_working %>% select(root_id, an_type)


# Remove root_id if it exists
collapsed_df <- collapsed_df %>% select(-any_of("root_id"))
library(tidyr)

type_by_type <- collapsed_df %>%
  # Convert to long format
  pivot_longer(cols = -an_type, names_to = "neuron_id", values_to = "value") %>%
  # Join with the neuron_to_type mapping
  left_join(neuron_to_type, by = c("neuron_id" = "root_id")) %>%
  # Group by row type and column type
  group_by(an_type.x, an_type.y) %>%
  # Sum the values
  summarise(total_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  # Spread back to wide format
  pivot_wider(names_from = an_type.y, values_from = total_value, values_fill = 0) %>%
  # Rename the row type column
  rename(an_type = an_type.x)


type_by_type <- subset(type_by_type, !is.na(an_type))

# Assuming type_by_type is your current dataframe
type_by_type <- type_by_type %>%
  column_to_rownames(var = "an_type")
type_by_type <- type_by_type[, !grepl("NA", colnames(type_by_type), ignore.case = TRUE)]


# now should make dendrogram

# dendrogram plots --------------------------------------------------------


# just dendrogram
# Perform clustering
dist_matrix <- dist(type_by_type)
hc <- hclust(dist_matrix, method = "complete")  # You can change the method as needed

# Plot the dendrogram
plot(hc, main = "Dendrogram of DN Types (NBLAST)", 
     xlab = "", sub = "", 
     hang = -1)  # hang = -1 makes the labels appear at the same level
library(pheatmap)

# Create a heatmap with dendrograms
pheatmap(type_by_type,
         main = "Clustered Heatmap of Type-by-Type Matrix",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize = 8,
         cellwidth = 10,
         cellheight = 10,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")


# Zaki csvs ---------------------------------------------------------------


# Set working directory
setwd("/Users/papers/BANC-project/data/influence/frankenbrain_v1.1/modalities")
# more expansive
setwd("/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/cell_function")


# Get list of all CSV files
csv_files <- list.files(pattern = "*.csv")

# Function to read a CSV and keep only 'id' and 'distance' columns
read_and_filter <- function(file) {
  dt <- fread(file)
  
  # Get column names containing 'influence'
  distance_cols <- grep("influence_count_unsigned", names(dt), value = TRUE, ignore.case = TRUE)
  
  # old more expansive version
  #distance_cols <- grep("distance", names(dt), value = TRUE, ignore.case = TRUE)
  # Keep only 'id' and 'distance' columns
  cols_to_keep <- c("id", distance_cols)
  dt <- dt[, ..cols_to_keep]
  
  return(dt)
}

# Read and filter all CSV files
data_list <- lapply(csv_files, read_and_filter)

# Merge all dataframes
influence_df <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), data_list)


# rename influence df influence_count_unsigned_forward_steady_state
influence_df <- influence_df %>%
  rename_with(~ gsub("_influence_count_unsigned_forward_steady_state", "", .x))

# Check the dimensions of the merged dataframe
print(dim(influence_df))

# View the first few rows
print(head(influence_df))
min(influence_df$chemosensory)

# now merging with frankenbrain
# Load libraries
library(bancr)
library(tidyverse)
library(arrow)

# frankenbrain ------------------------------------------------------------


# Get the meta data, cell types, etc.
franken.meta <- arrow::read_feather("/Volumes/Neurobio/wilsonlab/banc/connectivity/frankenbrain_v.1.5_meta.feather") %>%
  dplyr::group_by(cell_type) %>%
  #dplyr::arrange(dplyr::desc(top_p)) %>%
  dplyr::mutate(top_nt=top_nt[1]) %>%
  dplyr::ungroup()

# merge frankenbrain, influence metric, and an type -----------------------

franken_ans <- subset(franken.meta, cell_class == "ascending_neuron")
if (dn == "YES") {
  franken_ans <- subset(franken.meta, cell_class == "descending_neuron")
}
matched_franken_ans <- subset(franken_ans, banc_id != "")
influence_ans <- subset(influence_df, id %in% matched_franken_ans$id)

merging <- matched_franken_ans[, c("id", "banc_id", "cell_sub_class")]
# merging$banc_id <- banc_islatest(merging$banc_id)
influence_ans$id <- as.character(influence_ans$id)
influence_ans <- influence_ans %>% left_join(merging, "id")


# right now, only ~1100 are matched across datasets and in zaki's data
# need real AN types in data

# now merge an influence with working ans
influence_ans <- influence_ans %>% left_join(an_working, by = c("banc_id" = "root_id"))

type_by_type_influence <- influence_ans %>% group_by(an_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            id = first(id)) %>%
            filter(!is.na(an_type)) %>%
            column_to_rownames(var = "an_type") %>%
  ungroup()

type_by_type_influence <- type_by_type_influence %>% select(-id)
  

# 
influence_df_chemosensory <- type_by_type_influence %>%
  select(!contains("unknown") | contains("CX"))
type_by_type_influence <- influence_df_chemosensory


# plotting ----------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dendextend)

an_type_order <- hc$labels[hc$order]

# Reorder your data according to the dendrogram
type_by_type_influence_numeric_ordered <- type_by_type_influence[an_type_order, ]


# Convert data to matrix
mat <- as.matrix(type_by_type_influence_numeric_ordered)
# Store the original row names
original_rownames <- rownames(mat)

# Convert the matrix to numeric
mat <- apply(mat, 2, as.numeric)

# Reassign the original row names
rownames(mat) <- original_rownames

# Create color function
col_fun = colorRamp2(c(min(mat, na.rm = TRUE), 
                       max(mat, na.rm = TRUE)), 
                     c("white", "firebrick3"))
library(dendextend)

# Convert hc to a dendrogram if it's not already
dend <- as.dendrogram(hc)

# Remove the leaf labeled "100"
if (dn != "YES") {
  dend_pruned <- prune(dend, c("100", "65"))
} else {
  dend_pruned <- dend
}


# Get the labels in the order they appear in the dendrogram
correct_order <- labels(dend_pruned)

# Print this order to verify
print(correct_order)

mat <- mat[correct_order, ]

# Calculate the maximum width for row names
max_row_name_width <- max(strwidth(rownames(mat), units = "inches", cex = 0.8))

# Create the heatmap with adjusted parameters
ht <- Heatmap(mat,
              name = "Value",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_dend = TRUE,
              row_dend_width = unit(1.5, "cm"),  # Slightly larger dendrogram width
              column_names_rot = 45,  # 45-degree angle for column names
              column_names_gp = gpar(fontsize = 7),  # Slightly larger font for column names
              row_names_gp = gpar(fontsize = 7),  # Slightly larger font for row names
              row_names_max_width = unit(max_row_name_width, "inches"),  # Limit row name width
              column_names_max_height = unit(6, "cm"),  # More space for column names
              heatmap_legend_param = list(title = "Value", labels_gp = gpar(fontsize = 8)),
              width = ncol(mat) * unit(3, "mm"),  # Larger cell width
              height = nrow(mat) * unit(3, "mm"))  # Larger cell height

# Draw the heatmap
draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 20), "mm"))




plot(dend_pruned)

# plot normalized columns
# Function to compute Z-scores
zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Save row names
row_names <- rownames(mat)

# Apply Z-score normalization across columns
mat_zscored <- apply(mat, 2, zscore)

# Reapply row names
rownames(mat_zscored) <- row_names

# Calculate the maximum width for row names
max_row_name_width <- max(strwidth(rownames(mat_zscored), units = "inches", cex = 0.8))

# Create the heatmap with Z-scored data
ht <- Heatmap(mat_zscored,
              name = "Z-score",
              col = colorRamp2(c(min(mat_zscored), max(mat_zscored)), c("white", "red")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_dend = TRUE,
              row_dend_width = unit(1.5, "cm"),  # Slightly larger dendrogram width
              column_names_rot = 45,  # 45-degree angle for column names
              column_names_gp = gpar(fontsize = 7),  # Slightly larger font for column names
              row_names_gp = gpar(fontsize = 7),  # Slightly larger font for row names
              row_names_max_width = unit(max_row_name_width, "inches"),  # Limit row name width
              column_names_max_height = unit(6, "cm"),  # More space for column names
              heatmap_legend_param = list(title = "Z-score", labels_gp = gpar(fontsize = 8)),
              width = ncol(mat_zscored) * unit(3, "mm"),  # Larger cell width
              height = nrow(mat_zscored) * unit(3, "mm"))  # Larger cell height

# Draw the heatmap
draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 20), "mm"))


# look at ans
an <- subset(bc, an_type == "11")
an <- an[, c("root_id", "fafb_cell_type", "manc_cell_type", "an_type")]
table(an$fafb_cell_type)
table(an$manc_cell_type)

# now connectivity --------------------------------------------------------


# pull in data
# more expansive
setwd("/Users/papers/BANC-project/data/connectivity_embedding/cell_sub_class/ascending")
if (dn == "YES") {
  setwd("/Users/papers/BANC-project/data/connectivity_embedding/cell_sub_class/descending")
}



# Get list of all CSV files
csv_files_connectivity <- list.files(pattern = "*.csv")

# pick specifics
# Read the CSV without row names
an_type_connectivity_matrix <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/cell_sub_class/ascending/ascending_joint_cosine_distances.csv")
# an_type_connectivity_matrix <- read_csv("/Users/papers/BANC-project/data/connectivity_embedding/cell_sub_class/descending/descending_joint_cosine_distances.csv")
# Set the first column as row names and remove the first column from the data frame
an_type_connectivity_matrix <- an_type_connectivity_matrix %>%
  column_to_rownames(var = "post_cell_sub_class")

# View the modified data
View(an_type_connectivity_matrix)

# cluster
dist_matrix <- dist(an_type_connectivity_matrix)
hc <- hclust(dist_matrix, method = "complete")  # You can change the method as needed

# Plot the dendrogram
plot(hc, main = "Dendrogram of DN Types, Connectivity", 
     xlab = "", sub = "", 
     hang = -1)  # hang = -1 makes the labels appear at the same level
library(pheatmap)

# Create a heatmap with dendrograms
pheatmap(an_type_connectivity_matrix,
         main = "Clustered Heatmap of Type-by-Type Matrix",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize = 8,
         cellwidth = 10,
         cellheight = 10,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")