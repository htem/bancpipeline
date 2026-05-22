### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# load influence scores

######################
### load libraries ###
######################

library(tidyverse)
library(progress)
library(dplyr)
library(pheatmap)

#############################
### load influence scores ###
#############################

# Read influence data
banc.path <- "C:/Users/papers/BANC-project/"

# get all csv files from selected folders
csvs1 <- list.files(file.path(banc.path, "data/influence/frankenbrain_v1.5/Use_this/Modality/"),
                    pattern="csv",
                    full.names = TRUE)
csvs2 <- list.files(file.path(banc.path, "data/influence/frankenbrain_v1.5/Use_this/Rachel_request/"),
                    pattern="csv",
                    full.names = TRUE)
csvs3 <- list.files(file.path(banc.path, "data/influence/frankenbrain_v1.5/Use_this/Seed/"),
                    pattern="csv",
                    full.names = TRUE)

csvs <- c(csvs2,
          csvs1[grepl("olfactory|hygrosensory|thermosensory|visual|chemosensory|proprioceptive|tactile|visual_projection",csvs1)],
          csvs3[grepl("visual_projection",csvs3)])
csvs <- csvs[!grepl("motor|unknown|mixed|chemosensory_tactile|proprioceptive_tactile|chemosensory_proprioceptive",csvs)]

# Create progress bar
pb <- progress_bar$new(total = length(csvs), format = "[:bar] :percent ETA: :eta")

# extract scores from each csv
influence.list <- list()
for(csv in csvs){
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = col_character()))
  data <- data[, c(2, 3)]
  colnames(data) <- c("id", "influence") # influence_norm_unsigned_forward_steady_state
  if(all(data$influence == "0")){
    warning("No values for: ", csv)
    pb$tick()
    next
  }
  data$seed <- gsub(".*from_|\\.csv*|_influence.*","",basename(csv))
  data$resolution <- "forward_coarse"
  data$influence <- as.numeric(data$influence)
  const <- abs(min(data$influence[data$influence != 0], na.rm = T)/10)
  data$influence_log <- 1-log(const + data$influence)
  data$influence_scaled <- as.vector((data$influence_log - median(data$influence_log, na.rm = TRUE))/stats::mad(data$influence_log, na.rm = TRUE))
  influence.list[[length(influence.list) + 1]] <- data
  pb$tick()
}
influence.df <- do.call(rbind, influence.list)

##########################################################
### load metadata and add metadata to influence scores ###
##########################################################

# get franken metadata
franken.meta <- franken_meta()
franken.meta <- franken.meta %>%
  dplyr::mutate(id = neuron_id)
# fix region name
franken.meta$region[grepl("o,p,t,i,c,_,l,o,b,e,s", franken.meta$region)] <- "optic_lobes"
# fix cell_function name
franken.meta$cell_function[grepl("o,c,e,l,l,a,r", franken.meta$cell_function)] <- "ocellar"
franken.meta$cell_function[grepl("o,c,e,l,l,a,r", franken.meta$cell_function)] <- "ocellar"
franken.meta$cell_function[grepl("v,i,s,u,a,l,_,c,h,r,o,m,a,t,i,c", franken.meta$cell_function)] <- "visual_chromatic"
franken.meta$cell_function[grepl("v,i,s,u,a,l,_,a,c,h,r,o,m,a,t,i,c", franken.meta$cell_function)] <- "visual_achromatic"

# get banc metadata
bc.orig <- banctable_query()
banc.meta <- bc.orig %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON", status)) %>%
  dplyr::mutate(cell_type = ifelse(grepl("auto\\:", cell_type), NA, cell_type),
                fafb_cell_type = ifelse(grepl("auto\\:", fafb_cell_type), NA, cell_type),
                manc_cell_type = ifelse(grepl("auto\\:", manc_cell_type), NA, cell_type),
                id = root_id) %>%
  dplyr::distinct(id, .keep_all = TRUE)

# Make information vector
influence.meta <- left_join(influence.df, franken.meta, by = "id")

#############################
### load influence scores ###
#############################

cts <- c("oviDNa_a", "oviDNa_b", "oviDNb")

resolution_type = "forward_coarse"

influence_matrix <- influence.meta %>%
  dplyr::filter(cell_type %in% cts) %>%
  # choose particular resolution type
  dplyr::filter(resolution == resolution_type,
                !is.na(influence)) %>%
  # reshape dataframe
  reshape2::dcast(seed.x ~ cell_type,
                  fun.aggregate = median,
                  value.var = "influence_scaled",
                  fill = 10) %>%
  dplyr::rename(seed = seed.x)

# Set row names and remove the seed column
rownames(influence_matrix) <- influence_matrix$seed
influence_matrix$seed <- NULL

# Convert to matrix
influence_matrix <- as.matrix(influence_matrix)

# Region annotations
region_annotation <- franken.meta %>%
  dplyr::filter(!is.na(cell_function)) %>%
  dplyr::select(cell_function, region) %>%
  dplyr::distinct(cell_function, .keep_all = TRUE) %>%
  as.data.frame()
rownames(region_annotation) <- region_annotation$cell_function
region_annotation <- region_annotation[rownames(region_annotation) %in% rownames(influence_matrix),]
region_annotation$cell_function <- NULL

# Get cell types for each pre_root_id
col_annotation <- influence.meta %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::select(cell_type, super_class) %>%
  dplyr::filter(cell_type %in% cts) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  as.data.frame()
rownames(col_annotation) <- as.character(col_annotation$cell_type)
col_annotation <- col_annotation[rownames(col_annotation) %in% colnames(influence_matrix),]
col_annotation <- data.frame(cell_type = as.character(col_annotation$cell_type))
rownames(col_annotation) <- col_annotation$cell_type

# Create the heatmap
ph4 <- pheatmap(
  influence_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  #color = scaled_heatmap_palette4,
  #breaks = scaled_heatmap_breaks4,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "",
  fontsize = 8,
  cellwidth = 10,
  cellheight = 10,
  annotation_row = region_annotation,
  annotation_col = col_annotation,
  annotation_legend = FALSE,
  legend = TRUE,
  width = 6,
  height = 6,
  treeheight_row = 5, 
  treeheight_col = 5
)

# Plot
print(ph4)