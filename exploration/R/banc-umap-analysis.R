### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

library(bancr)
library(Matrix)
library(uwot)
library(tidyverse)

cts.chosen <- c("oviDNa_a","oviDNa_b")
nam <- paste(cts.chosen,collapse="_")

################
### GET DATA ###
################

el <- banc_edgelist()

franken.meta <- bancr::franken_meta()

bc.orig <- banctable_query() 
banc.meta <- bc.orig %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON",status)) %>%
  dplyr::mutate(cell_type = ifelse(grepl("auto\\:", cell_type),NA,cell_type),
                fafb_cell_type = ifelse(grepl("auto\\:",fafb_cell_type),NA,cell_type),
                manc_cell_type = ifelse(grepl("auto\\:",manc_cell_type),NA,cell_type),
                id=root_id) %>%
  dplyr::distinct(id, .keep_all = TRUE)

banc.edgelist.simple <- el %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id)) %>%
  dplyr::filter(pre!=post) %>%
  dplyr::rename(count = n) %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count)) %>%
  dplyr::mutate(norm = count/post_count) %>%
  dplyr::ungroup()

######################
### CHOOSE NEURONS ###
######################

# Get DN meta
banc.meta.chosen <- banc.meta %>%
  dplyr::filter(cell_type%in%cts.chosen)
ids.chosen <- unique(banc.meta.chosen$id)

#################
### INFLUENCE ###
#################

# Read influence data
csvs1 <- list.files(file.path("C:/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/Use_this/cell_function/"),
                    pattern="csv",
                    full.names = TRUE)
csvs2 <- list.files(file.path("C:/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/Use_this/Rachel_request/"),
                    pattern="csv",
                    full.names = TRUE)
csvs3 <- list.files(file.path("C:/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/Use_this/DN_type/"),
                    pattern="csv",
                    full.names = TRUE)
csvs4 <- list.files(file.path("C:/Users/papers/BANC-project/data/influence/frankenbrain_v1.5/Use_this/AN_type/"),
                    pattern="csv",
                    full.names = TRUE)
csvs <- c(csvs2,csvs4,csvs2,csvs1[grepl("olfactory|hygrosensory|thermosensory|visual|chemosensory|proprioceptive|tactile|visual_projection",csvs1)])
csvs <- csvs[!grepl("motor|unknown|mixed|chemosensory_tactile|proprioceptive_tactile|chemosensory_proprioceptive",csvs)]
influence.df <- data.frame()
for(csv in csvs){
  data <- readr::read_csv(csv, col_types = readr:::cols(.default = col_character()))
  data <- data[,c(2,3)]
  colnames(data) <- c("id","influence") # influence_norm_unsigned_forward_steady_state
  if(all(data$influence=="0")){
    warning("No values for: ", csv)
    next
  }
  data$seed <- gsub(".*from_|\\.csv*|_influence.*","",basename(csv))
  data$resolution <- "forward_coarse"
  data$influence <- as.numeric(data$influence)
  influence.df <- rbind(influence.df,data)
}

# Format
influence.meta <- influence.df %>%
  dplyr::group_by(seed, resolution) %>%
  dplyr::mutate(influence_log = 1-log(influence),
                influence_scaled = as.vector((influence - median(influence, na.rm = TRUE))/stats::mad(influence, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(banc.meta %>%
                     dplyr::select(id, side,
                                   region, super_class, hemilineage, 
                                   cell_function, nerve, cell_class, sez_class, cell_sub_class, 
                                   cell_type, top_nt, input_connections) %>%
                     dplyr::distinct(id, .keep_all = TRUE),
                   by = "id") %>%
  dplyr::filter(super_class!="glia", 
                !grepl("afferent|sensory",super_class),
                !is.na(super_class)) %>%
  dplyr::mutate(super_class = case_when(
    super_class == "brain_central_other" & !is.na(sez_class) ~ "sez",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(cell_type = ifelse(grepl("KCg-s",cell_type),"KCg-s",cell_type)) %>%
  dplyr::mutate(cell_sub_class = ifelse(is.na(cell_sub_class),gsub("_.*","",cell_type),cell_sub_class))

####################
### WRANGLE DATA ###
####################

# Join banc.meta data
banc.meta.pre <- banc.meta
colnames(banc.meta.pre) <- paste0("pre_",colnames(banc.meta.pre))
banc.meta.post <- banc.meta
colnames(banc.meta.post) <- paste0("post_",colnames(banc.meta.post))

# Make DN edgelist
elist.chosen.downstream <- banc.edgelist.simple %>%
  dplyr::filter(pre %in% ids.chosen) %>%
  dplyr::filter(norm >= 0.001, count > 10)
elist.chosen.upstream <- banc.edgelist.simple %>%
  dplyr::filter(post %in% ids.chosen) %>%
  dplyr::filter(norm >= 0.001, count > 10)

# Make DN edgelist
elist.chosen.pre <- banc.edgelist.simple %>%
  dplyr::filter(pre %in% elist.chosen.downstream$post) %>%
  dplyr::left_join(banc.meta.post %>% 
                     dplyr::select(post_id, 
                                   post_side, post_region, post_super_class, 
                                   post_hemilineage, post_cell_function, post_nerve, 
                                   post_cell_class, post_cell_sub_class, post_cell_type),
                   by = c("post"="post_id"))
elist.chosen.post <- banc.edgelist.simple %>%
  dplyr::filter(post %in% elist.chosen.upstream$pre) %>%
  dplyr::left_join(banc.meta.pre %>% 
                     dplyr::select(pre_id, 
                                   pre_side, pre_region, pre_super_class, 
                                   pre_hemilineage, pre_cell_function, pre_nerve, 
                                   pre_cell_class, pre_cell_sub_class, pre_cell_type),
                   by = c("pre"="pre_id"))

# Create a matrix of undirected connections
elist.chosen.cat <- rbind(elist.chosen.pre %>%
                        dplyr::mutate(id = pre, partner_id = paste0("post_",post)) %>%
                        dplyr::select(id, partner_id, count, norm),
                      elist.chosen.post %>%
                        dplyr::mutate(id = post, partner_id = paste0("pre_",pre)) %>%
                        dplyr::select(id, partner_id, count, norm)) %>%
  dplyr::select(id, partner_id, count, norm)

##################################################
### INPUT+OUTPUT INFLUENCE UMAP DATA BY NEURON ###
##################################################

# Function to calculate cosine similarity for sparse matrices
cosine_similarity_sparse <- function(mat) {
  # Calculate the norm of each column
  col_norms <- sqrt(colSums(mat^2))
  
  # Normalize the matrix
  mat_normalized <- mat %*% Diagonal(x = 1 / col_norms)
  
  # Calculate cosine similarity
  sim <- t(mat_normalized) %*% mat_normalized
  
  return(as.matrix(sim))
}

# Make matrix
inout_connection_matrix <- elist.chosen.cat %>%
  dplyr::mutate(weight = norm) %>%
  reshape2::dcast(partner_id ~ id,
                  fun.aggregate = mean,
                  value.var = "weight",
                  fill = 0)
rownames(inout_connection_matrix) <- inout_connection_matrix$partner_id
inout_connection_matrix$partner_id <- NULL

# Remove all-zero rows from the original matrix
non_zero_rows <- which(rowSums(abs(inout_connection_matrix)) > 0.0001)
inout_connection_matrix <- inout_connection_matrix[non_zero_rows, ]
non_zero_cols <- which(colSums(abs(inout_connection_matrix)) > 0.0001)
inout_connection_matrix <- inout_connection_matrix[,non_zero_cols]

# Calculate cosine similarity
sparsity <- sum(inout_connection_matrix == 0) / prod(dim(inout_connection_matrix))
print(paste("Sparsity:", sparsity))
sparse_matrix <- as(as.matrix(t(inout_connection_matrix)), "dgCMatrix")

# Calculate cosine similarity
undirected_cosine_sim_matrix <- cosine_similarity_sparse(t(sparse_matrix))
undirected_cosine_sim_matrix[is.infinite(undirected_cosine_sim_matrix)] <- 0
dimnames(undirected_cosine_sim_matrix) <- list(colnames(inout_connection_matrix),colnames(inout_connection_matrix))

# Represent as UMAP
set.seed(42)  
umap_result <- uwot::umap(t(undirected_cosine_sim_matrix),
                          metric = "cosine",
                          n_epochs = 500,
                          n_neighbors = 100, 
                          min_dist = 0,
                          n_trees = 100,
                          n_components = 2)

# Create a data frame with UMAP coordinates
umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  id = colnames(inout_connection_matrix)) %>% 
  dplyr::left_join(banc.meta %>%
                     dplyr::select(id, top_nt, 
                                   side, region, super_class, 
                                   hemilineage, cell_function, nerve, 
                                   cell_class, cell_sub_class, cell_type) %>%
                     dplyr::distinct(id, .keep_all = TRUE),
                   by = "id") %>%
  dplyr::left_join(influence.meta %>%
                      dplyr::select(cell_type,
                                    cell_function = cell_function) %>%
                      dplyr::distinct(cell_type, .keep_all = TRUE),
                    by = "cell_type") %>%
  dplyr::left_join(elist.chosen.downstream %>%
                     dplyr::select(post,output_count=count,input_norm=norm),
                   by=c("id"="post")) %>%
  dplyr::left_join(elist.chosen.upstream %>%
                     dplyr::select(pre,output_count=count,output_norm=norm),
                   by=c("id"="pre")) %>%
  dplyr::ungroup()

# Perform spectral clustering
sc <- kernlab::specc(x = as.matrix(umap_df[, c("UMAP1", "UMAP2")]), 
                     centers = 6)
umap_df$cluster <- factor(sc)

# Create a function to generate n colors
cerise_limon_base <- c("#EE5B32", "#F6B83C", "#4BA747", "#5BB6E4", "#7C378A")
cerise_limon_palette <- grDevices::colorRampPalette(cerise_limon_base)

# Ensure we have enough colors for all clusters
n_clusters <- length(unique(umap_df$cluster))
cluster_colors <- cerise_limon_palette(n_clusters)
names(cluster_colors) <- sort(unique(umap_df$cluster))
umap_df$colours <- cluster_colors[umap_df$cluster]

# Calculate cluster centroids
cluster_centroids <- umap_df %>%
  group_by(cluster) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1),
            UMAP2 = mean(UMAP2))

# If you want to add density contours
p_clusters <- ggplot(umap_df, aes(x = UMAP1, 
                                  y = UMAP2, 
                                  color = cluster)) +
  geom_density_2d(col="grey70", alpha = 0.5) +
  geom_point(aes(size = input_connections), alpha = 0.7, size = 2) +
  geom_text(data = cluster_centroids, 
            aes(label = cluster),
            colour = "black",
            size = 6, 
            fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  labs(title = "",
       x = "UMAP1",
       y = "UMAP2") +
  theme(
    legend.position = "bottom"
  )  +
  guides(color = "none",
         fill = guide_legend(title = "function"))

# If you want to add density contours
p_input <- ggplot(umap_df, aes(x = UMAP1, 
                                  y = UMAP2, 
                                  color = input_norm)) +
  geom_density_2d(col="grey70", alpha = 0.5) +
  geom_point(aes(size = input_connections), alpha = 0.7, size = 2) +
  geom_text(data = cluster_centroids, 
            aes(label = cluster),
            colour = "black",
            size = 6, 
            fontface = "bold") +
  viridis::scale_colour_viridis() +
  theme_void() +
  labs(title = "",
       x = "UMAP1",
       y = "UMAP2") +
  theme(
    legend.position = "bottom"
  )

# If you want to add density contours
p_output <- ggplot(umap_df, aes(x = UMAP1, 
                               y = UMAP2, 
                               color = output_norm)) +
  geom_density_2d(col="grey70", alpha = 0.5) +
  geom_point(aes(size = input_connections), alpha = 0.7, size = 2) +
  geom_text(data = cluster_centroids, 
            aes(label = cluster),
            colour = "black",
            size = 6, 
            fontface = "bold") +
  viridis::scale_colour_viridis() +
  theme_void() +
  labs(title = "",
       x = "UMAP1",
       y = "UMAP2") +
  theme(
    legend.position = "bottom"
  )

# Save
ggsave(plot = p_clusters,
       filename = file.path("images",paste0(nam,"_inout_connectivity_umap_density.png")),
       width = 10, height = 10, dpi = 300)