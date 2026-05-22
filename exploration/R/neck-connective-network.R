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
library(ggraph)
library(igraph)
library(dplyr)
library(tidygraph)
library(plotly)

# Make sure all functions query BANC and not FAFB
choose_banc()

################
### get data ###
################ # can skip, example data saved as .csv in repo

# # Get meta data
# bc <- banctable_query()
# meta <- bc %>%
#   dplyr::filter(region%in%c("neck_connective")) %>%
#   dplyr::mutate(cell_type = dplyr::case_when(
#     is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
#     is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
#     TRUE ~ cell_type
#   ),
#   top_nt = gsub("auto\\:","",top_nt),
#   cell_sub_class = gsub("auto\\:","",cell_sub_class),
#   cell_class = gsub("auto\\:","",cell_class),
#   super_class = gsub("auto\\:","",super_class)) %>%
#   dplyr::filter(!is.na(cell_type),
#                 super_class %in% c("ascending","descending","afferent, ascending"))
#   
# # Update IDs    
# meta$root_id <- banc_latestid(meta$root_id)
# ids <- unique(meta$root_id)
# ids <- ids[!is.na(ids)]
# ids <- ids[ids!="0"]
# ids <- unique(meta$root_id)
# 
# # Get connectivity
# up.list <- pbapply::pblapply(ids, function(id) {
#   try({
#     banc_partner_summary(id, partners = "inputs")
#   })
# })
# up <- do.call(rbind,up.list)
# 
# # Build edgelist
# neuron_edgelist <- up %>%
#   dplyr::rename(post_id = query) %>%
#   dplyr::filter(pre_id %in% ids, 
#                 post_id %in% ids) %>%
#   dplyr::left_join(meta[,c("root_id","cell_type","super_class","side","top_nt")],by = c("post_id"="root_id")) %>%
#   dplyr::rename(post_cell_type = cell_type, post_super_class = super_class,
#                 post_side = side, post_top_nt = top_nt) %>%
#   dplyr::left_join(meta[,c("root_id","cell_type","super_class","side","top_nt")],by = c("pre_id"="root_id")) %>%
#   dplyr::rename(pre_cell_type = cell_type, pre_super_class = super_class,
#                 pre_side = side, pre_top_nt = top_nt) %>%
#   dplyr::filter(!is.na(post_cell_type)|!is.na(pre_cell_type)) %>%
#   dplyr::group_by(post_cell_type,pre_super_class) %>%
#   dplyr::mutate(weight = sum(as.numeric(weight),na.rm=TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(pre_super_class %in% c("ascending","descending","afferent, ascending"),
#                 post_super_class %in% c("ascending","descending","afferent, ascending")) %>%
#   dplyr::filter(post_cell_type!=pre_cell_type)
# write_csv(neuron_edgelist,file="data/banc_neck_connectivity_edgelist.csv")

#####################
### organise data ###
#####################

# Read fetched data
neuron_edgelist <- readr::read_csv(file="data/banc_neck_connectivity_edgelist.csv",
                                   col_types = hemibrainr:::sql_col_types )

# Collapse the edgelist by cell type
collapsed_edgelist <- neuron_edgelist %>%
  group_by(pre_cell_type, post_cell_type, pre_super_class, post_super_class) %>%
  summarise(weight = sum(weight), .groups = 'drop') %>%
  filter(pre_cell_type != post_cell_type)  # Remove self-loops if desired

# Create a unique vertex list
vertices <- bind_rows(
  collapsed_edgelist %>% select(cell_type = pre_cell_type, super_class = pre_super_class),
  collapsed_edgelist %>% select(cell_type = post_cell_type, super_class = post_super_class)
) %>% 
  distinct(cell_type, .keep_all = TRUE)

##########################
### community overview ###
##########################

# Create the graph object
g <- graph_from_data_frame(
  d = collapsed_edgelist %>% select(pre_cell_type, post_cell_type, weight),
  directed = TRUE,
  vertices = vertices
)

# Perform community detection using Infomap
communities <- cluster_infomap(g)

# Convert to tbl_graph, calculate degree, and add module information
g <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(mode = "total"),
    module = as.factor(communities$membership)
  )

# Calculate node layout using a hierarchical layout algorithm
layout <- create_layout(g, layout = "sugiyama", circular = T)

# Create the plot with hierarchical layout
p <- ggraph(layout) +
  geom_edge_link(aes(width = weight, alpha = weight),
                 arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(3, 'mm'),
                 edge_colour = "grey50") +
  geom_node_point(aes(color = super_class, size = degree), alpha = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 1, max.overlaps = 15) +
  scale_edge_width(range = c(0.1, 2)) +
  scale_edge_alpha(range = c(0.1, 0.1)) +
  scale_size_continuous(range = c(0.5, 5)) +
  scale_color_brewer(palette = "Set2") +
  scale_y_reverse() +
  facet_nodes(~module, scales = "free_x", ncol = 4) +  # Adjust ncol as needed
  labs(title = "Hierarchical Cell Type Connectivity by Infomap Module",
       color = "Super Class",
       size = "Degree",
       edge_width = "Total Weight",
       edge_alpha = "Total Weight") +
  theme_graph() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.spacing = unit(0.5, "cm")
  ) +
  coord_flip()

# Show the plot
print(p)

# save
ggsave(plot = p,
       filename = "inst/images/_banc_neck_connective_graph.png", 
       width = 8, height = 8, dpi = 300)

##############################
### sub-community overview ###
##############################

# Create the graph object
g <- graph_from_data_frame(
  d = collapsed_edgelist %>% select(pre_cell_type, post_cell_type, weight),
  directed = TRUE,
  vertices = vertices
)

# First level of community detection using Infomap
communities_level1 <- cluster_infomap(g)

# Function to perform second level of community detection
detect_subcommunities <- function(graph, community) {
  subgraph <- induced_subgraph(graph, which(communities_level1$membership == community))
  if (vcount(subgraph) > 3) {  # Only detect subcommunities for communities with more than 3 nodes
    subcommunities <- cluster_infomap(subgraph)
    return(subcommunities$membership)
  } else {
    return(rep(1, vcount(subgraph)))
  }
}

# Perform second level of community detection
subcommunities <- lapply(unique(communities_level1$membership), function(comm) {
  detect_subcommunities(g, comm)
})

# Combine community and subcommunity information
node_communities <- data.frame(
  name = V(g)$name,
  community = communities_level1$membership,
  subcommunity = unlist(subcommunities)
)

# Convert to tbl_graph, add community information
g <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  left_join(node_communities, by = "name") %>%
  mutate(
    degree = centrality_degree(mode = "total"),
    community = as.factor(community),
    subcommunity = as.factor(paste(community, subcommunity, sep = "_"))
  )

# # Plot only the most prominent communities
# top_communities <- as.numeric(names(sort(table(V(g)$community), decreasing = TRUE)[1:5]))
# g <- g %>%
#   activate(nodes) %>%
#   filter(as.numeric(as.character(community)) %in% top_communities)

# Calculate node layout using a hierarchical layout algorithm
layout <- create_layout(g, layout = "sugiyama")

# Create the plot with hierarchical layout
p <- ggraph(layout) +
  geom_edge_link(aes(width = weight, alpha = weight),
                 arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(3, 'mm'),
                 edge_colour = "grey50") +
  geom_node_point(aes(color = super_class, size = degree)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2, max.overlaps = 10) +
  scale_edge_width(range = c(0.1, 1)) +
  scale_edge_alpha(range = c(0.1, 0.3)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_brewer(palette = "Set2") +
  scale_y_reverse() +
  facet_nodes(community ~ subcommunity, scales = "free") +
  labs(title = "Hierarchical Cell Type Connectivity by Infomap Communities",
       subtitle = "Main communities (rows) and sub-communities (columns)",
       color = "Super Class",
       size = "Degree",
       edge_width = "Total Weight",
       edge_alpha = "Total Weight") +
  theme_graph() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.spacing = unit(0.2, "cm"),
    strip.text = element_text(size = 8)
  ) +
  coord_flip()

# Show the plot
print(p)

# Function to create plot for a single community
plot_community <- function(comm_graph, comm_number) {
  ggraph(comm_graph, layout = "sugiyama") +
    geom_edge_link(aes(width = weight, alpha = weight),
                   arrow = arrow(length = unit(2, 'mm')), 
                   end_cap = circle(3, 'mm'),
                   edge_colour = "grey50") +
    geom_node_point(aes(color = super_class, size = degree)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 15) +
    scale_edge_width(range = c(0.1, 1)) +
    scale_edge_alpha(range = c(0.1, 0.3)) +
    scale_size_continuous(range = c(2, 7)) +
    scale_color_brewer(palette = "Set2") +
    facet_nodes(~subcommunity, scales = "free") +
    labs(title = paste("Community", comm_number),
         color = "Super Class",
         size = "Degree") +
    theme_graph() +
    theme(legend.position = "right")
}

# Create and save plots for each community
comdir <- "inst/images/communities"
dir.create(comdir, showWarnings = FALSE)
for (comm in unique(V(g)$community)) {
  comm_graph <- g %>% 
    activate(nodes) %>% 
    filter(community == comm) %>%
    mutate(subcommunity = as.factor(subcommunity))
  p <- plot_community(comm_graph, comm)
  plot(p)
  ggsave(file.path(comdir,paste0("community_", comm, ".png")), p, width = 12, height = 8)
}

########################
### interactive plot ###
########################