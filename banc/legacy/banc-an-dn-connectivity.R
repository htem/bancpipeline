######################
### load libraries ###
######################

# Load packages
source("banc/banc-startup.R")

# Make sure all functions query BANC and not FAFB
choose_banc()

################
### get data ###
################

# Get meta data
bc <- banctable_query() %>%
  dplyr::select(-starts_with("_"))
meta <- bc %>%
  dplyr::filter(!grepl("DELETE",status),!grepl("NOT_A_NEURON",status),!grepl("GLIA",status)) %>%
  dplyr::filter(grepl("neck",region)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
    is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
    TRUE ~ cell_type
  ),
  top_nt = gsub("auto\\:","",top_nt),
  cell_sub_class = gsub("auto\\:","",cell_sub_class),
  cell_class = gsub("auto\\:","",cell_class),
  super_class = gsub("auto\\:","",super_class)) %>%
  dplyr::filter(super_class %in% c("ascending","descending","afferent, ascending","sensory_ascending","sensory_descending")|
                  grepl("ascending, sensory|ascending|descending|sensory_ascending",cell_class)|
                  super_class %in% c("AN","DN")|grepl("AN|DN",cell_class)) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    !is.na(an_type) ~ paste0("AN_",an_type),
    !is.na(dn_type) ~ paste0("DN_",dn_type),
    TRUE ~ cell_sub_class
  )
  ) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    grepl("sensory",super_class)|grepl("sensory",cell_class) ~ "sensory_ascending",
    grepl("ascending|^AN",super_class)|grepl("ascending",cell_class) ~ "ascending",
    grepl("descending|^DN",super_class)|grepl("descending|^DN",cell_class) ~ "descending",
    TRUE ~ super_class
  )
  ) %>%
  dplyr::arrange(cell_type) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Add in some other data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
top.nts <- fw.meta %>%
  dplyr::select(cell_type, top_nt) %>%
  rbind(mc.meta %>% dplyr::select(cell_type, top_nt)) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::count(top_nt) %>%
  dplyr::slice_max(n, with_ties = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup()
meta <- meta %>%
  dplyr::select(-top_nt) %>%
  dplyr::left_join(top.nts, by = "cell_type")

# Update IDs
meta <- banc_updateids(meta)
ids <- unique(meta$root_id)
ids <- ids[!is.na(ids)]
ids <- ids[ids!="0"]
ids <- unique(meta$root_id)

# Save folders
split.master.folder <- banc.l2split.save.path
split.folder <- file.path(split.master.folder,"swc")
images.folder <- file.path(split.master.folder,"images")
synapses.folder <- file.path(split.master.folder,"synapses")
dones.syns <- list.files(synapses.folder)
done.ids <- gsub("\\.csv","",dones.syns)
missing.ids <- setdiff(ids,done.ids)
message("Missing synapse data for ", length(missing.ids))

# # Get connectivity directly
# up.list <- pbapply::pblapply(ids, function(id) {
#    tryCatch({
#      df <- banc_partner_summary(id, partners = "inputs")
#      df$query <- id
#      df
#    }, error = function(e){
#      Sys.sleep(1)
#      try({
#        df <- banc_partner_summary(id, partners = "inputs")
#        df$query <- id
#        df
#      })
#    } 
#    ) 
#  })
# up <- do.call(rbind,up.list)

# # Build edgelist
# neuron_edgelist <- up %>%
#   dplyr::mutate(weight = as.numeric(weight)) %>%
#   dplyr::rename(post_id = query) %>%
#   dplyr::group_by(post_id) %>%
#   dplyr::mutate(post_count = sum(weight, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(pre_id %in% ids,
#                 post_id %in% ids) %>%
#   dplyr::left_join(meta[,c("root_id","cell_type","super_class","side","top_nt")],by = c("post_id"="root_id")) %>%
#   dplyr::rename(post_cell_type = cell_type, 
#                 post_super_class = super_class,
#                 post_side = side, 
#                 post_top_nt = top_nt) %>%
#   dplyr::left_join(meta[,c("root_id","cell_type","super_class","side","top_nt")],by = c("pre_id"="root_id")) %>%
#   dplyr::rename(pre_cell_type = cell_type, 
#                 pre_super_class = super_class,
#                 pre_side = side, 
#                 pre_top_nt = top_nt) %>%
#   dplyr::filter(!is.na(post_cell_type)|!is.na(pre_cell_type)) %>%
#   dplyr::group_by(post_id, pre_id) %>%
#   dplyr::mutate(count = sum(as.numeric(weight),na.rm=TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(norm = count/post_count) %>%
#   dplyr::filter(pre_super_class %in% c("ascending","descending","afferent, ascending"),
#                 post_super_class %in% c("ascending","descending","afferent, ascending")) %>%
#   dplyr::filter(post_cell_type!=pre_cell_type) %>%
#   dplyr::select(-count)

# Assemble edgelist
split.master.folder <- banc.l2split.save.path
synapses.folder <- file.path(split.master.folder,"synapses")
csvs.all <- list.files(synapses.folder, full.names = TRUE)
csvs.ids <- gsub("\\.csv","",basename(csvs.all))
# csvs.ids <- banc_updateids(csvs.ids)
csvs <- csvs.all[csvs.ids%in%ids]

# Get synapses
synapses <- read_synapse_csvs(csvs, ids = NULL)

# Label lookup
lookup <- synapses %>%
  dplyr::filter(prepost==0)  %>%
  dplyr::distinct(connector_id, label = hemibrainr:::standard_compartments(Label)) %>%
  dplyr::rename(pre_label = label) %>% 
  dplyr::ungroup()

# Update IDs
post.synapses <- synapses %>%
  dplyr::filter(prepost==1) %>%
  dplyr::filter(post_id %in% ids) %>%
  dplyr::mutate(post_label = hemibrainr:::standard_compartments(Label)) %>%
  dplyr::left_join(lookup, by = "connector_id") %>%
  dplyr::select(-Label)
post.synapses$pre_id <- bancr::banc_rootid(post.synapses$pre_svid)
post.synapses <- post.synapses %>%
  dplyr::filter(pre_id %in% ids, 
                post_id %in% ids)

# Add neuropils
chunk_size <- 100000
synapses.np <- post.synapses %>%
  dplyr::mutate(neuropil = NA,
                region = NA,
                side = NA,
                chunk = ceiling(dplyr::row_number() / chunk_size)) %>%
  dplyr::group_by(chunk) %>%
  dplyr::group_split()
pb <- progress::progress_bar$new(
  format = "  Processing [:bar] :percent | Elapsed: :elapsed | ETA: :eta | :current/:total chunks",
  total = length(synapses.np),
  clear = FALSE,
  width = 100
)
synapses.np <- purrr::map(synapses.np, function(chunk) {
  result <- pointsinside_banc(chunk)
  result$chunk <- NULL
  pb$tick()
  return(result)
})
synapses.np <- do.call(rbind,synapses.np)

# Run
elist.raw <- read_elist_csvs(csvs, 
                             split = FALSE, 
                             ids = NULL)
elist.raw$pre <- bancr::banc_rootid(elist.raw$pre_svid)
elist <- elist.raw %>%
  dplyr::filter(pre %in% ids, 
                post %in% ids) %>%
  dplyr::arrange(dplyr::desc(count), 
                 dplyr::desc(norm)) %>%
  dplyr::left_join(meta[,c("root_id","top_nt", "side","cell_type","cell_sub_class","cell_class","super_class")], 
                   by = c("post"="root_id")) %>%
  dplyr::rename(post_top_nt = top_nt,
                post_side = side,
                post_cell_type = cell_type,
                post_sub_class = cell_sub_class,
                post_cell_class = cell_class,
                post_super_class = super_class) %>%
  dplyr::left_join(meta[,c("root_id","top_nt", "side", "cell_type","cell_sub_class","cell_class","super_class")], 
                   by = c("pre"="root_id")) %>%
  dplyr::rename(pre_top_nt = top_nt,
                pre_side = side,
                pre_cell_type = cell_type,
                pre_sub_class = cell_sub_class,
                pre_cell_class = cell_class,
                pre_super_class = super_class)

# Split elist
elist.ad.raw <- read_elist_csvs(csvs, 
                                split = TRUE, 
                                lookup = lookup,
                                ids = NULL)
elist.raw$pre <- bancr::banc_rootid(elist.raw$pre_svid)
elist.ad <- elist.ad.raw %>%
  dplyr::filter(pre %in% ids, 
                post %in% ids) %>%
  dplyr::arrange(dplyr::desc(count), 
                 dplyr::desc(norm)) %>%
  dplyr::left_join(meta[,c("root_id","top_nt", "side","cell_type","cell_sub_class","cell_class","super_class")], by = c("post"="root_id")) %>%
  dplyr::rename(post_top_nt = top_nt,
                post_side = side,
                post_cell_type = cell_type,
                post_sub_class = cell_sub_class,
                post_cell_class = cell_class,
                post_super_class = super_class) %>%
  dplyr::left_join(meta[,c("root_id","top_nt", "side", "cell_type","cell_sub_class","cell_class","super_class")], by = c("pre"="root_id")) %>%
  dplyr::rename(pre_top_nt = top_nt,
                pre_side = side,
                pre_cell_type = cell_type,
                pre_sub_class = cell_sub_class,
                pre_cell_class = cell_class,
                pre_super_class = super_class)

# Write
arrow::write_feather(elist,
                     file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_edgelist_simple.feather"))
arrow::write_feather(elist.ad,
                     file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_edgelist.feather"))
write_connectome_data(synapses.np,
                      file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_synapses.parquet"),
                      format = "parquet")
arrow::write_feather(meta,
                     file.path(banc.connectivity.save.path, "banc_v1.2_meta.feather"))

#######################
#### Synapses check ###
#######################

# Get the data
syns <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_synapses.parquet"))

syns.missed <- syns %>%
  dplyr::group_by(connector_id) %>%
  dplyr::mutate(connector_id_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(connector_id_count==1)
connector_id.missed <- syns$connector_id[!duplicated(syns$connector_id)]

##########
### QC ###
##########

# Get connectivity from CAVE
banc.el.orig <- banc_edgelist(fetch_all_rows = TRUE)
banc.el.1 <- banc.el.orig %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id))  %>%
  dplyr::filter(pre %in% ids, 
                post %in% ids) %>%
  dplyr::rename(count = n) %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(norm = round(count/post_count,6)) %>%
  dplyr::filter(pre_pt_root_id!=post_pt_root_id) %>%
  dplyr::arrange(dplyr::desc(count))
shared.ids <- intersect(ids, c(banc.el.1$pre,banc.el.1$post))

banc.el.2 <- elist %>%
  dplyr::filter(pre_pt_root_id %in% shared.ids, 
                post_pt_root_id %in% shared.ids) %>%
  dplyr::arrange(dplyr::desc(count))

banc.el.comp <- banc.el.2 %>%
  dplyr::left_join(banc.el.1 %>%
                     dplyr::rename(count = count_1, norm = norm_1) %>%
                     dplyr::select(pre, post, count_1, norm_1),
                   by = c("pre","post"))

################
### Analysis ###
################
library(ggraph)
library(igraph)
library(tidygraph)

# Get the data
elist <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_edgelist_simple.feather")) %>%
  dplyr::filter(count >= 5)

# Collapse the edgelist by cell type
collapsed_edgelist <- elist %>%
  dplyr::group_by(pre_cell_type, post_cell_type, pre_super_class, post_super_class, post_top_nt, pre_top_nt, pre_side, post_side) %>%
  dplyr::summarise(weight = sum(count), .groups = 'drop') %>%
  dplyr::filter(pre_cell_type != post_cell_type,
                weight >= 25)

# Create a unique vertex list
vertices <- bind_rows(
  collapsed_edgelist %>% select(cell_type = pre_cell_type, super_class = pre_super_class, top_nt = pre_top_nt),
  collapsed_edgelist %>% select(cell_type = post_cell_type, super_class = post_super_class, top_nt = post_top_nt)) %>% 
  dplyr::distinct(cell_type, .keep_all = TRUE)

# Create the graph object
g <- graph_from_data_frame(
  d = collapsed_edgelist %>% dplyr::select(pre_cell_type, post_cell_type, weight),
  directed = TRUE,
  vertices = vertices
)

# Perform community detection using Infomap
communities <- cluster_infomap(g)

# Convert to tbl_graph, calculate degree, and add module information
g <- as_tbl_graph(g) %>%
  tidygraph::activate(nodes) %>%
  dplyr::mutate(
    degree = centrality_degree(mode = "total"),
    module = as.factor(communities$membership)
  )

# Calculate node layout using a hierarchical layout algorithm
layout <- create_layout(g, layout = "sugiyama", circular = T)

# Create the plot with hierarchical layout
p <- ggraph(layout) +
  geom_edge_link(aes(width = log10(weight)),
                 arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(3, 'mm'),
                 edge_colour = "grey70") +
  geom_node_point(aes(fill = super_class, color = top_nt, size = degree), shape = 21, alpha = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 1, max.overlaps = 15) +
  scale_edge_width(range = c(0.1, 2)) +
  scale_edge_alpha(range = c(0.1, 0.1)) +
  scale_size_continuous(range = c(0.5, 5)) +
  scale_color_manual(values = paper.cols) +
  scale_fill_manual(values = paper.cols) +
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
       filename = "inst/images/banc_neck_connective_graph.png", 
       width = 8, 
       height = 8, 
       dpi = 300)

##############################
### sub-community overview ###
##############################

# Create the graph object
g <- graph_from_data_frame(
  d = collapsed_edgelist %>% dplyr::select(pre_cell_type, post_cell_type, weight),
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
  dplyr::left_join(node_communities, by = "name") %>%
  dplyr::mutate(
    degree = centrality_degree(mode = "total"),
    community = as.factor(community),
    subcommunity = as.factor(paste(community, subcommunity, sep = "_"))
  )

# Calculate node layout using a hierarchical layout algorithm
layout <- create_layout(g, layout = "sugiyama")

# # Create the plot with hierarchical layout
# p <- ggraph(layout) +
#   geom_edge_link(aes(width = log10(weight), 
#                      alpha = weight),
#                  arrow = arrow(length = unit(2, 'mm')), 
#                  end_cap = circle(3, 'mm'),
#                  edge_colour = "grey70") +
#   geom_node_point(aes(fill = super_class, color = top_nt, size = degree), shape = 21, alpha = 0.75) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 2, max.overlaps = 10) +
#   scale_edge_width(range = c(0.1, 1)) +
#   scale_edge_alpha(range = c(0.1, 0.3)) +
#   scale_size_continuous(range = c(1, 5)) +
#   scale_color_manual(values = paper.cols) +
#   scale_fill_manual(values = paper.cols) +
#   scale_y_reverse() +
#   facet_nodes(community ~ subcommunity, scales = "free") +
#   labs(title = "Hierarchical Cell Type Connectivity by Infomap Communities",
#        subtitle = "Main communities (rows) and sub-communities (columns)",
#        color = "Super Class",
#        size = "Degree",
#        edge_width = "Total Weight",
#        edge_alpha = "Total Weight") +
#   theme_graph() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, size = 12),
#     panel.spacing = unit(0.2, "cm"),
#     strip.text = element_text(size = 8)
#   ) +
#   coord_flip()
# 
# # Show the plot
# print(p)

# Function to create plot for a single community
plot_community <- function(layout.sc, comm_number) {
  ggraph(layout.sc) +
    geom_edge_link(aes(width = log10(weight), 
                       alpha = 0.95),
                   arrow = arrow(length = unit(2, 'mm')), 
                   end_cap = circle(3, 'mm'),
                   edge_colour = "grey70") +
    geom_node_point(aes(fill = super_class, color = top_nt, size = degree), shape = 21, alpha = 0.75) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 15) +
    scale_edge_width(range = c(0.1, 1)) +
    scale_edge_alpha(range = c(0.1, 0.3)) +
    scale_size_continuous(range = c(2, 7)) +
    scale_color_manual(values = paper.cols) +
    scale_fill_manual(values = paper.cols) +
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

# Get communities with more than one node
communities_to_plot <- unique(V(g)$community[tabulate(V(g)$community) > 1])
for (comm in unique(communities_to_plot)) {
  comm_graph <- g %>% 
    activate(nodes) %>% 
    filter(community == comm) %>%
    mutate(subcommunity = as.factor(subcommunity))
  
  # Check if the community has more than one node
  if (gorder(comm_graph) > 1) {
    layout.sc <- create_layout(comm_graph, layout = "sugiyama")
    layout.sc$top_nt[is.na(layout.sc$top_nt)] <- "unclear"
    p <- plot_community(layout.sc, comm)
    plot(p)
    ggsave(file.path(comdir, paste0("community_", comm, ".png")), p, width = 12, height = 8)
  }
}

############################
#### Axon-dendrite split ###
############################
library(ggsankey)

# Get the data
elist.ad <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "banc_v1.2_an_dn_edgelist.feather")) %>%
  dplyr::filter(count >= 5)

# Create separate dataframes for each step
collapsed_edgelist_broken <- elist.ad %>%
  dplyr::group_by(pre_cell_type, post_cell_type, pre_super_class, post_super_class, post_top_nt, pre_top_nt, post_label, pre_label, pre_side, post_side) %>%
  dplyr::summarise(weight = sum(count), .groups = 'drop') %>%
  dplyr::filter(pre_cell_type != post_cell_type,
                weight >= 15) %>% # Remove self-loops if desired
  dplyr::filter(post_label %in% c("axon","dendrite"),
                pre_label %in% c("axon","dendrite")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(connection = paste0(pre_label,"_",post_label)) %>%
  dplyr::filter(!is.na(pre_top_nt),
                !is.na(post_top_nt),
                pre_top_nt%in%c("acetylcholine","glutamate","gaba"),
                post_top_nt%in%c("acetylcholine","glutamate","gaba")) %>%
  dplyr::mutate(sidedness = dplyr::case_when(
    pre_side=="right"&post_side=="right" ~ "ipsilateral",
    pre_side=="left"&post_side=="left" ~ "ipsilateral",
    pre_side=="right"&post_side=="left" ~ "contralateral",
    pre_side=="left"&post_side=="right" ~ "contralateral",
    TRUE ~ "ipsilateral"
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(top_nt = paste0(pre_super_class,"_",pre_top_nt),
                connection_type = paste0(post_super_class,"_",connection),#paste0(connection,"_",pre_top_nt),
                target = paste0(sidedness,"_",post_super_class))#paste0(connection,"_",post_top_nt))

# Create separate dataframes for each step
df1 <- collapsed_edgelist_broken %>%
  select(pre_super_class, weight) %>%
  rename(node = pre_super_class)
df2 <- collapsed_edgelist_broken %>%
  select(pre_super_class, top_nt, weight) %>%
  rename(node = top_nt, next_node = pre_super_class)
df3 <- collapsed_edgelist_broken %>%
  select(top_nt, connection_type, weight) %>%
  rename(node = connection_type, next_node = top_nt)
df4 <- collapsed_edgelist_broken %>%
  select(connection_type, target, weight) %>%
  rename(node = target, next_node = connection_type)

# Combine the dataframes
sankey_data <- bind_rows(
  mutate(df1, step = "pre_super_class", next_step = "top_nt"),
  mutate(df2, step = "top_nt", next_step = "pre_super_class"),
  mutate(df3, step = "connection_type", next_step = "top_nt"),
  mutate(df4, step = "target", next_step = "connection_type")
) %>%
  mutate(step = factor(step, levels = c("pre_super_class", "top_nt", "connection_type", "target")),
         next_step = factor(next_step, levels = c("pre_super_class", "top_nt", "connection_type", "target"))) %>%
  mutate(colour = gsub(".*_","",node))

# Create the plot
ggplot(sankey_data, aes(x = step, 
                        next_x = next_step, 
                        node = node, 
                        next_node = next_node,
                        fill = colour,
                        label = node,
                        value = weight)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black") +
  geom_sankey_label(size = 3, color = "white", fill = "black") +
  scale_fill_manual(values = paper.cols) +
  theme_sankey(base_size = 18) +
  labs(title = "neck connective network sankey diagram",
       x = NULL) +
  theme(legend.position = "none")

# Save the plot
ggsave("sankey_plot_axon_dendrite.png", width = 16, height = 10, dpi = 300)
