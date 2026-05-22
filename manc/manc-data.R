#############################################
### Calculate MANC flow centrality splits ###
#############################################
source("banc/banc-startup.R")
manc.confidence.threshold <- 0.05
nt.cols <- c("gaba", "acetylcholine", "glutamate", "octopamine", "serotonin", "dopamine", "neither")

# File paths
split.master.folder <- banc.nblast.manc.split.save.path
split.folder <- file.path(split.master.folder,"swc")
synapses.folder <- file.path(split.master.folder,"synapses")
metrics.folder <- file.path(split.master.folder,"metrics")

# # Write meta to .sqlite
# mc.meta <- read_metrics_csvs(metrics.folder)
# mc.nt <- mc.meta[,c("bodyid","top_nt","top_p")] %>%
#   dplyr::rename(top_nt = top_nt, top_nt_p = top_p) %>%
#   dplyr::arrange(dplyr::desc(top_nt_p)) %>%
#   dplyr::distinct(bodyid, .keep_all = TRUE)
# mc.ids <- unique(mc.meta$bodyid)
# DBI::dbWriteTable(con, 
#                   name = "meta", 
#                   value = mc.nt, 
#                   overwrite = TRUE)

# # Simple synapses
# message("##### MANC reading synapses ...")
# mc.synapses <- read_synapse_csvs(synapses.folder)
# colnames(mc.synapses) <- snakecase::to_snake_case(colnames(mc.synapses))
# mc.synapses <- mc.synapses[,!duplicated(colnames(mc.synapses))]
# mc.synapses <- mc.synapses %>%
#   dplyr::select(-root_id)

# # Extract all the synapses from complete/adequate neurons
# message("##### MANC cleaning synapses ...")
# mc.synapses <- mc.synapses %>%
#   dplyr::filter(confidence>manc.confidence.threshold, 
#                 pre_id != post_id) %>%
#   dplyr::select(c(connector_id,
#                   x, y, z, 
#                   confidence, 
#                   prepost, 
#                   pre_id, post_id, 
#                   inside, 
#                   status,
#                   #strahler_order, 
#                   label)) %>%
#   dplyr::left_join(mc.nt, by = c("pre_id"="bodyid")) %>%
#   dplyr::mutate(label = hemibrainr:::standard_compartments(label)) %>% 
#   dplyr::group_by(pre_id, post_id, prepost, label) %>% 
#   dplyr::mutate(count = dplyr::n()) %>% 
#   dplyr::ungroup() %>%
#   as.data.frame(stringsAsFactors = FALSE)
# mc.synapses <- round_dataframe(mc.synapses, digits = 4)

# # Write to .sqlite
# gc()
# DBI::dbWriteTable(con, 
#                   name = "synapses", 
#                   value = mc.synapses, 
#                   overwrite = TRUE)
# Read data
mc.meta <- mc.nt <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "manc_1.2.1_meta.feather"))
mc.synapses <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, "manc_1.2.1_synapses.parquet"))
mc.nt <- mc.meta[,c("bodyid","top_nt","top_p")] %>%
  dplyr::rename(top_nt = top_nt, top_nt_p = top_p) %>%
  dplyr::arrange(dplyr::desc(top_nt_p)) %>%
  dplyr::distinct(bodyid, .keep_all = TRUE)

# Save presynapses separately
message("### MANC separating pre-synapses ...")
mc.synapses.presynapses.label <- mc.synapses %>%
  dplyr::filter(prepost==0)  %>%
  dplyr::distinct(connector_id, label) %>%
  dplyr::rename(pre_label = label) %>% 
  dplyr::ungroup()

# Save postsynapses separately
message("### MANC separating post-synapses ...")
mc.synapses.postsynapses <- mc.synapses %>%
  dplyr::filter(prepost==1) %>%
  dplyr::ungroup()
rm('mc.synapses')

# Do not separate connections by compartment
mc.neurons.syns.elist <- mc.synapses.postsynapses %>%
  dplyr::filter(post_id %in% mc.ids) %>% 
  dplyr::filter(pre_id %in% mc.ids) %>% 
  dplyr::rename(post = post_id) %>% 
  dplyr::rename(pre = pre_id) %>% 
  dplyr::group_by(post, pre) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(post) %>% 
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::mutate(norm = count/post_count) %>% 
  dplyr::distinct(post, 
                  pre, 
                  count, 
                  norm,
                  post_count,
                  .keep_all = FALSE) %>% 
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count)) %>%
  dplyr::filter(count > 0) %>% 
  dplyr::left_join(mc.nt, by = c("post"="bodyid")) %>%
  dplyr::mutate(top_nt=ifelse(is.na(top_nt),"unknown",top_nt)) %>%
  dplyr::rename(post_top_nt = top_nt,
                post_top_nt_p = top_nt_p) %>%
  dplyr::left_join(mc.nt, by = c("pre"="bodyid")) %>%
  dplyr::mutate(top_nt=ifelse(is.na(top_nt),"unknown",top_nt)) %>%
  dplyr::rename(pre_top_nt = top_nt,
                pre_top_nt_p = top_nt_p) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Fix some columns
mc.neurons.syns.elist$norm = round(mc.neurons.syns.elist$norm, digits = 6)
colnames(mc.neurons.syns.elist) = snakecase::to_snake_case(colnames(mc.neurons.syns.elist))

# Write edgelist_simple
arrow::write_feather(mc.neurons.syns.elist,
                     file.path(banc.connectivity.save.path, "manc_1.2.1_edgelist_simple.feather"))
all.edgelist=nrow(mc.neurons.syns.elist)
message("Saved edges: ", all.edgelist, ", targets: ", length(unique(mc.neurons.syns.elist$post)), " sources: ", length(unique(mc.neurons.syns.elist$pre)))

# Get edgelist
message("### MANC getting edge list ...")
mc.neurons.syns.ac.elist <- mc.synapses.postsynapses %>%
  dplyr::filter(prepost == 1,
                post_id %in% mc.ids & pre_id %in% mc.ids) %>% 
  dplyr::rename(post = post_id) %>% 
  dplyr::rename(pre = pre_id) %>% 
  dplyr::rename(post_label = label) %>% 
  dplyr::rowwise() %>%
  dplyr::left_join(mc.synapses.presynapses.label, by = "connector_id") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pre_label = ifelse(is.na(pre_label),"unknown", pre_label)) %>% 
  dplyr::group_by(post, pre, post_label, pre_label) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(post, post_label) %>% 
  dplyr::mutate(post_label_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(post) %>% 
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(norm = count/post_count) %>% 
  dplyr::mutate(norm_label = count/post_label_count) %>% 
  dplyr::distinct(post, 
                  pre, 
                  post_label, 
                  pre_label, 
                  count, 
                  norm,
                  norm_label,
                  post_label_count,
                  post_count,
                  .keep_all = FALSE) %>% 
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post) %>%
  dplyr::filter(!is.na(pre_label), count > 0) %>% 
  dplyr::mutate(connection = paste(pre_label,post_label, sep="-"))  %>%
  dplyr::mutate(pre_label = hemibrainr:::standard_compartments(pre_label)) %>%
  dplyr::mutate(post_label = hemibrainr:::standard_compartments(post_label)) %>%
  dplyr::left_join(mc.nt, by = c("post"="bodyid")) %>%
  dplyr::mutate(top_nt=ifelse(is.na(top_nt),"unknown",top_nt)) %>%
  dplyr::rename(post_top_nt = top_nt,
                post_top_nt_p = top_nt_p) %>%
  dplyr::left_join(mc.nt, by = c("pre"="bodyid")) %>%
  dplyr::mutate(top_nt=ifelse(is.na(top_nt),"unknown",top_nt)) %>%
  dplyr::rename(pre_top_nt = top_nt,
                pre_top_nt_p = top_nt_p) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Fix some columns
mc.neurons.syns.ac.elist$norm = round(mc.neurons.syns.ac.elist$norm, digits = 6)
mc.neurons.syns.elist$norm_label = round(mc.neurons.syns.elist$norm_label, digits = 6)
colnames(mc.neurons.syns.ac.elist) = snakecase::to_snake_case(colnames(mc.neurons.syns.ac.elist))

# Write edgelist
arrow::write_feather(mc.neurons.syns.ac.elist,
                     file.path(banc.connectivity.save.path, "manc_1.2.1_edgelist.feather"))
all.edgelist=nrow(mc.neurons.syns.ac.elist)
message("Saved split edges: ", all.edgelist, ", targets: ", length(unique(mc.neurons.syns.ac.elist$post)), " sources: ", length(unique(mc.neurons.syns.ac.elist$pre)))

