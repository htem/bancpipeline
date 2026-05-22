#############################################
### Calculate BANC flow centrality splits ###
#############################################
source("banc/banc-startup.R")
redo <- FALSE
banc.size.threshold <- 5
nt.cols <- c("gaba", "acetylcholine", "glutamate", "octopamine", "serotonin", "dopamine", "neither")

# File paths
split.master.folder <- banc.l2split.save.path
split.folder <- file.path(split.master.folder,"swc")
synapses.folder <- file.path(split.master.folder,"synapses")
metrics.folder <- file.path(split.master.folder,"metrics")

# Compile meta data
bc <- banctable_query()
bc.nt <- bc[,c("root_id","top_nt","top_nt_conf")] %>%
  dplyr::rename(conf_nt = top_nt, conf_nt_p = top_nt_conf)
bc.ids <- unique(bc$root_id)
bc.meta <- read_metrics_csvs(metrics.folder)
bc.meta <- bc.meta[!duplicated(bc.meta$root_id),]
bc.meta$root_id <- banc_updateids(bc.meta$root_id)
bc$root_id <- banc_updateids(bc$root_id)
bc.meta <- bc.meta[!duplicated(bc.meta$root_id),]
bc.meta <- bc.meta[,c(setdiff(colnames(bc.meta),colnames(bc)),"root_id")]
banc.meta <- dplyr::left_join(bc,bc.meta,by="root_id")

# Write meta
arrow::write_feather(bc.meta, file.path(banc.connectivity.save.path, "banc_meta.feather"))

# Simple synapses
message("##### banc reading synapses ...")
bc.synapses <- read_synapse_csvs(synapses.folder)
colnames(bc.synapses) <- snakecase::to_snake_case(colnames(bc.synapses))
bc.synapses <- bc.synapses[,!duplicated(colnames(bc.synapses))]
bc.synapses <- bc.synapses %>%
  dplyr::select(-root_id)

# Extract all the synapses from complete/adequate neurons
message("##### banc cleaning synapses ...")
bc.synapses <- bc.synapses %>%
  dplyr::filter(size>banc.size.threshold, 
                pre_id != post_id) %>%
  dplyr::select(c(connector_id,
                  x, y, z, 
                  size, 
                  prepost, 
                  pre_id, post_id, 
                  inside, 
                  status,
                  #strahler_order, 
                  label)) %>%
  dplyr::left_join(bc.nt, by = c("pre_id"="root_id")) %>%
  dplyr::mutate(label = hemibrainr:::standard_compartments(label)) %>% 
  dplyr::group_by(pre_id, post_id, prepost, label) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE)
bc.synapses <- round_dataframe(bc.synapses, digits = 4)

# Write synapses
gc()
write_connectome_data(bc.synapses,
                      file.path(banc.connectivity.save.path, "banc_synapses.parquet"),
                      format = "parquet")

# Save presynapses separately
message("### banc separating pre-synapses ...")
bc.synapses.presynapses <- bc.synapses %>%
  dplyr::filter(prepost==0)  %>%
  dplyr::ungroup()

# Save postsynapses separately
message("### banc separating post-synapses ...")
bc.synapses.postsynapses <- bc.synapses %>%
  dplyr::filter(prepost==1) %>%
  dplyr::ungroup()
rm('bc.synapses')

# Get edgelist
message("### banc getting edge list ...")
bc.neurons.syns.ac.elist <- bc.synapses.postsynapses %>%
  dplyr::filter(prepost == 1,
                post_id %in% bc.ids| pre_id %in% bc.ids) %>% 
  dplyr::rename(post = post_id) %>% 
  dplyr::rename(pre = pre_id) %>% 
  dplyr::rename(post_label = label) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(pre_label = bc.synapses.presynapses[match(as.character(connector_id),
                                                          as.character(bc.synapses.presynapses$connector_id))
                                                    ,"label"]) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(pre_label = ifelse(is.na(pre_label),"unknown", pre_label)) %>% 
  dplyr::group_by(post, pre, post_label,pre_label) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(post, post_label) %>% 
  dplyr::mutate(post_label_count = sum(count, na.rm = TRUE)) %>%
  dplyr::mutate(norm = count/post_label_count) %>% 
  dplyr::distinct(post, 
                  pre, 
                  post_label, 
                  pre_label, 
                  count, 
                  norm,
                  post_label_count,
                  .keep_all = FALSE) %>% 
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count)) %>%
  dplyr::filter(!is.na(pre_label), count > 0) %>% 
  dplyr::mutate(connection = paste(pre_label,post_label, sep="-"))  %>%
  dplyr::mutate(pre_label = hemibrainr:::standard_compartments(pre_label)) %>%
  dplyr::mutate(post_label = hemibrainr:::standard_compartments(post_label)) %>%
  dplyr::left_join(bc.nt, by = c("post"="root_id")) %>%
  dplyr::mutate(conf_nt=ifelse(is.na(conf_nt),"unknown",conf_nt)) %>%
  dplyr::rename(post_conf_nt = conf_nt,
                post_conf_nt_p = conf_nt_p) %>%
  dplyr::left_join(bc.nt, by = c("pre"="root_id")) %>%
  dplyr::mutate(conf_nt=ifelse(is.na(conf_nt),"unknown",conf_nt)) %>%
  dplyr::rename(pre_conf_nt = conf_nt,
                pre_conf_nt_p = conf_nt_p) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Fix some columns
bc.neurons.syns.ac.elist$norm = round(bc.neurons.syns.ac.elist$norm , digits = 4)
colnames(bc.neurons.syns.ac.elist) = snakecase::to_snake_case(colnames(bc.neurons.syns.ac.elist))

# Write edgelist
arrow::write_feather(bc.neurons.syns.ac.elist,
                     file.path(banc.connectivity.save.path, "banc_edgelist.feather"))
all.edgelist=nrow(bc.neurons.syns.ac.elist)
message("Saved edges: ", all.edgelist, ", targets: ", length(unique(bc.neurons.syns.ac.elist$post)), " sources: ", length(unique(bc.neurons.syns.ac.elist$pre)))