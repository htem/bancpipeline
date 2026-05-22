#################################
### Choose neurons to process ###
#################################
source("banc/banc-startup.R")
library(fancr)

# FANC cave tables
## https://cave.fanc-fly.com/annotation/views/aligned_volume/fanc_v4?middle_auth_url=global.daf-apis.com%2Fsticky_auth

# Define the mapping for each specific column you want to create
soma_side_tags <- c("soma side", "soma segment")
additional_fields <- c("publication") #"fast neurotransmitter"
cell_class_tags <- c("motor neuron", 
                     "sensory neuron",
                     "chordotonal neuron",
                     "UM neuron")
cell_sub_class_tags <- c("motor neuron primary neurite bundle",
                         "left-right projection pattern", 
                         "leg motor neuron",
                         "innervates leg", 
                         "innervates T1 leg", 
                         "innervates T2 leg", 
                         "innervates T3 leg",
                         "efferent non-motor neuron",
                         "body part innervated",
                         "muscle innervated")
primary_class_tags <- "primary class"
neuron_identity_tags <- "neuron identity"

# Get cell info table
## peripheral_nerves=fanc_cave_query('peripheral_nerves')
## cell_ids_v2=fanc_cave_query('cell_ids_v2')
neuron_information <- fanc_cave_query('neuron_information')
neuron_information$pt_position <- sapply(neuron_information$pt_position, paste, collapse=", ")
fanc.info <- neuron_information %>%
  dplyr::filter(tag2 %in% c(soma_side_tags, additional_fields, cell_class_tags, 
                            cell_sub_class_tags, primary_class_tags, neuron_identity_tags)) %>%
  dplyr::mutate(tag2 = dplyr::recode(tag2,
                       !!!setNames(rep("soma_side", length(soma_side_tags)), soma_side_tags),
                       !!!setNames(additional_fields, additional_fields),
                       !!!setNames(rep("cell_class", length(cell_class_tags)), cell_class_tags),
                       !!!setNames(rep("cell_sub_class", length(cell_sub_class_tags)), cell_sub_class_tags),
                       !!primary_class_tags := "super_class",
                       !!neuron_identity_tags := "cell_type")) %>%
  dplyr::group_by(pt_root_id, tag2) %>%
  dplyr::summarize(tag = paste(sort(unique(tag)), collapse = ", "), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = tag2, 
                   values_from = tag,
                   values_fill = list(tag = NA)) %>%
  dplyr::distinct() 
colnames(fanc.info) <- snakecase::to_snake_case(colnames(fanc.info))
fanc.info <-dplyr::left_join(fanc.info, 
                             neuron_information[,c("pt_root_id",
                                        "pt_supervoxel_id",
                                        "pt_position")] %>%
                               dplyr::distinct(pt_root_id, 
                                               .keep_all = TRUE), 
                   by = "pt_root_id")

# Join with proofread table
proofreading_status <- fanc_cave_query('proofreading_status_table_v0') %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  rbind(fanc_cave_query('proofread_first_pass') %>%
          dplyr::rowwise() %>%
          dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
          dplyr::ungroup()) %>%
  rbind(fanc_cave_query('proofread_second_pass') %>%
          dplyr::rowwise() %>%
          dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
          dplyr::ungroup()) %>%
  dplyr::filter(valid == "t") %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::mutate(proofread = TRUE) %>%
  dplyr::select(pt_root_id, proofread)

# Add cell ID information  
cell_ids_v2 <- fanc_cave_query('cell_ids_v2') %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pt_root_id,.keep_all = TRUE) %>%
  dplyr::distinct(cell_id = id, pt_root_id, root_position = pt_position)
xyz <- nat::xyzmatrix(cell_ids_v2$root_position)
fanc.voxdims <- fanc_voxdims()
xyz.nm <- fanc_raw2nm(xyz)
xyz.manc <- transform_fanc2manc(xyz.nm)/1000
xyz.t <- mirror_manc(xyz.manc)
lrdiffs <- xyz.t[, 1] - xyz.manc[, 1]
sides <- ifelse(lrdiffs>0,"right","left")
cell_ids_v2$side <- sides

# Add nerve information  
peripheral_nerves <- fanc_cave_query('peripheral_nerves') %>%
  dplyr::distinct(pt_root_id,.keep_all = TRUE) %>%
  dplyr::distinct(pt_root_id, nerve = tag)
  
# Add transmitter information  
nt.table <- fanc_cave_query('neurotransmitter_hemilineage_table') %>%
  dplyr::mutate(known_transmitter = dplyr::case_when(
    classification_system=="Glutamatergic" ~ "glutamate",
    classification_system=="Cholinergic" ~ "acetylcholine",
    classification_system=="GABAergic" ~ "gaba",
    TRUE ~ NA
  )) %>%
  dplyr::distinct(pt_root_id,.keep_all = TRUE) %>%
  dplyr::distinct(pt_root_id, hemilineage = cell_type, known_transmitter) %>%
  dplyr::mutate(known_nt_source = "Lacin et al., 2019 (immuno)")

# Perform joins
fanc.meta <- dplyr::left_join(fanc.info, proofreading_status, by = "pt_root_id") %>%
  dplyr::mutate(proofread = ifelse(is.na(proofread),FALSE,proofread)) %>%
  dplyr::left_join(nt.table, by = "pt_root_id") %>%
  dplyr::left_join(peripheral_nerves, by = "pt_root_id") %>%
  dplyr::left_join(cell_ids_v2, by = "pt_root_id") %>%
  dplyr::rename(root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position) %>%
  dplyr::mutate(side = dplyr::case_when(
    grepl("right",soma_side) ~ "right",
    grepl("left",soma_side) ~ "left",
    TRUE ~ side
  ))

# Save
readr::write_csv(fanc.meta, file.path(banc.meta.save.path,"fanc_meta.csv"))

# Announce
message("##### BANCpipeline: fanc meta updated #####")
message(sprintf("##### we have meta for : %s neurons", nrow(fanc.meta)))
  
