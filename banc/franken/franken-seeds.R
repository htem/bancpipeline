#' franken-seeds — Organise seed labels in franken-brain (descending / ascending hierarchical clusters).
#'
#' Reads the franken-brain meta and the cluster CSVs from BANC-project, then
#' aligns naming conventions (e.g. `SD_*` / `ED_*` → `DN_*`, `SA_*` / `EA_*`
#' → `AN_*`), summarises cluster composition, and emits seed worksheets used
#' to seed downstream hierarchical clustering.
#'
#' @section Reads:
#'   - `franken_meta()`
#'   - cluster CSVs under `data/banc_annotations/v888/`
#'
#' @section Writes:
#'   - seed worksheet CSVs under `<banc.meta.save.path>/`

#######################################
### ORGANISE SEEDS IN FRANKEN BRAIN ###
#######################################
source("banc/banc-startup.R")
library(knitr)

# # Get UMAP clusters
# umap.dn.df <- read_csv(file = "/home/ab714/BANC-project/data/banc_neck_functional_classes.csv", col_types = banc.col.types) %>%
#   dplyr::mutate(cluster = gsub("SD_|ED_","DN_",cluster),
#                 cluster = gsub("SA_|EA_","AN_",cluster)) %>%
#   dplyr::mutate(clusterno = gsub(".*_","",cluster)) %>%
#   dplyr::mutate(id = ifelse(is.na(id),composite_cell_type,id))

# Helper functions
extract_three_letters <- function(text) {
  sapply(text, function(t) {
    three_letters <- stringr::str_extract(t, "^[A-Za-z]{3}")
    if (!is.na(three_letters)) {
      return(three_letters)
    }
    two_letters <- stringr::str_extract(t, "^[A-Za-z]{2}")
    if (!is.na(two_letters)) {
      return(two_letters)
    }
    one_letter <- stringr::str_extract(t, "^[A-Za-z]{1}")
    return(one_letter)
  })
}
summary_table <- function(meta, 
                          col1 = "region",
                          col2 = "seed_01",
                          copy = TRUE
){
  # Create the table
  meta <- meta[!is.na(meta[[col2]]),]
  summary_table <- table(meta[[col1]], meta[[col2]])
  
  # Convert to a data frame
  summary_df <- as.data.frame.matrix(summary_table)
  summary_df <- summary_df[order(rownames(summary_df)),order(colnames(summary_df))]
  
  # Add row names as a column (for better markdown formatting)
  summary_df <- summary_df %>%
    tibble::rownames_to_column(var = col1)
  
  # Generate markdown table
  markdown_table <- kable(summary_df, format = "markdown")
  
  # Copy to clipboard
  #clipr::write_clip(markdown_table)
  
  # Print the markdown table (you can copy this output to your README.md)
  if(copy){
    clipr::write_clip(markdown_table)
  }else{
    print(markdown_table)
  }
}

# # Move optic sensories forward
# banctable_move_to_bigdata(where = "`super_class` = 'sensory'",
#                           table = "franken_meta", 
#                           base = "cns_meta",
#                           invert = TRUE,
#                           row_ids = b$`_id`)

# Get singular functions
cns.functions <- bancr::banctable_query(sql = "SELECT * FROM functions") %>%
  dplyr::select(-starts_with("_")) %>%
  dplyr::filter(!is.na(cell_type))
cns.functions.singular <- cns.functions %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    !is.na(response)&response!="" ~ response,
    !is.na(behaviour)&behaviour!="" ~ behaviour,
    !is.na(valence)&valence!="" ~ valence,
    !is.na(modality)&modality!="" ~ modality,
    TRUE ~ NA
  )) %>%
  dplyr::distinct(cell_type,.keep_all = TRUE) %>%
  dplyr::distinct(cell_type, cell_function)
  

# Get franken meta
franken.orig <- bancr::franken_meta()
cns.meta <- franken.orig %>%
  dplyr::left_join(cns.functions.singular,
                  by = "cell_type") %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    !is.na(cell_function.x) ~ cell_function.x,
    TRUE ~ cell_function.y
  )) %>%
  dplyr::select(-cell_function.x,-cell_function.y) %>%
  dplyr::mutate(id = neuron_id) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(multi_neuromere = length(unique(neuromere))>1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(composite_cell_type = dplyr::case_when(
    is.na(neuromere)|neuromere=="" ~ cell_type,
    multi_neuromere ~ paste0(cell_type,"_",neuromere),
    TRUE ~ cell_type
  )) %>%
  dplyr::ungroup() %>%
  # dplyr::left_join(umap.dn.df %>%
  #                    dplyr::select(cell_type, dn_cluster = cluster) %>%
  #                    dplyr::distinct(cell_type, .keep_all = TRUE),
  #                  by = c('cell_type')) %>%
  dplyr::mutate(cell_function = gsub("\\/|\\,","_",cell_function))

# APPLY SEED LOGIC
cns.meta.new <- cns.meta %>%
  dplyr::mutate(neuromere = dplyr::case_when(
    is.na(neuromere) ~ "",
    TRUE ~ neuromere
  )) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(repeated = length(unique(neuromere))>1) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !repeated ~ cell_type,
    TRUE ~ paste0(cell_type,"_",neuromere)
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cell_function = dplyr::case_when(
      is.na(cell_function)|cell_function=="" ~ "unknown",
      TRUE ~ cell_function), 
    body_part_sensory = dplyr::case_when(
      grepl("head|^frontal|^frontoorbital|^orbital|^interocellar|^vibrissa|^interommatidial|^occipital_dorsal|^occipital_ventral|^postorbital_dorsal|^postorbital_ventral|^vertical|^postocellar|^supracervical",body_part_sensory) ~ "head",
      is.na(body_part_sensory)|body_part_sensory=="" ~ "unknown",
      TRUE ~ body_part_sensory), 
    cell_function_detailed = dplyr::case_when(
      !is.na(cell_function_detailed)&!is.na(cell_function) ~ paste0(cell_function,"_",cell_function_detailed),
      !is.na(cell_function_detailed)&!is.na(cell_function) ~ cell_function,
      is.na(cell_function) ~ "unknown",
      TRUE ~ cell_function_detailed), 
    cell_function_nerve = dplyr::case_when(
      !is.na(cell_function_detailed) ~ cell_function_detailed,
      !is.na(nerve) ~ gsub("_r$|_l$|_left$|_right$|_R$|_L$|^right_|^left_","",nerve),
      TRUE ~ cell_function), 
    seed_00 = dplyr::case_when(
      is.na(side)| !side %in% c("left","right") ~ NA,
      grepl("sensory",super_class) ~ cell_type,
      grepl("central_complex_output|mushroom_body_output_neuron",cell_class) ~ cell_type,
      grepl("visual_projection",super_class) ~ cell_type,
      grepl("ascending|descending",super_class) ~ cell_type,
      TRUE ~ NA), 
    seed_01 = dplyr::case_when(
      grepl("efferent",flow) ~ NA,
      is.na(cell_function)|cell_function=="unknown" ~ NA,
      TRUE ~ paste0(super_class,"_",cell_function)), 
    seed_02 = dplyr::case_when(
      is.na(side)|!side %in% c("left","right")|is.na(cell_sub_class) ~ NA,
      grepl("sensory",super_class) ~ cell_sub_class,
      TRUE ~ NA),
    seed_03 = dplyr::case_when(
      !grepl("sensory|visual_projection",super_class)|is.na(cell_function_detailed) ~ NA,
      grepl("visual_projection",super_class) ~  paste0("visual_projection_",cell_function_detailed),
      grepl("sensory",super_class) ~ paste0(body_part_sensory,"_",cell_function_detailed),
      TRUE ~ NA),
    seed_04 = dplyr::case_when(
      is.na(side)| !side %in% c("left","right")|is.na(cell_sub_class) ~ NA,
      grepl("sensory",super_class) ~ paste0(cell_sub_class,"_",side),
      TRUE ~ NA),
    seed_05 = dplyr::case_when(
      grepl("central_complex_output",cell_class) & !is.na(cell_sub_class) & cell_sub_class !="" ~ cell_sub_class,
      grepl("ascending|descending",super_class) & !is.na(cluster) & cluster !="" ~ cluster,
      grepl("mushroom_body_output_neuron",cell_class) & !is.na(cell_class) & cell_class !="" ~ cell_class,
      grepl("visual_projection",super_class) & !is.na(cell_type) & cell_type !="" ~ extract_three_letters(cell_type),
      TRUE ~ NA),
    seed_06 = dplyr::case_when(
      is.na(side)|!side %in% c("left","right","midline","center") ~ NA,
      grepl("central_complex_output",cell_class) & !is.na(cell_sub_class) & cell_sub_class !="" ~ paste0(cell_sub_class,"_",side),
      grepl("ascending|descending",super_class) & !is.na(cluster) & cluster !="" ~ paste0(cluster,"_",side),
      grepl("mushroom_body_output_neuron",cell_class) & !is.na(cell_class) & cell_class !="" ~ paste0(cell_class,"_",side),
      grepl("visual_projection",super_class) & !is.na(cell_type) & cell_type !="" ~ paste0(extract_three_letters(cell_type),"_",side),
      is.na(nerve)|nerve=="" ~ NA,
      TRUE ~ NA),
    seed_07 = dplyr::case_when(
      grepl("central_complex_output|mushroom_body_output_neuron",cell_class) ~ cell_type,
      grepl("^EPG|^EL",cell_type) ~ cell_type,
      grepl("sensory_ascending",super_class) & grepl("SA",composite_cell_type) ~ composite_cell_type,
      grepl("ascending|descending",cell_class) & !is.na(composite_cell_type) & composite_cell_type !="" ~ composite_cell_type,
      grepl("visual_projection",super_class) ~ cell_type,
      TRUE ~ NA),
    seed_08 = dplyr::case_when(
      !grepl("sensory",super_class) ~ NA,
      grepl("proprio|tactile|contract|vib",cell_function) ~ paste0(cell_function,"_",body_part_sensory,"_",cell_function_nerve,"_",side),
      TRUE ~ NA),
    seed_09 = dplyr::case_when(
      grepl("efferent",flow) ~ NA,
      is.na(peripheral_target_type) ~ NA,
      TRUE ~ peripheral_target_type),
    seed_10 = dplyr::case_when(
      !grepl("sensory",super_class)|is.na(peripheral_target_type) ~ NA,
      TRUE ~ paste0(body_part_sensory,"_",peripheral_target_type)),
    seed_11 = dplyr::case_when(
      !is.na(cluster) ~ cluster,
      TRUE ~ NA),
    seed_12 = dplyr::case_when(
      grepl("ascending|descending|visual_projection",super_class) ~ paste0(cell_type,"_",id),
      TRUE ~ NA)
  ) %>%
  dplyr::mutate(
    seed_00 = gsub(";","_",seed_00),
    seed_01 = gsub(";","_",seed_01),
    seed_02 = gsub(";","_",seed_02),
    seed_03 = gsub(";","_",seed_03),
    seed_04 = gsub(";","_",seed_04),
    seed_05 = gsub(";","_",seed_05),
    seed_06 = gsub(";","_",seed_06),
    seed_07 = gsub(";","_",seed_07),
    seed_08 = gsub(";","_",seed_08),
    seed_09 = gsub(";","_",seed_09),
    seed_10 = gsub(";","_",seed_10),
    seed_11 = gsub(";","_",seed_11),
    seed_12 = gsub(";","_",seed_12)
  ) %>%
  as.data.frame()

# Copy to .Md file
#summary_table(cns.meta.new, col1 = "super_class", col2 = "seed_00")
summary_table(cns.meta.new, col1 = "body_part_sensory", col2 = "seed_01", copy = TRUE)
summary_table(cns.meta.new, col1 = "region", col2 = "seed_02")
summary_table(cns.meta.new, col1 = "cell_function", col2 = "seed_03")
summary_table(cns.meta.new, col1 = "region", col2 = "seed_04")
summary_table(cns.meta.new, col1 = "super_class", col2 = "seed_05")
summary_table(cns.meta.new, col1 = "super_class", col2 = "seed_06")
summary_table(cns.meta.new, col1 = "super_class", col2 = "seed_07")
summary_table(cns.meta.new, col1 = "body_part_sensory", col2 = "seed_08")
summary_table(cns.meta.new, col1 = "body_part_sensory", col2 = "seed_09")
summary_table(cns.meta.new, col1 = "body_part_sensory", col2 = "seed_10")
summary_table(cns.meta.new, col1 = "top_nt", col2 = "seed_11")

# UPDATE
franken.meta.upload <- cns.meta.new %>%
  dplyr::select(`_id`, seed_00, seed_01, seed_02, seed_03, seed_04, seed_05, seed_06, seed_07, seed_08, seed_09, seed_10, seed_11, seed_12) %>%
  dplyr::filter(if_any(starts_with("seed_0"), ~!is.na(.)))
franken.meta.upload[is.na(franken.meta.upload)] <- ""
banctable_update_rows(base='cns_meta', 
                      table = 'franken_meta', 
                      df = franken.meta.upload, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

################
### CLUSTERS ###
################
cluster.update <- left_join(franken.orig %>%
                              dplyr::select(`_id`,cell_type, super_class),
                            umap.dn.df %>%
                              dplyr::select(cell_type, dn_cluster = cluster) %>%
                              dplyr::distinct(cell_type, .keep_all = TRUE),
                            by = c('cell_type')) %>%
  # left_join(umap.an.df %>%
  #             dplyr::select(cell_type, an_cluster = cluster) %>%
  #             dplyr::distinct(cell_type, .keep_all = TRUE),
  #           by = c('cell_type'))  %>%
  # left_join(umap.eff.df %>%
  #             dplyr::select(cell_type, eff_cluster = cluster) %>%
  #             dplyr::distinct(cell_type, .keep_all = TRUE),
  #           by = c('cell_type'))  %>%
  # left_join(umap.sez.df %>%
  #             dplyr::select(cell_type, sez_cluster = cluster) %>%
  #             dplyr::distinct(cell_type, .keep_all = TRUE),
  #           by = c('cell_type'))  %>%
  dplyr::mutate(cluster = dplyr::case_when(
    !is.na(dn_cluster)&grepl("descending|ascending",super_class) ~ dn_cluster,
    #!is.na(eff_cluster)&grepl("efferent|motor|endocrine",super_class) ~ eff_cluster,
    #!is.na(sez_cluster) ~ sez_cluster,
    TRUE ~ NA
  )) %>%
  dplyr::select(-dn_cluster) %>%
  dplyr::filter(!is.na(cluster))

# Update!
banctable_update_rows(base='cns_meta', 
                      table = 'franken_meta', 
                      df = cluster.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

#############
### FIXES ###
#############

# # FAFB
# fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
# fw.update <- fw.meta %>%
#   dplyr::mutate(total_inputs = as.numeric(total_inputs),
#                 cell_dcv_density = as.numeric(cell_dcv_density),
#                 soma_dcv_density = as.numeric(soma_dcv_density)) %>%
#   dplyr::select(fafb_id = root_783,
#                 input_connections = total_inputs,
#                 cell_dcv_density, 
#                 soma_dcv_density,
#                 hemibrain_type,
#                 known_np) %>%
#   dplyr::distinct(fafb_id, .keep_all = TRUE)
# fw.update <- franken.orig %>%
#   dplyr::select(`_id`, fafb_id) %>%
#   dplyr::left_join(fw.update,
#                    by = "fafb_id")  %>%
#   dplyr::distinct(`_id`, .keep_all = TRUE) %>%
#   dplyr::filter(!is.na(fafb_id), !is.na(`_id`)) %>%
#   dplyr::anti_join(franken.orig %>%
#                      dplyr::distinct(`_id`, .keep_all = TRUE), 
#                    by = c("_id","fafb_id","input_connections","cell_dcv_density","soma_dcv_density","known_np","hemibrain_type")) %>%
#   dplyr::arrange(input_connections) %>%
#   as.data.frame()
# fw.update$input_connections[is.na(fw.update$input_connections)] <- -1
# fw.update$cell_dcv_density[is.na(fw.update$cell_dcv_density)] <- -1
# fw.update$soma_dcv_density[is.na(fw.update$soma_dcv_density)] <- -1
# fw.update$known_np[is.na(fw.update$known_np)] <- ''
# fw.update$hemibrain_type[is.na(fw.update$hemibrain_type)] <- ''
# banctable_update_rows(base='cns_meta', 
#                       table = 'franken_meta', 
#                       df = fw.update[,c("_id", "input_connections", "cell_dcv_density","soma_dcv_density", "known_np", "hemibrain_type")], 
#                       append_allowed = FALSE, 
#                       chunksize = 1000)
# 
# # MANC
# mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
# mc.update <- mc.meta %>%
#   dplyr::mutate(post = as.numeric(post)) %>%
#   dplyr::select(manc_id = bodyid,
#                 input_connections = post) %>%
#   dplyr::distinct(manc_id, .keep_all = TRUE)
# mc.update <- franken.orig %>%
#   dplyr::select(`_id`, manc_id) %>%
#   dplyr::left_join(mc.update,
#                    by = "manc_id") %>%
#   dplyr::distinct(`_id`, .keep_all = TRUE) %>%
#   dplyr::filter(!is.na(manc_id), !is.na(`_id`))
# banctable_update_rows(base='cns_meta', 
#                       table = 'franken_meta', 
#                       df = mc.update, 
#                       append_allowed = FALSE, 
#                       chunksize = 1000)

# ## Lateral horn neuron updates
# franken.meta <- franken_meta()
# franken.meta.update <- franken.meta.orig <- franken.meta %>%
#   dplyr::filter(cell_type %in% na.omit(ton.info$cell.type) | hemibrain_type %in% na.omit(ton.info$cell.type)) %>%
#   dplyr::select(`_id`, neuron_id, cell_type, hemibrain_type, super_class, cell_class, cell_sub_class)
# for(i in 1:nrow(ton.info)){
#   ct <- ton.info[i,"cell.type"]
#   cc <- ton.info[i,"class"]
#   bi <- ton.info[i,"bodyid"]
#   ct.simp <- gsub("_[a-zA-Z].*", "", ct)
#   if(is.na(ct)){
#     next
#   }
#   if(cc %in% c("LHLN","LHN","LHON","TOON","WEDPN","LHCENT")){
#     if(cc=="LHN"){
#       cc="LHON"
#     }
#     fw <- franken.meta.update %>%
#       dplyr::filter(cell_type == ct | hemibrain_type == ct)
#     if(nrow(fw)){
#       choose <- (franken.meta.update$cell_type == ct | franken.meta.update$hemibrain_type == ct | franken.meta.update$cell_type == ct.simp | franken.meta.update$hemibrain_type == ct.simp)
#       choose[is.na(choose)] <- FALSE
#       for(chosen in which(choose)){
#         if(is.na(franken.meta.update[chosen,"cell_class"])|franken.meta.update[chosen,"cell_class"]%in%c("LHON","LHLN","LHN")){
#           franken.meta.update[chosen,"cell_class"] <- cc
#         }
#         if(cc %in% c("LHLN","LHN","LHON","LHCENT")){
#           if(is.na(franken.meta.update[chosen,"super_class"])|franken.meta.update[chosen,"super_class"]=="brain_central_other"){
#             franken.meta.update[chosen,"super_class"] <- 'lateral_horn'
#           }
#         } 
#       }
#     }
#   }
# }
# banctable_update_rows(base='cns_meta',
#                       table = 'franken_meta',
#                       df = franken.meta.update,
#                       append_allowed = FALSE,
#                       chunksize = 1000)
# f = franken_meta()
# a = subset(f, super_class == "lateral_horn")
# b = a %>% 
#   select(root_783 = fafb_id, cell_type, hemibrain_type, cell_sub_class, cell_class, super_class, hemilineage, top_nt, known_nt, known_nt_source)
# write_csv(franken.meta.update,"/Users/abates/Downloads/fafb_783_lh_neurons.csv")
# write_csv(ton.info,"/Users/abates/Downloads/hemibrain_1_2_1_toon_neurons.csv")


