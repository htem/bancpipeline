##################
### BUILD BANC ###
##################
source("banc/banc-startup.R")

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

#####################
### GET META DATA ###
#####################

# Get BANC meta data
franken.meta <- franken_meta()
banc.meta <- banctable_query() %>%
  dplyr::filter(!super_class %in% c("glia",
                                    "trachea",
                                    "not_a_neuron"))

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

# Update functions for banc.meta
cns.meta <- banc.meta %>%
  dplyr::left_join(cns.functions.singular,
                   by = "cell_type") %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    !is.na(cell_function.x) ~ cell_function.x,
    TRUE ~ cell_function.y
  )) %>%
  dplyr::select(-cell_function.x,-cell_function.y) %>%
  dplyr::mutate(id = root_id,
                cell_function = gsub("\\/|\\,","_",cell_function)) %>%
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
      grepl("sensory_ascending",super_class) & grepl("SA",cell_type) ~ cell_type,
      grepl("ascending|descending",cell_class) & !is.na(cell_type) & cell_type !="" ~ cell_type,
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
        !is.na(cns_network) ~ cns_network,
        TRUE ~ NA),
      seed_12 = dplyr::case_when(
        grepl("ascending|descending",super_class) ~ paste0(cell_type,"_",root_id),
        TRUE ~ NA),
      seed_13 = root_id,
      seed_14 = cns_network
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
  dplyr::distinct(root_id,
                  root_626,
                  supervoxel_id,
                  position,
                  side,
                  root_region,
                  region,
                  hemilineage,
                  nerve,
                  super_class,
                  cell_class,
                  cell_sub_class,
                  cell_type,
                  body_part_sensory,
                  body_part_effector,
                  peripheral_target_type,
                  cell_function,
                  cell_function_detailed,
                  neurotransmitter_predicted, 
                  neurotransmitter_score, 
                  seed_00,
                  seed_01,
                  seed_02,
                  seed_03,
                  seed_04,
                  seed_05,
                  seed_06,
                  seed_07,
                  seed_08,
                  seed_09,
                  seed_10,
                  seed_11,
                  seed_12, 
                  seed_13,
                  seed_14) %>%
  as.data.frame()

#####################
### GET EDGELIST  ###
#####################

# Get edgelist
banc.view <- FALSE
if(banc.view){
  banc.el <- data.frame()
  while(nrow(banc.el)<=600000){
    banc.el <- banc_edgelist()
    Sys.sleep(60*15)
  }
}else{
  # # Or read from bucket
  banc.version <- banc_version()
  system2("gsutil cp gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v%s/synapses_250226_human_readable.csv.gz /n/data1/hms/neurobio/wilson/banc/exports/synapses_250226_human_readable.csv.gz")
  system2("gunzip /n/data1/hms/neurobio/wilson/banc/exports/synapses_250226_human_readable.csv.gz")
  
  # csv info
  banc.ids <- unique(cns.meta$root_id)
  desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
  col_types <- cols(
    id = col_character(),
    size = col_double(),
    pre_root_id = col_character(),
    post_root_id = col_character(),
    ctr_x = col_double(),
    ctr_y = col_double(),
    ctr_z = col_double(),
    .default = col_double()
  )
  column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                    'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                    'pre_root_id', 'post_supervoxel_id', 'post_root_id')
  bancsynapses <-  "/n/data1/hms/neurobio/wilson/banc/exports/synapses_250226_human_readable.csv.gz"
  banc.syns <- vroom::vroom(bancsynapses,
                            col_names = column_names,
                            col_select = dplyr::all_of(desired_columns),
                            col_types = col_types,
                            skip = 1) %>%
    dplyr::filter(pre_root_id %in% banc.ids,
                  post_root_id %in% banc.ids) %>%
    dplyr::rename(X=ctr_x,
                  Y=ctr_y,
                  Z=ctr_z) %>%
    tibble::as_tibble()
  
  # Calculate edgelist
  banc.el <- banc.syns %>%
    dplyr::group_by(pre_root_id, post_root_id) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pre_pt_root_id = pre_root_id, 
                    post_pt_root_id = post_root_id, 
                    n)
}

# Simplify
banc.elist.simp <- banc.el %>%
  dplyr::rename(post = post_pt_root_id, 
                pre = pre_pt_root_id, 
                count = n) %>%
  dplyr::mutate(pre = as.character(pre),
                post = as.character(post)) %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre, post, .keep_all = TRUE) %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(norm = round(count/post_count,6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(pre, post, count, norm, post_count, pre_count)

# Save versioned data
cns.meta[cns.meta==""] <- NA
arrow::write_feather(banc.elist.simp,
                     file.path(banc.connectivity.save.path,
                               sprintf("banc_%s_edgelist_simple.feather", banc.version)))
arrow::write_feather(cns.meta,
                     file.path(banc.connectivity.save.path,
                               sprintf("banc_%s_meta.feather", banc.version)))

# Test read
banc.meta.test <- arrow::read_feather(
  file.path(banc.connectivity.save.path, sprintf("banc_%s_meta.feather", banc.version)))
nrow(banc.meta.test)

# Test read
elist <- arrow::read_feather(
  file.path(banc.connectivity.save.path, sprintf("banc_%s_edgelist_simple.feather", banc.version)))
nrow(elist)

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'
sync_files(path(A, "connectivity/"), path(B, "connectivity/"), extensions = c("csv", "feather", "parquet"), move.old = FALSE)
remove_empty_dirs(path(A, "connectivity/"))

# Define the remote name
remote_name <- "hms"
local_files <- "/n/data1/hms/neurobio/wilson/banc/connectivity/"
local_files <- list.files(local_files, pattern = "feather$|parquet$", full.names = TRUE)
for(local_file in local_files){
  message("Working on: ", local_file)
  remote_path <- "neuroanat/connectomes/" 
  system(paste("rclone copy", local_file, paste0(remote_name, ":", remote_path)))
  cat("File transfer complete.\n")  
}

