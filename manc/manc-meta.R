##########################
### GET MANC META DATA ###
##########################
source("banc/banc-startup.R")
drosophila_neurotransmitters <- "/home/ab714/drosophila_neurotransmitters/"
drosophila_neuropeptides <- "/home/ab714/drosophila_neuropeptides/"

# Run locally
banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
banc.connectivity.save.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity"
drosophila_neurotransmitters <- "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters"
drosophila_neuropeptides <- "/Users/GD/LMBD/Papers/synister/drosophila_neuropeptides/"

# Get the meta data
malevnc:::choose_malevnc_dataset('MANC')
mc.find <- neuprintr::neuprint_search("Traced",field="status",dataset="manc:v1.2.1")
mc.ids <- unique(mc.find$bodyid)
mc.meta.orig <- manc_neuprint_meta(mc.ids)

# Organise
mc.meta <- mc.meta.orig %>%
  dplyr::rename(top_nt = predictedNt,
                top_p = predictedNtProb,
                cell_class = class,
                cell_sub_class = subclass,
                other_names = synonyms,
                neuromere = somaNeuromere,
                vfb_id = vfbId) %>%
  dplyr::mutate(cell_type = type,
                entryNerve = ifelse(is.na(entryNerve),NA,entryNerve),
                exitNerve = ifelse(is.na(exitNerve),NA,exitNerve),
                top_p = round(as.numeric(top_p),3),
                top_nt = ifelse(is.na(top_nt),"unknown",top_nt),
                top_nt = ifelse(top_nt%in%"neither","unknown",top_nt)) %>%
  dplyr::mutate(side = dplyr::case_when(
    somaSide=="LHS" ~ "left",
    somaSide=="RHS" ~ "right",
    somaSide=="Midline" ~ "midline",
    rootSide=="LHS" ~ "left",
    rootSide=="RHS" ~ "right",
    rootSide=="Midline" ~ "midline",
    TRUE ~ rootSide 
  )) %>%
  dplyr::mutate(flow = dplyr::case_when(
    grepl("sensory",cell_class) ~ "vnc_afferent",
    grepl("motor",cell_class) ~ "vnc_efferent",
    grepl("ascending_neuron",cell_class) ~ "vnc_brain_transfer",
    grepl("descending_neuron",cell_class) ~ "brain_vnc_transfer",
    TRUE ~ "intrinsic" 
  )) %>%
  dplyr::mutate(target = gsub(" ","_",target)) %>%
  dplyr::mutate(origin = gsub(" ","_",origin)) %>%
  dplyr::mutate(cell_class = gsub(" ","_",cell_class)) %>%
  dplyr::mutate(super_class = cell_class) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    super_class=="motor_neuron" ~ "vnc_motor_neuron",
    cell_sub_class=="EN" ~ "endocrine",
    super_class=="EN" ~ "endocrine",
    cell_sub_class=="efferent_ascending" ~ "endocrine",
    super_class%in%c("TBD","Interneuron_TBD") ~ "",
    grepl("ascending_neuron",super_class) ~ "ascending",
    grepl("descending_neuron",super_class) ~ "descending",
    grepl("sensory_ascending",super_class) ~ "sensory_ascending",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(region = dplyr::case_when(
    super_class%in%c("ascending","descending","sensory_ascending","sensory_descending","efferent_ascending","efferent_descending")~ "neck_connective",
    TRUE ~ "vnc"
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(entryNerve = ifelse(entryNerve%in%c("None","none"),NA,entryNerve),
                exitNerve = ifelse(exitNerve%in%c("None","none"),NA,exitNerve)) %>%
  dplyr::mutate(nerve = dplyr::case_when(
    !is.na(entryNerve) & !is.na(exitNerve) ~ entryNerve,
    !is.na(entryNerve) ~ entryNerve,
    !is.na(exitNerve) ~ exitNerve,
    TRUE ~ NA
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(nerve = gsub("CvC","CV",nerve)) %>%
  dplyr::select(bodyid, 
                post,
                pre,
                voxels,
                top_nt,
                top_nt_p = top_p,
                region,
                side,
                nerve,
                neuromere,
                hemilineage,
                flow,
                super_class,
                cell_class,
                cell_sub_class,
                cell_type,
                type,
                name,
                other_names,
                origin,
                serialMotif,
                cell_function,
                receptorType,
                somaLocation,
                tosomaLocation,
                rootLocation,
                subclassabbr,
                soma,
                target
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    super_class=="vnc_motor_neuron" ~ paste0(target,"_motor_neuron"),
    cell_class=="efferent_ascending" ~ "endocrine_ascending",
    super_class=="endocrine" ~ paste0(target,"_endocrine"),
    super_class=="sensory_neuron" ~ paste0(cell_function,"_sensory_neuron"),
    super_class=="sensory_ascending" ~ paste0(cell_function,"_sensory_ascending"),
    TRUE ~ cell_class
  )) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    super_class=="vnc_motor_neuron" ~ paste0(cell_sub_class,"_motor_neuron"),
    cell_class=="endocrine" ~ paste0(cell_sub_class,"_endocrine"),
    flow=="vnc_afferent"&cell_function=="unknown"&receptorType=="mechanosensory bristle" ~ "tactile",
    TRUE ~ cell_sub_class
  )) %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    cell_class=="vnc_motor_neuron" ~ cell_sub_class,
    cell_class=="putative_efferent_neuron" ~ cell_sub_class,
    cell_class=="ascending" ~ "ascending",
    cell_class=="descending" ~ "descending",
    cell_sub_class=="EN" ~ "endocrine",
    #flow=="vnc_afferent"&receptorType=="campaniform sensilla" ~ "campaniform_sensilla",
    #flow=="vnc_afferent"&receptorType=="chordotonal organ" ~ "chordotonal_organ",
    #flow=="vnc_afferent"&receptorType=="hair plate" ~ "hair_plate",
    flow=="vnc_afferent"&cell_function=="unknown"&receptorType=="mechanosensory bristle" ~ "tactile",
    flow=="vnc_afferent"&cell_function=="unknown"&receptorType=="taste bristle" ~ "chemosensory",
    !is.na(super_class) & (super_class%in%c("ascending","descending","descending_sensory","ascending_sensory","sensory_ascending","sensory_descending")) ~ super_class,
    #grepl("unknown", cell_function) ~ NA,
    TRUE ~ cell_function
  )) %>%
  dplyr::mutate(seed = dplyr::case_when(
    (flow%in%"vnc_afferent") & !is.na(cell_function) & !is.na(nerve) ~ paste(cell_function,snakecase::to_snake_case(nerve),collapse="_",sep="_"),
    (flow%in%"vnc_efferent") & !is.na(cell_function) & !is.na(target) ~ paste(cell_function,snakecase::to_snake_case(target),side,collapse="_",sep="_"),
    (flow%in%"vnc_efferent") & !is.na(cell_function) & is.na(target) ~ paste(cell_function,side,collapse="_",sep="_"),
    (super_class%in%c("ascending","descending","descending_sensory","ascending_sensory","sensroy_ascending","sensory_descending")) & !is.na(side) ~ paste(super_class,side,collapse="_",sep="_"),
    TRUE ~ NA
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(seed = gsub("_R$","_right",seed),
                seed = gsub("_L$","_left",seed),
                seed = gsub("_R_","_right_",seed),
                seed = gsub("_L_","_left_",seed)) %>%
  dplyr::mutate(flow = dplyr::case_when(
    cell_class%in%c("descending_neuron") ~ "brain_vnc_transfer",
    cell_class%in%c("ascending_neuron") ~ "vnc_brain_transfer",
    flow=="intrinsic" ~ "vnc_intrinsic",
    flow=="afferent" ~ "vnc_afferent",
    flow=="transfer" ~ "vnc_afferent",
    flow=="efferent" ~ "vnc_efferent",
    TRUE ~ flow
  )) %>%
  dplyr::mutate(dataset="MANC")
colnames(mc.meta) <- snakecase::to_snake_case(colnames(mc.meta))

# Get known transmitters
nt.data <- readr::read_csv(file.path(drosophila_neurotransmitters,"gt_data.csv"),
                           col_types = banc.col.types)
np.data <- readr::read_csv(file.path(drosophila_neuropeptides,"gt_np_data.csv"),
                           col_types = banc.col.types) %>%
  dplyr::rename(known_nt_source=known_np_source,
                known_nt_evidence=known_np_evidence,
                known_nt_confidence=known_np_confidence)
data <- plyr::rbind.fill(nt.data,np.data) %>%
  dplyr::filter(region%in%c("vnc","neck","neck_connective","VNC"),
                species=="adult_drosophila_melanogaster") 
data[is.na(data)] <- 0
nt.data.processed <- data %>%
  # Select only numeric columns
  dplyr::select(cell_type, where(is.numeric)) %>%
  # Reshape the data to long format
  tidyr::pivot_longer(cols = -cell_type, names_to = "transmitter", values_to = "value") %>%
  # Filter for values >= 1
  dplyr::filter(value >= 1) %>%
  # Join back with original data
  dplyr::right_join(data, by = "cell_type", relationship = "many-to-many") %>%
  # Group by cell_type and concatenate transmitters
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt_source = sprintf("%s (%s)",known_nt_source, known_nt_evidence)) %>%
  dplyr::mutate(known_nt = paste(unique(transmitter), collapse = ","), .groups = "drop") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(known_nt = paste(unique(known_nt), collapse = ";"), .groups = "drop") %>%
  dplyr::mutate(known_nt_source = paste(unique(known_nt_source), collapse = "; "), .groups = "drop") %>%
  # Replace NA in known_nt with empty string
  dplyr::mutate(known_nt = replace_na(known_nt, "")) %>%
  # Filter out rows where known_nt is empty
  dplyr::filter(known_nt != "") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(known_nt_source = paste(unique(known_nt_source),collapse=", ")) %>%
  dplyr::distinct(cell_type, known_nt, known_nt_source)

# Add
mc.meta <- dplyr::left_join(mc.meta, nt.data.processed, by = c("cell_type"))

# Add metrics
metrics.folder <- file.path(banc.nblast.manc.split.save.path,"metrics")
mc.metrics <- read_metrics_csvs(metrics.folder)%>%
  dplyr::select(-neuromore)
keep <- setdiff(colnames(mc.metrics),colnames(mc.meta))
mc.meta <- mc.meta %>%
  dplyr::left_join(mc.metrics[,c("bodyid",keep)], by = "bodyid")

# Save
readr::write_csv(mc.meta, file.path(banc.meta.save.path,"manc_meta.csv"))

# Write meta
arrow::write_feather(mc.meta,
                     file.path(banc.connectivity.save.path, "manc_1.2.1_meta.feather"))

# Announce
message("##### BANCpipeline: manc meta updated #####")
message(sprintf("##### we have meta for : %s neurons", nrow(mc.meta)))

