##########################
### GET FAFB META DATA ###
##########################
source("banc/banc-startup.R")

# Run locally
drosophila_neurotransmitters <- "/home/ab714/drosophila_neurotransmitters/"
banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
#banc.connectivity.save.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity"
fafb.orig.save.path <- banc.meta.save.path
fafb.orig.save.path <- '/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/783/'

# Get fafbpipeline meta data product
# source("~/fafbpipeline/flywire/flywire-ids.R")
# source("~/fafbpipeline/flywire/flywire-metrics.R")
fw.meta.orig <- readr::read_csv(file.path(fafb.orig.save.path,'flywire_meta.csv'), 
                                col_types = hemibrainr:::sql_col_types)
fw.meta.orig$root_id <- fw.meta.orig$root_783

# Order statuses correctly for duplicate removal
statuses <- c("adequate", "complete", "proofread",
              "incomplete", "unassessed", "unknown", "to_review",
              "wrong_hemilineage", "wrong_side", "outlier_bio", "bad_nucleus",
              "NA", "outlier_seg",
              "needs_extending", "merge_error", "large_fragment",
              "fragment", "tiny", "hard", "tadpole",
              "duplicate", "not_neuron", "not_a_neuron")

# Get meta data
ft <- fafbseg::flytable_query("select _id, root_id, root_630, root_783, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, soma_x, soma_y, soma_z, side, ito_lee_hemilineage, hartenstein_hemilineage, nerve, top_nt, flow, super_class, cell_class, cell_sub_class, cell_type, hemibrain_match, hemibrain_type, cb_type, root_duplicated, morphology_group, known_nt, known_nt_source, notes from info")
ft$region <- 'midbrain'
ft.optic <- fafbseg::flytable_query("select * from optic")
ft.optic <- ft.optic[, intersect(colnames(ft.optic),colnames(ft))]
ft.optic <- ft.optic[!ft.optic$root_id %in% ft$root_id,]
ft.optic$region <- 'optic_lobes'
statuses <- c(statuses, setdiff(unique(ft$status),statuses))
fw.meta <- plyr::rbind.fill(ft, ft.optic) %>% 
  dplyr::filter(!duplicated(root_id)) %>%
  dplyr::mutate(hemilineage = ito_lee_hemilineage,
                cell_type = dplyr::case_when(
                  !is.na(cell_type) ~ cell_type,
                  !is.na(hemibrain_type) ~ hemibrain_type,
                  !is.na(cb_type) ~ cb_type,
                  !is.na(cell_sub_class) ~ cell_sub_class,
                  !is.na(morphology_group) ~ morphology_group,
                  TRUE ~ cell_type)
  ) %>%
  dplyr::filter(
    !super_class %in% c('not_a_neuron'),
    !cell_class %in% c('glia', 'trachea','tadpole', "not_a_neuron", "putative_glia", "fragment"), 
    !status %in% c('duplicate', 'not_a_neuron', 'tiny', 'merge_error', 'tadpole', 'large_fragment')
  ) %>%
  dplyr::mutate(ito_lee_hemilineage = ifelse(grepl("prim|Prim",ito_lee_hemilineage),"putative_primary",ito_lee_hemilineage)) %>%
  dplyr::mutate(cell_class = ifelse(grepl("cell_class==clock|cell_class=clock",notes),'clock',cell_class),
                super_class = ifelse(grepl("cell_class==clock|cell_class=clock",notes),'clock',super_class)) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    grepl("^LH",cell_type)&flow=="intrinsic"&is.na(cell_class)~ "LHON",
    TRUE ~ cell_class
  )) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    cell_class %in% c('glia','putative_glia') ~ "glia",
    cell_class%in%c("ALIN","ALLN","ALON","ALPN") ~ "antennal_lobe",
    cell_class%in%c("LHLN","LHON","TOON","LHCENT") ~ "lateral_horn",
    cell_class%in%c("DAN","Kenyon_Cell","MBIN","MBON","MB") ~ "mushroom_body",
    cell_type%in%c("APL","DPM") ~ "mushroom_body",
    cell_class%in%c("CX") ~ "central_complex",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    cell_class=="CX" & grepl("hDelta",cell_type) ~ "hDelta",
    cell_class=="CX" & grepl("vDelta",cell_type) ~ "vDelta",
    cell_class=="CX" & grepl("^FB|^CB.FB|^SA1|^SA2|^SA3",cell_type) ~ "tangential",
    cell_class=="CX" & grepl("^ExR",cell_type) ~ "ExR",
    cell_class=="CX" & grepl("^FC",cell_type) ~ "FS",
    cell_class=="CX" & grepl("^FS",cell_type) ~ "FS",
    cell_class=="CX" & grepl("^PFL",cell_type) ~ "PFL",
    cell_class=="CX" & grepl("^PFN",cell_type) ~ "PFN",
    cell_class=="CX" & grepl("^ER",cell_type) ~ "ER",
    cell_class=="CX" & grepl("^FR1|^FR2",cell_type) ~ "FR",
    cell_class=="CX" & grepl("^PFR$",cell_type) ~ "PFR",
    cell_class=="CX" & grepl("^PFGs$",cell_type) ~ "PFG",
    cell_class=="CX" & grepl("^PEG",cell_type) ~ "PEG",
    cell_class=="CX" & grepl("^PEN_",cell_type) ~ "PEN",
    cell_class=="CX" & grepl("^GLNO$|^LCNOp$|^LCNOpm$|^LNO1$|^LNO1$|^LNO2$",cell_type) ~ "NO_input",
    cell_class=="CX" & grepl("^IbSpsP|^LPsP|^SpsP",cell_type) ~ "PB_input",
    TRUE ~ cell_sub_class
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    cell_class=="CX" & grepl("hDelta|vDelta|FC|columnar|^PEN$",cell_sub_class) ~ "CX_intrinsic",
    cell_class=="CX" & grepl("^EL$|^EPG|^Delta7|^P6\\-8P9|^P1\\-9",cell_type) ~ "CX_intrinsic",
    cell_class=="CX" & grepl("^NO$|PB_input|ER|tangential|ExR|PFN|NO_input",cell_sub_class) ~ "CX_input",
    cell_class=="CX" & grepl("^FS|^PFL|^PFR|^FR|^PFG|^PEG",cell_sub_class) ~ "CX_output",
    cell_class=="CX" ~ "CX_intrinsic",
    grepl("motor",super_class) & grepl("salivary|crop|haustellum|ciberial",cell_sub_class) ~ "ingestion_motor_neuron",
    grepl("motor",super_class) ~ cell_sub_class,
    super_class=="ascending" ~ "ascending_neuron",
    super_class=="descending" ~ "descending_neuron",
    super_class=="sensory_ascending" ~ "sensory_ascending",
    TRUE ~ cell_class
  )) %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    grepl("sensory",super_class) & cell_class=='unknown_sensory' ~ NA,
    grepl("sensory",super_class) & grepl("^AN_4|^SA_DMT",cell_type) ~ "haltere",
    grepl("sensory",super_class) & grepl("^SA_VTV_",cell_type) ~ "gustatory",
    grepl("sensory",super_class) & super_class=='sensory_ascending' ~ "sensory_ascending",
    grepl("sensory",super_class) & cell_class=='olfactory' ~ "olfactory",
    grepl("sensory",super_class) & cell_sub_class=='ocellar' ~ cell_sub_class,
    grepl("sensory",super_class) & cell_class=='ocellar' ~ "ocellar",
    grepl("sensory",super_class) & grepl("^JO\\-",cell_type) & !is.na(cell_sub_class) ~ cell_sub_class,
    grepl("sensory",super_class) & grepl("^JO\\-",cell_type) & is.na(cell_sub_class) ~ "johnstons_organ_other",
    # grepl("sensory",super_class) & cell_sub_class=='auditory' ~ cell_sub_class,
    # grepl("sensory",super_class) & cell_sub_class=='eye bristle' ~ "eye_bristle",
    # grepl("sensory",super_class) & cell_sub_class=='head bristle' ~ "head_bristle",
    # grepl("sensory",super_class) & cell_type=='TPMN' ~ "proboscis_bristle",
    grepl("sensory",super_class) & cell_class=='visual' & cell_type %in% "R1-6" ~ "visual_achromatic",
    grepl("sensory",super_class) & cell_class=='visual' & cell_type %in% c("R7","R8") ~ "visual_chromatic",
    grepl("sensory",super_class) & cell_class=='thermosensory' ~ "thermosensory",
    grepl("sensory",super_class) & cell_class=='mechanosensory' ~ "mechanosensory",
    grepl("sensory",super_class) & cell_class=='hygrosensory' ~ "hygrosensory",
    grepl("sensory",super_class) & cell_class=='gustatory' ~ "gustatory",
    grepl("sensory",super_class) & !is.na(cell_class) ~ cell_class,
    grepl("endocrine",super_class) ~ "endocrine",
    grepl("motor",super_class) & grepl("salivary|crop|haustellum|ciberial",cell_sub_class) ~ "ingestion_motor_neuron",
    grepl("motor",super_class) ~ cell_sub_class,
    grepl("motor",cell_class) ~ cell_sub_class,
    grepl("efferent",super_class) & cell_class=="pars_intercerebralis" ~ "endocrine",
    grepl("efferent", super_class) & cell_class=="pars_lateralis" ~ "endocrine",
    grepl("efferent", super_class) & cell_class=="ocellar" ~ "ocellar",
    cell_class=="ascending" ~ "ascending",
    cell_class=="descending" ~ "descending",
    TRUE ~ NA
  )) %>%
  dplyr::arrange(proofread, cell_type, cell_class, hemilineage) %>%
  dplyr::mutate(status = factor(status, levels = statuses),
                side = ifelse(side=="na",NA,side)) %>%
  dplyr::arrange(status) %>%
  dplyr::filter(!duplicated(root_783)) %>%
  dplyr::select(-`_id`, 
                -ito_lee_hemilineage, 
                -hartenstein_hemilineage,
                -root_duplicated,
                -cb_type) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(seed = dplyr::case_when(
    grepl("motor",cell_class) ~ cell_sub_class,
    grepl("sensory",super_class) & cell_sub_class=='eye bristle' ~ "eye_bristle",
    grepl("sensory",super_class) & cell_sub_class=='head bristle' ~ "head_bristle",
    grepl("sensory",super_class) & cell_type=='TPMN' ~ "proboscis_bristle",
    grepl("sensory",super_class) & cell_class=='visual' & cell_type %in% "R1-6" ~ "visual_achromatic",
    grepl("sensory",super_class) & cell_class=='visual' & cell_type %in% c("R7","R8") ~ "visual_chromatic",
    grepl("sensory",super_class) & grepl("^JO\\-",cell_type) ~ cell_type,
    grepl("sensory",super_class) & cell_function=="gustatory" ~ cell_sub_class,
    grepl("sensory",super_class) & cell_function=="thermosensory" ~ cell_sub_class,
    grepl("sensory",super_class) & cell_function=="hygrosensory" ~ cell_sub_class,
    grepl("sensory",super_class) & cell_function=="olfactory" & !is.na(cell_sub_class) ~ "olfactory",
    grepl("sensory",super_class) & cell_function=="olfactory" & !is.na(cell_sub_class) ~ cell_sub_class,
    cell_class == "MBON" ~ "MBON",
    cell_class == "CX_output" ~  "CX_output",
    super_class == "visual_centrifugal" ~ "visual_centrifugal",
    super_class == "visual_projection" ~ "visual_projection",
    super_class == "descending" ~ "descending",
    super_class == "ascending" ~ "ascending",
    grepl("efferent",super_class) & cell_class=="pars_intercerebralis" ~ "pars_intercerebralis",
    grepl("efferent", super_class) & cell_class=="pars_lateralis" ~ "pars_lateralis",
    super_class=="endocrine" ~ "endocrine",
    cell_class == "clock" ~ "clock",
    TRUE ~ cell_function
  )) %>%
  dplyr::mutate(seed = gsub(" |\\/|\\-","_",seed)) %>%
  dplyr::mutate(seed = dplyr::case_when(
    !is.na(seed) ~ paste(seed,side,collapse="_",sep="_"),
    TRUE ~ NA
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(region = dplyr::case_when(
    grepl("neck",region) ~ region,
    grepl("^LB|^MX|^MD|^TRd|^WEDa1|^WEDd1|^PSp2|^FLA",hemilineage) ~ "sez",
    grepl("^aPhN|^MxLbN|^MD|^ON$",nerve) ~ "sez",
    super_class == "visual_projection" ~ "optic_transfer",
    super_class == "visual_centrifugal" ~ "optic_transfer",
    TRUE ~ region
  )) %>%
  dplyr::mutate(flow = dplyr::case_when(
    super_class == "ascending" ~ "vnc_brain_transfer",
    super_class == "descending" ~ "brain_vnc_transfer",
    grepl("motor",super_class) ~ "brain_efferent",
    grepl("sensory",cell_class) ~ "brain_afferent",
    cell_class%in%c("descending_neuron") ~ "brain_vnc_transfer",
    cell_class%in%c("ascending_neuron") ~ "vnc_brain_transfer",
    flow=="intrinsic" ~ "brain_intrinsic",
    flow=="afferent" ~ "brain_afferent",
    flow=="efferent" ~ "brain_efferent",
    TRUE ~ flow
  ))

# Combine with extra columns from the fafbpipeline project
keep.cols <- c("root_id",setdiff(colnames(fw.meta.orig),colnames(fw.meta)))
keep.cols <- c("root_id",setdiff(keep.cols,colnames(ft)))
fw.meta.comb <- left_join(fw.meta,fw.meta.orig[,keep.cols], by = c("root_783"="root_id"))

# Choose potential GNG classes
presynapses.sez <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, "flywire_783_presynapses.parquet")) %>%
  dplyr::select(offset, inside, cleft_scores, pre_id, post_id, prepost) %>%
  dplyr::filter(cleft_scores >= 50)
postsynapses.sez <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, "flywire_783_postsynapses.parquet")) %>%
  dplyr::select(offset, inside, cleft_scores, pre_id, post_id, prepost) %>%
  dplyr::filter(cleft_scores >= 50)

# Determine classes
sez.out <- presynapses.sez %>%
  dplyr::group_by(pre_id) %>%
  dplyr::mutate(total_output = sum(prepost==0)) %>%
  dplyr::mutate(sez_output = sum(prepost==0&grepl("SAD|GNG|PRW|FLA|AMMC|CANT|WED|IPS",inside))) %>%
  dplyr::mutate(sez_output_prop = sez_output/total_output) %>%
  dplyr::mutate(superior_output = sum(prepost==0&grepl("LH|SMP|SLP|FB|EB|CRE|SCL|SIP|ICL",inside))) %>%
  dplyr::mutate(superior_output_prop = superior_output/total_output) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(root_783 = pre_id, total_output, sez_output, sez_output_prop, superior_output, superior_output_prop)
sez.in <- postsynapses.sez %>%
  dplyr::group_by(post_id) %>%
  dplyr::mutate(total_input = sum(prepost==1)) %>%
  dplyr::mutate(sez_input = sum(prepost==1&grepl("SAD|GNG|PRW|FLA|AMMC|CANT|WED|IPS",inside))) %>%
  dplyr::mutate(sez_input_prop = sez_input/total_input) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(root_783 = pre_id, total_input, sez_input, sez_input_prop)
sez.classes <- dplyr::full_join(sez.in,sez.out,by="root_783") %>%
  dplyr::left_join(fw.meta[,c("side","region","flow","super_class","cell_class","cell_sub_class","cell_type","root_783")], 
                   by = "root_783") %>%
  dplyr::filter(sez_output_prop>=0.1, sez_input_prop>=0.1) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(sez_class = dplyr::case_when(
    # GNG specific
    super_class=="ascending"&any(sez_output_prop>=0.5) ~ "sez_ascending",
    super_class=="motor"&any(sez_input_prop>=0.5) ~ "sez_motor",
    super_class=="endocrine"&(any(sez_output_prop>=0.5)|any(sez_input_prop>=0.5)) ~ "sez_endocrine",
    super_class=="descending"&any(sez_input_prop>=0.5) ~ "sez_descending",
    super_class=="sensory"&any(sez_output_prop>=0.5) ~ "sez_sensory",
    super_class=="sensory_ascending"&any(sez_output_prop>=0.5) ~ "sez_sensory_ascending",
    # GNG less specific
    super_class=="ascending" ~ "sez_midbrain_ascending",
    super_class=="motor" ~ "sez_midbrain_motor",
    super_class=="endocrine" ~ "sez_midbrain_endocrine",
    super_class=="descending" ~ "sez_midbrain_descending",
    super_class=="sensory" ~ "sez_midbrain_sensory",
    super_class=="sensory_ascending" ~ "sez_midbrain_sensory_ascending",
    super_class=="visual_centrifugal" ~ "sez_optic",
    # GNG general
    any(sez_output_prop>=0.5)&any(sez_input_prop>=0.5) ~ "sez_intrinsic",
    any(superior_output_prop>=0.5)&any(sez_input_prop>=0.5) ~ "sez_superior_projection",
    any(sez_output_prop<=0.5)&any(sez_input_prop>=0.5) ~ "sez_projection",
    any(sez_output_prop>=0.5)&any(sez_input_prop<=0.5) ~ "sez_centrifugal",
    any(sez_output_prop>=0.5)&any(sez_input_prop>=0.1) ~ "sez_peripheral_projection",
    any(sez_output_prop>=0.5)&any(sez_input_prop<=0.1) ~ "sez_peripheral_centrifugal",
    TRUE ~ "sez_peripheral"
  )) %>%
  dplyr::ungroup()

# Write
fw.meta.comb <- dplyr::left_join(fw.meta.comb,sez.classes[,c("root_783","sez_class")], by = "root_783")
readr::write_csv(fw.meta.comb, file.path(banc.meta.save.path,"flywire_meta.csv"))

# Write to feather
arrow::write_feather(fw.meta.comb,
                     file.path(banc.connectivity.save.path, "flywire_783_meta.feather"))
# sql_file <- file.path(banc.connectivity.save.path,"flywire_783_data.sqlite")
# con <- dbConnect(RSQLite::SQLite(), sql_file)
# DBI::dbWriteTable(con,
#                   name = "meta",
#                   value = fw.meta,
#                   overwrite = TRUE)
# dbDisconnect(con)

# Update franken meta seatable

# Update!
ft <- franken_meta("SELECT _id, fafb_id from franken_meta")
cluster.update <- sez.classes %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(sez_class = dplyr::case_when(
    all(is.na(sez_class)) ~ '',
    TRUE ~ sort(na.omit(sez_class))[1]
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(fafb_id = root_783, sez_group = sez_class) %>%
  dplyr::left_join(ft, by = "fafb_id") %>%
  dplyr::filter(!is.na(`_id`)) %>%
  dplyr::distinct() %>%
  as.data.frame()
banctable_update_rows(base='cns_meta', 
                      table = 'franken_meta', 
                      df = cluster.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

# Announce
message("##### BANCpipeline: FAFB-flywire meta updated #####")
message(sprintf("##### we have meta data for : %s neurons", nrow(fw.meta)))



