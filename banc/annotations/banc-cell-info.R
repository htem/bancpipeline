#' banc-cell-info — Compile and publish BANC cell-info reference tables to GCS.
#'
#' Brain hemilineage table, NT / neuropeptide controlled vocabularies, plus
#' other annotation reference CSVs.
#'
#' @section Reads:
#'   - FAFB `flytable_query` (hemilineage source)
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `banc_brain_hemilineages.csv` (local + GCS `meta/`)

############################
### ORGANISE BANC LABELS ###
############################
source("banc/banc-startup.R")

# Direct us to the BANC dataset
bancr::choose_banc()

########################################
### CREATE CELL TYPE REFERENCE TABLE ###
########################################

hl.orig <- flytable_query("SELECT ito_lee_hemilineage, hartenstein_hemilineage from info")
hl <- hl.orig %>%
  dplyr::filter(!is.na(ito_lee_hemilineage)|!is.na(hartenstein_hemilineage)) %>%
  dplyr::distinct(ito_lee_hemilineage, hartenstein_hemilineage)
readr::write_csv(hl, file.path(banc.meta.save.path,"banc_brain_hemilineages.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "banc_brain_hemilineages.csv"))

# Transmitters we care about
fast.nts <- c("acetylcholine", "glutamate",  "gaba", "glycine",
              "dopamine", "serotonin", "octopamine", "tyramine",
              "histamine", "nitric oxide")
neg.fast.nts <- c("acetylcholine-negative", "gaba-negative", "glutamate-negative",
                  "dopamine-negative", "serotonin-negative", "octopamine-negative",
                  "nitric oxide-negative", "histamine-negative", "tyramine-negative", "glycine-negative")
all.fast.nts <- c(fast.nts, neg.fast.nts)

# Function to process known_nt column
filter_words <- function(input_string, words_to_keep, invert = FALSE){
  words <- unlist(strsplit(input_string, ",|, |;|; "))
  words <- gsub("^ | $","",words)
  if (invert){
    filtered_words <- words[! words %in% words_to_keep]
  }else{
    filtered_words <- words[words %in% words_to_keep]
  }
  paste(unique(filtered_words), collapse = ", ")
}

# Get our franken brain data
franken.meta.orig <- franken_meta()

# Build our cell type reference table
ct.ref <- franken.meta.orig %>%
  dplyr::rename(manc_cell_type = MANC_type, 
                fafb_cell_type = FAFB_cell_type) %>%
  dplyr::select(-starts_with("_FAFB|_MANC")) %>%
  dplyr::mutate(neuromere = gsub(";.*","",neuromere),
                cell_function = gsub(",","",cell_function),
                super_class = gsub(",","",super_class),
                region = gsub(",","",region),
                manc_cell_type = gsub(";.*","",manc_cell_type),
                fafb_cell_type = gsub(";.*","",fafb_cell_type)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type)|cell_type=="" ~ cell_type,
    (grepl("ascending",super_class))&!is.na(manc_cell_type) ~ gsub(";.*","",manc_cell_type),
    (grepl("descending",super_class))&!is.na(fafb_cell_type) ~ gsub(";.*","",fafb_cell_type),
    TRUE ~ gsub(";.*","",cell_type)
  )) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(multi_neuromere = length(unique(neuromere))>1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(composite_cell_type = dplyr::case_when(
    is.na(neuromere)|neuromere=="" ~ cell_type,
    multi_neuromere ~ paste0(cell_type,"_",neuromere),
    TRUE ~ cell_type
  )) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cell_type = composite_cell_type,
                  parent_dataset = dataset,
                  region,
                  flow,
                  hemilineage,
                  nerve,
                  super_class,
                  cell_class,
                  cell_sub_class,
                  body_part_sensory,
                  body_part_effector,
                  neuromere,
                  cell_function, 
                  cell_function_detailed,
                  known_nt) %>%
  group_by(cell_type) %>%
  dplyr::summarize(
    parent_dataset = first(na.omit(parent_dataset)),
    region = first(na.omit(region)),
    flow = first(na.omit(flow)),
    hemilineage = first(na.omit(hemilineage)),
    nerve = first(na.omit(nerve)),
    super_class = first(na.omit(super_class)),
    cell_class = first(na.omit(cell_class)),
    cell_sub_class = first(na.omit(cell_sub_class)),
    body_part_sensory = first(na.omit(body_part_sensory)),
    body_part_effector = first(na.omit(body_part_effector)),
    neuromere = first(na.omit(neuromere)),
    cell_function = first(na.omit(cell_function)), 
    cell_function_detailed = first(na.omit(cell_function_detailed)),
    known_nt = first(na.omit(known_nt))
  ) %>%
  ungroup() %>%
  tidyr::separate_longer_delim(known_nt, delim = ";") %>%
  group_by(cell_type) %>%
  dplyr::mutate(known_nt = paste0(unique(na.omit(known_nt)),collapse=", ")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt = filter_words(known_nt, all.fast.nts)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(known_nt = ifelse(known_nt=="",NA,known_nt)) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::arrange(cell_type)
 
# Save in CAVE table format
ct.ref.m <- reshape2::melt(ct.ref, id = "cell_type", value.name = "tag3") %>%
  dplyr::rename(tag = cell_type, 
                tag2 = variable)

# Save
# banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
readr::write_csv(ct.ref, file.path(banc.meta.save.path,"banc_cell_type_reference_table_matrix.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "banc_cell_type_reference_table_matrix.csv"))
readr::write_csv(ct.ref.m, file.path(banc.meta.save.path,"banc_cell_type_reference_table_list.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "banc_cell_type_reference_table_list.csv"))

#################################
### MODIFY FOR CELL INFO TAGS ###
#################################

cell.info.updates <- franken.meta.orig %>%
  dplyr::rename(manc_cell_type = MANC_type, 
                fafb_cell_type = FAFB_cell_type) %>%
  dplyr::mutate(neuromere = gsub(";.*","",neuromere),
                cell_function = gsub(",","",cell_function),
                super_class = gsub(",","",super_class),
                region = gsub(",","",region),
                manc_cell_type = gsub(";.*","",manc_cell_type),
                fafb_cell_type = gsub(";.*","",fafb_cell_type)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type)|cell_type=="" ~ cell_type,
    (grepl("ascending",super_class))&!is.na(manc_cell_type) ~ gsub(";.*","",manc_cell_type),
    (grepl("descending",super_class))&!is.na(fafb_cell_type) ~ gsub(";.*","",fafb_cell_type),
    TRUE ~ gsub(";.*","",cell_type)
  )) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(multi_neuromere = length(unique(neuromere))>1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(composite_cell_type = dplyr::case_when(
    is.na(neuromere)|neuromere=="" ~ cell_type,
    multi_neuromere ~ paste0(cell_type,"_",neuromere),
    TRUE ~ cell_type
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(primary_class = dplyr::case_when(
    super_class=="sensory"&grepl("JO-",cell_type) ~ "Johnston's organ neuron",
    super_class=="sensory"&grepl("claw$",cell_function_detailed) ~ "claw chordotonal neuron",
    super_class=="sensory"&grepl("club$",cell_function_detailed) ~ "club chordotonal neuron",
    super_class=="sensory"&grepl("hook$",cell_function_detailed) ~ "hook chordotonal neuron",
    super_class=="sensory"&grepl("bristle$",cell_function_detailed) ~ "bristle mechanosensory neuron",
    super_class=="sensory"&grepl("peg$",cell_function_detailed) ~ "peg neuron",
    super_class=="sensory"&grepl("hair_plate$",cell_function_detailed) ~ "hair plate neuron",
    super_class=="sensory"&grepl("campaniform$",cell_function_detailed) ~ "campaniform sensillum neuron",
    super_class=="sensory"&grepl("hair_plate$",cell_function_detailed) ~ "campaniform sensillum neuron",
    super_class=="sensory"&grepl("stretch$",cell_function_detailed) ~ "campaniform sensillum neuron",
    super_class=="sensory"&grepl("olfactory$",cell_function) ~ "olfactory receptor neuron",
    super_class=="sensory"&grepl("gustatory|chemosensory|pheromone",cell_function) ~ "gustatory neuron",
    super_class=="sensory"&grepl("thermo$",cell_function) ~ "thermosensory neuron",
    super_class=="sensory"&grepl("hygrosensory$",cell_function) ~ "hygrosensory neuron",
    super_class=="sensory"&grepl("ocell$",cell_function) ~ "ocellar neuron",
    super_class=="sensory"&grepl("visual$",cell_function) ~ "photoreceptor neuron",
    super_class=="sensory"&grepl("nocicept$",cell_function) ~ "nociceptive neuron",
    super_class=='visual_centrifugal' ~ 'visual projection',
    super_class=="visual_centrifugal" ~  'visual centrifugal', 
    grepl("ascending",super_class) ~ 'ascending',
    grepl("descending",super_class) ~ 'descending',
    grepl("endocrine",super_class) ~ 'endocrine',
    grepl("motor",super_class) ~ 'motor neuron',
    grepl("brain_intrinsic",super_class) ~ 'central brain intrinsic',
    grepl("vnc_intrinsic",super_class) ~ 'VNC intrinsic',
    TRUE ~ NA
  )) %>%
  dplyr::mutate(body_part_innervated = dplyr::case_when(
    !grepl('endocrine',super_class) ~ NA,
    !is.na(body_part_sensory) ~ paste0("innervates ", gsub("_"," ",body_part_sensory)),
    !is.na(body_part_effector) ~ paste0("innervates ", gsub("_"," ",body_part_effector)),
    TRUE ~ NA
  )) %>%
  tidyr::separate_longer_delim(known_nt, delim = ";") %>%
  group_by(cell_type) %>%
  dplyr::mutate(known_nt = paste0(unique(na.omit(known_nt)),collapse=", ")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt = filter_words(known_nt, fast.nts)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fast_neurotransmitter = dplyr::case_when(
    grepl("acetylcholine",known_nt) ~ 'cholinergic',
    grepl("GABA",known_nt) ~ 'GABAergic',
    grepl("glutamate",known_nt) ~ 'glutamatergic',
    TRUE ~ NA)) %>%
  dplyr::mutate(other_neurotransmitter = dplyr::case_when(
    grepl("octopamine",known_nt) ~ 'octopaminergic',
    grepl("serotonergic",known_nt) ~ 'serotonergic',
    grepl("dopamine",known_nt) ~ 'dopaminergic',
    grepl("tyramine",known_nt) ~ 'tyraminergic',
    grepl("histamine",known_nt) ~ 'histaminergic',
    grepl("glycine",known_nt) ~ 'glycinergic',
    TRUE ~ NA)) %>%
  dplyr::distinct(cell_type,
                  primary_class,
                  fast_neurotransmitter,
                  other_neurotransmitter,
                  body_part_innervated,
                  hemilineage) 

# Reformat
cell.info.updates.for.cave <- reshape2::melt(cell.info.updates, id = "cell_type", variable = "tag2") %>%
  dplyr::rename(tag = value) %>%
  dplyr::filter(!is.na(tag),
                !grepl("Unknown|unknown|NA|None|none",tag)) %>%
  dplyr::mutate(user_id = 355)

# Save
#banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
readr::write_csv(cell.info.updates.for.cave, file.path(banc.meta.save.path,"banc_cell_info_updates.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "banc_cell_info_updates.csv"))

##########################
### COMPILE MATCH DATA ###
##########################

# Compile all match data
# Non-mirror feathers now use CAVE column names (validation, query_root_id);
# convert back to internal names (valid, query_id) for compatibility with CSV data
.cave_to_internal <- function(df) {
  if ("query_root_id" %in% names(df)) df <- dplyr::rename(df, query_id = query_root_id)
  if ("validation" %in% names(df)) {
    df$valid <- ifelse(df$validation, 't', 'f')
    df$validation <- NULL
  }
  df
}
banc.meta.fafb.nb <- .cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather")))
banc.meta.manc.nb <- .cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather")))
banc.meta.hemibrain.nb <- .cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather")))
banc.meta.fanc.nb <- .cave_to_internal(arrow::read_feather(file.path(banc.meta.save.path, "banc_fanc_1116_nblast.feather")))
banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_mirror_nblast.feather"))

######################################
### BANC seatable labels from CAVE ###
######################################

# Sort out data to share, human verified:
bc <- banctable_query()
bc.cts <- bc %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON",status)) %>%
  dplyr::filter(!is.na(root_id), 
                !is.na(supervoxel_id), 
                !is.na(position), 
                root_id!="", 
                root_id!="0") %>%
  dplyr::arrange(
    fafb_cell_type,
    manc_cell_type,
    hemibrain_cell_type,
    cell_type
  ) %>%
  dplyr::filter(!duplicated(root_id)) %>%
  dplyr::mutate(
    fafb_cell_type = dplyr::case_when(
      grepl("auto",fafb_cell_type) ~ NA,
      is.na(fafb_match) ~ NA,
      TRUE ~ fafb_cell_type
    ),
    manc_cell_type = dplyr::case_when(
      grepl("auto",manc_cell_type) ~ NA,
      is.na(manc_match) ~ NA,
      TRUE ~ manc_cell_type
    ),
    hemibrain_cell_type = dplyr::case_when(
      grepl("auto",hemibrain_cell_type) ~ NA,
      is.na(hemibrain_match) ~ NA,
      TRUE ~ hemibrain_cell_type
    )
  ) %>%
  dplyr::filter(!is.na(fafb_match)|!is.na(manc_match)|!is.na(hemibrain_match)) %>%
  dplyr::select(root_id, supervoxel_id, position, fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::arrange(fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::distinct(root_id, fafb_cell_type, manc_cell_type, hemibrain_cell_type, .keep_all = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, fafb_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, manc_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, hemibrain_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!duplicate_flag) %>%
  dplyr::rename(pt_root_id = root_id, 
                pt_supervoxel_id = supervoxel_id, 
                pt_position = position) %>%
  dplyr::select(-duplicate_flag) %>%
  reshape2::melt(id = c("pt_root_id",  "pt_supervoxel_id", "pt_position"),
                 value.name = "tag",
                 variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA","unknown","Uknown","fragment","Fragment","None","none","no_match"),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity",
                valid = "t") %>%
  dplyr::mutate(tag2 = gsub("_"," ",tag2)) %>%
  dplyr::filter(!grepl("auto",tag), !grepl("no_match",tag), 
                !grepl("auto",tag2), !grepl("no_match",tag2), 
                !grepl("None|^BM",tag),
                !is.na(tag), !is.na(tag2)) %>%
  dplyr::distinct(pt_root_id, pt_supervoxel_id, pt_position, tag2, tag, user_id, valid)

####################################
### High NBLAST labels from CAVE ###
####################################

# Sort out data to share, NBLAST high:
nb.thresh <- 0.5
fw.nblast.scores <- banc.meta.fafb.nb %>% # banc.meta.fafb.nb
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh, !is.na(match_cell_type)) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(fafb_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
mc.nblast.scores <- banc.meta.manc.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh, !is.na(match_cell_type)) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(manc_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
hb.nblast.scores <- banc.meta.hemibrain.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh, !is.na(match_cell_type)) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(hemibrain_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
fc.nblast.scores <- banc.meta.fanc.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh, !is.na(match_cell_type)) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(fanc_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
nblast.scores <- fw.nblast.scores %>%
  dplyr::full_join(mc.nblast.scores, by = c("pt_root_id","pt_supervoxel_id","pt_position")) %>%
  dplyr::full_join(hb.nblast.scores, by = c("pt_root_id","pt_supervoxel_id","pt_position")) %>%
  dplyr::full_join(fc.nblast.scores, by = c("pt_root_id","pt_supervoxel_id","pt_position")) %>%
  dplyr::filter(!is.na(pt_root_id), !is.na(pt_supervoxel_id), !is.na(pt_position), pt_root_id!="", pt_root_id!="0") %>%
  dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, 
                fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  reshape2::melt(id = c("pt_root_id",  "pt_supervoxel_id", "pt_position"),
                 value.name = "tag",
                 variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA","unknown","Uknown","fragment","Fragment","None","none","no_match"),
                !grepl("None",tag),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity",
                tag = paste0(tag,"?"),
                user_id = 355,
                valid = "t") %>% 
  dplyr::filter(!is.na(tag), !is.na(tag2))
nblast.scores <- bancr::banc_updateids(nblast.scores, root.column = "pt_root_id",
                                       supervoxel.column = "pt_supervoxel_id",
                                       position.column = "pt_position")
nblast.scores <- nblast.scores %>% 
  dplyr::group_by(pt_root_id, tag2) %>%
  dplyr::filter(!duplicated(pt_position), 
                !duplicated(pt_root_id),
                pt_root_id!="0") %>%
  dplyr::ungroup()

#############################
### Read labels from CAVE ###
#############################

# Read from CAVE
banc.cell.info <- banc_cell_info(rawcoords = TRUE) %>%
  dplyr::filter(user_id %in% c(355), 
                tag2 == "neuron identity",
                valid == "t") %>%
  dplyr::mutate(pt_supervoxel_id = as.character(pt_supervoxel_id),
                pt_root_id = as.character(pt_root_id)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = hemibrainr:::paste_coords(pt_position),
                pt_position = gsub("\\(|\\)","",pt_position)) %>%
  dplyr::select(id, valid, tag, tag2, user_id, pt_supervoxel_id, pt_root_id, pt_position)

#######################################
### Make modifications to my labels ###
#######################################

# Get the CAVE rows we need to flip to FALSE for the high NBLAST matches
banc.cell.info.nblast.false <- banc.cell.info %>%
  dplyr::filter(grepl("\\?",tag)) %>%
  dplyr::anti_join(nblast.scores, by = c("pt_root_id","tag")) %>%
  dplyr::mutate(valid = "f")
high.nblast <- rbind.fill(banc.cell.info.nblast.false,
                          banc.cell.info %>%
                            dplyr::filter(grepl("\\?",tag)) %>%
                            dplyr::semi_join(nblast.scores, by = c("pt_root_id","tag")),
                          nblast.scores %>%
                            dplyr::filter(grepl("\\?",tag)) %>%
                            dplyr::anti_join(banc.cell.info, by = c("pt_root_id","tag")))

# Get the CAVE rows we need to flip to FALSE for the manual matches
banc.cell.info.verified.false <- banc.cell.info %>%
  dplyr::filter(!grepl("\\?",tag)) %>%
  dplyr::anti_join(bc.cts, by = c("pt_root_id","tag")) %>%
  dplyr::mutate(valid = "f")
matches.verified <- rbind.fill(banc.cell.info.verified.false,
                               banc.cell.info %>%
                                 dplyr::filter(!grepl("\\?",tag)) %>%
                                 dplyr::semi_join(bc.cts, by = c("pt_root_id","tag")),
                               bc.cts %>%
                                 dplyr::filter(!grepl("\\?",tag)) %>%
                                 dplyr::anti_join(banc.cell.info, by = c("pt_root_id","tag")))

# Decommissioned 2026-05-13: banc_cell_type_verified.csv and
# banc_cell_type_high_nblast_match.csv have no readers (bancpipeline or
# BANC-project); the GCS copies were 15-month-old curation artifacts. Local
# writes and GCS pushes removed. `matches.verified` and `high.nblast` are
# still computed above in case they're needed in an interactive session.

#########################################
### Update cell type labels from CAVE ###
#########################################

# Updated cell_type_source column based on CAVE
bc.all <- banctable_query("SELECT _id, root_id, cell_type, other_names, super_class, cell_class, cell_type_source from banc_meta")
banc.cell.info <- banc_cave_cell_types(cave_id = c("392","103", "104", "7", "5476","13","17","116","1152","847","60","1234","3153","188","96","27"), 
                                       invert = FALSE)
banc.cell.info.mod <- banc.cell.info %>%
  dplyr::filter(!user_id %in% c("355","52"), 
                cell_type!="",
                !is.na(cell_type_source),
                cell_type_source!="611") %>%
  dplyr::mutate(root_id=as.character(root_id)) %>%
  dplyr::distinct(root_id, cell_type, cell_type_source,user_id)

# Process
bc.ct <- bc.all %>%
  dplyr::filter(root_id %in% banc.cell.info.mod$root_id) %>%
  dplyr::left_join(banc.cell.info.mod,
                   by = "root_id") %>%
  dplyr::mutate(
    other_names = ifelse(is.na(other_names),'',other_names),
    cell_type_source.y = gsub("Rachel Wilson Lab", "Wilson lab", cell_type_source.y),
    cell_type_source.y = ifelse(is.na(cell_type_source.y),NA,tolower(cell_type_source.y)),
    cell_type_source.x = ifelse(is.na(cell_type_source.x),NA,tolower(cell_type_source.x)),
    cell_type_source.x = ifelse(grepl("NA|na|None",cell_type_source.x),NA,cell_type_source.x),
    cell_type_source.x = ifelse(cell_type_source.x%in%c("","NA"),NA,cell_type_source.x),
    cell_type_source.y = ifelse(cell_type_source.y%in%c("","NA"),NA,cell_type_source.y)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type.x)|grepl("auto",cell_type.x) ~ cell_type.y,
    is.na(cell_type.y)|grepl("\\?",cell_type.x) ~ cell_type.x,
    TRUE ~ cell_type.x),
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(other_names = dplyr::case_when(
    !is.na(cell_type.y) & (cell_type.y!= cell_type) ~ paste(setdiff(sort(unique(c(unlist(strsplit(other_names,split=", "),cell_type.y)))),cell_type),collapse=", "),
    !is.na(cell_type.x) & (cell_type.x!= cell_type) ~ paste(setdiff(sort(unique(c(unlist(strsplit(other_names,split=", "),cell_type.x)))),cell_type),collapse=", "),
    TRUE ~ paste(setdiff(sort(unique(c(unlist(strsplit(other_names,split=", "))))),cell_type),collapse=", ")
  )) %>%
  dplyr::mutate(
    cell_type_source = dplyr::case_when(
      is.na(cell_type_source.x) ~ cell_type_source.y,
      is.na(cell_type_source.y) ~ cell_type_source.x,
      cell_type_source.x=="NA" ~ cell_type_source.y,
      cell_type_source.y=="NA" ~ cell_type_source.x,
      cell_type_source.x=="cave"&!is.na(cell_type_source.y) ~ cell_type_source.y,
      cell_type_source.x=="community"&!is.na(cell_type_source.y) ~ cell_type_source.y,
      cell_type_source.x==""&!is.na(cell_type_source.y) ~ cell_type_source.y,
      !is.na(cell_type_source.x)&!is.na(cell_type_source.y) ~ paste(tolower(sort(unique(unlist(c(strsplit(cell_type_source.x,split=",")),cell_type_source.y)),
                                                                         decreasing=TRUE)),
                                                                    collapse=","),
      TRUE ~ cell_type_source.x
    )) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("\\?|intrinsic|sensory|glia|ascending|descending|central|optic|,|GABA|Ach|Glut",cell_type)) %>%
  dplyr::filter(!is.na(cell_type_source), cell_type_source!="") %>%
  dplyr::distinct(`_id`, root_id, .keep_all = TRUE) %>%
  dplyr::select(`_id`, root_id, cell_type, 
                other_names, cell_type_source,
                super_class, cell_class) %>%
  dplyr::mutate(other_names = gsub("^,|^ ,|^ ","",other_names),
                cell_type_source = ifelse(cell_type_source=='151184',
                                          NA,
                                          cell_type_source))
bc.new.cts <- bc.ct %>%
  dplyr::anti_join(bc.all,by=c("_id","root_id","cell_type"))

# Add cell type source labels
if(nrow(bc.new.cts)){
  bc.update <- as.data.frame(bc.new.cts)
  bc.update$cell_type_source <- tolower(bc.update$cell_type_source)
  bc.update <- subset(bc.update, !is.na(super_class))
  bc.update[is.na(bc.update)] <- ''
  banctable_update_rows(base='banc_meta',
                        table = "banc_meta",
                        df = bc.update[,c("_id","cell_type", "other_names", "cell_type_source")],
                        append_allowed = FALSE,
                        chunksize = 1000)
}

############################
### Update super classes ###
############################

# banc.cell.class.info <- banc_cave_cell_types()
# bc.sp <- bc.all %>%
#   dplyr::filter(root_id %in% banc.cell.class.info$root_id) %>%
#   dplyr::left_join(banc.cell.class.info %>%
#                      dplyr::mutate(root_id = as.character(root_id)) %>%
#                      dplyr::select(root_id, super_class) %>%
#                      dplyr::distinct(root_id, .keep_all = TRUE),
#                    by = "root_id") %>%
#   dplyr::mutate(super_class = dplyr::case_when(
#     !is.na(super_class.x)&super_class.x!="" ~ super_class.x,
#     grepl("sensory",super_class.y) ~ "sensory",
#     grepl("motor",super_class.y) ~ "motor",
#     grepl("endocrine",super_class.y) ~ "endocrine",
#     grepl("glia",super_class.y) ~ "glia",
#     grepl("optic",super_class.y) ~ "optic",
#     TRUE ~ super_class.x
#   )) %>%
#   dplyr::select(`_id`,root_id,super_class)
# bc.sp <- bc.sp %>%
#   dplyr::anti_join(bc.all,by=c("_id","super_class"))
# 
# if(nrow(bc.sp)){
#   bc.update <- as.data.frame(bc.sp)
#   bc.update[is.na(bc.update)] <- ''
#   banctable_update_rows(base='banc_meta',
#                         table = "banc_meta",
#                         df = bc.update,
#                         append_allowed = FALSE,
#                         chunksize = 1000)
# }

###################################
### Update our verified matches ###
###################################

# Older results
fafb.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"),
                                         col_types = banc.col.types, 
                                         show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id))
manc.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"),
                                         col_types = banc.col.types, 
                                         show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))
hemibrain.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"),
                                              col_types = banc.col.types,
                                              show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))

# Update fafb
match.ids <- na.omit(unique(banc.meta.fafb.nb$match_cell_type))
banc.meta.fafb.nb.2 <- banc.meta.fafb.nb %>%
  dplyr::left_join(bc.ct %>% 
                     dplyr::filter(cell_type %in% match.ids,
                                   !is.na(cell_type)) %>%
                     dplyr::select(cell_type, root_id) %>% 
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::mutate(in_df2 = TRUE),
                   by = c("match_cell_type"="cell_type", "pt_root_id"="root_id")) %>%
  dplyr::mutate(valid = ifelse(!is.na(in_df2), 't', valid)) %>%
  dplyr::select(-in_df2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(valid = case_when(
    valid == 't' & cumsum(valid == 't') > 1 ~ 'f',
    TRUE ~ valid
  )) %>%
  dplyr::ungroup()
fafb.matches.df.valid.2 <- banc.meta.fafb.nb.2 %>%
  dplyr::filter(valid=='t', score>0.3) %>%
  dplyr::anti_join(fafb.matches.df.valid,
                   by = c("pt_root_id","match_cell_type")) %>%
  dplyr::select(-score) %>%
  rbind(fafb.matches.df.valid) %>%
  dplyr::distinct()

# Update manc
match.ids <- na.omit(unique(banc.meta.manc.nb$match_cell_type))
banc.meta.manc.nb.2 <- banc.meta.manc.nb %>%
  dplyr::left_join(bc.ct %>% 
                     dplyr::filter(cell_type %in% match.ids,
                                   !is.na(cell_type)) %>%
                     dplyr::select(cell_type, root_id) %>% 
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::mutate(in_df2 = TRUE),
                   by = c("match_cell_type"="cell_type", "pt_root_id"="root_id")) %>%
  dplyr::mutate(valid = ifelse(!is.na(in_df2), 't', valid)) %>%
  dplyr::select(-in_df2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(valid = case_when(
    valid == 't' & cumsum(valid == 't') > 1 ~ 'f',
    TRUE ~ valid
  )) %>%
  dplyr::ungroup()
manc.matches.df.valid.2 <- banc.meta.manc.nb.2 %>%
  dplyr::filter(valid=='t', score>0.3) %>%
  dplyr::anti_join(manc.matches.df.valid,
                   by = c("pt_root_id","match_cell_type")) %>%
  dplyr::select(-score) %>%
  rbind(manc.matches.df.valid) %>%
  dplyr::distinct()

# Update hemibrain
match.ids <- na.omit(unique(banc.meta.hemibrain.nb$match_cell_type))
banc.meta.hemibrain.nb.2 <- banc.meta.hemibrain.nb %>%
  dplyr::left_join(bc.ct %>% 
                     dplyr::filter(cell_type %in% match.ids,
                                   !is.na(cell_type)) %>%
                     dplyr::select(cell_type, root_id) %>% 
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::mutate(in_df2 = TRUE),
                   by = c("match_cell_type"="cell_type", "pt_root_id"="root_id")) %>%
  dplyr::mutate(valid = ifelse(!is.na(in_df2), 't', valid)) %>%
  dplyr::select(-in_df2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(valid = case_when(
    valid == 't' & cumsum(valid == 't') > 1 ~ 'f',
    TRUE ~ valid
  )) %>%
  dplyr::ungroup()
hemibrain.matches.df.valid.2 <- banc.meta.hemibrain.nb.2 %>%
  dplyr::filter(valid=='t', score>0.3) %>%
  dplyr::anti_join(hemibrain.matches.df.valid,
                   by = c("pt_root_id","match_cell_type")) %>%
  dplyr::select(-score) %>%
  rbind(hemibrain.matches.df.valid) %>%
  dplyr::distinct()

# Write
readr::write_csv(fafb.matches.df.valid.2, file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"))
readr::write_csv(manc.matches.df.valid.2, file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"))
readr::write_csv(hemibrain.matches.df.valid.2, file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"))





