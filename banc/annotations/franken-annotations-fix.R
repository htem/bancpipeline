#' franken-annotations-fix — Standardise franken-brain annotations against published labels.
#'
#' Diff `franken_meta()` cell_types against the published FAFB Supplemental
#' table and apply targeted corrections to the franken-brain CSV (column
#' harmonisation, dataset-specific fix blocks).
#'
#' @section Reads:
#'   - `franken_meta()`
#'   - `BANC-project/data/banc_annotations/Supplemental_file1_neuron_annotations.tsv`
#'
#' @section Writes:
#'   - `<banc.meta.save.path>/frankenbrain_v*_meta.csv` (per-fix-block)
#'
#' @section Notes:
#'   - Manual; many blocks are currently commented out — uncomment to apply.

############################################
### STANDARDISE THE ANNOTATIONS FOR BANC ###
############################################
source("banc/banc-startup.R")

# get current table
franken.meta <- franken_meta()

#######################################
## DIFFERENCE WITH PUBLISHED LABELS ###
#######################################

# # read published data
# fafb.pub <- read_tsv("/Users/papers/BANC-project/data/banc_annotations/Supplemental_file1_neuron_annotations.tsv",
#                      col_types = banc.col.types)
# 
# # cell type difference
# franken.cts <- franken.meta %>%
#   dplyr::filter(dataset=="FAFB", region!="optic_lobe") %>%
#   dplyr::distinct(cell_type) %>%
#   dplyr::pull(cell_type)
# fafb.cts <- fafb.pub %>%
#   dplyr::mutate(cell_type = dplyr::case_when(
#     is.na(cell_type)&!is.na(hemibrain_type) ~ hemibrain_type,
#     is.na(cell_type)&!is.na(morphology_group) ~ morphology_group,
#     TRUE ~ cell_type
#   )) %>%
#   dplyr::distinct(cell_type) %>%
#   dplyr::pull(cell_type)
# 
# # diff
# ct.diff <- setdiff(franken.cts, fafb.cts)

######################################
## ADD MISSING NEURONS SINCE PRINT ###
######################################
# ft.midbrain <- fafbseg::flytable_query("select root_id, root_783, supervoxel_id, nucleus_id, super_class, cell_class, cell_sub_class, top_nt, side from info")
# ft.optic <- fafbseg::flytable_query("select root_id, root_783, supervoxel_id, nucleus_id, super_class, cell_class, cell_sub_class, top_nt, side from optic")
# ft.glia <- rbind(ft.midbrain,ft.optic) %>%
#   dplyr::filter(grepl("glia",super_class)|grepl("glia",cell_class)|grepl("glia",cell_sub_class)) %>%
#   dplyr::distinct(root_id, .keep_all = TRUE)
# readr::write_csv(ft.glia, file = "/Users/abates/Downloads/fafb_flywire_glia.csv")

# ft.midbrain <- fafbseg::flytable_query("select _id, root_783, supervoxel_id, nucleus_id, flow, super_class, cell_class, cell_sub_class, cell_type, hemibrain_type, ito_lee_hemilineage, hartenstein_hemilineage, morphology_group, top_nt, top_nt_conf, side, nerve, vfb_id, fbbt_id, status, soma_x from info")
# ft.midbrain$region <- "midbrain"
# ft.optic <- fafbseg::flytable_query("select * from optic")
# ft.optic$region <- "optic_lobe"
# ft.update <- ft.midbrain %>%
#   rbind.fill(ft.optic.update) %>%
#   dplyr::filter(!root_783 %in% na.omit(unique(franken.meta$fafb_id))) %>%
#   dplyr::filter(!cell_class %in% c("fragment","trachea","large_fragment","ECM","tadpole")) %>%
#   dplyr::filter(!super_class %in% c("not_a_neuron")) %>%
#   dplyr::filter(!grepl("glia",cell_class)) %>%
#   dplyr::rename(neuron_id = root_783,
#                 FAFB_supervoxel_id=supervoxel_id, 
#                 FAFB_nucleus_id=nucleus_id, 
#                 FAFB_flow=flow, 
#                 FAFB_super_class=super_class, 
#                 FAFB_cell_class=cell_class, 
#                 FAFB_cell_sub_class=cell_sub_class, 
#                 FAFB_cell_type=cell_type, 
#                 FAFB_hemibrain_type=hemibrain_type, 
#                 FAFB_ito_lee_hemilineage=ito_lee_hemilineage, 
#                 FAFB_hartenstein_hemilineage=hartenstein_hemilineage, 
#                 FAFB_morphology_group=morphology_group, 
#                 FAFB_top_nt=top_nt, 
#                 FAFB_top_nt_conf=top_nt_conf, 
#                 FAFB_side=side, 
#                 FAFB_nerve=nerve, 
#                 FAFB_vfb_id=vfb_id, 
#                 FAFB_fbbt_id=fbbt_id, 
#                 FAFB_status=status) %>%
#   dplyr::mutate(fafb_id = neuron_id,
#                 flow=FAFB_flow, 
#                 super_class=FAFB_super_class, 
#                 cell_type=dplyr::case_when(
#                   !is.na(matsliah_type) ~ matsliah_type,
#                   TRUE ~ FAFB_cell_type
#                 ), 
#                 cell_class = FAFB_cell_class,
#                 cell_sub_class=FAFB_cell_sub_class, 
#                 hemilineage=dplyr::case_when(
#                   !is.na(FAFB_ito_lee_hemilineage) ~ FAFB_ito_lee_hemilineage,
#                   TRUE ~ FAFB_hartenstein_hemilineage
#                 ), 
#                 top_nt=FAFB_top_nt, 
#                 side=FAFB_side, 
#                 nerve=FAFB_nerve) %>%
#   dplyr::select(neuron_id,
#                 fafb_id,
#                 region,
#                 flow,
#                 super_class,
#                 cell_class,
#                 cell_sub_class,
#                 cell_type,
#                 top_nt,
#                 side,
#                 nerve,
#                 hemilineage,
#                 FAFB_supervoxel_id, 
#                 FAFB_nucleus_id, 
#                 FAFB_flow, 
#                 FAFB_super_class, 
#                 FAFB_cell_class, 
#                 FAFB_cell_sub_class, 
#                 FAFB_cell_type, 
#                 FAFB_hemibrain_type, 
#                 FAFB_ito_lee_hemilineage, 
#                 FAFB_hartenstein_hemilineage, 
#                 FAFB_morphology_group, 
#                 FAFB_top_nt, 
#                 FAFB_top_nt_conf, 
#                 FAFB_side, 
#                 FAFB_nerve, 
#                 FAFB_vfb_id, 
#                 FAFB_fbbt_id, 
#                 FAFB_status) %>%
#   dplyr::filter(!is.na(cell_type))
# ft.update[is.na(ft.update)] <- ''
# banctable_append_rows(base='cns_meta',
#                      table = "franken_meta",
#                      df = as.data.frame(ft.update),
#                      chunksize = 1000)

##################################
## USE ARIE'S OPTIC LOBE TYPES ###
##################################
# Change to matisliah types
# ft.optic <- fafbseg::flytable_query("select root_783, matsliah_type from optic")

# allowed neurotransmitters (for neurotransmitters_verified column)
allowed_neurotransmitters <- c(
  "acetylcholine", "gaba", "glutamate", "glycine", "dopamine", 
  "serotonin", "histamine", "tyramine", "octopamine", "nitric_oxide"
)

# execute changes to annotations
# documented here: https://github.com/wilson-lab/bancpipeline/blob/main/annotations/annotations.md
franken.meta.update <- franken.meta %>%
  ### CELL TYPE CHANGES ####
# # cell type left join
# dplyr::left_join(ft.optic %>% 
#                    dplyr::distinct(root_783, 
#                                  matsliah_type), 
#                  by = c("fafb_id"="root_783")) %>%
# # cell_type 
# dplyr::mutate(cell_type = dplyr::case_when(
#   cell_type=="ocellar retinula cell" ~ "ocellar_retinula_cell",
#   !is.na(matsliah_type) ~ matsliah_type,
#   TRUE ~ cell_type
# )) %>%
# dplyr::select(-matsliah_type) %>%
### MOVING TO NEW POLICY ####
# citation 
dplyr::mutate(citation_cell_type = dplyr::case_when(
  !is.na(FAFB_cell_type)&cell_type==FAFB_cell_type ~ "Schlegel et al., 2024 (cell type)",
  !is.na(MANC_type)&cell_type==MANC_type ~ "Marin et al., 2024 (cell type)",
  TRUE ~ NA
) ) %>%
  # flow
  dplyr::mutate(flow = dplyr::case_when(
    cell_type %in%  c("AN_4_22ac4f12", "AN_4_63", "AN_4_None", "AN_4_c964e2c2", "AN_GNG_69", # corrections
                      "SA_MDA_2", "SA_MDA_3", "Dm9", "L4", "C3", "Dm11", "L3", "Lawf1", "Lawf2") ~ "intrinsic", 
    grepl("afferent",flow) ~ "afferent",
    grepl("sensory",super_class)&!grepl("efferent|intrinsic",flow) ~ "afferent",
    grepl("^R1-6$|^R7$|^R8$", cell_type) ~ "afferent",
    grepl("motor|endocrine|efferent",super_class) ~ "efferent",
    grepl("efferent",flow) ~ "efferent",
    grepl("^ascending$",super_class) ~ "intrinsic",
    grepl("^descending$",super_class) ~ "intrinsic",
    grepl("ascending|descending",super_class) ~ "intrinsic",
    grepl("intrinsic",flow) ~ "intrinsic",
    is.na(nerve) ~ "intrinsic",
    # none get NA, as desired
    TRUE ~ NA 
  ) ) %>%
  # side
  dplyr::mutate(side = dplyr::case_when(
    grepl("right",side) ~ "right",
    grepl("left",side) ~ "left",
    grepl("center|midline|unpaired|both",side) ~ "center", # midline and unpaired unused
    # small number get NA, some sensory and others with missing class info
    TRUE ~ NA
  ) ) %>%
  # region — "neck_connective" disabled as a region value. Ascending/descending
  # neurons are identified via super_class; cervical_connective is a tract entry.
  dplyr::mutate(region = dplyr::case_when(
    grepl("visual_centrifugal",super_class) ~ "central_brain",
    (grepl("visual_projection|optic",super_class) &
       !grepl("ocellar_interneuron",super_class)) ~ "optic_lobe",
    grepl("ocellar",super_class) ~ "central_brain",
    grepl("ocellar",cell_class) ~ "central_brain",
    grepl("vnc|VNC",region) ~ "ventral_nerve_cord",
    grepl("midbrain|central_brain|sez",region) ~ "central_brain",
    grepl("OL|optic",region) ~ "optic_lobe",
    TRUE ~ NA
  ) ) %>%
  dplyr::mutate(body_part_sensory = dplyr::case_when(
    grepl("^R1-6$|^R7$|^R8$",cell_type) ~ "retina",
    cell_type %in% "R1-6" ~ "retina",
    cell_type %in% c("R7","R8") ~ "retina",
    !grepl("sensory|visceral",super_class) ~ NA,
    grepl("^ISN$",cell_type) ~ "hemolymph_sensory",
    # head
    ## head bristles
    grepl("BM_Ant",cell_type) ~ "antenna",
    grepl("BM_dOcci",cell_type) ~ "occipital_dorsal",
    grepl("BM_dPoOr",cell_type) ~ "postorbital_dorsal",
    grepl("BM_FrOr", cell_type) ~ "frontoorbital",
    grepl("BM_Fr",cell_type) ~ "frontal",
    grepl("BM_Hau",cell_type) ~ "haustellum",
    grepl("BM_InOc",cell_type) ~ "interocellar",
    grepl("BM_InOm",cell_type) ~ "interommatidial",
    grepl("BM_MaPa",cell_type) ~ "maxillary_palp",
    grepl("BM_Oc",cell_type) ~ "ocellar",
    grepl("BM_Or",cell_type) ~ "orbital",
    grepl("BM_Taste",cell_type) ~ "labellum",
    grepl("BM_Vib",cell_type) ~ "vibrissa",
    grepl("BM_vOcci_vPoOr",cell_type) ~ "postorbital_ventral",
    grepl("BM_Vt_PoOc",cell_type) ~ "postocellar",
    ## by prior annotation
    grepl("^eye$",body_part_sensory) ~ "eye",
    grepl("head$",body_part_sensory) ~ "head",
    grepl("antenna",nerve) ~ "antenna",
    grepl("antenna",body_part_sensory) ~ "antenna",
    grepl("aorta",body_part_sensory) ~ "aorta",
    grepl("CB0991",cell_type) ~ "aorta",
    grepl("crop",body_part_sensory) ~ "crop",
    grepl("cibarium",body_part_sensory) ~ "cibarium",
    grepl("retina",body_part_sensory) ~ "retina",
    cell_type %in% "HBeyelet" ~ "eyelet",
    grepl("ocellar_retinula_cell",cell_type) ~ "ocellus",
    grepl("ocellar|ocelli",body_part_sensory)&peripheral_target_type!="bristle" ~ "ocellus",
    grepl("ocellar|ocelli",cell_function)&peripheral_target_type!="bristle" ~ "ocellus",
    grepl("ocellar|ocelli",cell_function_detailed)&peripheral_target_type!="bristle" ~ "ocellus",
    grepl("vibrissa",body_part_sensory) ~ "vibrissa",
    grepl("labellum",body_part_sensory) ~ "labellum",
    grepl("proboscis",cell_function)&flow=="efferent" ~ "proboscis",
    grepl("pharynx",cell_function)&flow=="efferent" ~ "pharynx",
    grepl("proboscis-pharynx",body_part_sensory) ~ "proboscis-pharynx",
    grepl("proboscis",body_part_sensory) ~ "proboscis",
    grepl("pharynx",body_part_sensory) ~ "pharynx",
    grepl("palps",body_part_sensory) ~ "maxillary_palp",
    grepl("eye",body_part_sensory) ~ "eye",
    # CNS
    grepl("pars_intercerebralis",body_part_sensory) ~ "pars_intercerebralis",
    grepl("pars_lateralis",body_part_sensory) ~ "pars_lateralis",
    grepl("sez|subesophageal_zone",body_part_sensory) ~ "subesophageal_zone",
    grepl("ventral_nerve_cord|vnc",body_part_sensory) ~ "ventral_nerve_cord",
    # blood
    grepl("haemolymph|hemolymph",body_part_sensory) ~ "hemolymph",
    # thorax
    grepl("thorax",body_part_sensory) ~ "thorax",
    #grepl("neck",body_part_sensory) ~ "neck",
    # Legs
    grepl("front_leg",body_part_sensory) ~ "front_leg",
    grepl("middle_leg",body_part_sensory) ~ "middle_leg",
    grepl("hind_leg",body_part_sensory) ~ "hind_leg",
    grepl("leg",body_part_sensory) ~ "leg",
    # Flight
    grepl("haltere",body_part_sensory) ~ "haltere",
    grepl("notum|thorax",body_part_sensory) ~ "thorax",
    grepl("wing_base",body_part_sensory) ~ "wing_base",
    grepl("wing_margin",body_part_sensory) ~ "wing_margin",
    grepl("wing_tegula",body_part_sensory) ~ "wing_tegula",
    grepl("wing",body_part_sensory) ~ "wing",
    grepl("wing_endocrine",cell_function) ~ "wing",
    # Abdomen
    grepl("abdomen",body_part_sensory) ~ "abdomen",
    grepl("abdominal_wall",body_part_sensory) ~ "abdominal_wall",
    # Chordotonal organs
    grepl("metathoracic_chordotonal_organ",body_part_sensory) ~ "metathoracic_chordotonal_organ",
    grepl("prothoracic_chordotonal_organ",body_part_sensory) ~ "prothoracic_chordotonal_organ",
    grepl("prosternal_organ",body_part_sensory) ~ "prosternal_organ",
    grepl("wheelers_organ",body_part_sensory) ~ "wheelers_organ",
    # other
    !is.na(body_part_sensory) ~ body_part_sensory,
    ## by nerves
    grepl("^AN$",FAFB_nerve) ~ "antenna",
    grepl("^aPhN$",FAFB_nerve) ~ "pharynx",
    grepl("^PhN$",FAFB_nerve) ~ "pharynx",
    grepl("AbN2|AbN3|AbN4|AbNT",MANC_entryNerve) ~ "abdomen",
    grepl("DMetaN",MANC_entryNerve) ~ "haltere",
    grepl("DProN",MANC_entryNerve) ~ "front_leg",
    grepl("ProLN",MANC_entryNerve) ~ "front_leg",
    grepl("VProN",MANC_entryNerve) ~ "front_leg",
    grepl("MesoLN",MANC_entryNerve) ~ "middle_leg",
    grepl("MetaLN", MANC_entryNerve) ~ "hind_leg",
    grepl("PDMN",MANC_entryNerve) ~ "thorax",
    grepl("ADMN",MANC_entryNerve) ~ "wing",
    grepl("PrN",MANC_entryNerve) ~ "prosternal_organ",
    grepl("ProCN",MANC_entryNerve) ~ "prothoracic_chordotonal_organ",
    # remainder
    TRUE ~ "unknown"
  ) ) %>%
  dplyr::mutate(body_part_effector = dplyr::case_when(
    !grepl("efferent",flow) ~ NA,
    grepl("retrocerebral_complex",body_part_effector) ~ "retrocerebral_complex",
    grepl("enteric_complex",body_part_effector) ~ "digestive_tract",
    grepl("corpus_allatum",body_part_effector) ~ "corpus_allatum",
    grepl("salivary_gland",body_part_effector) ~ "salivary_gland",
    grepl("neurohemal_complex",body_part_effector) ~ "neurohemal_complex",
    grepl("antenna",body_part_effector) ~ "antenna",
    grepl("crop",body_part_effector) ~ "crop",
    grepl("eye",body_part_effector) ~ "eye",
    grepl("front_leg",body_part_effector) ~ "front_leg",
    grepl("middle_leg",body_part_effector) ~ "middle_leg",
    grepl("hind_leg",body_part_effector) ~ "hind_leg",
    grepl("haltere",body_part_effector) ~ "haltere",
    grepl("neck",body_part_effector) ~ "neck",
    grepl("proboscis",body_part_effector) ~ "proboscis",
    grepl("pharynx",body_part_effector) ~ "pharynx",
    grepl("wing",body_part_effector) ~ "wing",
    grepl("abdomen",body_part_effector) ~ "abdomen",
    grepl("prothorax",body_part_effector) ~ "prothorax",
    ## by nerves
    !is.na(body_part_effector) ~ body_part_effector,
    #    grepl("^AN$",FAFB_nerve) ~ "antenna",
    #    grepl("^aPhN$",FAFB_nerve) ~ "pharynx",
    #    grepl("^PhN$",FAFB_nerve) ~ "pharynx",
    #    grepl("^CV$|CvN|CVC",FAFB_nerve) ~ "neck",
    #    grepl("AbN2|AbN3|AbN4|AbNT",MANC_exitNerve) ~ "abdomen",
    #    grepl("DMetaN",MANC_exitNerve) ~ "haltere",
    #    grepl("DProN",MANC_exitNerve) ~ "prothorax",
    #    grepl("ProAN|ProLN",MANC_exitNerve) ~ "front_leg",
    #    grepl("VProN",MANC_exitNerve) ~ "front_leg",
    #    grepl("MesoLN",MANC_exitNerve) ~ "middle_leg",
    #    grepl("MetaLN",MANC_exitNerve) ~ "hind_leg",
    #    grepl("PDMN",MANC_exitNerve) ~ "wing",
    #    grepl("ADMN|MesoAN",MANC_exitNerve) ~ "wing",
    #    grepl("CvN|CVC|^CV$",MANC_exitNerve) ~ "neck",
    TRUE ~ "unknown"
  ) ) %>%
  dplyr::mutate(nerve = dplyr::case_when(
    
    grepl("NCC",nerve) & side=="right" ~ "right_corpus_cardiacum_nerve",
    grepl("NCC",nerve) &side=="left" ~ "left_corpus_cardiacum_nerve",
    grepl("NCC",nerve) ~ "left_corpus_cardiacum_nerve",
    
    grepl("^AN",nerve) & side=="right" ~ "right_antennal_nerve",
    grepl("^AN",nerve) &side=="left" ~ "left_antennal_nerve",
    grepl("^AN",nerve) ~ "antennal_nerve",
    
    grepl("MxLbN",nerve) & side=="right" ~ "right_maxillary-labial_nerve",
    grepl("MxLbN",nerve) &side=="left" ~ "left_maxillary-labial_nerve",
    grepl("MxLbN",nerve) ~ "maxillary-labial_nerve",
    
    grepl("PhN",nerve) & side=="right" ~ "right_pharyngeal_nerve",
    grepl("PhN",nerve) &side=="left" ~ "left_pharyngeal_nerve",
    grepl("PhN",nerve) ~ "pharyngeal_nerve",
    
    grepl("aPhN",nerve) & side=="right" ~ "right_accessory_pharyngeal_nerve",
    grepl("aPhN",nerve) &side=="left" ~ "left_accessory_pharyngeal_nerve",
    grepl("aPhN",nerve) ~ "accessory_pharyngeal_nerve",
    
    grepl("OCN",nerve) & side=="right" ~ "right_ocellar_nerve",
    grepl("OCN",nerve) &side=="left" ~ "left_ocellar_nerve",
    grepl("OCN",nerve) ~ "ocellar_nerve",
    
    grepl("^ON",nerve) & side=="right" ~ "right_optic_nerve",
    grepl("^ON",nerve) &side=="left" ~ "left_optic_nerve",
    grepl("^ON",nerve) ~ "optic_nerve",
    
    grepl("CvN_R",nerve)  ~ "right_cervical_nerve",
    grepl("CvN_L",nerve) ~ "left_cervical_nerve",
    grepl("CvN",nerve) ~ "cervical_nerve",
    grepl("CV",nerve) & side=="right" ~ "right_cervical_nerve",
    grepl("CV",nerve) & side=="left" ~ "left_cervical_nerve",
    grepl("CV",nerve) & side=="left" ~ "cervical_nerve",
    
    grepl("DProN_R",nerve)  ~ "right_dorsal_prothoracic_nerve",
    grepl("DProN_L",nerve)  ~ "left_dorsal_prothoracic_nerve",
    grepl("DProN",nerve) ~ "dorsal_prothoracic_nerve",
    
    grepl("ProLN_R",nerve)  ~ "right_prothoracic_leg_nerve",
    grepl("ProLN_L",nerve)  ~ "left_prothoracic_leg_nerve",
    grepl("ProLN",nerve) ~ "prothoracic_leg_nerve",
    
    grepl("PrN_R",nerve) ~ "right_prosternal_nerve",
    grepl("PrN_L",nerve)  ~ "left_prosternal_nerve",
    grepl("PrN",nerve) ~ "prosternal_nerve",
    
    grepl("ProAN_R",nerve)  ~ "right_prothoracic_accessory_nerve",
    grepl("ProAN_L",nerve) ~ "left_prothoracic_accessory_nerve",
    grepl("ProAN",nerve) ~ "prothoracic_accessory_nerve",
    
    grepl("ProCN_R",nerve)  ~ "right_prothoracic_chordotonal_nerve",
    grepl("ProCN_L",nerve)  ~ "left_prothoracic_chordotonal_nerve",
    grepl("ProCN",nerve) ~ "prothoracic_chordotonal_nerve",
    
    grepl("VProN_R",nerve)  ~ "right_ventral_prothoracic_nerve",
    grepl("VProN_L",nerve)  ~ "left_ventral_prothoracic_nerve",
    grepl("VProN",nerve) ~ "ventral_prothoracic_nerve",
    
    grepl("ADMN_R",nerve) ~ "right_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN_L",nerve)  ~ "left_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN",nerve) ~ "anterior_dorsal_mesothoracic_nerve",
    
    grepl("PDMN_R",nerve)  ~ "right_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN_L",nerve)  ~ "left_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN",nerve) ~ "posterior_dorsal_mesothoracic_nerve",
    
    grepl("MesoLN_R",nerve)  ~ "right_mesothoracic_leg_nerve",
    grepl("MesoLN_L",nerve)  ~ "left_mesothoracic_leg_nerve",
    grepl("MesoLN",nerve) ~ "mesothoracic_leg_nerve",
    
    grepl("MesoAN_R",nerve)  ~ "right_mesothoracic_accessory_nerve",
    grepl("MesoAN_L",nerve)  ~ "left_mesothoracic_accessory_nerve",
    grepl("MesoAN",nerve) ~ "mesothoracic_accessory_nerve",
    
    grepl("DMetaN_R",nerve)  ~ "right_dorsal_metathoracic_nerve",
    grepl("DMetaN_L",nerve)  ~ "left_dorsal_metathoracic_nerve",
    grepl("DMetaN",nerve) ~ "dorsal_metathoracic_nerve",
    
    grepl("MetaLN_R",nerve)  ~ "right_metathoracic_leg_nerve",
    grepl("MetaLN_L",nerve)  ~ "left_metathoracic_leg_nerve",
    grepl("MetaLN",nerve) ~ "metathoracic_leg_nerve",
    
    grepl("AbN1_R",nerve)  ~ "right_first_abdominal_nerve",
    grepl("AbN1_L",nerve)  ~ "left_first_abdominal_nerve",
    grepl("AbN1",nerve) ~ "first_abdominal_nerve",
    
    grepl("AbN2_R",nerve)  ~ "right_second_abdominal_nerve",
    grepl("AbN2_L",nerve)  ~ "left_second_abdominal_nerve",
    grepl("AbN2",nerve) ~ "second_abdominal_nerve",
    
    grepl("AbN3_R",nerve)  ~ "right_third_abdominal_nerve",
    grepl("AbN3_L",nerve)  ~ "left_third_abdominal_nerve",
    grepl("AbN3",nerve) ~ "third_abdominal_nerve",
    
    grepl("AbN4_R",nerve)  ~ "right_fourth_abdominal_nerve",
    grepl("AbN4_L",nerve)  ~ "left_fourth_abdominal_nerve",
    grepl("AbN4",nerve) ~ "fourth_abdominal_nerve",
    
    grepl("AbNT_R",nerve)  ~ "right_abdominal_nerve_trunk",
    grepl("AbNT_L",nerve)  ~ "left_abdominal_nerve_trunk",
    grepl("AbNT",nerve) ~ "abdominal_nerve_trunk",
    
    grepl("AbNX_R",nerve)  ~ "right_abdominal_nerve_other",
    grepl("AbNX_L",nerve)  ~ "left_abdominal_nerve_other",
    grepl("AbNX",nerve) ~ "abdominal_nerve_other",
    
    TRUE ~ NA
  ) ) %>%
  # translate neuropeptide names to FlyBase ones in neuropeptide_verified column
  dplyr::mutate(neuropeptide_verified = gsub("allatostatin-a","AstA",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("allatostatin-c","AstC",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("amnesiac","amn",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("capability","Capa",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("ccap","CCAP",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("ccha1","CCHa1",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("ccha2","CCHa2",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("cnma","CNMa",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("corazonin","Crz",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("dARC1","dARC1",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("dh31|Dh331","Dh31",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("dh44","Dh44",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("DILP","Ilp",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("dsk","Dsk",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("dnpf|^NPF","NPF",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("drosulfakinin","Dsk",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("eclosion hormone","Eh",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("fmrfa","FMRFa",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("gpa2","Gpa2",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("gpb5","Gpb5",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("hugin","Hug",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("itp","ITP",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("leucokinin","Lk",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("MIP","Mip",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("myosuppressin|myosupressin","Ms",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("natalisin","Natalisin",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("Nplp1","Nplp1",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("Nplp2","Nplp2",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("Nplp3","Nplp3",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("Nplp4","Nplp4",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("orcokinin","Orcokinin",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("burisconin","Pburs",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("pdf","Pdf",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("proctolin","Proc",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("ptth","Ptth",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("ryanimide","RYa",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("^sp$","SP",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("SIFamide","SIFa",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("sNPF","sNPF",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("space blanket","Sb",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("tachykinin","Tk",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = gsub("trissin","Trissin",neuropeptide_verified)) %>%
  dplyr::mutate(neuropeptide_verified = 
           ifelse(grepl("neuropeptide-negative|negative", neuropeptide_verified, ignore.case = TRUE), 
                  NA, neuropeptide_verified)) %>%
  # translate neuropeptide names to FlyBase ones in neurotransmitter_verified column
  dplyr::mutate(neurotransmitter_verified = gsub("allatostatin-a","AstA",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("allatostatin-c","AstC",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("amnesiac","amn",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("capability","Capa",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("ccap","CCAP",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("ccha1","CCHa1",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("ccha2","CCHa2",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("cnma","CNMa",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("corazonin","Crz",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("dARC1","dARC1",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("dh31|Dh331","Dh31",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("dh44","Dh44",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("DILP","Ilp",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("dsk","Dsk",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("dnpf|^NPF","NPF",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("drosulfakinin","Dsk",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("eclosion hormone","Eh",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("fmrfa","FMRFa",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("gpa2","Gpa2",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("gpb5","Gpb5",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("hugin","Hug",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("itp","ITP",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("leucokinin","Lk",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("MIP","Mip",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("myosuppressin|myosupressin","Ms",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("natalisin","Natalisin",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("Nplp1","Nplp1",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("Nplp2","Nplp2",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("Nplp3","Nplp3",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("Nplp4","Nplp4",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("orcokinin","Orcokinin",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("burisconin","Pburs",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("pdf","Pdf",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("proctolin","Proc",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("ptth","Ptth",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("ryanimide","RYa",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("^sp$","SP",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("SIFamide","SIFa",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("sNPF","sNPF",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("space blanket","Sb",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("tachykinin","Tk",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("trissin","Trissin",neurotransmitter_verified)) %>%
  dplyr::mutate(neurotransmitter_verified = 
                  ifelse(grepl("neuropeptide-negative|negative", neurotransmitter_verified, ignore.case = TRUE), 
                         NA, neurotransmitter_verified)) %>%
  # change nitric oxide to include _
  dplyr::mutate(neurotransmitter_verified = gsub("nitric oxide","nitric_oxide",neurotransmitter_verified)) %>%
  # save existing neuropeptide_verified, neurotransmitter_verified, 
  #  make additional updates to new columns banc_neuropeptide_verified and
  #  banc_neurotransmitter_verified
  dplyr::mutate(banc_neuropeptide_verified = neuropeptide_verified) %>%
  # join with semicolons, remove duplicates
  dplyr::mutate(banc_neuropeptide_verified = 
           sapply(banc_neuropeptide_verified, function(entry) {
             # Check if entry is empty or NA
             if (is.na(entry) || entry == "") {
               return(NA)  # Return NA if the entry is empty or NA
             }
             
             # Otherwise, proceed to split the entry
             split_entries <- unlist(str_split(entry, "[ ,;]+")) %>%
               unique()  # Remove duplicates
             
             # Join the unique entries with "; "
             paste(split_entries, collapse = "; ")
           })
  ) %>%
  dplyr::mutate(banc_neurotransmitter_verified = neurotransmitter_verified) %>%
  # clean up banc_neurotransmitter_verified column to remove neuropeptides and
  #  notes about being negative for specific neurotransmitter
  dplyr::mutate(banc_neurotransmitter_verified = map_chr(
    str_split(banc_neurotransmitter_verified, "[,;]\\s*"), 
    function(values) {
      # Filter to keep only the allowed neurotransmitters
      valid_values <- values[values %in% allowed_neurotransmitters]
      
      # Remove any empty values that might remain
      valid_values <- valid_values[valid_values != ""]
      
      # Remove duplicates
      unique_values <- unique(valid_values)
      
      # Join with semicolons
      if(length(unique_values) > 0) {
        paste(unique_values, collapse = "; ")
      } else {
        NA  # Return NA if no valid neurotransmitters
      }
    }
  )) %>%
  #dplyr::rowwise() %>%
  #dplyr::mutate(neurotransmitter_verified = strsplit(neurotransmitter_verified,split=",")) %>%
  #dplyr::ungroup() %>%
  dplyr::mutate(super_class = dplyr::case_when(
    # not neurons
    cell_class %in% c('glia','putative_glia') ~ "glia",
    grepl("trachea",super_class) ~ "trachea",
    grepl("not_a_neuron|NOT_A_NEURON",super_class) ~ "not_a_neuron",
    # Neck
    cell_type == "SNpp54" ~"ventral_nerve_cord_sensory",
    grepl("^ISN$",cell_type) ~ "central_brain_sensory",
    super_class %in% c("descending") ~ "descending",
    super_class %in% c("ascending") ~ "ascending",
    super_class %in% c("sensory_ascending") ~ "sensory_ascending",
    grepl("ascending|descending", super_class) & grepl('sensory_descending', super_class) ~ "sensory_descending",
    grepl("ascending|descending", super_class) & grepl('sensory', super_class) ~ "sensory_ascending",
    grepl("ascending|descending", super_class) & grepl('efferent|motor|endocrine', super_class) ~ "ascending_visceral_circulatory",
    grepl("ascending|descending", super_class) & grepl('afferent', flow) ~ "sensory_ascending",
    super_class %in% c("sensory_descending") ~ "sensory_descending",
    super_class %in% c("visceral_circulatory_ascending") ~ "ascending_motor",
    super_class %in% c("efferent_ascending") ~ "efferent_ascending",
    super_class %in% c("efferent_descending") ~ "efferent_descending",
    # Afferent
    grepl("afferent",flow) ~ paste0(region,"_sensory"),
    grepl("sensory",flow) ~ paste0(region,"_sensory"),
    cell_type %in% "HBeyelet" ~ "optic_lobe_sensory",
    grepl("ocellar_retinula", cell_type) ~ "central_brain_sensory",
    grepl("^R1-6$|^R7$|^R8$", cell_type) ~ "optic_lobe_sensory",
    # efferent
    grepl("endocrine|visceral_circulatory",cell_class) ~ paste0(region,"_visceral_circulatory"),
    grepl("endocrine|visceral_circulatory",super_class) ~ paste0(region,"_visceral_circulatory"),
    grepl("endocrine|visceral_circulatory",cell_function) ~ paste0(region,"_visceral_circulatory"),
    grepl("motor",cell_class) ~ paste0(region,"_motor"),
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ paste0(region,"_motor"),
    grepl("motor",cell_function) ~ paste0(region,"_motor"),
    # optic transfer
    super_class %in% c("visual_centrifugal") ~ "visual_centrifugal",
    super_class %in% c("visual_projection") ~ "visual_projection",
    # intrinsic
    grepl("ocellar",cell_class) ~ "central_brain_intrinsic",
    grepl("midbrain|central_brain",region) ~ "central_brain_intrinsic",
    grepl("vnc|ventral_nerve_cord",region) ~ "ventral_nerve_cord_intrinsic",
    grepl("vnc|ventral_nerve_cord",super_class) ~ "ventral_nerve_cord_intrinsic",
    grepl("vnc|ventral_nerve_cord",cell_class) ~ "ventral_nerve_cord_intrinsic",
    grepl("optic",region) ~ "optic_lobe_intrinsic",
    super_class %in% c("optic") ~ "optic_lobe_intrinsic",
    grepl("central",flow) ~ "central_brain_intrinsic",
    grepl("central",super_class) ~ "central_brain_intrinsic",
    TRUE ~ "unknown"
  )) %>%
  # flow again
  dplyr::mutate(flow = dplyr::case_when(
    super_class=="ascending_sensory"&grepl("^SA",cell_type) ~ "afferent",
    super_class=="ascending$"&!grepl("^SA",cell_type) ~ "intrinsic",
    cell_type %in% "HBeyelet" ~ "afferent",
    grepl("ocellar_retinula", cell_type) ~ "afferent",
    # none get NA, as desired
    TRUE ~ flow 
  ) ) %>%
  # peripheral_target_type
  dplyr::mutate(peripheral_target_type = dplyr::case_when(
    # pain
    grepl("SNch01|SNxx19|SNxx20|SNxx21",MANC_type) ~ "multidendritic", # class_iv
    grepl("nociception",cell_function) ~ "multidendritic", 
    grepl("pharynx|cibarium",body_part_sensory)&(grepl("mechanosensory|tactile|pumping|ciberial",cell_function)|grepl("mechanosensory|tactile|pumping|ciberial",cell_function_detailed)) ~ "multidendritic",
    # taste
    grepl("taste peg",cell_sub_class) ~ "taste_peg",
    grepl("taste_peg",cell_sub_class) ~ "taste_peg",
    grepl("taste bristle",cell_sub_class) ~ "taste_peg",
    grepl("taste_bristle",cell_sub_class) ~ "taste_peg",
    grepl("taste_peg",cell_function_detailed) ~ "taste_peg",
    grepl("hair plate",cell_sub_class) ~ "hair_plate",
    grepl("eye_bristle",cell_sub_class) ~ "bristle",
    grepl("head_bristle",cell_sub_class) ~ "bristle",
    grepl("bristle",cell_function_detailed) ~ "bristle",
    grepl("bristle",cell_function_detailed) ~ "bristle",
    grepl("mechanosensory_bristle",cell_sub_class) ~ "bristle",
    grepl("bristle",body_part_sensory) ~ "bristle",
    grepl("^JO-",cell_type) ~ "chordotonal_organ",
    grepl("chordotonal",cell_type) ~ "chordotonal_organ",
    grepl("chordotonal",cell_sub_class) ~ "chordotonal_organ",
    grepl("chordotonal",cell_function) ~ "chordotonal_organ",
    grepl("DNx01",cell_type) ~ "campaniform_sensillum",
    grepl("campaniform",cell_sub_class) ~ "campaniform_sensillum",
    grepl("campaniform",cell_function_detailed) ~ "campaniform_sensillum",
    grepl("campaniform",cell_function) ~ "campaniform_sensillum",
    # hemolymph
    grepl("haemolymph|hemolymph",body_part_sensory) ~ "nutrient_receptor",
    # both
    grepl("taste peg",cell_function_detailed) ~ "taste_peg",
    grepl("taste_peg",cell_function_detailed) ~ "taste_peg",
    grepl("taste_peg",cell_sub_class) ~ "taste_peg",
    grepl("bristle",cell_function_detailed) ~ "bristle",
    grepl("bristle",cell_sub_class) ~ "bristle",
    grepl("hair_plate",cell_function_detailed) ~ "hair_plate",
    grepl("hair_plate",cell_sub_class) ~ "hair_plate",
    grepl("taste peg",cell_function) ~ "taste_peg",
    grepl("taste_peg",cell_function) ~ "taste_peg",
    grepl("bristle",cell_function) ~ "bristle",
    grepl("bristle",cell_sub_class) ~ "bristle",
    grepl("hair_plate",cell_function) ~ "hair_plate",
    # from MANC
    grepl("campaniform sensilla",MANC_receptorType) ~ "campaniform_sensillum",
    grepl("chordotonal",MANC_receptorType) ~ "chordotonal_organ",
    grepl("chordotonal",cell_function_detailed) ~ "chordotonal_organ",
    grepl("hair plate",MANC_receptorType) ~ "hair_plate",
    grepl("putative sweet taste bristle",MANC_receptorType) ~ "taste_peg",
    grepl("sweet taste bristle",MANC_receptorType) ~ "taste_peg",
    grepl("taste bristle",MANC_receptorType) ~ "taste_peg",
    grepl("strand receptor",MANC_receptorType) ~ "strand",
    grepl("taste bristle",MANC_receptorType) ~ "taste_peg",
    # visual
    grepl("ocellar retinula",FAFB_cell_type) ~ "photoreceptor",
    grepl("ocellar_retinula",cell_type) ~ "photoreceptor",
    grepl("^R1-6|^R1-6|^R7|^R8",cell_type) ~ "photoreceptor",
    grepl("^R1-6|^R1-6|^R7|^R8",cell_type) ~ "photoreceptor",
    grepl("^HBeyelet",cell_type) ~ "photoreceptor",
    grepl("^HBeyelet",cell_type) ~ "photoreceptor",
    # from FAFB
    grepl("^ORN",cell_type)&cell_function=="olfactory" ~ "olfactory_sensillum",
    grepl("^TRN",cell_type)&cell_function=="thermosensory" ~ "thermosensory_sensillum",
    grepl("^HRN",cell_type)&cell_function=="hygrosensory" ~ "hygrosensory_sensillum",
    grepl("AC neuron|AC_neuron|ac_neuron",cell_type) ~ "internal_thermosensory_receptor",
    grepl("leg",body_part_sensory)&grepl("chemosensory|gustatory",cell_function) ~ "taste_peg",
    grepl("labellum",body_part_sensory)&grepl("chemosensory|gustatory",cell_function) ~ "external_taste_sensillum",
    grepl("abdomen",body_part_sensory)&grepl("chemosensory|gustatory",cell_function) ~ "external_taste_sensillum",
    grepl("pharynx|crop",body_part_sensory)&grepl("chemosensory|gustatory",cell_function) ~ "internal_taste_sensillum",
    grepl("wing_margin",body_part_sensory)&grepl("chemosensory|gustatory",cell_function) ~ "taste_sensillum",
    grepl("sensory",super_class) ~ "orphan",
    # motor
    # leg muscles - from MANC (note, subtle differences from FANC)
    grepl("Acc._ti_flexor", cell_class) ~ "accessory_tibia_flexor_muscle",
    grepl("Fe_reductor", cell_class) ~ "femur_reductor_muscle",
    grepl("ltm1-tibia", cell_class) ~ "long_tendon_muscle_1",
    grepl("ltm2-femur", cell_class) ~ "long_tendon_muscle_2",
    grepl("ltm", cell_class) ~ "long_tendon_muscle",
    grepl("Pleural_remotor/abductor", cell_class) ~ "pleural_remotor_and_abductor_muscle",
    grepl("Sternal_adductor", cell_class) ~ "sternal_adductor_muscle",
    grepl("Sternal_anterior_rotator", cell_class) ~ "sternal_anterior_rotator_muscle",
    grepl("Sternal_posterior_rotator", cell_class) ~ "sternal_posterior_rotator_muscle",
    grepl("Sternotrochanter", cell_class) ~ "sternotrochanter_extensor_muscle",
    grepl("Ta_depressor", cell_class) ~ "tarsus_depressor_muscle",
    grepl("Ta_levator", cell_class) ~ "tarsus_levator_muscle",
    grepl("Tergopleural/Pleural_promotor", cell_class) ~ "tergopleural_or_pleural_promotor_muscle",
    grepl("Tergotr.", cell_class) ~ "tergotrochanter_extensor_muscle",
    grepl("Ti_extensor", cell_class) ~ "tibia_extensor_muscle",
    grepl("Ti_flexor", cell_class) ~ "tibia_flexor_muscle",
    grepl("Tr_extensor", cell_class) ~ "trochanter_extensor_muscle",
    grepl("Tr_flexor", cell_class) ~ "trochanter_flexor_muscle",
    grepl("leg_motor", cell_function) ~ "leg_muscle",
    # wing muscles - from MANC (note, subtle differences from FANC)
    grepl("b1_motor_neuron", cell_class) ~ "b1_muscle",
    grepl("b2_motor_neuron", cell_class) ~ "b2_muscle",
    grepl("b3_motor_neuron", cell_class) ~ "b3_muscle",
    grepl("^DLM", cell_class) ~ "dorsal_longitudinal_muscle",
    grepl("^DVM", cell_class) ~ "dorsoventral_muscle",
    grepl("^hg1_motor_neuron", cell_class) ~ "iv1_muscle",
    grepl("^hg2_motor_neuron", cell_class) ~ "iv2_muscle",
    grepl("^hg3_motor_neuron", cell_class) ~ "iv3_muscle",
    grepl("^hg4_motor_neuron", cell_class) ~ "iv4_muscle",
    grepl("^i1_motor_neuron", cell_class) ~ "i1_muscle",
    grepl("^i2_motor_neuron", cell_class) ~ "i2_muscle",
    grepl("^iii1_motor_neuron", cell_class) ~ "iii1_muscle",
    grepl("^iii3_motor_neuron", cell_class) ~ "iii2_muscle",
    grepl("^ps1_motor_neuron", cell_class) ~ "ps1_muscle",
    grepl("^ps2_motor_neuron", cell_class) ~ "ps2_muscle",
    grepl("^tp1_motor_neuron", cell_class) ~ "tp1_muscle",
    grepl("^tp2_motor_neuron", cell_class) ~ "tp2_muscle",
    grepl("^tp_motor_neuron", cell_class) ~ "tergopleural_muscle",
    grepl("^TTM_motor_neuron", cell_class) ~ "trochanter_extensor_muscle",
    ((grepl("wing", body_part_effector) & grepl("motor", super_class)) | 
       (grepl("wing_motor_neuron", cell_class))) ~ "wing_muscle",
    # haltere muscles - from MANC
    grepl("^hDVM_motor_neuron", cell_class) ~ "haltere_dorsoventral_muscle",
    grepl("^hi1_motor_neuron", cell_class) ~ "hi1_muscle",
    grepl("^hi2_motor_neuron", cell_class) ~ "hi2_muscle",
    grepl("^hiii2_motor_neuron", cell_class) ~ "hiii2_muscle",
    grepl("haltere_motor_neuron", cell_class) ~ "haltere_direct_control_muscle",
    # neck muscles - from MANC
    grepl("^CvN_A1$", cell_type) ~ "transverse_horizontal_1_muscle",
    grepl("^CvN_A2$", cell_type) ~ "transverse_horizontal_2_muscle",
    (grepl("^FNM2$", cell_type) & grepl("AD_motor_neuron", cell_class)) ~ "adductor_muscle",
    grepl("neck_motor_neuron", cell_class) ~ "neck_muscle",
    # antenna muscles - from FAFB (specific muscle targets unknown)
    grepl("antennal_motor_neuron", cell_class) ~ "antenna_muscle",
    # retina muscles - from FAFB (specific muscle targets unknown)
    grepl("eye_motor_neuron", cell_class) ~ "retina_muscle",
    # proboscis muscles (specific muscle targets unknown)
    (grepl("proboscis_motor_neuron", cell_class) | grepl("proboscis", body_part_effector)) ~ "proboscis_muscle",
    # pharynx muscles (specific muscle targets unknown)
    grepl("pharynx", body_part_effector) ~ "pharynx_muscle",
    # proboscis muscles (specific muscle targets unknown)
    grepl("proboscis_motor_neuron", cell_class) ~ "proboscis_muscle",
    # salivary muscle
    grepl("salivary_motor_neuron", cell_class) ~ "m13_muscle",
    # crop
    grepl("crop_motor_neuron", cell_class) ~ "crop_muscle",
    # remainder
    TRUE ~ NA
  ) ) %>%
  dplyr::mutate(cell_function_detailed = dplyr::case_when(
    # gustatory
    grepl("sugar_or_water",cell_function)|grepl("sugar_or_water",cell_function_detailed) ~ "sweet_or_water",
    grepl("putative sweet taste bristle",MANC_receptorType) ~ "sweet",
    grepl("sweet taste bristle",MANC_receptorType) ~ "sweet",
    grepl("sweet taste bristle",MANC_receptorType) ~ "sweet",
    grepl("sugar",cell_function)|grepl("sugar",cell_function_detailed) ~ "sweet",
    grepl("sweet",cell_function)|grepl("sweet",cell_function_detailed) ~ "sweet",
    grepl("water",cell_function)|grepl("water",cell_function_detailed)|grepl("water",cell_sub_class) ~ "water",
    grepl("bitter",cell_function)|grepl("bitter",cell_function_detailed)|grepl("bitter",cell_sub_class) ~ "bitter",
    grepl("low_salt",cell_function)|grepl("low_salt",cell_function_detailed)|grepl("low_salt",cell_sub_class) ~ "low_salt",
    grepl("pheromone_contact",cell_function)|grepl("pheromone_contact",cell_function_detailed) ~ "pheromone_contact",
    grepl("gustatory",cell_function)&grepl("pheromone",cell_function_detailed) ~ "pheromone_contact",
    grepl("chemosensory",cell_function)&grepl("pheromone",cell_function_detailed) ~ "pheromone_contact",
    # olfactory --add
    grepl("^ORN_VL2p$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_DL5$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_DM2$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_VM7d$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_VM3$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_VA2$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_DP1m$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_DP1l$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_VC4$",cell_type)&cell_function=="olfactory" ~ "alcoholic_fermentation_volatile",
    grepl("^ORN_V$",cell_type)&cell_function=="olfactory" ~ "carbon_dioxide_volatile",
    grepl("^ORN_DM1$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_DM4$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VM7d$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VM5v$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_DM3$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VA3$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VA4$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VC2$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_DM1$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VA3$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_VA6$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_DM5$",cell_type)&cell_function=="olfactory" ~ "yeasty_volatile",
    grepl("^ORN_DM3$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_VA3$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_VM5v$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_DM2$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_VM2$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_VM5d$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_DL1$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_DA3$",cell_type)&cell_function=="olfactory" ~ "fruity_volatile",
    grepl("^ORN_DL2d$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_DL2v$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VC5$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VL2a$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VM4$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VL1$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VM1$",cell_type)&cell_function=="olfactory" ~ "decaying_fruit_volatile",
    grepl("^ORN_VA1v$",cell_type)&cell_function=="olfactory" ~ "pheromone_volatile",
    grepl("^ORN_VA1d$",cell_type)&cell_function=="olfactory" ~ "pheromone_volatile",
    grepl("^ORN_DA1$",cell_type)&cell_function=="olfactory" ~ "pheromone_volatile",
    grepl("^ORN_DL3$",cell_type)&cell_function=="olfactory" ~ "pheromone_volatile",
    grepl("^ORN_DC3$",cell_type)&cell_function=="olfactory" ~ "pheromone_volatile",
    grepl("^ORN_DA2$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DM3$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DC1$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DL5$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DL4$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DC4$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_VC3$",cell_type)&cell_function=="olfactory" ~ "aversive_volatile",
    grepl("^ORN_DC2$",cell_type)&cell_function=="olfactory" ~ "plant_matter_volatile",
    grepl("^ORN_VC1$",cell_type)&cell_function=="olfactory" ~ "plant_matter_volatile",
    grepl("^ORN_VA5$",cell_type)&cell_function=="olfactory" ~ "animal_matter_volatile",
    grepl("^ORN_VA7l$",cell_type)&cell_function=="olfactory" ~ "animal_matter_volatile",
    grepl("co2|carbon_dioxide",cell_function)|grepl("co2|carbon_dioxide",cell_function_detailed) ~ "carbon_dioxide",
    grepl("pheromone_volatile",cell_function)|grepl("pheromone_volatile",cell_function_detailed) ~ "pheromone_volatile",
    grepl("olfactory",cell_function)&grepl("pheromone",cell_function_detailed) ~ "pheromone_volatile",
    # thermosensory
    grepl("TRN_VP1",cell_type) ~ "evaporation",
    grepl("TRN_VP2",cell_type) ~ "heating",
    grepl("TRN_VP3",cell_type) ~ "cooling",
    grepl("HRN_VP4",cell_type) ~ "dry",
    grepl("HRN_VP5",cell_type) ~ "humid",
    grepl("heating",cell_sub_class) ~ "heating",
    grepl("cooling",cell_sub_class) ~ "cooling",
    grepl("humid",cell_sub_class) ~ "humid",
    grepl("dry",cell_sub_class) ~ "dry",
    # mechanosensory
    ### chordotonal
    grepl("JO-A",cell_type) ~ "auditory_low_frequency",
    grepl("JO-B",cell_type) ~ "auditory_high_frequency",
    grepl("JO-C",cell_type) ~ "direction",
    grepl("JO-D",cell_type) ~ "position",
    grepl("JO-E",cell_type) ~ "position",
    grepl("JO-F",cell_type) ~ "position",
    grepl("chordotonal_club|club_chordotonal",cell_function_detailed)|grepl("club",cell_sub_class) ~ "vibro_tactile",
    grepl("chordotonal_claw|claw_chordotonal",cell_function_detailed)|grepl("claw",cell_sub_class)  ~ "position",
    grepl("chordotonal_hook|hook_chordotonal",cell_function_detailed)|grepl("hook",cell_sub_class)  ~ "direction",
    grepl("club",cell_sub_class) ~ "vibro_tactile",
    grepl("claw",cell_sub_class) ~ "position",
    grepl("hook",cell_sub_class) ~ "direction",
    grepl("vibro-tactile|vibro_tactile",cell_function)|grepl("vibro-tactile|vibro_tactile",cell_function_detailed) ~ "vibro_tactile",
    grepl("vibration",cell_function)|grepl("vibration",cell_function_detailed) ~ "vibration",
    grepl("position",cell_function)|grepl("position",cell_function_detailed) ~ "position",
    grepl("stretch",cell_function)|grepl("stretch",cell_function_detailed) ~ "stretch",
    grepl("deflection",cell_function)|grepl("deflection",cell_function_detailed) ~ "deflection",
    ### other
    grepl("_bristle",cell_function_detailed) ~ cell_function_detailed,
    grepl("campaniform",peripheral_target_type) ~ "mechanical_strain",
    grepl("hair_plate",peripheral_target_type) ~ "joint_angle",
    grepl("chordotonal",peripheral_target_type) ~ "vibro_position",
    ### putative pump
    grepl("ciberial",cell_function)|grepl("ciberial",cell_function_detailed) ~ "contractile",
    grepl("pumping",cell_function)|grepl("pumping",cell_function_detailed) ~ "contractile",
    grepl("ciberial",cell_function)|grepl("ciberial",cell_function_detailed) ~ "contractile",
    grepl("CB0991",cell_type) ~ "contractile",
    # respiratory
    grepl("oxygenation",cell_function)|grepl("oxygenation",cell_function_detailed) ~ "oxygenation",
    grepl("respiratory",cell_function)|grepl("respiratory",cell_function_detailed) ~ "respiratory",
    # interoceptive
    grepl("interoceptive",cell_function)|grepl("interoceptive",cell_function_detailed)|grepl("^ISN",cell_type) ~ "interoceptive",
    grepl("AC neuron|AC_neuron|ac_neuron",cell_type) ~ "internal_thermosensory",
    # motor related
    # leg 
    grepl("Acc._ti_flexor", cell_class) ~ "flex_femur_tibia_joint",
    grepl("Fe_reductor", cell_class) ~ "unknown_leg_movement",
    grepl("ltm1-tibia", cell_class) ~ "pull_long_tendon",
    grepl("ltm2-femur", cell_class) ~ "pull_long_tendon",
    grepl("ltm", cell_class) ~ "pull_long_tendon",
    grepl("Pleural_remotor/abductor", cell_class) ~ "move_coxa_posterior_lateral",
    grepl("Sternal_adductor", cell_class) ~ "move_coxa_medial",
    grepl("Sternal_anterior_rotator", cell_class) ~ "move_coxa_anterior",
    grepl("Sternal_posterior_rotator", cell_class) ~ "move_coxa_posterior",
    grepl("Sternotrochanter", cell_class) ~ "extend_coxa_trochanter_joint",
    grepl("Ta_depressor", cell_class) ~ "extend_tibia_tarsus_joint",
    grepl("Ta_levator", cell_class) ~ "flex_tibia_tarsus_joint",
    grepl("Tergopleural/Pleural_promotor", cell_class) ~ "move_coxa_anterior",
    grepl("Tergotr.", cell_class) ~ "extend_coxa_trochanter_joint",
    grepl("Ti_extensor", cell_class) ~ "extend_femur_tibia_joint",
    grepl("Ti_flexor", cell_class) ~ "flex_femur_tibia_joint",
    grepl("Tr_extensor", cell_class) ~ "extend_coxa_trochanter_joint",
    grepl("Tr_flexor", cell_class) ~ "flex_coxa_trochanter_joint",
    grepl("leg_motor", cell_function) ~ "unknown_leg_movement",
    # wing - from MANC (note, subtle differences from FANC)
    grepl("^b1_motor_neuron", cell_class) ~ "tonic_wing_steering",
    grepl("^b2_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("^b3_motor_neuron", cell_class) ~ "tonic_wing_steering",
    grepl("^DLM", cell_class) ~ "wing_power",
    grepl("^DVM", cell_class) ~ "wing_power",
    grepl("hg1_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("hg2_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("hg3_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("hg4_motor_neuron", cell_class) ~ "tonic_wing_steering",
    grepl("^i1_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("^i2_motor_neuron", cell_class) ~ "tonic_wing_steering",
    grepl("^iii1_motor_neuron", cell_class) ~ "phasic_wing_steering",
    grepl("^iii3_motor_neuron", cell_class) ~ "unknown_wing_steering",
    grepl("^ps1_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^ps2_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp1_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp2_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^TTM_motor_neuron", cell_class) ~ "middle_leg_extension",
    # haltere muscles - from MANC
    grepl("hDVM_motor_neuron", cell_class) ~ "haltere_power",
    grepl("hi1_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("hi2_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("hiii2_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("haltere_motor_neuron", cell_class) ~ "unknown_haltere",
    # visceral_circulatory - neuropeptide(s) for identified neurons (by Meet)
    grepl("^m_NSC_DILP$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^m_NSC_DH44$", cell_type) ~ banc_neuropeptide_verified,
    (grepl("^ventral_nerve_cord_visceral_circulatory$", super_class) &
       grepl("Dh44", banc_neuropeptide_verified) &
       grepl("Lk", banc_neuropeptide_verified)) ~ banc_neuropeptide_verified,
    grepl("^m_NSC_DMS$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^l_NSC_ITP$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^l_NSC_CRZ$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^l_NSC_DH31$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^SEZ_NSC_Hugin$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^SEZ_NSC_CAPA$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^ENXXX012$", cell_type) ~ banc_neuropeptide_verified,
    grepl("^EN27X010$", cell_type) ~ banc_neuropeptide_verified,
    # for remainder, wipe cell_function_detailed from visceral_circulatory neurons
    grepl("visceral_circulatory", super_class) ~ NA, 
    # remainder
    cell_function_detailed=="NA" ~ NA,
    is.na(cell_function_detailed) ~ cell_function_detailed,
    TRUE ~ cell_function_detailed
  ) ) %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    grepl("CB0991",cell_type) ~ "aorta",
    # gustatory
    grepl("taste_peg",peripheral_target_type) ~ "gustatory_tactile",
    grepl("sensory",super_class) & cell_type=='TPMN' ~ "gustatory_tactile",
    grepl("bitter",cell_function)|grepl("bitter",cell_function_detailed) ~ "gustatory",
    grepl("sugar_or_water",cell_function)|grepl("sugar_or_water",cell_function_detailed) ~ "gustatory",
    grepl("sugar",cell_function)|grepl("sugar",cell_function_detailed) ~ "gustatory",
    grepl("water",cell_function)|grepl("water",cell_function_detailed) ~ "gustatory",
    grepl("pheromone_contact",cell_function)|grepl("pheromone_contact",cell_function_detailed) ~ "gustatory",
    grepl("pheromone_contact",cell_function)|grepl("low_salt",cell_function_detailed) ~ "gustatory",
    grepl("chemosensory|gustatory",cell_function) ~ "gustatory",
    grepl("chemosensory|gustatory",cell_function_detailed) ~ "gustatory",
    grepl("sensory",super_class) & grepl("^SA_VTV_",cell_type) ~ "gustatory",
    # olfactory
    grepl("co2|carbon_dioxide",cell_function)|grepl("co2|carbon_dioxide",cell_function_detailed) ~ "olfactory",
    grepl("pheromone_volatile",cell_function)|grepl("pheromone_volatile",cell_function_detailed) ~ "olfactory",
    grepl("olfactory",cell_function) ~ "olfactory",
    # thermosensory
    grepl("thermosensory",cell_function) ~ "thermosensory",
    grepl("^TRN",cell_type) ~ "thermosensory",
    grepl("hygrosensory",cell_function) ~ "hygrosensory",
    grepl("^HRN",cell_type) ~ "hygrosensory",    
    # mechanosensory
    grepl("JO-",cell_type) ~ "proprioception",
    grepl("vibro-tactile|vibro_tactile",cell_function)|grepl("vibro-tactile|vibro_tactile",cell_function_detailed) ~ "proprioception",
    grepl("proprioception|proprioceptive",cell_function)|grepl("proprioception|proprioceptive",cell_function_detailed) ~ "proprioception",
    grepl("vibration",cell_function)|grepl("vibration",cell_function_detailed) ~ "proprioception",
    grepl("position",cell_function)|grepl("position",cell_function_detailed) ~ "proprioception",
    grepl("stretch",cell_function)|grepl("stretch",cell_function_detailed) ~ "proprioception",
    grepl("campaniform|hair_plate|chordotonal",peripheral_target_type) ~ "proprioception",
    ### bristle
    grepl("^tactile$",cell_function)|grepl("^tactile$",cell_function_detailed) ~ "tactile",
    grepl("mechanosensory|bristle",cell_function) ~ "tactile",
    grepl("bristle",peripheral_target_type) ~ "tactile",
    ### nociception
    grepl("nociception",cell_function)|grepl("nociception",cell_function_detailed) ~ "putative_nociception",
    ### contractile
    grepl("ciberial",cell_function)|grepl("ciberial",cell_function_detailed) ~ "contractile",
    grepl("pumping",cell_function)|grepl("pumping",cell_function_detailed) ~ "contractile",
    # visual
    cell_type %in% "R1-6" ~ "visual_achromatic",
    cell_type %in% c("L1","L2","L3","L4","L5") ~ "visual_achromatic",
    cell_type %in% "HBeyelet" ~ "visual_achromatic",
    cell_type %in% c("R7","R8") ~ "visual_chromatic",
    (grepl("^DmDRA", cell_type)) ~ "visual_polarized_light",
    grepl("ocellar",cell_function)|grepl("ocellar",cell_function_detailed) ~ "visual_ocellar",
    grepl("ocellar_retinula_cell",cell_type) ~ "visual_ocellar",
    grepl("^OCG",cell_type) ~ "visual_ocellar",
    grepl("visual_chromatic",cell_function)|grepl("visual_chromatic",cell_function_detailed) ~ "visual_chromatic",
    grepl("visual_achromatic",cell_function)|grepl("visual_achromatic",cell_function_detailed) ~ "visual_achromatic",
    grepl("ocellar",cell_function)|grepl("ocellar",cell_function_detailed) ~ "visual_ocellar",
    # respiratory
    grepl("oxygenation",cell_function)|grepl("oxygenation",cell_function_detailed) ~ "oxygenation",
    grepl("respiratory",cell_function)|grepl("respiratory",cell_function_detailed) ~ "respiratory",
    # interoceptive
    grepl("interoceptive",cell_function)|grepl("interoceptive",cell_function_detailed)|grepl("^ISN",cell_type) ~ "interoceptive",
    grepl("AC neuron",cell_type) ~ "thermosensory",
    # motor
    # leg
    grepl("leg_motor", cell_function) ~ "leg_motor",
    # wing muscles
    grepl("b1_motor_neuron", cell_class) ~ "wing_steering",
    grepl("b2_motor_neuron", cell_class) ~ "wing_steering",
    grepl("b3_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^DLM", cell_class) ~ "wing_power",
    grepl("^DVM", cell_class) ~ "wing_power",
    grepl("^hg1_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^hg2_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^hg3_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^hg4_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^i1_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^i2_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^iii1_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^iii3_motor_neuron", cell_class) ~ "wing_steering",
    grepl("^ps1_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^ps2_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp1_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp2_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^tp_motor_neuron", cell_class) ~ "wing_tension",
    grepl("^TTM_motor_neuron", cell_class) ~ "jump_escape",
    ((grepl("wing", body_part_effector) & grepl("motor", super_class)) | 
      (grepl("wing_motor_neuron", cell_class))) ~ "unknown_wing_motor",
    # haltere
    grepl("^hDVM_motor_neuron", cell_class) ~ "haltere_power",
    grepl("^hi1_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("^hi2_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("^hiii2_motor_neuron", cell_class) ~ "haltere_steering",
    grepl("haltere_motor_neuron", cell_class) ~ "unknown_haltere_motor",
    # neck - from MANC
    grepl("neck_motor_neuron", cell_class) ~ "neck_motor",
    grepl("neck_motor_neuron", cell_sub_class) ~ "neck_motor",
    # antenna 
    grepl("antennal_motor_neuron", cell_class) ~ "antenna_motor",
    # retina
    grepl("eye_motor_neuron", cell_class) ~ "retina_motor",
    # proboscis
    (grepl("proboscis_motor_neuron", cell_class) | grepl("proboscis", body_part_effector)) ~ "proboscis_motor",
    # pharynx
    grepl("pharynx", body_part_effector) ~ "pharynx_motor",
    # proboscis muscles (specific muscle targets unknown)
    grepl("proboscis_motor_neuron", cell_class) ~ "proboscis_motor",
    # salivary 
    grepl("salivary_motor_neuron", cell_class) ~ "salivary_motor",
    # crop
    grepl("crop_motor_neuron", cell_class) ~ "crop_motor",
    # visceral_circulatory - functions defined by Meet
    grepl("^m_NSC_DILP$", cell_type) ~ "energy_metabolism; food_intake",
    grepl("^m_NSC_DH44$", cell_type) ~ 
      "diuresis; food_intake; energy_metabolism; post_mating_responses",
    (grepl("^ventral_nerve_cord_visceral_circulatory$", super_class) &
       grepl("Dh44", banc_neuropeptide_verified) &
       grepl("Lk", banc_neuropeptide_verified)) ~ 
      "diuresis; stress_tolerance; pain_threshold",
    grepl("^m_NSC_DMS$", cell_type) ~ "food_intake",
    grepl("^l_NSC_ITP$", cell_type) ~ 
      "hormone_regulation; stress_tolerance; energy_metabolism; osmotic_balance",
    grepl("^l_NSC_CRZ$", cell_type) ~ 
      "cardiostimulation; energy_metabolism; stress_tolerance; hormone_regulation",
    grepl("^l_NSC_DH31$", cell_type) ~ 
      "diuresis; gut_motility; hormone_regulation",
    grepl("^SEZ_NSC_CAPA$", cell_type) ~ "anti_diuresis; stress_tolerance",
    grepl("^ENXXX012$", cell_type) ~ 
      "anti_diuresis; stress_tolerance; energy_metabolism",
    grepl("^EN27X010$", cell_type) ~ "myomodulation; energy_metabolism",
    # remainder
    TRUE ~ NA
  ) ) %>%
  dplyr::mutate(cell_function_detailed = dplyr::case_when(
    cell_function == cell_function_detailed ~ NA,
    TRUE ~ cell_function_detailed
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    # 'endocrine'
    grepl("^ISN$",cell_type) ~ "nutrient_sensory_neuron",
    grepl("^PSI$",cell_type) ~ "peripheral_intrinsic_neuron",
    grepl("^ascending_visceral_circulatory",super_class) ~  "ventral_nerve_cord_neurosecretory_cell",
    grepl("^m_NSC_unknown$",cell_type) ~  "pars_intercerebralis_neurosecretory_cell",
    grepl("^DNg28$",cell_type) ~  "subesophageal_zone_neurosecretory_cell",
    grepl("^m-NSC$|medial_NSC|medial_neurosecretory_cell",cell_sub_class) ~ "pars_intercerebralis_neurosecretory_cell",
    grepl("^l-NSC$|lateral_NSC|lateral_neurosecretory_cell",cell_sub_class) ~ "pars_lateralis_neurosecretory_cell",
    grepl("SEZ-NSC|subesophageal_zone_neurosecretory_cell",cell_sub_class) ~ "subesophageal_zone_neurosecretory_cell",
    grepl("^PI$",FAFB_cell_type) ~  "pars_intercerebralis_efferent_neuron",
    grepl("FMRFa|CAPA",neuropeptide_verified)&grepl("efferent",flow)&grepl("ventral_nerve_cord",region)&grepl("ANm",MANC_origin) ~  "ventral_nerve_cord_neurosecretory_cell",
    grepl("endocrine",cell_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector)&grepl("ANm",MANC_origin) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$|ANm_endocrine",cell_class)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^EN$",cell_sub_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector)&grepl("ANm",MANC_origin) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("FMRFa|CAPA",neuropeptide_verified)&grepl("efferent",flow)&grepl("ventral_nerve_cord",region) ~  "ventral_nerve_cord_neurosecretory_cell",
    grepl("endocrine",cell_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^EN$",cell_sub_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$|ANm_endocrine",cell_class)&!grepl("hemal",body_part_effector)&grepl("ANm",MANC_origin) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$",cell_sub_class)&grepl("ANm",MANC_origin) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("endocrine",cell_class)&grepl("ventral_nerve_cord",region)&!grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^EN$",cell_sub_class)&grepl("ventral_nerve_cord",region)&!grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^multi_endocrine$",cell_class) ~ "ventral_nerve_cord_efferent_neuron",
    grepl("^EN$",cell_sub_class)&grepl("brain",region) ~ "central_brain_efferent_neuron",
    grepl("endocrine",cell_class)&grepl("brain",region) ~ "central_brain_efferent_neuron",
    cell_sub_class=="NSC" ~ paste0(region,"_neurosecretory_cell"),
    # sensory
    ## classics
    grepl("^ORN",cell_type) ~ "olfactory_receptor_neuron",
    grepl("^GRN",cell_type) ~ "gustatory_receptor_neuron",
    grepl("^HRN",cell_type) ~ "hygrosensory_receptor_neuron",
    grepl("^TRN",cell_type) ~ "thermosensory_receptor_neuron",
    grepl("AC neuron|AC_neuron|ac_neuron",cell_type) ~ "thermosensory_receptor_neuron",
    ### conjunctions
    grepl("sensory",super_class)&!is.na(peripheral_target_type) ~ paste0(peripheral_target_type,"_neuron"),
    grepl("sensory",super_class)&is.na(peripheral_target_type) ~ "orphan_sensory_neuron",
    # grepl("sensory",super_class)&grepl("head|^frontal|^frontoorbital|^orbital|^interocellar|^vibrissa|^interommatidial|^occipital_dorsal|^occipital_ventral|^postorbital_dorsal|^postorbital_ventral|^vertical|^postocellar|^supracervical|^haustellum|^ocellar",body_part_sensory)&!is.na(peripheral_target_type) ~ paste0("head_",peripheral_target_type,"_neuron"),
    # central complex
    grepl("CX|central_complex",cell_class) & grepl("^hDelta|^vDelta|^FC|columnar|^PEN$",cell_sub_class) ~ "central_complex_intrinsic_neuron",
    grepl("CX|central_complex",cell_class) & grepl("^EL$|^EPG|^Delta7|^P6\\-8P9|^P1\\-9|^PFN|^FC",cell_type) ~ "central_complex_intrinsic_neuron",
    grepl("CX|central_complex",cell_class) & grepl("^NO$|PB_input|ER|fb_tangential|FB_tangential|ExR|NO_input|CB.FB|^SA1|^SA2|^SA3|^SAF|^FB",cell_sub_class) ~ "central_complex_input_neuron",
    grepl("CX|central_complex",cell_class) & grepl("^NO$|PB_input|ER|fb_tangential|FB_tangential|ExR|NO_input|CB.FB|^SA1|^SA2|^SA3|^SAF|^FB",cell_type) ~ "central_complex_input_neuron",
    grepl("CX|central_complex",cell_class) & grepl("^FS|^PFL|^PFR|^FR|^PFG|^PEG",cell_sub_class) ~ "central_complex_output_neuron",
    # motor
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & is.na(body_part_effector) ~ "orphan_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & grepl("leg",body_part_effector) ~ "leg_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & grepl("spiracle", cell_class) ~ "spiracle_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & grepl("wing",body_part_effector) ~ "wing_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & grepl("neck",body_part_effector) ~ "neck_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) & grepl("salivary",body_part_effector) ~ "salivary_motor_neuron",
    grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ paste0(body_part_effector,"_motor_neuron"),
    # neck
    super_class=="ascending" ~ "ascending_neuron",
    super_class=="descending" ~ "descending_neuron",
    grepl("sensory_ascending",cell_class) ~ "sensory_ascending_neuron",
    grepl("sensory_descending",cell_class) ~ "sensory_descending_neuron",
    grepl("efferent_ascending",cell_class) ~ "efferent_ascending_neuron",
    grepl("efferent_descending",cell_class) ~ "efferent_descending_neuron",
    # circadian
    grepl("^clock$|circadian",cell_class) ~ "circadian_neuron",
    cell_type %in% "HBeyelet" ~ "circadian_neuron",
    grepl("APDN3|LTe71|SLP1500|LMTe01|SLP249|PLP080",cell_type)~ "circadian_neuron",
    grepl("^LNd|^LNv|^DN1p|^DN1a|^DN1b|^DN3|^LPN$|CB1215|PV7c11",cell_type)~ "circadian_neuron",
    grepl("SLPpl1",FAFB_ito_lee_hemilineage)&grepl("clock|circadian",cell_class)~ "circadian_neuron",
    # optic 
    (grepl("^C2$|^C3$", cell_type)) ~ "optic_lobe_intrinsic_centrifugal",
    (grepl("^Dm\\d+", cell_type)) ~ "distal_medulla",
    (grepl("^CB3849$", cell_type)) ~ "distal_medulla",
    (grepl("^DmDRA", cell_type)) ~ "distal_medulla_dorsal_rim_area",
    (grepl("^Lai$", cell_type)) ~ "lamina_intrinsic",
    (grepl("^L\\d$", cell_type)) ~ "lamina_monopolar",
    (grepl("^L\\d-\\d$", cell_type)) ~ "lamina_monopolar",
    (grepl("^Lat$", cell_type)) ~ "lamina_tangential",
    (grepl("Lawf", cell_type)) ~ "lamina_wide_field",
    (grepl("^Li\\d+$", cell_type)) ~ "lobula_intrinsic",
    (grepl("^CB3815$", cell_type)) ~ "lobula_intrinsic",
    (grepl("^LLPt$", cell_type)) ~ "lobula_lobula_plate_tangential",
    (grepl("^CT1$|^LMa\\d$", cell_type)) ~ "lobula_medulla_amacrine",
    (grepl("^CB3818$", cell_type)) ~ "lobula_medulla_amacrine",
    (grepl("^LMt\\d$", cell_type)) ~ "lobula_medulla_tangential",
    (grepl("CB3820", cell_type)) ~ "lobula_medulla_tangential",
    (grepl("^LPi\\d+$", cell_type)) ~ "lobula_plate_intrinsic",
    (grepl("^CB3826$", cell_type)) ~ "lobula_plate_intrinsic",
    (grepl("^Mi\\d+$", cell_type)) ~ "medulla_intrinsic",
    (grepl("^Am1$", cell_type)) ~ "medulla_lobula_lobula_plate_amacrine",
    (grepl("^MLt\\d$", cell_type)) ~ 'medulla_lobula_tangential',
    (grepl("^CB3833$", cell_type)) ~ 'medulla_lobula_tangential',
    (grepl("^PDt$", cell_type)) ~ "proximal_distal_medulla_tangential",
    (grepl("^Pm\\d+$", cell_type)) ~ "proximal_medulla",
    (grepl("^Sm\\d+$", cell_type)) ~ "serpentine_medulla",
    (grepl("^CB3825$|^CB3832$", cell_type)) ~ "serpentine_medulla",
    (grepl("^T\\d", cell_type)) ~ 'transverse_neuron',
    (grepl("^Tlp\\d+$", cell_type)) ~ "translobula_plate",
    (grepl("^Tm\\d+", cell_type)) ~ "transmedullary",
    (grepl("^CB3851$|^CB3864$", cell_type)) ~ "transmedullary", 
    (grepl("^TmY\\d+", cell_type)) ~ "transmedullary_y",
    (grepl("^CB3816$", cell_type)) ~ "transmedullary_y",
    (grepl("^Y\\d+$", cell_type)) ~ "y_neuron",
    (grepl("^CB3846$", cell_type)) ~ "y_neuron",
    (grepl("^LC", cell_type)) ~ "lobula_columnar",
    (grepl("^LT", cell_type)) ~ "lobula_tangential",
    (grepl("^LPC\\d+$", cell_type)) ~ "lobula_plate_columnar",
    (grepl("^LPLC\\d+$", cell_type)) ~ "lobula_plate_lobula_columnar",
    (grepl("^LLPC\\d+$", cell_type)) ~ "lobula_lobula_plate_columnar",
    (grepl("^LPT|^Nod|^VS|^HS|^H1$|^H2$|^DCH$|^FD1$|^FD3$|^V1$|^vCal1$|^VCH$", cell_type)) 
    ~ "lobula_plate_tangential_cell",
    (grepl("^MeMe", cell_type)) ~  "medulla_medulla",
    (grepl("^MeTu", cell_type)) ~ "medulla_tubercle",
    (grepl("^MeLp", cell_type)) ~ "medulla_lobula_plate",
    (grepl("^aMe", cell_type)) ~ "amacrine_medulla",
    (grepl("^MTe|MC65", cell_type)) ~ "medulla_tangential",
    (grepl("^mALC\\d+", cell_type)) ~ "medial_antennal_lobula",
    (grepl("^cM\\d+", cell_type)) ~ "centrifugal_medulla",
    (grepl("^cML\\d+", cell_type)) ~ "centrifugal_medulla_lobula",
    (grepl("^cMLLP\\d+", cell_type)) ~ "centrifugal_medulla_lobula_lobula_plate",
    (grepl("^cL\\d+", cell_type)) ~ "centrifugal_lobula",
    (grepl("^cLM\\d+", cell_type)) ~ "centrifugal_lobula_medulla",
    (grepl("^cLP", cell_type)) ~ "centrifugal_lobula_plate",
    (grepl("^cLLP\\d+", cell_type)) ~ "centrifugal_lobula_lobula_plate",
    (grepl("cLLPM\\d+", cell_type)) ~ "centrifugal_lobula_lobula_plate_medulla",
    (grepl("^OA-A", cell_type)) ~ "optic_anterior",
    (grepl("^R\\d$|^R1-6$", cell_type)) ~ "photoreceptor",
    grepl("^TuBu$",cell_class)|grepl("^TuBu$",cell_type) ~ "tubercular_bulbar_neuron",
    grepl("ocellar_projection",cell_class)|grepl("ocellar_projection",cell_sub_class) ~ "ocellar_projection_neuron",
    grepl("^OCC",cell_type)|grepl("ocellar_centrifugal",cell_sub_class) ~ "ocellar_centrifugal_neuron",
    grepl("^OCI|^OCL",cell_type)|grepl("ocellar_interneuron",cell_sub_class) ~ "ocellar_intrinsic_neuron",
    grepl("ocellar",cell_sub_class) ~ "ocellar_intrinsic_neuron",
    # classics
    grepl("^ALIN$",cell_class) ~ "antennal_lobe_centrifugal_neuron",
    grepl("^ALON$",cell_class) ~ "antennal_lobe_output_neuron",
    grepl("^ALPN$",cell_class) ~ "antennal_lobe_projection_neuron",
    grepl("^ALLN$",cell_class) ~ "antennal_lobe_local_neuron",
    grepl("^TPN$",cell_class) ~ "subesophageal_zone_projection_neuron",
    grepl("^water_PN|^mAL",cell_sub_class) ~ "subesophageal_zone_projection_neuron",
    grepl("^water_PN|^mAL$",cell_class)|grepl("^mAL$",cell_sub_class)  ~ "subesophageal_zone_projection_neuron",
    grepl("^WEDPN$",cell_class) ~ "wedge_projection_neuron",
    grepl("^MBON$",cell_class) ~ "mushroom_body_output_neuron",
    grepl("^MBIN$",cell_class) ~ "mushroom_body_extrinsic_neuron",
    cell_type %in% c("APL","DPM") ~ "mushroom_intrinsic_body",
    grepl("^DAN$",cell_class)&grepl("^PAM|^PPL",cell_type) ~ "mushroom_body_dopaminergic_neuron",
    grepl("^DAN$",cell_class)&grepl("^PPM",cell_type) ~ "PPM_dopaminergic_neuron",
    #grepl("^TOON$",cell_class) ~ "third_order_olfactory_neuron",
    grepl("^LHON$",cell_class) ~ "lateral_horn_output_neuron",
    grepl("^LHLN$",cell_class) ~ "lateral_horn_local_neuron",
    grepl("^LHCENT$",cell_class) ~ "lateral_horn_centrifugal_neuron",
    grepl("^KC$|^Kenyon_cell|Kenyon_Cell",cell_class) ~ "kenyon_cell",
    grepl("^bilateral$",cell_class) ~ "bilateral_neuron",
    # VNC intrinsic
    grepl("independent leg",MANC_serialMotif) ~ "single_leg_neuromere",
    grepl("sequential",MANC_serialMotif) ~ "sequential_leg_neuromeres",
    grepl("centrifugal",MANC_serialMotif) ~ "ventral_nerve_cord_centrifugal",
    grepl("dorsal",MANC_serialMotif) ~ "ventral_nerve_cord_dorsal",
    grepl("centripetal",MANC_serialMotif) ~ "ventral_nerve_cord_centripetal",
    grepl("convergent",MANC_serialMotif) ~ "ventral_nerve_cord_serially_convergent",
    # bad
    grepl("^Interneuron_TBD$",cell_class) ~ NA,
    grepl("^TBD$",FAFB_cell_class) ~ NA,
    grepl("^TBD$",cell_class) ~ NA,
    TRUE ~ NA
  )) %>%
  #dplyr::rowwise() %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    grepl("^ISN$",cell_type) ~ "hemolymph_sensory_neuron",
    # chordotonals
    grepl("chordotonal_club|club_chordotonal",cell_function_detailed) ~ paste0(body_part_sensory,"_club_chordotonal_organ_neuron"),
    grepl("chordotonal_claw|claw_chordotonal",cell_function_detailed) ~ paste0(body_part_sensory,"_claw_chordotonal_organ_neuron"),
    grepl("chordotonal_hook|hook_chordotonal",cell_function_detailed) ~ paste0(body_part_sensory,"_hook_chordotonal_organ_neuron"),
    grepl("chordotonal_club|club_chordotonal",cell_sub_class) ~ paste0(body_part_sensory,"_club_chordotonal_organ_neuron"),
    grepl("chordotonal_claw|claw_chordotonal",cell_sub_class) ~ paste0(body_part_sensory,"_claw_chordotonal_organ_neuron"),
    grepl("chordotonal_hook|hook_chordotonal",cell_sub_class) ~ paste0(body_part_sensory,"_hook_chordotonal_organ_neuron"),
    grepl("chordotonal_claw",cell_function_detailed) ~ paste0(body_part_sensory,"_claw_chordotonal_organ_neuron"),
    grepl("chordotonal_club",cell_function_detailed) ~ paste0(body_part_sensory,"_club_chordotonal_organ_neuron"),
    grepl("chordotonal_hook",cell_function_detailed) ~ paste0(body_part_sensory,"_hook_chordotonal_organ_neuron"),
    grepl("^JO-A",cell_type) ~ "johnstons_organ_A_neuron",
    grepl("^JO-B",cell_type) ~ "johnstons_organ_B_neuron",
    grepl("^JO-C",cell_type) ~ "johnstons_organ_C_neuron",
    grepl("^JO-D",cell_type) ~ "johnstons_organ_D_neuron",
    grepl("^JO-E",cell_type) ~ "johnstons_organ_E_neuron",
    grepl("^JO-F",cell_type) ~ "johnstons_organ_F_neuron",
    grepl("^JO-",cell_type) ~ "johnstons_organ_other_neuron",
    grepl("chordotonal organ, JO",cell_class) ~ "johnstons_organ_neuron",
    # sensory
    grepl("AC neuron|AC_neuron|ac_neuron",cell_type) ~ "internal_thermosensory_receptor_neuron",
    grepl("sensory",super_class)&!is.na(body_part_sensory)&!grepl("metathoracic|wheelers|prothoracic",cell_class) ~ paste0(gsub("_chordotonal_organ|_organ","",body_part_sensory),"_",cell_class),
    grepl("sensory",super_class)&!is.na(body_part_sensory) ~ paste0(body_part_sensory,"_",cell_class),
    # ALPNS
    grepl("multiglomerular",cell_sub_class) ~ "multiglomerular_projection_neuron",
    grepl("uniglomerular",cell_sub_class) ~ "uniglomerular_projection_neuron",
    # bristles
    grepl("_bristle",cell_sub_class) ~ cell_sub_class,
    # central complex
    FAFB_cell_class=="CX" & grepl("hDelta",cell_type) ~ "hDelta",
    FAFB_cell_class=="CX" & grepl("vDelta",cell_type) ~ "vDelta",
    FAFB_cell_class=="CX" & grepl("^FB|^CB.FB|^SA1|^SA2|^SA3",cell_type) ~ "FB_tangential",
    FAFB_cell_class=="CX" & grepl("^ExR",cell_type) ~ "ExR",
    FAFB_cell_class=="CX" & grepl("^FC",cell_type) ~ "FC",
    FAFB_cell_class=="CX" & grepl("^FS",cell_type) ~ "FS",
    FAFB_cell_class=="CX" & grepl("^PFL",cell_type) ~ "PFL",
    FAFB_cell_class=="CX" & grepl("^PFN",cell_type) ~ "PFN",
    FAFB_cell_class=="CX" & grepl("^ER",cell_type) ~ "EB_input",
    FAFB_cell_class=="CX" & grepl("^FR1|^FR2",cell_type) ~ "FR",
    FAFB_cell_class=="CX" & grepl("^PFR$",cell_type) ~ "PFR",
    FAFB_cell_class=="CX" & grepl("^PFGs$",cell_type) ~ "PFG",
    FAFB_cell_class=="CX" & grepl("^PEG",cell_type) ~ "PEG",
    FAFB_cell_class=="CX" & grepl("^PEN_",cell_type) ~ "PEN",
    FAFB_cell_class=="CX" & grepl("^GLNO$|^LCNOp$|^LCNOpm$|^LNO1$|^LNO1$|^LNO2$",cell_type) ~ "NO_input",
    FAFB_cell_class=="CX" & grepl("^IbSpsP|^LPsP|^SpsP",cell_type) ~ "PB_input",
    # kenyon cells
    grepl("^KCab",cell_type)~ "KCab",
    grepl("^KCapbp",cell_type)~ "KCapbp",
    grepl("^KCg",cell_type)~ "KCg",
    # wing
    #grepl("^w-cHIN",cell_type)~ "w-cHIN",
    #grepl("^n-cHIN",cell_type)~ "n-cHIN",
    # dopamine
    grepl("^PAM",cell_type) ~ "PAM_dopaminergic_neuron",
    grepl("^CB2730",cell_type) ~ "PAM_dopaminergic_neuron",
    grepl("^PPL1",cell_type)~ "PPL1_dopaminergic_neuron",
    grepl("^PPL2",cell_type)~ "PPL2_dopaminergic_neuron",
    grepl("^PPM",cell_type)|grepl("^PPM",cell_type)~ "PPM_dopaminergic_neuron",
    # circadian
    grepl("APDN3|LTe71|SLP1500|LMTe01|SLP249|PLP080",cell_type)~ "APDN3",
    grepl("^LNd",cell_type)~ "LNd",
    grepl("^LNv",cell_type)~ "LNv",
    grepl("^DN1p",cell_type)~ "DN1p",
    grepl("^DN1a",cell_type)~ "DN1a",
    grepl("^DN1b",cell_type)~ "DN1b",
    grepl("^DN3",cell_type)~ "DN3",
    grepl("^DN1p",cell_sub_class)~ "DN1p",
    grepl("^DN1a",cell_sub_class)~ "DN1a",
    grepl("^DN1b",cell_sub_class)~ "DN1b",
    grepl("^DN3",cell_sub_class)~ "DN3",
    grepl("^LPN$|CB1215|PV7c11",cell_type)~ "LPN",
    grepl("SLPpl1",FAFB_ito_lee_hemilineage)&grepl("clock|circadian",cell_class)~ "s-CPDN",
    cell_type %in% "HBeyelet" ~ "hofbauer_buchner_eyelet_neuron",
    # ventral nerve cord
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="BR" ~ "ventral_nerve_cord_bilateral_restricted", 
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="IR" ~ "ventral_nerve_cord_ipsilateral_restricted", 
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="CR" ~ "ventral_nerve_cord_contralateral_restricted", 
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="BI" ~ "ventral_nerve_cord_bilateral_interconnecting", 
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="II" ~ "ventral_nerve_cord_ipsilateral_interconnecting", 
    grepl("ventral_nerve_cord_intrinsic",super_class)&MANC_subclass=="CI" ~ "ventral_nerve_cord_contralateral_interconnecting",
    # motor
    grepl("front_leg",body_part_effector)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "front_leg_motor_neuron",
    grepl("haltere_power",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "haltere_power_neuron",
    grepl("haltere_steering",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "haltere_steering_neuron",
    grepl("hind_leg",body_part_effector)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "hind_leg_motor_neuron",
    grepl("jump_escape",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "middle_leg_extension_jump_motor_neuron",
    grepl("middle_leg",body_part_effector)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "middle_leg_motor_neuron",
    grepl("neck_pitch",cell_function_detailed)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "neck_pitch_motor_neuron",
    grepl("neck_roll",cell_function_detailed)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "neck_roll_motor_neuron",
    grepl("neck_yaw",cell_function_detailed)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "neck_yaw_motor_neuron",
    grepl("wing_power",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "wing_power_motor_neuron",
    grepl("wing_steering",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "wing_steering_motor_neuron",
    grepl("wing_tension",cell_function)&grepl("central_brain_motor|ventral_nerve_cord_motor",super_class) ~ "wing_tension_motor_neuron",
    # endocrine
    grepl("^PSI$",cell_type) ~ paste0(body_part_effector, "_peripheral_intrinsic_neuron"),
    grepl("visceral",super_class) ~ paste0(body_part_effector,gsub("pars_intercerebralis_|pars_lateralis_|ventral_nerve_cord_|central_brain_|subesophageal_zone_|abdominal_neuromere_","_",cell_class)),
    #grepl("^l_NSC_unknown|l_NSC_DH31",cell_type)~ "circadian_neuroendocrine_neuron", # DNES
    #grepl("abdomen_endocrine",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "abdomen_endocrine",
    # grepl("pars_intercerebralis_endocrine_enteric",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "pars_intercerebralis_enteric_endocrine_neuron",
    # grepl("pars_intercerebralis_endocrine_unknown",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "pars_intercerebralis_unknown_endocrine_neuron",
    # grepl("pars_lateralis_endocrine_corpus_allatum",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "pars_lateralis_corpus_allatum_endocrine_neuron",
    # grepl("pars_lateralis_endocrine_retrocerebral_complex",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "pars_lateralis_retrocerebral_complex_endocrine_neuron",
    # grepl("pars_lateralis_endocrine_unknown",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "pars_lateralis_unknown_endocrine_neuron",
    # !is.na(body_part_effector)&grepl("endocrine|visceral_circulatory",super_class) ~ paste0(body_part_effector,"_endocrine_neuron"),
    # !is.na(body_part_effector)&grepl("^EN$",cell_class)|grepl("^EN$",cell_sub_class) ~ paste(body_part_effector,"_endocrine"),
    # grepl("SEZ_endocrine",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "subesophageal_zone_endocrine",
    # grepl("wing_endocrine",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "wing_endocrine_neuron",
    # grepl("neurohemal_complex_endocrine",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "neurohemal_complex_endocrine_neuron",
    # grepl("ovulation",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "ovulation_endocrine_neuron",
    # grepl("peripheral_serotonin",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "peripheral_serotonin_neuron",
    # grepl("^endocrine$",cell_function)&grepl("endocrine|visceral_circulatory",super_class) ~ "endocrine_neuron",
    # other
    grepl("^water_PN",cell_sub_class) ~ "BiT",
    grepl("^mAL$",cell_class)|grepl("^mAL$",cell_sub_class) ~ "mAL",
    TRUE ~ NA
  )) 
#dplyr::mutate(cell_class = gsub("\\>|\\.","-",cell_class)) %>%


# select changed columns
franken.meta.update.toPush <- franken.meta.update %>%
  select(flow, super_class, cell_class, cell_sub_class, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, nerve, banc_neuropeptide_verified, 
         banc_neurotransmitter_verified, side, region, `_id`)

# add old franken_meta columns to updated version, for relevant columns, as 
#  franken_[column name]
# columns that have changed with above code
franken.meta.keep <- franken.meta %>%
  # rename columns
  mutate(franken_flow = flow, 
         franken_super_class = super_class,
         franken_cell_class = cell_class,
         franken_cell_sub_class = cell_sub_class,
         franken_cell_function = cell_function,
         franken_cell_function_detailed = cell_function_detailed,
         franken_peripheral_target_type = peripheral_target_type,
         franken_body_part_sensory = body_part_sensory,
         franken_body_part_effector = body_part_effector,
         franken_nerve = nerve,
         franken_side = side,
         franken_region = region,) %>%
  select(franken_flow, franken_super_class, franken_cell_class,
         franken_cell_sub_class, franken_cell_function, 
         franken_cell_function_detailed, franken_peripheral_target_type,
         franken_body_part_sensory, franken_body_part_effector, 
         franken_nerve, franken_side, franken_region,
         neurotransmitter_verified, neuropeptide_verified, `_id`)

# join with franken.meta.update.toPush
franken.meta.joined <- franken.meta.update.toPush %>%
  inner_join(franken.meta.keep, by = "_id") 

# write to franken_meta SeaTable
banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.joined),
                      append_allowed = FALSE,
                      chunksize = 1000)


# 250619 - policy update to super_class - implement here, given previous update
#  has already been pushed to the table
franken.meta.update <- franken.meta %>%
  mutate(super_class = case_when(
    grepl("central_brain_motor|ventral_nerve_cord_motor", super_class) ~ "motor",
    grepl("central_brain_sensory|optic_lobe_sensory|ventral_nerve_cord_sensory", 
          super_class) ~ "sensory",
    grepl("central_brain_visceral_circulatory|ventral_nerve_cord_visceral_circulatory",
          super_class) ~ "visceral_circulatory",
    TRUE ~ super_class
  ))

franken.meta.update.toPush <- franken.meta.update %>%
  select(super_class, banc_id, `_id`)

banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.update.toPush),
                      append_allowed = FALSE,
                      chunksize = 1000)

# 250626 - update to nerve - right and left cervical nerve will only apply to 
#  neurons that leave the side of the neck and not regular ascending and 
#  descending neurons

franken.meta.update <- franken.meta %>%
  filter((grepl("cervical_nerve", nerve) & !grepl("ventral", nerve) & 
            !grepl("motor", super_class))) %>%
  mutate(nerve = case_when(
    (grepl("cervical_nerve", nerve) & !grepl("ventral", nerve) & 
      !grepl("motor", super_class)) ~ NA,
    TRUE ~ nerve
  ))

franken.meta.update.toPush <- franken.meta.update %>%
  select(nerve, sexual_dimorphism, `_id`)

banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.update.toPush),
                      append_allowed = FALSE,
                      chunksize = 1000)

# 250627 - super_class for TmY14 is wrong
franken.meta.update <- franken_meta %>%
  filter(grepl("^TmY14$", cell_type)) %>%
  mutate(super_class = case_when(
    grepl("^TmY14$", cell_type) ~ "optic_lobe_intrinsic",
    TRUE ~ super_class
    )) %>%
      select(super_class, sexual_dimorphism, `_id`)

banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.update),
                      append_allowed = FALSE,
                      chunksize = 1000)

# 250627 - there are a subset of ascending neurons with flow afferent, which is wrong
franken.meta.update <- franken_meta %>%
  filter(grepl("afferent", flow) & grepl("^ascending$", super_class)) %>%
  mutate(flow = case_when(
    (grepl("afferent", flow) & grepl("^ascending$", super_class)) ~ "intrinsic",
    TRUE ~ flow
  )) %>%
  select(flow, sexual_dimorphism, `_id`)

banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.update),
                      append_allowed = FALSE,
                      chunksize = 1000)

# 6/29/25 - missing organ_neuron after hook_chordotonal for small subset
franken.meta.update <- franken.meta %>%
  filter(grepl("hook_chordotonal", cell_sub_class) & !grepl("neuron", cell_sub_class)) %>%
  mutate(cell_sub_class = gsub("hook_chordotonal", "hook_chordotonal_organ_neuron",
                               cell_sub_class)) %>%
  select(`_id`, cell_sub_class, sexual_dimorphism)
    
banctable_update_rows(base='cns_meta',
                      table = "franken_meta",
                      df = as.data.frame(franken.meta.update),
                      append_allowed = FALSE,
                      chunksize = 1000)




# View and check a chosen snippet
snippet <- franken.meta.update %>% 
  #  dplyr::filter(grepl("visceral",super_class)) %>%
  dplyr::group_by(flow,super_class, cell_class, cell_sub_class,cell_function, cell_function_detailed) %>%
  dplyr::summarize(
    count = dplyr::n()
  ) %>%
  dplyr::arrange(flow,super_class, cell_class, cell_sub_class, cell_function, cell_function_detailed)
knitr::kable(snippet)
clipr::write_clip(knitr::kable(snippet))

### STANDARDIZATION TESTS FOR BANC ANNOTATIONS ###

# Helper function to compare original and updated values for a column
compare_and_summarize <- function(original_df, updated_df, column_name) {
  if(!(column_name %in% colnames(original_df) && column_name %in% colnames(updated_df))) {
    return(paste("Column", column_name, "not found in one or both dataframes"))
  }
  
  # Get original and updated values
  orig_values <- original_df[[column_name]]
  updt_values <- updated_df[[column_name]]
  
  # Create mapping table
  mapping_df <- data.frame(
    original = orig_values,
    updated = updt_values,
    row_id = 1:length(orig_values)
  ) %>%
    filter(original != updated | (is.na(original) & !is.na(updated)) | (!is.na(original) & is.na(updated)))
  
  # Create unique mapping combinations with counts
  unique_mapping_df <- mapping_df %>%
    dplyr::group_by(original, updated) %>%
    dplyr::summarize(
      count = dplyr::n(), 
      example_row_id = first(row_id), 
      .groups = 'drop'
    ) %>%
    dplyr::arrange(desc(count))
  
  # Generate summary
  orig_counts <- table(orig_values, useNA = "ifany")
  updt_counts <- table(updt_values, useNA = "ifany")
  
  orig_summary <- data.frame(
    value = names(orig_counts),
    count = as.numeric(orig_counts),
    dataset = "original"
  )
  
  updt_summary <- data.frame(
    value = names(updt_counts),
    count = as.numeric(updt_counts),
    dataset = "updated"
  )
  
  summary_df <- rbind(orig_summary, updt_summary)
  
  return(list(
    mapping = mapping_df,
    unique_mapping = unique_mapping_df,
    summary = summary_df,
    changed_count = nrow(mapping_df),
    total_count = length(orig_values),
    change_percent = round(nrow(mapping_df) / length(orig_values) * 100, 2)
  ))
}

# Function to print a nicely formatted report for a column
print_column_report <- function(comparison_result, column_name) {
  cat("\n\n======================================\n")
  cat("ANALYSIS FOR COLUMN:", column_name, "\n")
  cat("======================================\n\n")
  
  cat("CHANGE STATISTICS:\n")
  cat("Total rows:", comparison_result$total_count, "\n")
  cat("Changed rows:", comparison_result$changed_count, "\n")
  cat("Change percentage:", comparison_result$change_percent, "%\n\n")
  
  if(nrow(comparison_result$unique_mapping) > 0) {
    cat("UNIQUE MAPPING COMBINATIONS (top 150):\n")
    unique_mappings_to_show <- head(comparison_result$unique_mapping, 150)
    
    # Format NA values for better display
    unique_mappings_to_show$original <- ifelse(is.na(unique_mappings_to_show$original), 
                                               "NA", 
                                               as.character(unique_mappings_to_show$original))
    unique_mappings_to_show$updated <- ifelse(is.na(unique_mappings_to_show$updated), 
                                              "NA", 
                                              as.character(unique_mappings_to_show$updated))
    
    # Calculate percentage of total changes
    unique_mappings_to_show$percent <- round(unique_mappings_to_show$count / comparison_result$changed_count * 100, 1)
    
    # Create a nice display format
    for(i in 1:nrow(unique_mappings_to_show)) {
      cat(sprintf("'%s' → '%s': %d occurrences (%.1f%% of changes, example row: %d)\n", 
                  unique_mappings_to_show$original[i],
                  unique_mappings_to_show$updated[i], 
                  unique_mappings_to_show$count[i],
                  unique_mappings_to_show$percent[i],
                  unique_mappings_to_show$example_row_id[i]))
    }
    
    if(nrow(comparison_result$unique_mapping) > 150) {
      cat("... and", nrow(comparison_result$unique_mapping) - 150, "more combinations\n")
    }
    cat("\n")
  }
  
  cat("VALUE INVENTORY (top 150 values before and after):\n")
  orig_top <- comparison_result$summary %>% 
    dplyr::filter(dataset == "original") %>%
    dplyr::arrange(desc(count)) %>%
    head(150)
  
  updt_top <- comparison_result$summary %>% 
    dplyr::filter(dataset == "updated") %>%
    dplyr::arrange(desc(count)) %>%
    head(150)
  
  cat("Original top values:\n")
  print(orig_top)
  cat("\nUpdated top values:\n")
  print(updt_top)
  
  # Check if any values disappeared or appeared
  orig_values <- comparison_result$summary$value[comparison_result$summary$dataset == "original"]
  updt_values <- comparison_result$summary$value[comparison_result$summary$dataset == "updated"]
  
  disappeared <- setdiff(orig_values, updt_values)
  appeared <- setdiff(updt_values, orig_values)
  
  if(length(disappeared) > 0) {
    cat("\nValues that disappeared after standardization:\n")
    print(disappeared)
  }
  
  if(length(appeared) > 0) {
    cat("\nNew values that appeared after standardization:\n")
    print(appeared)
  }
}

# List of columns to check
columns_to_check <- c(
  "citation_cell_type", "flow", "side", "region", "peripheral_target_type",
  "cell_function_detailed", "cell_function", "nerve", "body_part_sensory",
  "body_part_effector", "super_class", "cell_sub_class", "cell_class"
)

# Run tests for all columns
test_results <- list()
for(column in columns_to_check) {
  test_results[[column]] <- compare_and_summarize(franken.meta, franken.meta.update, column)
  print_column_report(test_results[[column]], column)
}

# Additional validation tests
cat("\n\n======================================\n")
cat("ADDITIONAL VALIDATION CHECKS\n")
cat("======================================\n\n")

# Check for any NAs introduced in critical columns
for(column in columns_to_check) {
  orig_na_count <- sum(is.na(franken.meta[[column]]))
  updt_na_count <- sum(is.na(franken.meta.update[[column]]))
  
  if(updt_na_count > orig_na_count) {
    cat("WARNING: Column", column, "has more NAs after standardization.\n")
    cat("  Original NA count:", orig_na_count, "\n")
    cat("  Updated NA count:", updt_na_count, "\n")
    cat("  Difference:", updt_na_count - orig_na_count, "\n\n")
  }
}

# Check for standardization consistency
for(column in columns_to_check) {
  mapping_df <- test_results[[column]]$mapping
  if(nrow(mapping_df) > 0) {
    # Check if same original value maps to different updated values
    inconsistencies <- mapping_df %>%
      dplyr::filter(!is.na(original)) %>%
      dplyr::group_by(original) %>%
      dplyr::summarize(unique_updates = dplyr::n_distinct(updated, na.rm = TRUE)) %>%
      dplyr::filter(unique_updates > 1)
    
    if(nrow(inconsistencies) > 0) {
      cat("WARNING: Inconsistent mappings found in column", column, "\n")
      for(i in 1:nrow(inconsistencies)) {
        orig_val <- inconsistencies$original[i]
        cat("  Original value '", orig_val, "' maps to multiple updated values:\n", sep="")
        
        examples <- mapping_df %>%
          dplyr::filter(original == orig_val) %>%
          dplyr::group_by(updated) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup()
        
        for(j in 1:nrow(examples)) {
          cat("    -> '", examples$updated[j], "' (example row: ", examples$row_id[j], ")\n", sep="")
        }
      }
      cat("\n")
    }
  }
}

# Print summary statistics for all columns
cat("\n\n======================================\n")
cat("OVERALL STANDARDIZATION SUMMARY\n")
cat("======================================\n\n")

summary_df <- data.frame(
  column = character(),
  total_rows = integer(),
  changed_rows = integer(),
  change_percent = numeric(),
  original_distinct = integer(),
  updated_distinct = integer(),
  distinct_change = integer(),
  stringsAsFactors = FALSE
)

for(column in columns_to_check) {
  result <- test_results[[column]]
  orig_distinct <- length(unique(franken.meta[[column]]))
  updt_distinct <- length(unique(franken.meta.update[[column]]))
  
  summary_df <- rbind(summary_df, data.frame(
    column = column,
    total_rows = result$total_count,
    changed_rows = result$changed_count,
    change_percent = result$change_percent,
    original_distinct = orig_distinct,
    updated_distinct = updt_distinct,
    distinct_change = updt_distinct - orig_distinct
  ))
}

summary_df <- summary_df %>% dplyr::arrange(desc(change_percent))
print(summary_df)

# Make full changelog
franken.meta.changelog <- data.frame(
  timestamp = Sys.time(),
  columns_changed = length(columns_to_check),
  total_rows = nrow(franken.meta),
  total_changes = sum(sapply(test_results, function(x) x$changed_count))
)

# Add a detailed description of changes
changelog_text <- paste0(
  "BANC METADATA STANDARDIZATION CHANGELOG\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Total rows processed: ", nrow(franken.meta), "\n",
  "Total columns modified: ", length(columns_to_check), "\n\n",
  "SUMMARY OF CHANGES BY COLUMN:\n"
)

# Add summary stats for each column
for (column in columns_to_check) {
  result <- test_results[[column]]
  orig_distinct <- length(unique(franken.meta[[column]]))
  updt_distinct <- length(unique(franken.meta.update[[column]]))
  
  changelog_text <- paste0(
    changelog_text,
    "Column: ", column, "\n",
    "  - Changed values: ", result$changed_count, " (", result$change_percent, "% of data)\n",
    "  - Original distinct values: ", orig_distinct, "\n",
    "  - Updated distinct values: ", updt_distinct, "\n",
    "  - Net change in distinct values: ", updt_distinct - orig_distinct, "\n\n"
  )
  
  # Add top 150 most frequent mappings
  if (nrow(result$unique_mapping) > 0) {
    changelog_text <- paste0(
      changelog_text,
      "  Top 150 most common transformations:\n"
    )
    
    top_mappings <- head(result$unique_mapping, 150)
    
    for (i in 1:nrow(top_mappings)) {
      orig_val <- if(is.na(top_mappings$original[i])) "NA" else as.character(top_mappings$original[i])
      updt_val <- if(is.na(top_mappings$updated[i])) "NA" else as.character(top_mappings$updated[i])
      
      changelog_text <- paste0(
        changelog_text,
        "    '", orig_val, "' → '", updt_val, "': ", 
        top_mappings$count[i], " occurrences (", 
        round(top_mappings$count[i] / result$changed_count * 100, 1), "% of changes)\n"
      )
    }
    changelog_text <- paste0(changelog_text, "\n")
  }
  
  # Add values that disappeared
  orig_values <- result$summary$value[result$summary$dataset == "original"]
  updt_values <- result$summary$value[result$summary$dataset == "updated"]
  
  disappeared <- setdiff(orig_values, updt_values)
  appeared <- setdiff(updt_values, orig_values)
  
  if (length(disappeared) > 0) {
    changelog_text <- paste0(
      changelog_text,
      "  Values that were removed (up to 150):\n    ",
      paste(head(disappeared, 150), collapse = ", "),
      if(length(disappeared) > 150) paste0(" (and ", length(disappeared) - 150, " more)") else "",
      "\n\n"
    )
  }
  
  if (length(appeared) > 0) {
    changelog_text <- paste0(
      changelog_text,
      "  New values that were added (up to 150):\n    ",
      paste(head(appeared, 150), collapse = ", "),
      if(length(appeared) > 150) paste0(" (and ", length(appeared) - 150, " more)") else "",
      "\n\n"
    )
  }
  
  # Check for inconsistent mappings
  inconsistencies <- result$mapping %>%
    dplyr::filter(!is.na(original)) %>%
    dplyr::group_by(original) %>%
    dplyr::summarize(unique_updates = dplyr::n_distinct(updated, na.rm = TRUE)) %>%
    dplyr::filter(unique_updates > 1)
  
  if (nrow(inconsistencies) > 0) {
    changelog_text <- paste0(
      changelog_text,
      "  WARNING: Inconsistent mappings found for ", nrow(inconsistencies), " original value(s):\n"
    )
    
    for (i in 1:min(nrow(inconsistencies), 5)) {
      orig_val <- inconsistencies$original[i]
      
      examples <- result$mapping %>%
        dplyr::filter(original == orig_val) %>%
        dplyr::group_by(updated) %>%
        dplyr::slice(1)
      
      changelog_text <- paste0(
        changelog_text,
        "    '", orig_val, "' maps to: ", 
        paste(sapply(examples$updated, function(x) ifelse(is.na(x), "NA", as.character(x))), collapse = ", "),
        "\n"
      )
    }
    
    if (nrow(inconsistencies) > 5) {
      changelog_text <- paste0(
        changelog_text,
        "    (and ", nrow(inconsistencies) - 5, " more inconsistent mappings)\n"
      )
    }
    changelog_text <- paste0(changelog_text, "\n")
  }
}

# Complete the object with the text
franken.meta.changelog$changes_detail <- changelog_text

# Save original
franken.meta.concise <- franken.meta.update %>%
  dplyr::distinct(neuron_id,
                  dataset,
                  side, 
                  nerve,
                  hemilineage,
                  flow,
                  super_class,
                  cell_class,
                  cell_sub_class,
                  cell_type,
                  body_part_sensory,
                  body_part_effector,
                  peripheral_target_type,
                  cell_function,
                  cell_function_detailed,
                  top_nt)
save.path <- "/Users/abates/projects/flyconnectome/bancpipeline/data/meta"
readr::write_csv(franken.meta, file.path(save.path,"franken_meta_v1.csv"))
readr::write_csv(franken.meta.update, file.path(save.path,"franken_meta_update_for_v2.csv"))
readr::write_csv(franken.meta.concise, file.path(save.path,"franken_meta_concise_update_for_v2.csv"))
writeLines(changelog_text, file.path(save.path,"franken_meta_v1_to_v2_changelog.txt"))

# # update  
# banctable_update_rows(base='cns_meta',
#                      table = "franken_meta",
#                      df = as.data.frame(franken.meta.update),
#                      append_allowed = FALSE,
#                      chunksize = 1000)

#######################
### Plots to assess ###
#######################

# Create a heatmap to visualize the complete hierarchy
# This helps show which combinations exist and their frequencies
hierarchy_heatmap <- franken.meta.update %>%
  dplyr::filter(!is.na(flow) & !is.na(super_class) & !is.na(cell_class)) %>%
  dplyr::count(flow, super_class, cell_class) %>%
  # Transform to log scale to better show differences
  dplyr::mutate(log_n = log10(n + 1))

# Get top cell classes by count for visualization clarity
top_cell_classes <- hierarchy_heatmap %>%
  dplyr::group_by(cell_class) %>%
  dplyr::summarize(total = sum(n)) %>%
  dplyr::arrange(desc(total)) %>%
  utils::head(40) %>%
  dplyr::pull(cell_class)

p1 <- hierarchy_heatmap %>%
  dplyr::filter(cell_class %in% top_cell_classes) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(cell_class, log_n), 
                               y = reorder(super_class, log_n), 
                               fill = log_n)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~flow, ncol = 1, scales = "free_y") +
  ggplot2::scale_fill_viridis_c(name = "log10(count+1)") +
  ggplot2::labs(title = "Hierarchical Relationship Matrix",
                subtitle = "flow > super_class > cell_class",
                x = "cell_class", y = "super_class") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
                 axis.text.y = ggplot2::element_text(size = 7),
                 panel.grid = ggplot2::element_blank())


# Create a heatmap to visualize the complete hierarchy
# This helps show which combinations exist and their frequencies
hierarchy_heatmap <- franken.meta.update %>%
  dplyr::filter(!is.na(flow) & !is.na(super_class) & !is.na(cell_class),
                super_class == "sensory") %>%
  dplyr::count(cell_class, cell_sub_class) %>%
  # Transform to log scale to better show differences
  dplyr::mutate(log_n = log10(n + 1))

# Get top cell classes by count for visualization clarity
top_cell_classes <- hierarchy_heatmap %>%
  dplyr::group_by(cell_sub_class) %>%
  dplyr::summarize(total = sum(n)) %>%
  dplyr::arrange(desc(total)) %>%
  utils::head(40) %>%
  dplyr::pull(cell_sub_class)

p2 <- hierarchy_heatmap %>%
  dplyr::filter(cell_sub_class %in% top_cell_classes) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(cell_sub_class, log_n), 
                               y = reorder(cell_class, log_n), 
                               fill = log_n)) +
  ggplot2::geom_tile() +
  #ggplot2::facet_wrap(~flow, ncol = 1, scales = "free_y") +
  ggplot2::scale_fill_viridis_c(name = "log10(Count+1)") +
  ggplot2::labs(title = "Hierarchical Relationship Matrix",
                subtitle = " cell_class > cell_sub_class",
                x = "cell_sub_class", y = "cell_class") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
                 axis.text.y = ggplot2::element_text(size = 7),
                 panel.grid = ggplot2::element_blank())








