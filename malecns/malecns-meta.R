#############################
### GET MALECNS META DATA ###
#############################
source("banc/banc-startup.R")
library(malecns)

# Franken meta authority
franken.meta <- franken_meta() %>%
  dplyr::distinct(neuron_id, .keep_all = TRUE) %>%
  dplyr::select(neuron_id,
                side,
                region,
                hemilineage,
                nerve,
                flow,
                super_class,
                cell_class,
                cell_sub_class,
                cell_type,
                manc_cell_type = MANC_type,
                fafb_cell_type = FAFB_cell_type,
                hemibrain_cell_type = hemibrain_type,
                cell_function,
                cell_function_detailed,
                body_part_sensory,
                body_part_effector,
                peripheral_target_type)

franken_annot <- franken.meta %>%
  dplyr::filter(!is.na(cell_type),cell_type!="unknown") %>%
  dplyr::select(
    cell_type_fr = cell_type,
    manc_cell_type_fr = manc_cell_type,
    fafb_cell_type_fr = fafb_cell_type,
    hemibrain_cell_type_fr = hemibrain_cell_type,
    nerve_fr = nerve,
    super_class_fr = super_class,
    body_part_sensory_fr       = body_part_sensory,
    body_part_effector_fr      = body_part_effector,
    cell_function_fr           = cell_function,
    cell_function_detailed_fr  = cell_function_detailed,
    peripheral_target_type_fr  = peripheral_target_type,
    cell_class_fr  = cell_class,
    cell_sub_class_fr  = cell_sub_class
  ) %>%
  dplyr::distinct()

franken_nerve <- franken.meta %>%
  dplyr::filter(!is.na(cell_type),cell_type!="unknown") %>%
  dplyr::select(
    cell_type_fr = cell_type,
    side_fr = side,
    manc_cell_type_fr = manc_cell_type,
    fafb_cell_type_fr = fafb_cell_type,
    hemibrain_cell_type_fr = hemibrain_cell_type,
    body_part_sensory_fr  = body_part_sensory,
    body_part_effector_fr = body_part_effector,
    nerve_fr              = nerve
  ) %>%
  dplyr::filter(!is.na(nerve_fr)) %>%
  dplyr::distinct()

overwrite_from_ref <- function(df,
                               ref,
                               by,
                               cols,
                               suffix = "_fr",
                               clean_cell_type_lhs = TRUE,
                               only_if_na = FALSE,
                               distinct_keys = NULL,
                               drop_if_all_ref_na = FALSE 
) {
  lhs_cols <- names(by)
  rhs_cols <- unname(by)
  
  # 0) Strip any *_fr columns already present on the LHS (from previous calls)
  df_clean <- df
  fr_pattern <- paste0(suffix, "$")
  fr_cols_lhs <- grep(fr_pattern, names(df_clean), value = TRUE)
  if (length(fr_cols_lhs)) {
    df_clean <- df_clean[, setdiff(names(df_clean), fr_cols_lhs), drop = FALSE]
  }
  
  # 1) Enforce uniqueness on ref, and **drop rows with NA in any join key on RHS**
  if (is.null(distinct_keys)) {
    distinct_keys <- rhs_cols
  }
  ref_clean <- ref %>%
    # remove any Franken rows where a join key is NA
    dplyr::filter(
      !dplyr::if_any(
        dplyr::all_of(rhs_cols),
        ~ is.na(.x)
      )
    ) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(distinct_keys)), .keep_all = TRUE)
  
  # 2) Optionally clean *_cell_type on LHS (only if those cols exist)
  if (clean_cell_type_lhs) {
    ct_lhs <- lhs_cols[grepl("cell_type$", lhs_cols)]
    ct_lhs <- intersect(ct_lhs, names(df_clean))
    if (length(ct_lhs)) {
      df_clean <- df_clean %>%
        dplyr::mutate(
          dplyr::across(
            dplyr::all_of(ct_lhs),
            ~ sub(",.*$", "", .x)
          )
        )
    }
  }
  
  # 3) Join (no rows dropped; NA keys on LHS will now never find a match)
  df_joined <- df_clean %>%
    dplyr::left_join(ref_clean, by = by)
  
  # 4) Overwrite requested columns from *_fr
  for (nm in cols) {
    ref_nm <- paste0(nm, suffix)
    if (!ref_nm %in% names(df_joined)) next
    
    if (only_if_na) {
      idx <- is.na(df_joined[[nm]]) & !is.na(df_joined[[ref_nm]])
      df_joined[[nm]][idx] <- df_joined[[ref_nm]][idx]
    } else {
      df_joined[[nm]] <- dplyr::coalesce(df_joined[[ref_nm]], df_joined[[nm]])
    }
  }
  
  # 5) Optionally: drop rows where all Franken cols for `cols` are NA
  fr_cols_rhs <- grep(fr_pattern, names(df_joined), value = TRUE)
  fr_target_cols <- intersect(paste0(cols, suffix), fr_cols_rhs)
  
  if (drop_if_all_ref_na && length(fr_target_cols)) {
    df_joined <- df_joined %>%
      dplyr::filter(
        dplyr::if_any(
          dplyr::all_of(fr_target_cols),
          ~ !is.na(.x)
        )
      )
  }
  
  # 6) Drop all *_fr columns from the result (they are only temporary)
  if (length(fr_cols_rhs)) {
    df_joined <- df_joined[, setdiff(names(df_joined), fr_cols_rhs), drop = FALSE]
  }
  
  df_joined
}

# Get the meta data
choose_mcns_dataset("male-cns:v0.9")
mcns.find <- neuprintr::neuprint_search("Traced",field="status",dataset="male-cns:v0.9")
mcns.ids <- unique(mcns.find$bodyid)
mcns.meta.orig <- mcns_neuprint_meta(mcns.ids)
mcns.meta <- mcns.meta.orig %>%
  dplyr::mutate(side = ifelse(is.na(rootSide),somaSide,rootSide)) %>%
  dplyr::mutate(hemilineage = ifelse(is.na(itoleeHl),trumanHl,itoleeHl)) %>%
  dplyr::mutate(nerve = ifelse(is.na(entryNerve),exitNerve,entryNerve)) %>%
  dplyr::mutate(status = ifelse(is.na(statusLabel),status,statusLabel),
                mancBodyid = as.character(mancBodyid)) %>%
  dplyr::select(malecns_09_id = bodyid,
                fafb_cell_type = flywireType,
                hemibrain_cell_type = hemibrainType,
                manc_cell_type = mancType,
                manc_match = mancBodyid,
                instance = name,
                hemilineage,
                nerve,
                exitNerve,
                entryNerve,
                side,
                status = status,
                super_class = superclass,
                cell_type = type,
                cell_class = class,
                cell_sub_class = subclass,
                serial_motif= serialMotif,
                dimorphism,
                synonyms,
                cell_function = class,
                cell_function_detailed = subclass,
                peripheral_target_type = receptorType,
                neurotransmitter_predicted = predictedNt,
                neurotransmitter_score = predictedNtConfidence,
                optic_lobe_hex_1 = assignedOlHex1,
                optic_lobe_hex_2 = assignedOlHex2
  ) %>%
  dplyr::mutate(side = dplyr::case_when(
    side=="R" ~ "right",
    side=="L" ~ "left",
    side=="M" ~ "center",
    TRUE ~ NA
  )) %>%
  # region
  dplyr::mutate(region = dplyr::case_when(
    grepl("descending|ascending",super_class) ~ "neck_connective",
    grepl("visual_centrifugal",super_class) ~ "central_brain",
    (grepl("visual_projection|optic",super_class) & 
       !grepl("ocellar_interneuron",super_class)) ~ "optic_lobe",
    grepl("ocellar",super_class) ~ "central_brain",
    grepl("ocellar",cell_class) ~ "central_brain",
    grepl("vnc|VNC",super_class) ~ "ventral_nerve_cord",
    grepl("midbrain|central_brain|sez|^cb",super_class) ~ "central_brain",
    grepl("OL|optic|^ol",super_class) ~ "optic_lobe",
    # none get NA, as desired
    TRUE ~ "unknown"
  ) ) %>%
  # flow
  dplyr::mutate(flow = dplyr::case_when(
    grepl("motor|endocrine|efferent",super_class) ~ "efferent",
    grepl("sensory",super_class)&!grepl("efferent|intrinsic",super_class) ~ "afferent",
    super_class=="ascending_sensory"&grepl("^SA",cell_type) ~ "afferent",
    super_class=="ascending$"&!grepl("^SA",cell_type) ~ "intrinsic",
    cell_type %in% "HBeyelet" ~ "afferent",
    grepl("ocellar_retinula", cell_type) ~ "afferent",
    grepl("^R1-6$|^R1-R6$|^R7$|^R8$", cell_type) ~ "afferent",
    grepl("^ascending$",super_class) ~ "intrinsic",
    grepl("^descending$",super_class) ~ "intrinsic",
    grepl("ascending|descending",super_class) ~ "intrinsic",
    grepl("intrinsic",super_class) ~ "intrinsic",
    is.na(nerve) ~ "intrinsic",
    # none get NA, as desired
    TRUE ~ "unknown" 
  ) ) %>%
  # super_class
  dplyr::mutate(super_class = dplyr::case_when(
    grepl("afferent",flow) ~ "sensory",
    grepl("sensory",super_class) ~ "sensory",
    grepl("motor",super_class) ~ "motor",
    cell_type == "SNpp54" ~"sensory",
    grepl("^ISN$|^CB0991",cell_type) ~ "sensory",
    grepl("sensory_descending",super_class) ~ "sensory_descending",
    grepl("sensory_ascending",super_class) ~ "sensory_ascending",
    grepl("efferent_descending",super_class) ~ "central_brain_visceral_circulatory",
    grepl("efferent_ascending",super_class) ~ "ascending_visceral_circulatory",
    grepl("ascending_neuron",super_class) ~ "ascending",
    grepl("descending_neuron",super_class) ~ "descending",
    grepl("cb_motor",super_class) ~ "motor",
    grepl("cb_sensory",super_class) ~ "sensory",
    grepl("cb_efferent",super_class) ~ "central_brain_visceral_circulatory",
    grepl("ol_intrinsic",super_class) ~ "optic_lobe_intrinsic",
    grepl("cb_intrinsic",super_class) ~ "central_brain_intrinsic",
    grepl("ol_sensory",super_class) ~ "optic_lobe_sensory",
    grepl("vnc_efferent",cell_class) ~ "ventral_nerve_cord_visceral_circulatory",
    grepl("vnc_motor",cell_class) ~ "motor",
    grepl("vnc_sensory",cell_class) ~ "sensory",
    grepl("vnc_intrinsic",cell_class) ~ "ventral_nerve_cord_intrinsic",
    grepl("vnc_tbc",cell_class) ~ "ventral_nerve_cord_intrinsic",
    # not neurons
    cell_class %in% c('glia','putative_glia') ~ "glia",
    grepl("trachea",super_class) ~ "trachea",
    grepl("not_a_neuron|NOT_A_NEURON",super_class) ~ "not_a_neuron",
    # general
    grepl("ol_intrinsic",super_class) ~ "optic_lobe_intrinsic",
    grepl("cb_intrinsic",super_class) ~ "central_brain_intrinsic",
    # Neck
    super_class %in% c("descending") ~ "descending",
    super_class %in% c("ascending") ~ "ascending",
    super_class %in% c("sensory_ascending") ~ "sensory_ascending",
    region %in% c("neck_connective")&grepl('sensory_descending',super_class) ~ "sensory_descending",
    region %in% c("neck_connective")&grepl('sensory',super_class) ~ "sensory_ascending",
    region %in% c("neck_connective")&grepl('efferent|motor|endocrine',super_class) ~ "ascending_visceral_circulatory",
    super_class %in% c("sensory_descending") ~ "sensory_descending",
    super_class %in% c("visceral_circulatory_ascending") ~ "motor",
    super_class %in% c("efferent_ascending") ~ "efferent_ascending",
    super_class %in% c("efferent_descending") ~ "efferent_descending",
    # Afferent
    cell_type %in% "HBeyelet" ~ "sensory",
    grepl("ocellar_retinula", cell_type) ~ "sensory",
    grepl("^R1-6$|^R7$|^R8$", cell_type) ~ "sensory",
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
  # side
  dplyr::mutate(side = dplyr::case_when(
    grepl("right",side) ~ "right",
    grepl("left",side) ~ "left",
    grepl("center|midline|unpaired|both",side) ~ "center", # midline and unpaired unused
    # small number get NA, some sensory and others with missing class info
    TRUE ~ NA
  ) ) %>%
  # body_part_sensory
  dplyr::mutate(body_part_sensory = dplyr::case_when(
    # specific cell type fixes
    cell_type %in% "HBeyelet" ~ "eyelet",
    grepl("^R1-6$|^R1-R6|^R7$|^R8$",cell_type) ~ "retina",
    cell_type %in% "R1-6" ~ "retina",
    cell_type %in% c("R7","R8") ~ "retina",
    !grepl("sensory|visceral",super_class) ~ NA,
    grepl("^ISN$",cell_type) ~ "hemolymph_sensory",
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
    grepl("labellar",cell_sub_class) ~ "labellum",
    grepl("BM_Taste",cell_type) ~ "labellum",
    grepl("BM_Vib",cell_type) ~ "vibrissa",
    grepl("BM_vOcci_vPoOr",cell_type) ~ "postorbital_ventral",
    grepl("BM_Vt_PoOc",cell_type) ~ "postocellar",
    ## by prior annotation
    grepl("CB0991",cell_type) ~ "aorta",
    grepl("CB0991",fafb_cell_type) ~ "aorta",
    cell_type %in% "HBeyelet" ~ "eyelet",
    grepl("ocellar_retinula_cell",cell_type) ~ "ocellus",
    grepl("ocellar|ocelli",cell_function)&peripheral_target_type!="bristle" ~ "ocellus",
    grepl("ocellar|ocelli",cell_function_detailed)&peripheral_target_type!="bristle" ~ "ocellus",
    grepl("proboscis",cell_function)&flow=="efferent" ~ "proboscis",
    grepl("pharyn",cell_function)&flow=="efferent" ~ "pharynx",
    grepl("^abdomen$",cell_sub_class) ~ "abdomen",
    ## by nerves
    grepl("^AN$",entryNerve) ~ "antenna",
    grepl("^aPhN$",entryNerve) ~ "pharynx",
    grepl("^PhN$",entryNerve) ~ "pharynx",
    grepl("AbN2|AbN3|AbN4|AbNT",entryNerve) ~ "abdomen",
    grepl("DMetaN",entryNerve) ~ "haltere",
    grepl("DProN",entryNerve) ~ "front_leg",
    grepl("ProLN",entryNerve) ~ "front_leg",
    grepl("VProN",entryNerve) ~ "front_leg",
    grepl("MesoLN",entryNerve) ~ "middle_leg",
    grepl("MetaLN",entryNerve) ~ "hind_leg",
    grepl("PDMN",entryNerve) ~ "thorax",
    grepl("ADMN",entryNerve) ~ "wing",
    grepl("PrN",entryNerve) ~ "prosternal_organ",
    grepl("ProCN",entryNerve) ~ "prothoracic_chordotonal_organ",
    grepl("^ad$",cell_sub_class) ~ "abdomen",
    grepl("^am$",cell_sub_class) ~ "antenna",
    grepl("^fl$",cell_sub_class) ~ "front_leg",
    grepl("^ml$",cell_sub_class) ~ "middle_leg",
    grepl("^hl$",cell_sub_class) ~ "hind_leg",
    grepl("^hm$",cell_sub_class) ~ "haltere",
    grepl("^nm$",cell_sub_class) ~ "neck",
    grepl("^pm$",cell_sub_class) ~ "proboscis",
    grepl("^rm$",cell_sub_class) ~ "eye",
    grepl("^wm$",cell_sub_class) ~ "wing",
    grepl("^xm$",cell_sub_class) ~ "unknown",
    grepl("^aPhN$",exitNerve) ~ "pharynx",
    grepl("^PhN$",exitNerve) ~ "pharynx",
    grepl("^abdomen$",cell_sub_class) ~ "abdomen",
    grepl("^ad$",cell_sub_class) ~ "abdomen",
    grepl("^am$",cell_sub_class) ~ "antenna",
    grepl("^fl$",cell_sub_class) ~ "front_leg",
    grepl("^ml$",cell_sub_class) ~ "middle_leg",
    grepl("^hl$",cell_sub_class) ~ "hind_leg",
    grepl("^hm$",cell_sub_class) ~ "haltere",
    grepl("^nm$",cell_sub_class) ~ "neck",
    grepl("^pm$",cell_sub_class) ~ "proboscis",
    grepl("^rm$",cell_sub_class) ~ "eye",
    grepl("^wm$",cell_sub_class) ~ "wing",
    grepl("^xm$",cell_sub_class) ~ "unknown",
    grepl("^AN$",exitNerve) ~ "antenna",
    grepl("^CV$|CvN|CVC",exitNerve) ~ "neck",
    grepl("AbN2|AbN3|AbN4|AbNT",exitNerve) ~ "abdomen",
    grepl("DMetaN",exitNerve) ~ "haltere",
    grepl("DProN",exitNerve) ~ "prothorax",
    grepl("ProAN|ProLN",exitNerve) ~ "front_leg",
    grepl("VProN",exitNerve) ~ "front_leg",
    grepl("MesoLN",exitNerve) ~ "middle_leg",
    grepl("MetaLN",exitNerve) ~ "hind_leg",
    grepl("PDMN",exitNerve) ~ "wing",
    grepl("ADMN|MesoAN",exitNerve) ~ "wing",
    grepl("CvN|CVC|^CV$",exitNerve) ~ "neck",
    TRUE ~ NA
  ) ) %>%
  dplyr::mutate(nerve = dplyr::case_when(
    
    grepl("^ENS",cell_type) & side=="right" ~ "right_stomodeal_nerve",
    grepl("^ENS",cell_type) &side=="left" ~ "left_stomodeal_nerve",
    grepl("NCC",cell_type) ~ "stomodeal_nerve",
    
    grepl("NCC",nerve) & side=="right" ~ "right_corpus_cardiacum_nerve",
    grepl("NCC",nerve) &side=="left" ~ "left_corpus_cardiacum_nerve",
    grepl("NCC",nerve) ~ "corpus_cardiacum_nerve",
    
    grepl("^AN",nerve) & side=="right" ~ "right_antennal_nerve",
    grepl("^AN",nerve) &side=="left" ~ "left_antennal_nerve",
    grepl("^AN",nerve) ~ "antennal_nerve",
    
    grepl("MxLbN",nerve) & side=="right" ~ "right_maxillary-labial_nerve",
    grepl("MxLbN",nerve) & side=="left" ~ "left_maxillary-labial_nerve",
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
    grepl("^ON",nerve) & side=="left" ~ "left_optic_nerve",
    grepl("^ON",nerve) ~ "optic_nerve",
    
    grepl("CvN_R",nerve)  ~ "right_cervical_nerve",
    grepl("CvN_L",nerve) ~ "left_cervical_nerve",
    grepl("CvN",nerve) ~ "cervical_nerve",
    grepl("CV",nerve) & side=="right" ~ "right_cervical_nerve",
    grepl("CV",nerve) & side=="left" ~ "left_cervical_nerve",
    grepl("CV",nerve) & side=="left" ~ "cervical_nerve",
    
    grepl("DProN_R",nerve)  ~ "right_dorsal_prothoracic_nerve",
    grepl("DProN_L",nerve)  ~ "left_dorsal_prothoracic_nerve",
    grepl("DProN",nerve)&side=="right"  ~ "right_dorsal_prothoracic_nerve",
    grepl("DProN",nerve)&side=="left"  ~ "left_dorsal_prothoracic_nerve",
    grepl("DProN",nerve) ~ "dorsal_prothoracic_nerve",
    
    # Prothoracic leg nerve
    grepl("ProLN_R", nerve)              ~ "right_prothoracic_leg_nerve",
    grepl("ProLN_L", nerve)              ~ "left_prothoracic_leg_nerve",
    grepl("ProLN",   nerve) & side=="right" ~ "right_prothoracic_leg_nerve",
    grepl("ProLN",   nerve) & side=="left"  ~ "left_prothoracic_leg_nerve",
    grepl("ProLN",   nerve)              ~ "prothoracic_leg_nerve",
    
    # Prosternal nerve
    grepl("PrN_R", nerve)                ~ "right_prosternal_nerve",
    grepl("PrN_L", nerve)                ~ "left_prosternal_nerve",
    grepl("PrN",   nerve) & side=="right"   ~ "right_prosternal_nerve",
    grepl("PrN",   nerve) & side=="left"    ~ "left_prosternal_nerve",
    grepl("PrN",   nerve)                ~ "prosternal_nerve",
    
    # Prothoracic accessory nerve
    grepl("ProAN_R", nerve)              ~ "right_prothoracic_accessory_nerve",
    grepl("ProAN_L", nerve)              ~ "left_prothoracic_accessory_nerve",
    grepl("ProAN",   nerve) & side=="right" ~ "right_prothoracic_accessory_nerve",
    grepl("ProAN",   nerve) & side=="left"  ~ "left_prothoracic_accessory_nerve",
    grepl("ProAN",   nerve)              ~ "prothoracic_accessory_nerve",
    
    # Prothoracic chordotonal nerve
    grepl("ProCN_R", nerve)              ~ "right_prothoracic_chordotonal_nerve",
    grepl("ProCN_L", nerve)              ~ "left_prothoracic_chordotonal_nerve",
    grepl("ProCN",   nerve) & side=="right" ~ "right_prothoracic_chordotonal_nerve",
    grepl("ProCN",   nerve) & side=="left"  ~ "left_prothoracic_chordotonal_nerve",
    grepl("ProCN",   nerve)              ~ "prothoracic_chordotonal_nerve",
    
    # Ventral prothoracic nerve
    grepl("VProN_R", nerve)              ~ "right_ventral_prothoracic_nerve",
    grepl("VProN_L", nerve)              ~ "left_ventral_prothoracic_nerve",
    grepl("VProN",   nerve) & side=="right" ~ "right_ventral_prothoracic_nerve",
    grepl("VProN",   nerve) & side=="left"  ~ "left_ventral_prothoracic_nerve",
    grepl("VProN",   nerve)              ~ "ventral_prothoracic_nerve",
    
    # Anterior dorsal mesothoracic nerve
    grepl("ADMN_R", nerve)              ~ "right_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN_L", nerve)              ~ "left_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN",   nerve) & side=="right" ~ "right_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN",   nerve) & side=="left"  ~ "left_anterior_dorsal_mesothoracic_nerve",
    grepl("ADMN",   nerve)              ~ "anterior_dorsal_mesothoracic_nerve",
    
    # Posterior dorsal mesothoracic nerve
    grepl("PDMN_R", nerve)              ~ "right_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN_L", nerve)              ~ "left_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN",   nerve) & side=="right" ~ "right_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN",   nerve) & side=="left"  ~ "left_posterior_dorsal_mesothoracic_nerve",
    grepl("PDMN",   nerve)              ~ "posterior_dorsal_mesothoracic_nerve",
    
    # Mesothoracic leg nerve
    grepl("MesoLN_R", nerve)            ~ "right_mesothoracic_leg_nerve",
    grepl("MesoLN_L", nerve)            ~ "left_mesothoracic_leg_nerve",
    grepl("MesoLN",   nerve) & side=="right" ~ "right_mesothoracic_leg_nerve",
    grepl("MesoLN",   nerve) & side=="left"  ~ "left_mesothoracic_leg_nerve",
    grepl("MesoLN",   nerve)            ~ "mesothoracic_leg_nerve",
    
    # Mesothoracic accessory nerve
    grepl("MesoAN_R", nerve)            ~ "right_mesothoracic_accessory_nerve",
    grepl("MesoAN_L", nerve)            ~ "left_mesothoracic_accessory_nerve",
    grepl("MesoAN",   nerve) & side=="right" ~ "right_mesothoracic_accessory_nerve",
    grepl("MesoAN",   nerve) & side=="left"  ~ "left_mesothoracic_accessory_nerve",
    grepl("MesoAN",   nerve)            ~ "mesothoracic_accessory_nerve",
    
    # Dorsal metathoracic nerve
    grepl("DMetaN_R", nerve)            ~ "right_dorsal_metathoracic_nerve",
    grepl("DMetaN_L", nerve)            ~ "left_dorsal_metathoracic_nerve",
    grepl("DMetaN",   nerve) & side=="right" ~ "right_dorsal_metathoracic_nerve",
    grepl("DMetaN",   nerve) & side=="left"  ~ "left_dorsal_metathoracic_nerve",
    grepl("DMetaN",   nerve)            ~ "dorsal_metathoracic_nerve",
    
    # Metathoracic leg nerve
    grepl("MetaLN_R", nerve)            ~ "right_metathoracic_leg_nerve",
    grepl("MetaLN_L", nerve)            ~ "left_metathoracic_leg_nerve",
    grepl("MetaLN",   nerve) & side=="right" ~ "right_metathoracic_leg_nerve",
    grepl("MetaLN",   nerve) & side=="left"  ~ "left_metathoracic_leg_nerve",
    grepl("MetaLN",   nerve)            ~ "metathoracic_leg_nerve",
    
    # First abdominal nerve
    grepl("AbN1_R", nerve)              ~ "right_first_abdominal_nerve",
    grepl("AbN1_L", nerve)              ~ "left_first_abdominal_nerve",
    grepl("AbN1",   nerve) & side=="right" ~ "right_first_abdominal_nerve",
    grepl("AbN1",   nerve) & side=="left"  ~ "left_first_abdominal_nerve",
    grepl("AbN1",   nerve)              ~ "first_abdominal_nerve",
    
    # Second abdominal nerve
    grepl("AbN2_R", nerve)              ~ "right_second_abdominal_nerve",
    grepl("AbN2_L", nerve)              ~ "left_second_abdominal_nerve",
    grepl("AbN2",   nerve) & side=="right" ~ "right_second_abdominal_nerve",
    grepl("AbN2",   nerve) & side=="left"  ~ "left_second_abdominal_nerve",
    grepl("AbN2",   nerve)              ~ "second_abdominal_nerve",
    
    # Third abdominal nerve
    grepl("AbN3_R", nerve)              ~ "right_third_abdominal_nerve",
    grepl("AbN3_L", nerve)              ~ "left_third_abdominal_nerve",
    grepl("AbN3",   nerve) & side=="right" ~ "right_third_abdominal_nerve",
    grepl("AbN3",   nerve) & side=="left"  ~ "left_third_abdominal_nerve",
    grepl("AbN3",   nerve)              ~ "third_abdominal_nerve",
    
    # Fourth abdominal nerve
    grepl("AbN4_R", nerve)              ~ "right_fourth_abdominal_nerve",
    grepl("AbN4_L", nerve)              ~ "left_fourth_abdominal_nerve",
    grepl("AbN4",   nerve) & side=="right" ~ "right_fourth_abdominal_nerve",
    grepl("AbN4",   nerve) & side=="left"  ~ "left_fourth_abdominal_nerve",
    grepl("AbN4",   nerve)              ~ "fourth_abdominal_nerve",
    
    # Abdominal nerve trunk
    grepl("AbNT_R", nerve)              ~ "right_abdominal_nerve_trunk",
    grepl("AbNT_L", nerve)              ~ "left_abdominal_nerve_trunk",
    grepl("AbNT",   nerve) & side=="right" ~ "right_abdominal_nerve_trunk",
    grepl("AbNT",   nerve) & side=="left"  ~ "left_abdominal_nerve_trunk",
    grepl("AbNT",   nerve)              ~ "abdominal_nerve_trunk",
    
    # Abdominal nerve "other"
    grepl("AbNX_R", nerve)              ~ "right_abdominal_nerve_other",
    grepl("AbNX_L", nerve)              ~ "left_abdominal_nerve_other",
    grepl("AbNX",   nerve) & side=="right" ~ "right_abdominal_nerve_other",
    grepl("AbNX",   nerve) & side=="left"  ~ "left_abdominal_nerve_other",
    grepl("AbNX",   nerve)              ~ "abdominal_nerve_other",
    
    TRUE ~ nerve
  ) ) %>%
  # Fill missing nerve from franken using side + type matches (first type before comma)
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("fafb_cell_type"      = "cell_type_fr",
             "side"                = "side_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("hemibrain_cell_type" = "cell_type_fr",
             "side"                = "side_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  # body part effector
  dplyr::mutate(body_part_effector = dplyr::case_when(
    !grepl("efferent",flow) ~ NA,
    grepl("^aPhN$",exitNerve) ~ "pharynx",
    grepl("^PhN$",exitNerve) ~ "pharynx",
    grepl("^abdomen$",cell_sub_class) ~ "abdomen",
    grepl("^ad$",cell_sub_class) ~ "abdomen",
    grepl("^am$",cell_sub_class) ~ "antenna",
    grepl("^fl$",cell_sub_class) ~ "front_leg",
    grepl("^ml$",cell_sub_class) ~ "middle_leg",
    grepl("^hl$",cell_sub_class) ~ "hind_leg",
    grepl("^hm$",cell_sub_class) ~ "haltere",
    grepl("^nm$",cell_sub_class) ~ "neck",
    grepl("^pm$",cell_sub_class) ~ "proboscis",
    grepl("^rm$",cell_sub_class) ~ "eye",
    grepl("^wm$",cell_sub_class) ~ "wing",
    grepl("^xm$",cell_sub_class) ~ "unknown",
    grepl("^AN$",exitNerve) ~ "antenna",
    grepl("^CV$|CvN|CVC",exitNerve) ~ "neck",
    grepl("AbN2|AbN3|AbN4|AbNT",exitNerve) ~ "abdomen",
    grepl("DMetaN",exitNerve) ~ "haltere",
    grepl("DProN",exitNerve) ~ "prothorax",
    grepl("ProAN|ProLN",exitNerve) ~ "front_leg",
    grepl("VProN",exitNerve) ~ "front_leg",
    grepl("MesoLN",exitNerve) ~ "middle_leg",
    grepl("MetaLN",exitNerve) ~ "hind_leg",
    grepl("PDMN",exitNerve) ~ "wing",
    grepl("ADMN|MesoAN",exitNerve) ~ "wing",
    grepl("CvN|CVC|^CV$",exitNerve) ~ "neck",
    TRUE ~ NA
  ) ) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("body_part_sensory","body_part_effector")) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("fafb_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
    cols = c("body_part_sensory","body_part_effector")) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("body_part_sensory","body_part_effector")) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("manc_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
    cols = c("body_part_sensory","body_part_effector"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("body_part_sensory","body_part_effector"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("body_part_sensory","body_part_effector"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("body_part_sensory","body_part_effector"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("manc_cell_type" = "manc_cell_type_fr",
             "side" = "side_fr",
             "body_part_sensory" = "body_part_sensory_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("manc_cell_type" = "manc_cell_type_fr",
             "side" = "side_fr",
             "body_part_effector" = "body_part_effector_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("manc_cell_type" = "cell_type_fr",
             "side" = "side_fr",
             "body_part_sensory" = "body_part_sensory_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("manc_cell_type" = "cell_type_fr",
             "side" = "side_fr",
             "body_part_effector" = "body_part_effector_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("cell_type" = "cell_type_fr",
             "side" = "side_fr",
             "body_part_sensory" = "body_part_sensory_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  overwrite_from_ref(
    ref  = franken_nerve,
    by   = c("cell_type" = "cell_type_fr",
             "side" = "side_fr",
             "body_part_effector" = "body_part_effector_fr"),
    cols = "nerve",
    suffix = "_fr",
    clean_cell_type_lhs = TRUE,
    only_if_na = TRUE
  ) %>%
  # peripheral_target_type
  dplyr::mutate(peripheral_target_type = dplyr::case_when(
    # pain
    grepl("SNch01|SNxx19|SNxx20|SNxx21",manc_cell_type) ~ "multidendritic", # class_iv
    grepl("nociception",cell_function) ~ "multidendritic", 
    grepl("pharynx|cibarium",body_part_sensory)&(grepl("mechanosensory|tactile|pumping|ciberial",cell_function)|grepl("mechanosensory|tactile|pumping|ciberial",cell_function_detailed)) ~ "multidendritic",
    # taste
    grepl("pharyngeal sensillum",cell_sub_class) ~ "internal_taste_sensillum",
    grepl("taste peg",cell_sub_class) ~ "taste_peg",
    grepl("taste_peg",cell_sub_class) ~ "taste_peg",
    grepl("taste bristle",cell_sub_class) ~ "taste_bristle",
    grepl("labellar bristle",cell_sub_class) ~ "taste_bristle",
    grepl("taste_bristle",cell_sub_class) ~ "taste_bristle",
    grepl("taste_peg",cell_function_detailed) ~ "taste_peg",
    grepl("hair plate",cell_sub_class) ~ "hair_plate",
    grepl("eye_bristle",cell_sub_class) ~ "bristle",
    grepl("head_bristle",cell_sub_class) ~ "bristle",
    grepl("bristle",cell_function_detailed) ~ "bristle",
    grepl("bristle",cell_function_detailed) ~ "bristle",
    grepl("mechanosensory_bristle",cell_sub_class) ~ "bristle",
    grepl("mechanosensory bristle",cell_sub_class) ~ "bristle",
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
    grepl("leg bristle",cell_function_detailed) ~ "bristle",
    grepl("mechanosensory bristle",cell_function_detailed) ~ "bristle",
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
    grepl("campaniform sensilla",cell_function_detailed) ~ "campaniform_sensillum",
    grepl("chordotonal",cell_function_detailed) ~ "chordotonal_organ",
    grepl("chordotonal",cell_function_detailed) ~ "chordotonal_organ",
    grepl("hair plate",cell_function_detailed) ~ "hair_plate",
    grepl("putative sweet taste bristle",cell_function_detailed) ~ "taste_peg",
    grepl("sweet taste bristle",cell_function_detailed) ~ "taste_peg",
    grepl("taste bristle",cell_function_detailed) ~ "taste_peg",
    grepl("strand receptor",cell_function_detailed) ~ "strand",
    grepl("taste bristle",cell_function_detailed) ~ "taste_peg",
    # visual
    grepl("ocellar retinula",fafb_cell_type) ~ "photoreceptor",
    grepl("ocellar_retinula",cell_type) ~ "photoreceptor",
    grepl("^R1-6|^R1-R6|^R7|^R8",cell_type) ~ "photoreceptor",
    grepl("^R1-6|^R1-R6|^R7|^R8",cell_type) ~ "photoreceptor",
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
    TRUE ~ NA_character_
  ) ) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("peripheral_target_type")) %>%
  overwrite_from_ref(ref=franken_annot,
                    by   = c("fafb_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                    cols = c("peripheral_target_type")) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("peripheral_target_type"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                    by   = c("hemibrain_cell_type"= "cell_type_fr","nerve"= "nerve_fr"),
                    cols = c("peripheral_target_type"),
                    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("peripheral_target_type"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("peripheral_target_type"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr","nerve"= "nerve_fr"),
                     cols = c("peripheral_target_type"),
                     only_if_na = FALSE) %>%
  # cell_function_detailed
  dplyr::mutate(cell_function_detailed = dplyr::case_when(
    #!is.na(cell_function_detailed) ~ cell_function_detailed,
    # gustatory
    grepl("sugar_or_water",cell_function)|grepl("sugar_or_water",cell_function_detailed) ~ "sugar_or_water",
    grepl("putative sweet taste bristle",cell_function_detailed) ~ "sugar",
    grepl("sweet taste bristle",cell_function_detailed) ~ "sugar",
    grepl("sweet taste bristle",cell_function_detailed) ~ "sugar",
    grepl("sugar",cell_function)|grepl("sugar",cell_function_detailed) ~ "sugar",
    grepl("sweet",cell_function)|grepl("sweet",cell_function_detailed) ~ "sugar",
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
    grepl("JO-C/D/E",cell_type) ~ "position",
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
    # for remainder, wipe cell_function_detailed from visceral_circulatory neurons
    grepl("visceral_circulatory", super_class) ~ NA_character_, 
    # remainder
    grepl("WG1",cell_type) ~ "contact_pheromone, Ir52a",
    grepl("LgLG6",cell_type) ~ "contact_pheromone, ppk23a",
    grepl("LgLG5|LgLG8",cell_type) ~ "contact_pheromone, ppk23, ppk25",
    grepl("^water$", cell_function_detailed)&grepl("sensory",super_class) ~ "water, ppk28", 
    grepl("^sweet$", cell_function_detailed)&grepl("sensory",super_class) ~ "sugar, Gr33a", 
    grepl("^Gr66a$", cell_function_detailed) ~ "bitter, Gr66a", 
    grepl("^ppk23, ppk25", cell_function_detailed) ~ "contact_pheromone, ppk23, ppk25", 
    grepl("^ppk23", cell_function_detailed) ~ "contact_pheromone, ppk23", 
    grepl("bristle|chordotonal|campaniform|hair|^BI|^BR|neck|notum|sensillum|^ad$|haltere",
          cell_function_detailed) ~ NA_character_, 
    cell_function_detailed == "NA" ~ NA_character_,
    is.na(cell_function_detailed) ~ cell_function_detailed,
    grepl("sensory",super_class) & is.na(cell_function_detailed) ~ NA_character_,
    TRUE ~ NA_character_
  ) ) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("cell_type"= "cell_type_fr"),
                     cols = c("cell_function_detailed")) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("fafb_cell_type"= "cell_type_fr"),
    cols = c("cell_function_detailed")) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "cell_type_fr"),
                     cols = c("cell_function_detailed"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("hemibrain_cell_type"= "cell_type_fr"),
    cols = c("cell_function_detailed"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr"),
                     cols = c("cell_function_detailed")) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr"),
                     cols = c("cell_function_detailed")) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr"),
                     cols = c("cell_function_detailed")) %>%
  # cell function
  dplyr::mutate(cell_function = dplyr::case_when(
    grepl("SNch01|SAxx02|SNxx24|SNxx25|SNxx27|SNxx29",cell_type) ~ "putative_nociception",
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
    cell_type %in% "R1-R6" ~ "visual_achromatic",
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
    grepl("Kenyon_Cell",cell_type) ~ "associative_memory_storage",
    grepl("CX",cell_type) ~ "navigation",
    grepl("DAN",cell_type) ~ "associative_learning",
    grepl("ALLN|ALPN|ALIN|ALON",cell_type) ~ "olfactory_processing",
    grepl("sensory",super_class)&is.na(cell_function) ~ NA_character_,
    TRUE ~ cell_function
  ) ) %>%
  dplyr::mutate(cell_function_detailed = dplyr::case_when(
    cell_function == cell_function_detailed ~ NA_character_,
    TRUE ~ cell_function_detailed
  )) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("cell_type"= "cell_type_fr"),
                     cols = c("cell_function"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("fafb_cell_type"= "cell_type_fr"),
    cols = c("cell_function"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "cell_type_fr"),
                     cols = c("cell_function"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("hemibrain_cell_type"= "cell_type_fr"),
    cols = c("cell_function"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr"),
                     cols = c("cell_function"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr"),
                     cols = c("cell_function"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr"),
                     cols = c("cell_function"),
                     only_if_na = FALSE) %>%
  # cell_class
  dplyr::mutate(cell_class = dplyr::case_when(
    grepl("ol_bilateral",cell_type) ~ "optic_lobe_bilateral",
    # 'endocrine'
    grepl("^ISN$",cell_type) ~ "nutrient_sensory_neuron",
    grepl("^PSI$",cell_type) ~ "peripheral_intrinsic_neuron",
    grepl("^ascending_visceral_circulatory",super_class) ~  "ventral_nerve_cord_neurosecretory_cell",
    grepl("^m_NSC_unknown$",cell_type) ~  "pars_intercerebralis_neurosecretory_cell",
    grepl("^DNg28$",cell_type) ~  "subesophageal_zone_neurosecretory_cell",
    grepl("^m-NSC$|medial_NSC|medial_neurosecretory_cell",cell_sub_class) ~ "pars_intercerebralis_neurosecretory_cell",
    grepl("^l-NSC$|lateral_NSC|lateral_neurosecretory_cell",cell_sub_class) ~ "pars_lateralis_neurosecretory_cell",
    grepl("SEZ-NSC|subesophageal_zone_neurosecretory_cell",cell_sub_class) ~ "subesophageal_zone_neurosecretory_cell",
    grepl("^PI$",fafb_cell_type) ~  "pars_intercerebralis_efferent_neuron",
    grepl("endocrine",cell_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector)&grepl("ventral",region) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$|ANm_endocrine",cell_class)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^EN$",cell_sub_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector)&grepl("ventral",region) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("endocrine",cell_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^EN$",cell_sub_class)&grepl("ventral_nerve_cord",region)&grepl("hemal",body_part_effector) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$|ANm_endocrine",cell_class)&!grepl("hemal",body_part_effector)&grepl("ventral",region) ~ "ventral_nerve_cord_neurosecretory_cell",
    grepl("^ad_encodrine_neuron$",cell_sub_class)&grepl("ventral",region) ~ "ventral_nerve_cord_neurosecretory_cell",
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
    grepl("SLPpl1",hemilineage)&grepl("clock|circadian",cell_class)~ "circadian_neuron",
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
    grepl("^Kenyon_Cell$",cell_class) ~ "kenyon_cell",
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
    grepl("independent leg",serial_motif) ~ "single_leg_neuromere",
    grepl("sequential",serial_motif) ~ "sequential_leg_neuromeres",
    grepl("centrifugal",serial_motif) ~ "ventral_nerve_cord_centrifugal",
    grepl("dorsal",serial_motif) ~ "ventral_nerve_cord_dorsal",
    grepl("centripetal",serial_motif) ~ "ventral_nerve_cord_centripetal",
    grepl("convergent",serial_motif) ~ "ventral_nerve_cord_serially_convergent",
    # bad
    grepl("^Interneuron_TBD$",cell_class) ~ NA_character_,
    grepl("^TBD$",cell_class) ~ NA_character_,
    grepl("^TBD$",cell_class) ~ NA_character_,
    TRUE ~ cell_class
  )) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("fafb_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
    cols = c("cell_class"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_class"),
                     only_if_na = TRUE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("hemibrain_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
    cols = c("cell_class"),
    only_if_na = TRUE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_class"),
                     only_if_na = FALSE) %>%
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
    grepl("CX|central_complex",cell_class)  & grepl("hDelta",cell_type) ~ "hDelta",
    grepl("CX|central_complex",cell_class)  & grepl("vDelta",cell_type) ~ "vDelta",
    grepl("CX|central_complex",cell_class)  & grepl("^FB|^CB.FB|^SA1|^SA2|^SA3",cell_type) ~ "FB_tangential",
    grepl("CX|central_complex",cell_class)  & grepl("^ExR",cell_type) ~ "ExR",
    grepl("CX|central_complex",cell_class)  & grepl("^FC",cell_type) ~ "FC",
    grepl("CX|central_complex",cell_class)  & grepl("^FS",cell_type) ~ "FS",
    grepl("CX|central_complex",cell_class)  & grepl("^PFL",cell_type) ~ "PFL",
    grepl("CX|central_complex",cell_class)  & grepl("^PFN",cell_type) ~ "PFN",
    grepl("CX|central_complex",cell_class)  & grepl("^ER",cell_type) ~ "EB_input",
    grepl("CX|central_complex",cell_class)  & grepl("^FR1|^FR2",cell_type) ~ "FR",
    grepl("CX|central_complex",cell_class)  & grepl("^PFR$",cell_type) ~ "PFR",
    grepl("CX|central_complex",cell_class)  & grepl("^PFGs$",cell_type) ~ "PFG",
    grepl("CX|central_complex",cell_class)  & grepl("^PEG",cell_type) ~ "PEG",
    grepl("CX|central_complex",cell_class)  & grepl("^PEN_",cell_type) ~ "PEN",
    grepl("CX|central_complex",cell_class)  & grepl("^GLNO$|^LCNOp$|^LCNOpm$|^LNO1$|^LNO1$|^LNO2$",cell_type) ~ "NO_input",
    grepl("CX|central_complex",cell_class)  & grepl("^IbSpsP|^LPsP|^SpsP",cell_type) ~ "PB_input",
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
    grepl("SLPpl1",hemilineage)&grepl("clock|circadian",cell_class)~ "s-CPDN",
    cell_type %in% "HBeyelet" ~ "hofbauer_buchner_eyelet_neuron",
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
    !is.na(body_part_effector)&grepl("motor",super_class) ~ paste0(body_part_effector,"_motor_neuron"),
    !is.na(body_part_effector)&grepl("visceral",super_class) ~ paste0(body_part_effector,"_endocrine"),
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
    # ventral nerve cord
    cell_sub_class=="BA" ~ "bilateral_connections_ascending", 
    cell_sub_class=="BA" ~ "contralateral_connections_ascending",
    cell_sub_class=="XA" ~ "no_connections_ascending", 
    cell_sub_class=="BR" ~ "ventral_nerve_cord_bilateral_restricted", 
    cell_sub_class=="IR" ~ "ventral_nerve_cord_ipsilateral_restricted", 
    cell_sub_class=="CR" ~ "ventral_nerve_cord_contralateral_restricted", 
    cell_sub_class=="BI" ~ "ventral_nerve_cord_bilateral_interconnecting", 
    cell_sub_class=="II" ~ "ventral_nerve_cord_ipsilateral_interconnecting", 
    cell_sub_class=="CI" ~ "ventral_nerve_cord_contralateral_interconnecting",
    # other
    grepl("^water_PN",cell_sub_class) ~ "BiT",
    grepl("^mAL$",cell_class)|grepl("^mAL$",cell_sub_class) ~ "mAL",
    cell_sub_class=="NANA" ~ NA_character_,
    TRUE ~ cell_sub_class
  )) %>%
  overwrite_from_ref(ref=franken_annot,
                     by = c("cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_sub_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("fafb_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
    cols = c("cell_sub_class"),
    only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("hemibrain_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
    cols = c("cell_sub_class"),
    only_if_na = TRUE) %>%
  overwrite_from_ref(ref=franken_annot,
    by   = c("manc_cell_type"= "cell_type_fr","super_class"= "super_class_fr"),
    cols = c("cell_sub_class"),
    only_if_na = TRUE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("fafb_cell_type"= "fafb_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_sub_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("manc_cell_type"= "manc_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_sub_class"),
                     only_if_na = FALSE) %>%
  overwrite_from_ref(ref=franken_annot,
                     by   = c("hemibrain_cell_type"= "hemibrain_cell_type_fr","super_class"= "super_class_fr"),
                     cols = c("cell_sub_class"),
                     only_if_na = FALSE) %>%
  dplyr::select(-exitNerve,
                -entryNerve)

# Other meta sources
banc.meta <- banctable_query() %>%
  dplyr::distinct(cell_type, nerve, hemilineage, .keep_all = TRUE) %>%
  dplyr::select(region,
                hemilineage,
                nerve,
                flow,
                super_class,
                cell_class,
                cell_sub_class,
                cell_type,
                fafb_cell_type,
                hemibrain_cell_type,
                manc_cell_type,
                cell_function,
                cell_function_detailed,
                body_part_sensory,
                body_part_effector,
                peripheral_target_type) 

# combine
prefer_cols <- c(
  "hemilineage","nerve","flow","super_class","cell_class","cell_sub_class",
  "cell_function","cell_function_detailed","body_part_sensory","body_part_effector","peripheral_target_type"
)

# Treat "" like NA for coalescing
na_empty <- function(x) dplyr::na_if(x, "")

# Prefer RHS non-empty; otherwise keep LHS
coalesce_preferring_rhs <- function(lhs, rhs) {
  dplyr::coalesce(na_empty(rhs), na_empty(lhs))
}

# Safe left-join:
# - keep only join keys + prefer_cols from RHS
# - ensure RHS is unique on the join keys (prevents explosive joins)
# - coalesce preferred columns with RHS priority
do_one_join <- function(df, right, by) {
  left_keys  <- names(by)
  right_keys <- unname(by)
  
  right_min <- right %>%
    dplyr::select(dplyr::any_of(c(right_keys, prefer_cols))) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(right_keys)), .keep_all = TRUE)
  
  if (anyDuplicated(right_min[, right_keys, drop = FALSE])) {
    stop("Right table is not unique on keys: ", paste(right_keys, collapse = ", "))
  }
  
  out <- dplyr::left_join(df, right_min, by = by, suffix = c("", ".j"))
  
  # Resolve preferred columns
  for (nm in prefer_cols) {
    jnm <- paste0(nm, ".j")
    if (jnm %in% names(out)) {
      out[[nm]] <- coalesce_preferring_rhs(out[[nm]], out[[jnm]])
      out[[jnm]] <- NULL
    }
  }
  out
}

# --- Slim source tables (keep only what we need + ensure uniqueness) ----------
franken.slim <- franken.meta %>%
  dplyr::select(dplyr::any_of(c(
    "neuron_id","cell_type","hemilineage","nerve","super_class", prefer_cols
  ))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x)) %>%
  dplyr::distinct()

banc.slim <- banc.meta %>%
  dplyr::select(dplyr::any_of(c(
    "cell_type","fafb_cell_type","hemibrain_cell_type","manc_cell_type",
    "hemilineage","nerve","super_class", prefer_cols
  ))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x)) %>%
  dplyr::distinct()

mcns.enriched <- mcns.meta
gc()

# --- JOIN ORDER: lowest -> highest priority ----------------------------------
# BANC by cell_type = cell_type AND super_class = super_class
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(cell_type), !is.na(super_class)),
  by = c("cell_type" = "cell_type", "super_class" = "super_class")
); gc()

# BANC by manc_cell_type = manc_cell_type AND nerve = nerve
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(manc_cell_type), !is.na(nerve)),
  by = c("manc_cell_type" = "manc_cell_type", "nerve" = "nerve")
); gc()

# BANC by manc_cell_type = manc_cell_type AND hemilineage = hemilineage
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(manc_cell_type), !is.na(hemilineage)),
  by = c("manc_cell_type" = "manc_cell_type", "hemilineage" = "hemilineage")
); gc()

# BANC by hemibrain_cell_type = hemibrain_cell_type
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(hemibrain_cell_type)),
  by = c("hemibrain_cell_type" = "hemibrain_cell_type")
); gc()

# BANC by fafb_cell_type = fafb_cell_type
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(fafb_cell_type)),
  by = c("fafb_cell_type" = "fafb_cell_type")
); gc()

# BANC by fafb_cell_type = fafb_cell_type AND hemilineage = hemilineage
mcns.enriched <- do_one_join(
  mcns.enriched, banc.slim %>% dplyr::filter(!is.na(fafb_cell_type),!is.na(hemilineage)),
  by = c("fafb_cell_type" = "fafb_cell_type", "hemilineage" = "hemilineage")
); gc()

# FRANKEN by cell_type = cell_type AND super_class = super_class
mcns.enriched <- do_one_join(
  mcns.enriched, franken.slim %>% dplyr::filter(!is.na(cell_type),!is.na(super_class)),
  by = c("cell_type" = "cell_type", "super_class" = "super_class")
); gc()

# FRANKEN by cell_type = fafb_cell_type AND hemilineage = hemilineage
mcns.enriched <- do_one_join(
  mcns.enriched, franken.slim %>% dplyr::filter(!is.na(cell_type),!is.na(hemilineage)),
  by = c("fafb_cell_type" = "cell_type", "hemilineage" = "hemilineage")
); gc()

# FRANKEN by cell_type = manc_cell_type AND nerve = nerve
mcns.enriched <- do_one_join(
  mcns.enriched, franken.slim %>% dplyr::filter(!is.na(cell_type),!is.na(nerve)),
  by = c("manc_cell_type" = "cell_type", "nerve" = "nerve")
); gc()

# FRANKEN by cell_type = manc_cell_type AND hemilineage = hemilineage
mcns.enriched <- do_one_join(
  mcns.enriched, franken.slim %>% dplyr::filter(!is.na(cell_type),!is.na(hemilineage)),
  by = c("manc_cell_type" = "cell_type", "hemilineage" = "hemilineage")
); gc()

# FRANKEN by neuron_id = manc_match   (HIGHEST PRIORITY)
mcns.enriched <- do_one_join(
  mcns.enriched, franken.slim %>% 
    dplyr::rename(neuron_id_fr = neuron_id) %>% 
    dplyr::filter(!is.na(neuron_id_fr),
                  !neuron_id_fr%in%c("17796","11540","22631","23505","27083","17530")),
  by = c("manc_match" = "neuron_id_fr")
); gc()

# --- Canonical column order (same headers you’ve been using) ------------------
final_cols <- c(
  # IDs you have in malecns + links
  "malecns_09_id","manc_match","fafb_cell_type","hemibrain_cell_type","manc_cell_type",
  # anatomy & function (resolved)
  "region","side","hemilineage","nerve","flow","super_class","cell_class","cell_sub_class",
  "cell_type","cell_function","cell_function_detailed",
  "body_part_sensory","body_part_effector",
  # NT fields if present
  "neurotransmitter_predicted","neurotransmitter_score",
  # misc you kept earlier
  "instance","status","serial_motif","dimorphism","optic_lobe_hex_1","optic_lobe_hex_2","synonyms",
  "peripheral_target_type"
)

# Create any missing columns as NA to stabilize select()
missing_cols <- setdiff(final_cols, names(mcns.enriched))
if (length(missing_cols)) {
  mcns.enriched[missing_cols] <- NA
}

mcns.enriched <- mcns.enriched %>%
  dplyr::select(dplyr::any_of(final_cols)) %>%
  dplyr::distinct()

# Save
readr::write_csv(mcns.enriched, file.path(banc.meta.save.path,"malecns_09_meta.csv"))

# Announce
message("##### BANCpipeline: malecns meta updated #####")
message(sprintf("##### we have meta for : %s neurons", nrow(mcns.enriched)))

# mcns.enriched[is.na(mcns.enriched)] <- ""
#banctable_append_rows(df = mcns.enriched, table = "malecns", base = "cns_meta", bigdata = TRUE)

# mcns.enriched <- mcns.meta

# Build per-dataset minimal tables
franken_sens <- franken.meta %>%
  dplyr::transmute(
    dataset = "franken",
    super_class,
    cell_function_detailed,
    body_part_sensory
  )

malecns_sens <- mcns.enriched %>%
  dplyr::transmute(
    dataset = "maleCNS",
    super_class,
    cell_function_detailed,
    body_part_sensory
  )

# Combine, keep only "sensory" super_classes, count per dataset
sensory_long <- dplyr::bind_rows(franken_sens, malecns_sens) %>%
  dplyr::filter(grepl("sensory", super_class)) %>%
  dplyr::group_by(dataset, cell_function_detailed, body_part_sensory) %>%
  dplyr::summarise(n_neurons = dplyr::n(), .groups = "drop")

# Wide table + diff column, sorted by largest diff first
sensory_wide <- sensory_long %>%
  tidyr::pivot_wider(
    id_cols    = c(cell_function_detailed, body_part_sensory),
    names_from = dataset,
    values_from = n_neurons,
    names_glue  = "n_neurons_{dataset}",
    values_fill = list(n_neurons = 0)
  ) %>%
  dplyr::mutate(
    diff = n_neurons_maleCNS - n_neurons_franken
  ) %>%
  dplyr::arrange(dplyr::desc(diff))

as.data.frame(sensory_wide)


