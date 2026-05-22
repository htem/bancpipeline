# script to update metadata based only off cell_type entry, when the cell_type
#  matches one in franken_meta
# do for all cell_class contains central_complex
# do for all neurons with cell_type_source contains lee_lab
# cell_info_id tag is cell_type

source("banc/curation/banc_curation_init.R")

# get banc_meta
banc_meta <- banctable_query()

# get franken meta
franken_meta <- franken_meta()

cell_info <- banc_cave_query("cell_info", live = 2)


# replace semi-colon separated list with commas
franken_meta <- franken_meta %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, ";\\s*", ",")))

# get list of all cell_types 
franken_meta_ct <- franken_meta %>%
  distinct(cell_type)

# prep franken_meta for joining
# manc, for join based on cell_type
franken_meta_manc_joinCT <- franken_meta %>%
  filter(is_not_empty_na(manc_id)) %>%
  select(flow, super_class, cell_class, cell_sub_class, cell_type, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter,
         neuropeptide, sexually_dimorphic, dataset) %>%
  mutate(neurotransmitter_verified = neurotransmitter) %>%
  mutate(neuropeptide_verified = neuropeptide) %>%
  select(-neurotransmitter, -neuropeptide) %>%
  distinct(cell_type, .keep_all = TRUE)

# fafb, for join based on cell_type
franken_meta_fafb_joinCT <- franken_meta %>%
  filter(is_not_empty_na(fafb_id)) %>%
  select(flow, super_class, cell_class, cell_sub_class, cell_type, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter,
         neuropeptide, sexually_dimorphic, dataset) %>%
  mutate(neurotransmitter_verified = neurotransmitter) %>%
  mutate(neuropeptide_verified = neuropeptide) %>%
  select(-neurotransmitter, -neuropeptide) %>%
  distinct(cell_type, .keep_all = TRUE)

# combine for cell_type join
franken_meta_both_joinCT <- bind_rows(franken_meta_manc_joinCT, 
                                      franken_meta_fafb_joinCT)

# remove duplicates, prioritizing FAFB for descending neurons and MANC for 
#  ascending neurons
franken_meta_both_joinCT <- franken_meta_both_joinCT %>%
  group_by(cell_type) %>%
  arrange(
    cell_type,
    case_when(
      grepl("descending", super_class, ignore.case = TRUE) & grepl("FAFB", dataset, ignore.case = TRUE) ~ 1,
      grepl("ascending", super_class, ignore.case = TRUE) & grepl("MANC", dataset, ignore.case = TRUE) ~ 1,
      TRUE ~ 2
    )
  ) %>%
  slice(1) %>%
  ungroup()

###

# central_complex neurons - cell type in franken_meta

central_complex_ct_fm <- franken_meta_both_joinCT %>%
  filter(grepl("central_complex", cell_class))

banc_meta_central_complex <- banc_meta %>%
  filter(grepl("central_complex", cell_class) | (cell_type %in% central_complex_ct_fm$cell_type))

# modify banc_meta_central_complex for joining to franken_meta
banc_meta_central_complex_mod <- banc_meta_central_complex %>%
  select(`_id`, cell_type)

# modify franken_meta for join (remove hemilineage because it differs across 
#  members of same cell type)
franken_meta_both_joinCT_cxMod <- franken_meta_both_joinCT %>%
  select(-hemilineage, -dataset)

# join and modify for updating banc_meta
banc_meta_central_complex_join <- banc_meta_central_complex_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_cxMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_cxMod, by = "cell_type") %>%
  select(-cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_central_complex_join, 
                      # table='banc_meta',base='banc_meta')


###

# lee_lab cell_type_source
# and fafb_match gives different cell_type (match is wrong)

banc_meta_lee <- banc_meta %>%
  filter(grepl("lee_lab", cell_type_source)) %>%
  filter(cell_type != fafb_cell_type)
  
# modify banc_meta_lee for joining to franken_meta
banc_meta_lee_mod <- banc_meta_lee %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_leeMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
banc_meta_lee_join <- banc_meta_lee_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_leeMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_leeMod, by = "cell_type") %>%
  select(-cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_lee_join, 
                      # table='banc_meta',base='banc_meta')

###

# cell_type came from cell_info
# fafb or manc match give different cell type
# exclude ascending and descending neurons
# exclude agrawal_lab from modification

banc_meta_cell_info <- banc_meta %>%
  filter(is_not_empty_na(cell_info_id)) %>%
  filter(!grepl("ascending|descending", super_class)) %>%
  filter(!grepl("agrawal_lab", cell_type_source)) %>%
  mutate(cell_info_id = as.integer(cell_info_id))

# cell_info for join
cell_info_mod <- cell_info %>%
  filter(id %in% banc_meta_cell_info$cell_info_id) %>%
  rename(cell_type = tag) %>%
  rename(cell_info_id = id) %>%
  mutate(cell_info_id = as.integer(cell_info_id)) %>%
  select(cell_info_id, cell_type)

# those rows of banc_meta where cell_type pulled from cell_info
banc_meta_cell_info_mod <- banc_meta_cell_info %>%
  select(`_id`, root_id, cell_type, cell_info_id, fafb_match, manc_match, fafb_cell_type,
         manc_cell_type) %>%
  left_join(cell_info_mod, by = "cell_info_id", suffix = c("_ci", "_bm")) %>%
  filter(cell_type_bm == cell_type_ci) %>%
  mutate(cell_type = cell_type_bm) %>%
  select(-cell_type_bm, -cell_type_ci)


# if cell_type is not fafb_cell_type or manc_cell_type, then pull metadata from
#  franken_meta cell_type (only when cell_type exists there)
banc_meta_cell_info_match_wrong <- banc_meta_cell_info_mod %>%
  filter(cell_type != coalesce(fafb_cell_type, "")) %>%
  filter(cell_type != coalesce(manc_cell_type, ""))

# modify banc_meta_cell_info_match_wrong for joining to franken_meta
banc_meta_cell_info_match_wrong_mod <- banc_meta_cell_info_match_wrong %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_ciMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
banc_meta_cell_info_match_wrong_join <- banc_meta_cell_info_match_wrong_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_ciMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_ciMod, by = "cell_type") %>%
  select(-cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_cell_info_match_wrong_join, 
                      # table='banc_meta',base='banc_meta')



###

# update based off cell_type for ORN, HRN, TRN neurons
# after Alex pushed connectivity-based fix on 1/21/26

banc_meta_orn <- banc_meta %>%
  filter(grepl("^ORN|^HRN|^TRN", cell_type))

# modify banc_meta_orn for joining to franken_meta
banc_meta_orn_mod <- banc_meta_orn %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_ornMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
banc_meta_orn_join <- banc_meta_orn_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_ornMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_ornMod, by = "cell_type") %>%
  select(-cell_type)


# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_orn_join, 
                      # table='banc_meta',base='banc_meta')


###

# for those cells with placeholder names "ORN", "KC", "PAM", fill in rest of
#  metadata

# ORN
banc_meta_orn_ph <- banc_meta %>%
  filter(grepl("^ORN$", cell_type)) %>%
  mutate(flow = "afferent") %>%
  mutate(super_class = "sensory") %>%
  mutate(cell_class = "olfactory_receptor_neuron") %>%
  mutate(cell_sub_class = case_when(
    grepl("antennal_nerve", nerve) ~ "antenna_olfactory_receptor_neuron",
    grepl("maxillary-labial_nerve", nerve) ~ "maxillary_palp_olfactory_receptor_neuron",
    TRUE ~ NA
  )) %>%
  mutate(cell_function = "olfactory") %>%
  mutate(cell_function_detailed = NA) %>%
  mutate(peripheral_target_type = "olfactory_sensillum") %>%
  mutate(body_part_sensory = case_when(
    grepl("antennal_nerve", nerve) ~ "antenna",
    grepl("maxillary-labial_nerve", nerve) ~ "maxillary_palp",
    TRUE ~ NA
  )) %>%
  mutate(body_part_effector = NA) %>%
  mutate(hemilineage = NA) %>%
  mutate(neurotransmitter_verified = "acetylcholine") %>%
  mutate(neuropeptide_verified = NA) %>%
  mutate(manc_match = NA) %>%
  mutate(manc_cell_type = NA) %>%
  select(`_id`, flow, super_class, cell_class, cell_sub_class, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter_verified,
         neuropeptide_verified, manc_match, manc_cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_orn_ph, 
                      # table='banc_meta',base='banc_meta')
  
# KC
banc_meta_kc_ph <- banc_meta %>%
  filter(grepl("^KC$", cell_type)) %>%
  mutate(flow = "intrinsic") %>%
  mutate(super_class = "central_brain_intrinsic") %>%
  mutate(cell_class = "kenyon_cell") %>%
  mutate(cell_sub_class = NA) %>%
  mutate(cell_function = NA) %>%
  mutate(cell_function_detailed = NA) %>%
  mutate(peripheral_target_type = NA) %>%
  mutate(body_part_sensory = NA) %>%
  mutate(body_part_effector = NA) %>%
  mutate(nerve = NA) %>%
  mutate(hemilineage = NA) %>%
  mutate(neurotransmitter_verified = "acetylcholine") %>%
  mutate(neuropeptide_verified = NA) %>%
  mutate(manc_match = NA) %>%
  mutate(manc_cell_type = NA) %>%
  select(`_id`, flow, super_class, cell_class, cell_sub_class, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, nerve, hemilineage, neurotransmitter_verified,
         neuropeptide_verified, manc_match, manc_cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_kc_ph, 
                      # table='banc_meta',base='banc_meta')


# PAM
banc_meta_pam_ph <- banc_meta %>%
  filter(grepl("^PAM$", cell_type) & grepl("bates", cell_type_source)) %>%
  mutate(flow = "intrinsic") %>%
  mutate(super_class = "central_brain_intrinsic") %>%
  mutate(cell_class = "mushroom_body_dopaminergic_neuron") %>%
  mutate(cell_sub_class = "PAM_dopaminergic_neuron") %>%
  mutate(cell_function = NA) %>%
  mutate(cell_function_detailed = NA) %>%
  mutate(peripheral_target_type = NA) %>%
  mutate(body_part_sensory = NA) %>%
  mutate(body_part_effector = NA) %>%
  mutate(nerve = NA) %>%
  mutate(hemilineage = NA) %>%
  mutate(neurotransmitter_verified = "dopamine") %>%
  mutate(neuropeptide_verified = NA) %>%
  mutate(manc_match = NA) %>%
  mutate(manc_cell_type = NA) %>%  
  select(`_id`, flow, super_class, cell_class, cell_sub_class, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, nerve, hemilineage, neurotransmitter_verified,
         neuropeptide_verified, manc_match, manc_cell_type)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_pam_ph, 
                      # table='banc_meta',base='banc_meta')


### 

# cell_type_source has bates

bates_cell_types <- banc_meta %>%
  filter(grepl("bates", cell_type_source)) %>%
  # these labels dealt with above
  filter(!(grepl("^ORN$|^KC$|^PAM$", cell_type))) %>%
  filter(is_not_empty_na(cell_type))

# modify bates_cell_types for joining to franken_meta
bates_cell_types_mod <- bates_cell_types %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_batesMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
bates_cell_types_join <- bates_cell_types_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_batesMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_batesMod, by = "cell_type") %>%
  select(-cell_type)


# WARNING: script not live, uncomment to run
# banctable_update_rows(bates_cell_types_join, 
                      # table='banc_meta',base='banc_meta')

###

# cell type source has bates and cell type contains LB3

bates_cell_types <- banc_meta %>%
  filter(grepl("bates", cell_type_source) & grepl("LB3", cell_type))

# modify bates_cell_types for joining to franken_meta
bates_cell_types_mod <- bates_cell_types %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_batesMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
bates_cell_types_join <- bates_cell_types_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_batesMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_batesMod, by = "cell_type") %>%
  select(-cell_type)


# WARNING: script not live, uncomment to run
# banctable_update_rows(bates_cell_types_join, 
                      # table='banc_meta',base='banc_meta')

###

# fill in based off hemibrain match (hemibrain cell_type exists in franken_meta)
# filter for being empty in other major columns (flow, super_class, cell_class)
# filter for no cell_type_source other than wilson_lab
#  implies there was no previous match or cell_typing
#
# will be replaced with explict comparison to hemibrain tab of cns_meta

banc_meta_only_hemi_cell_type <- banc_meta %>%
  # filter(is_na_or_empty(fafb_match) & is_na_or_empty(manc_match) &
  #          is_na_or_empty(hemibrain_match) & is_na_or_empty(fanc_match)) %>%
  # filter(is_na_or_empty(cell_info_id)) %>%
  filter(is_not_empty_na(cell_type)) %>%
  filter(cell_type == hemibrain_cell_type) %>%
  filter(cell_type %in% franken_meta_both_joinCT$cell_type) %>%
  filter(is_na_or_empty(flow) & is_na_or_empty(super_class) & 
           is_na_or_empty(cell_class)) %>%
  filter(grepl("^wilson_lab$", cell_type_source) | is_na_or_empty(cell_type_source))

# modify banc_meta_only_hemi_cell_type for joining to franken_meta
banc_meta_only_hemi_cell_type_mod <- banc_meta_only_hemi_cell_type %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_hemiMod <- franken_meta_both_joinCT %>%
  select(-dataset)

# join and modify for updating banc_meta
banc_meta_only_hemi_cell_type_join <- banc_meta_only_hemi_cell_type_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_hemiMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_hemiMod, by = "cell_type") %>%
  select(-cell_type)


# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_only_hemi_cell_type_join, 
                      # table='banc_meta',base='banc_meta')


###
# 4/9/26
# update metadata for manually added optic lobe cell_types from Helen 
#  has yang in cell_type_source; after join, filter for optic lobe neurons specifically

banc_meta_yang <- banc_meta %>%
  filter(grepl("yang", cell_type_source))
  
# modify banc_meta_yang for joining to franken_meta
yang_cell_types_mod <- banc_meta_yang %>%
  select(`_id`, cell_type)

# modify franken_meta for join
franken_meta_both_joinCT_yangMod <- franken_meta_both_joinCT %>%
  select(-dataset)


# join and modify for updating banc_meta
yang_cell_types_join <- yang_cell_types_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_yangMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_yangMod, by = "cell_type") %>%
  filter((grepl("optic_lobe", super_class) | grepl("photoreceptor", cell_class) |
           grepl("visual", super_class))) %>%
  select(-cell_type)


# WARNING: script not live, uncomment to run
# banctable_update_rows(yang_cell_types_join, 
                      # table='banc_meta',base='banc_meta')

  


###
# 4/9/26
# update metadata for correct FAFB alignment cell typing for optic lobe neurons
#  will have fafb_alignment_decision == T

banc_meta_optic_align <- banc_meta %>%
  filter((fafb_alignment_decision == "T") | (fafb_alignment_decision == "t"))


# modify banc_meta_optic_align for joining to franken_meta
banc_meta_optic_align_mod <- banc_meta_optic_align %>%
  mutate(cell_type = fafb_alignment_cell_type) %>%
  mutate(fafb_match = fafb_alignment_match) %>%
  select(`_id`, cell_type, fafb_match)

# modify franken_meta for join
franken_meta_both_joinCT_opticMod <- franken_meta_both_joinCT %>%
  select(-dataset)


# join and modify for updating banc_meta
banc_meta_optic_align_join <- banc_meta_optic_align_mod %>%
  filter(cell_type %in% franken_meta_both_joinCT_opticMod$cell_type) %>%
  left_join(franken_meta_both_joinCT_opticMod, by = "cell_type")


# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_optic_align_join, 
                      # table='banc_meta',base='banc_meta')
  

