# script to update banc_meta annotation metadata, based on matching and franken_meta

source("banc/curation/banc_curation_init.R")

# position that separates brain from VNC in y-axis (~middle of neck)
brain_vnc_div_y <- 107000 

# get banc_meta
banc_meta <- banctable_query()

# get franken meta
franken_meta <- franken_meta()

# manc, for join based on ID
franken_meta_manc_joinID <- franken_meta %>%
  filter(is_not_empty_na(manc_id)) %>%
  filter(is_not_empty_na(cell_type)) %>%
  select(manc_id, flow, super_class, cell_class, cell_sub_class, cell_type, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, banc_neurotransmitter_verified,
         neuropeptide, sexually_dimorphic, dataset) %>%
  mutate(manc_cell_type = cell_type) %>%
  mutate(neurotransmitter_verified = neurotransmitter) %>%
  mutate(neuropeptide_verified = neuropeptide) %>%
  select(-neurotransmitter, -neuropeptide) %>%
  distinct(manc_id, .keep_all = TRUE)

# fafb, for join based on ID
franken_meta_fafb_joinID <- franken_meta %>%
  filter(is_not_empty_na(fafb_id)) %>%
  filter(is_not_empty_na(cell_type)) %>%
  select(fafb_id, flow, super_class, cell_class, cell_sub_class, cell_type, cell_function,
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter,
         neuropeptide, sexually_dimorphic, dataset) %>%
  mutate(fafb_cell_type = cell_type) %>%
  mutate(neurotransmitter_verified = neurotransmitter) %>%
  mutate(neuropeptide_verified = neuropeptide) %>%
  select(-neurotransmitter, -neuropeptide) %>%
  distinct(fafb_id, .keep_all = TRUE)


# separate root and nucleus position 
banc_meta <- banc_meta %>%
  separate(nucleus_position, 
           into = c("nucleus_pos_x", "nucleus_pos_y", "nucleus_pos_z"), 
           sep = ",\\s*", 
           remove = FALSE,  # Keep original column
           convert = TRUE) %>%
  separate(root_position, 
           into = c("root_pos_x", "root_pos_y", "root_pos_z"), 
           sep = ",\\s*", 
           remove = FALSE,  # Keep original column
           convert = TRUE) 

# select which rows of banc_meta are eligible for change
# exclude those where hemibrain match is preferred (status == HEMIBRAIN_MATCH_PREFERRED)
banc_meta_include <- banc_meta %>%
  # filter(!grepl("HAS_MANUAL_ANNOTATION", status)) %>%
  # filter(is_na_or_empty(cell_info_id)) %>%
  filter(is_not_empty_na(fafb_match) | is_not_empty_na(manc_match)) %>%
  filter(!grepl("HEMIBRAIN_MATCH_PREFERRED", status)) %>%
  filter(!grepl("FANC_MATCH_PREFERRED", STATUS))


banc_meta_include_fafb <- banc_meta_include %>%
  filter(
    # fafb_match exists but manc_match doesn't
    (((is_not_empty_na(fafb_match) & (is_na_or_empty(manc_match))) |
        # Both exist and super_class contains "descending"
        (is_not_empty_na(fafb_match) & is_not_empty_na(manc_match) & grepl("descending", super_class)) |
        # Both exist and nucleus_pos_y or root_pos_y < brain_vnc_div_y
        (is_not_empty_na(fafb_match) & is_not_empty_na(manc_match) & 
           (nucleus_pos_y < brain_vnc_div_y | root_pos_y < brain_vnc_div_y))) & 
       !(grepl("ascending", super_class)))
  ) %>%
  # select requisite rows
  select(root_id, `_id`, fafb_match, flow, super_class, cell_class,
         cell_sub_class, cell_type, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         hemilineage, neurotransmitter_verified, neuropeptide_verified,
         sexually_dimorphic, status, cell_type_source, cell_info_id)

banc_meta_include_manc <- banc_meta_include %>%
  filter(
    # manc_match exists but fafb_match doesn't
    (((is_not_empty_na(manc_match) & (is_na_or_empty(fafb_match))) |
        # Both exist and super_class contains "ascending"
        (is_not_empty_na(fafb_match) & is_not_empty_na(manc_match) & grepl("ascending", super_class)) |
        # Both exist and nucleus_pos_y or root_pos_y > brain_vnc_div_y
        (is_not_empty_na(fafb_match) & is_not_empty_na(manc_match) & 
           (nucleus_pos_y > brain_vnc_div_y | root_pos_y > brain_vnc_div_y))) &
       !(grepl("descending", super_class)))
  ) %>%
  # select requisite rows
  select(root_id, `_id`, manc_match, flow, super_class, cell_class,
         cell_sub_class, cell_type, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         hemilineage, neurotransmitter_verified, neuropeptide_verified, 
         sexually_dimorphic, status, cell_type_source, cell_info_id)

####
# JOIN
####

# flags whether there's an error, and also indicates which columns it occurs in
# also flags rows where banc_meta entry is empty but franken_meta one is not

# join fafb dataframe with franken_meta fafb subsets
franken_meta_fafb_joinID_mod <- franken_meta_fafb_joinID %>%
  mutate(fafb_match = fafb_id) %>%
  select(-fafb_id)

# Get shared columns (excluding the join key)
shared_cols <- intersect(names(banc_meta_include_fafb), names(franken_meta_fafb_joinID_mod))
shared_cols <- shared_cols[shared_cols != "fafb_match"]

# Get the actual matching fafb_match values between the two dataframes
matching_fafb_ids <- intersect(banc_meta_include_fafb$fafb_match, 
                               franken_meta_fafb_joinID_mod$fafb_match)

banc_meta_include_fafb_fm <- banc_meta_include_fafb %>%
  left_join(franken_meta_fafb_joinID_mod, by = "fafb_match", suffix = c("_banc", "_franken")) %>%
  # Create conflict detection - only for rows with matching fafb_match values
  rowwise() %>%
  mutate(
    conflict_info = list({
      if(!(fafb_match %in% matching_fafb_ids)) {
        list(has_conflict = FALSE, conflict_cols = NA_character_, has_bm_empty = FALSE)
      } else {
        conflict_cols <- c()
        bm_empty_cols <- c()
        
        for(col in shared_cols) {
          banc_col <- paste0(col, "_banc")
          franken_col <- paste0(col, "_franken")
          
          if(banc_col %in% names(cur_data()) && franken_col %in% names(cur_data())) {
            banc_val <- get(banc_col)
            franken_val <- get(franken_col)
            
            # Check for bm_empty condition: banc is empty but franken is not
            if(!is_not_empty_na(banc_val) & is_not_empty_na(franken_val)) {
              bm_empty_cols <- c(bm_empty_cols, col)
            }
            
            # Conflict if both exist and are different
            if(is_not_empty_na(banc_val) & is_not_empty_na(franken_val) & 
               (banc_val != franken_val)) {
              conflict_cols <- c(conflict_cols, col)
            }
          }
        }
        
        list(
          has_conflict = length(conflict_cols) > 0,
          conflict_cols = if(length(conflict_cols) > 0) paste(conflict_cols, collapse = ", ") else NA_character_,
          has_bm_empty = length(bm_empty_cols) > 0
        )
      }
    }),
    fafb_match_bm_fm_conflict = conflict_info$has_conflict,
    conflict_cols = conflict_info$conflict_cols,
    bm_empty = conflict_info$has_bm_empty
  ) %>%
  ungroup() %>%
  select(-conflict_info)

banc_meta_include_fafb_fm <- banc_meta_include_fafb_fm %>%
  distinct(`_id`, .keep_all = TRUE)



###

# when there is no conflict but there are NA/empty values in banc_meta, fill
#  those in from franken_meta
banc_meta_include_fafb_fm_noCon_empty <- banc_meta_include_fafb_fm %>%
  filter(!fafb_match_bm_fm_conflict & bm_empty)

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_fafb_fm_noCon_empty[[col]] <- 
    banc_meta_include_fafb_fm_noCon_empty[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_fafb_fm_noCon_empty <- banc_meta_include_fafb_fm_noCon_empty %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_fafb_fm_noCon_empty_update <- banc_meta_include_fafb_fm_noCon_empty %>%
  select(-root_id, -fafb_match, -status, -cell_type_source, -dataset,
         -fafb_cell_type, -fafb_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_fafb_fm_noCon_empty_update, 
                      # table='banc_meta',base='banc_meta')

###

# EXCLUDE optic_lobe_intrinsic and photoreceptors

# when there is no conflict but there are NA/empty values in banc_meta, fill
#  those in from franken_meta
banc_meta_include_fafb_fm_noCon_empty <- banc_meta_include_fafb_fm %>%
  filter(!fafb_match_bm_fm_conflict & bm_empty)

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_fafb_fm_noCon_empty[[col]] <- 
    banc_meta_include_fafb_fm_noCon_empty[[franken_col]]
}

# remove optic_lobe_intrinsic matches
banc_meta_include_fafb_fm_noCon_empty_noOptic <- banc_meta_include_fafb_fm_noCon_empty %>%
  filter(!(grepl("optic_lobe_intrinsic", super_class_franken))) %>%
  filter(!(grepl("photoreceptor", cell_class_franken)))

# Remove suffixed columns
banc_meta_include_fafb_fm_noCon_empty_noOptic <- banc_meta_include_fafb_fm_noCon_empty_noOptic %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_fafb_fm_noCon_empty_noOptic_update <- banc_meta_include_fafb_fm_noCon_empty_noOptic %>%
  select(-root_id, -fafb_match, -status, -cell_type_source, -dataset,
         -fafb_cell_type, -fafb_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_fafb_fm_noCon_empty_noOptic_update, 
                      # table='banc_meta',base='banc_meta')

# for optic_lobe_intrinsic or photoreceptor, clear fafb_match, fafb_cell_type, cell_type
# append status FAFB_PNG_MATCH_WRONG
optic_clear <- banc_meta_include_fafb_fm_noCon_empty %>%
  filter(grepl("optic_lobe_intrinsic", super_class_franken) |
           grepl("photoreceptor", cell_class_franken)) %>%
  mutate(super_class = NA) %>%
  mutate(cell_type = NA) %>%
  mutate(cell_class = NA) %>%
  mutate(cell_sub_class = NA) %>%
  mutate(hemilineage = NA) %>%
  mutate(cell_function = NA) %>%
  mutate(cell_function_detailed = NA) %>%
  mutate(peripheral_target_type = NA) %>%
  mutate(body_part_sensory = NA) %>%
  mutate(body_part_effector = NA) %>%
  mutate(sexually_dimorphic = NA) %>%
  mutate(neurotransmitter_verified = NA) %>%
  mutate(neuropeptide_verified = NA) %>%
  mutate(fafb_match = NA) %>%
  mutate(fafb_cell_type = NA) %>%
  rowwise() %>%
  mutate(status = append_status(status, "FAFB_PNG_MATCH_WRONG")) %>%
  ungroup() %>%
  mutate(status_text = status) %>%
  select(`_id`, super_class, cell_class, cell_sub_class, cell_type, 
         hemilineage, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         sexually_dimorphic, neurotransmitter_verified, neuropeptide_verified,
         fafb_match, fafb_cell_type, status, status_text)

# WARNING: script not live, uncomment to run
# banctable_update_rows(optic_clear, 
                      # table='banc_meta',base='banc_meta')
  


###

# only conflict in cell_type and neuron is not: 
# flow==efferent, cell_type==wCHIN or nCHIN, cell_class contains central_complex,
#  cell_type_source does not contain lee_lab, seeds_hampel_lab, bates, yang
# exclude ORN/TRN/HRN and JOs, exclude things with cell_info_id (cell type from cell_info)
banc_meta_include_fafb_fm_cellTypeCon <- banc_meta_include_fafb_fm %>%
  filter(fafb_match_bm_fm_conflict & grepl("^cell_type$", conflict_cols)) %>%
  filter(!grepl("efferent", flow_banc)) %>%
  filter(!grepl("wCHIN|nCHIN", cell_type_banc)) %>%
  filter(!grepl("central_complex", cell_class_banc)) %>%
  filter(!grepl("lee_lab", cell_type_source)) %>%
  filter(!grepl("seeds_hampel_lab", cell_type_source)) %>%
  filter(!grepl("bates|yang", cell_type_source)) %>%
  filter(!grepl("ORN|TRN|HRN", cell_type_banc)) %>%
  filter(!grepl("JO-", cell_type_banc)) %>%
  filter(is_na_or_empty(cell_info_id))

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_fafb_fm_cellTypeCon[[col]] <- 
    banc_meta_include_fafb_fm_cellTypeCon[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_fafb_fm_cellTypeCon <- banc_meta_include_fafb_fm_cellTypeCon%>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_fafb_fm_cellTypeCon_update <- banc_meta_include_fafb_fm_cellTypeCon %>%
  select(-root_id, -fafb_match, -status, -cell_type_source, -dataset,
         -fafb_cell_type, -fafb_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_fafb_fm_cellTypeCon_update, 
                      # table='banc_meta',base='banc_meta')


###
# has conflict yet,
# cell_type == fafb_cell_type, cell_type_source has nothing but wilson_lab
#  status doesn't have HAS_MANUAL_ANNOTATION, and conflict not in flow or
#  super_class
# also, yao_lab annotations where there's more cell_function_detailed info
banc_meta_include_fafb_fm_sameCellType <- banc_meta_include_fafb_fm %>%
  filter(fafb_match_bm_fm_conflict) %>%
  filter(cell_type_banc == fafb_cell_type) %>%
  filter(grepl("^wilson_lab$", cell_type_source) | is_na_or_empty(cell_type_source) |
         (grepl("yao_lab", cell_type_source) & grepl("^cell_function_detailed$", conflict_cols))) %>%
  filter(!grepl("HAS_MANUAL_ANNOTATION", status)) %>%
  filter(!grepl("flow|super_class", conflict_cols)) %>%
  filter(`_id` == "3ff1sKwCT-SVwMbWIVB39w" | `_id` == "dfyO_6ZaSEeUiXAm-7DNWg")

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_fafb_fm_sameCellType[[col]] <- 
    banc_meta_include_fafb_fm_sameCellType[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_fafb_fm_sameCellType <- banc_meta_include_fafb_fm_sameCellType %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_fafb_fm_sameCellType_update <- banc_meta_include_fafb_fm_sameCellType %>%
  select(-root_id, -fafb_match, -status, -cell_type_source, -dataset,
         -fafb_cell_type, -fafb_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_fafb_fm_sameCellType_update, 
                      # table='banc_meta',base='banc_meta')



####
# MANC
####

# join manc dataframe with franken_meta manc subsets
franken_meta_manc_joinID_mod <- franken_meta_manc_joinID %>%
  mutate(manc_match = manc_id) %>%
  select(-manc_id)

# Get shared columns (excluding the join key)
shared_cols <- intersect(names(banc_meta_include_manc), names(franken_meta_manc_joinID_mod))
shared_cols <- shared_cols[shared_cols != "manc_match"]

# Get the actual matching manc_match values between the two dataframes
matching_manc_ids <- intersect(banc_meta_include_manc$manc_match, 
                               franken_meta_manc_joinID_mod$manc_match)

# flags whether there's an error, also indicates which columns it occurs in
# also flags rows where banc_meta entry is empty but franken_meta one is not

banc_meta_include_manc_fm <- banc_meta_include_manc %>%
  left_join(franken_meta_manc_joinID_mod, by = "manc_match", suffix = c("_banc", "_franken")) %>%
  # Create conflict detection - only for rows with matching fafb_match values
  rowwise() %>%
  mutate(
    conflict_info = list({
      if(!(manc_match %in% matching_manc_ids)) {
        list(has_conflict = FALSE, conflict_cols = NA_character_, has_bm_empty = FALSE)
      } else {
        conflict_cols <- c()
        bm_empty_cols <- c()
        
        for(col in shared_cols) {
          banc_col <- paste0(col, "_banc")
          franken_col <- paste0(col, "_franken")
          
          if(banc_col %in% names(cur_data()) && franken_col %in% names(cur_data())) {
            banc_val <- get(banc_col)
            franken_val <- get(franken_col)
            
            # Check for bm_empty condition: banc is empty but franken is not
            if(!is_not_empty_na(banc_val) & is_not_empty_na(franken_val)) {
              bm_empty_cols <- c(bm_empty_cols, col)
            }
            
            # Conflict if both exist and are different
            if(is_not_empty_na(banc_val) & is_not_empty_na(franken_val) & 
               (banc_val != franken_val)) {
              conflict_cols <- c(conflict_cols, col)
            }
          }
        }
        
        list(
          has_conflict = length(conflict_cols) > 0,
          conflict_cols = if(length(conflict_cols) > 0) paste(conflict_cols, collapse = ", ") else NA_character_,
          has_bm_empty = length(bm_empty_cols) > 0
        )
      }
    }),
    manc_match_bm_fm_conflict = conflict_info$has_conflict,
    conflict_cols = conflict_info$conflict_cols,
    bm_empty = conflict_info$has_bm_empty
  ) %>%
  ungroup() %>%
  select(-conflict_info)

banc_meta_include_manc_fm <- banc_meta_include_manc_fm %>%
  distinct(`_id`, .keep_all = TRUE)


###

# when there is no conflict but there are NA/empty values in banc_meta, fill
#  those in from franken_meta
banc_meta_include_manc_fm_noCon_empty <- banc_meta_include_manc_fm %>%
  filter(!manc_match_bm_fm_conflict & bm_empty)

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_manc_fm_noCon_empty[[col]] <- 
    banc_meta_include_manc_fm_noCon_empty[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_manc_fm_noCon_empty <- banc_meta_include_manc_fm_noCon_empty %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_manc_fm_noCon_empty_update <- banc_meta_include_manc_fm_noCon_empty %>%
  select(-root_id, -manc_match, -status, -cell_type_source, -dataset,
         -manc_cell_type, -manc_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_manc_fm_noCon_empty_update, 
                      # table='banc_meta',base='banc_meta')



###

# only conflict in cell_type and neuron is not: 
#  flow==efferent, cell_type==wCHIN or nCHIN, status!=FANC_MATCH_PREFERRED,
#  not just because a _subtype was appended to the end of the MANC name (e.g. 
#  for ANs, see above section that deals with them)

banc_meta_include_manc_fm_cellTypeCon <- banc_meta_include_manc_fm %>%
  filter(manc_match_bm_fm_conflict & grepl("^cell_type$", conflict_cols)) %>%
  filter(!grepl("efferent", flow_banc)) %>%
  filter(!grepl("wCHIN|nCHIN", cell_type_banc)) %>%
  filter(!grepl("FANC_MATCH_PREFERRED", status)) %>%
  rowwise() %>%
  filter(!grepl(cell_type_franken, cell_type_banc)) %>%
  ungroup()

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_manc_fm_cellTypeCon[[col]] <- 
    banc_meta_include_manc_fm_cellTypeCon[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_manc_fm_cellTypeCon <- banc_meta_include_manc_fm_cellTypeCon %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_manc_fm_cellTypeCon_update <- banc_meta_include_manc_fm_cellTypeCon %>%
  select(-root_id, -manc_match, -status, -cell_type_source, -dataset,
         -manc_cell_type, -manc_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_manc_fm_cellTypeCon_update, 
                      # table='banc_meta',base='banc_meta')


###

# has conflict yet,
# cell_type == manc_cell_type, cell_type_source has nothing but wilson_lab
#  status doesn't have HAS_MANUAL_ANNOTATION, and conflict not in flow or
#  super_class
banc_meta_include_manc_fm_sameCellType <- banc_meta_include_manc_fm %>%
  filter(manc_match_bm_fm_conflict) %>%
  filter(cell_type_banc == manc_cell_type) %>%
  filter(grepl("^wilson_lab$", cell_type_source) | is_na_or_empty(cell_type_source)) %>%
  filter(!grepl("HAS_MANUAL_ANNOTATION", status)) %>%
  filter(!grepl("flow|super_class", conflict_cols)) %>%
  # sexually_dimorphic conflicts resolve later, when Princeton has finalized list
  filter(!grepl("sexually_dimorphic", conflict_cols))

for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  franken_col <- paste0(col, "_franken")
  
  banc_meta_include_manc_fm_sameCellType[[col]] <- 
    banc_meta_include_manc_fm_sameCellType[[franken_col]]
}

# Remove suffixed columns
banc_meta_include_manc_fm_sameCellType <- banc_meta_include_manc_fm_sameCellType %>%
  select(-ends_with("_banc"), -ends_with("_franken"))

# more cleaning 
banc_meta_include_manc_fm_sameCellType_update <- banc_meta_include_manc_fm_sameCellType %>%
  select(-root_id, -manc_match, -status, -cell_type_source, -dataset,
         -manc_cell_type, -manc_match_bm_fm_conflict, -conflict_cols, 
         -bm_empty, -sexually_dimorphic)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_include_manc_fm_sameCellType_update, 
                      # table='banc_meta',base='banc_meta')
