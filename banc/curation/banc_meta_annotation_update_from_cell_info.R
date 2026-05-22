# script to update banc_meta metadata from cell_info CAVE table cell type 
#  information, based on matching to franken_meta

source("banc/curation/banc_curation_init.R")

fcc <- banc_cave_client()
cave_versions <- fcc$materialize$get_versions()

# when joining for optic lobe annotations, take ones from these trusted sources
user_id_overwrite_optic <- c(392, 1081, 6818, 6830)

banc_meta <- banctable_query()

update_timestamps <- data.frame(
  datetime = as.POSIXct(readLines(banc_meta_cell_info_update_file), 
                        format = "%Y-%m-%d %H:%M:%S")
)

banc_meta_updatedID <- banc_updateids(banc_meta)

banc_version <- as.integer(890)

# get cell_info
cell_info <- banc_cave_query("cell_info", live = 2)
cell_info_read_datetime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Get the latest datetime from update_timestamps
latest_datetime <- max(update_timestamps$datetime, na.rm = TRUE)

# Filter cell_info for rows more recent than last read and update
recent_cell_info <- cell_info %>%
  filter(created > latest_datetime)

# get franken meta
franken_meta <- franken_meta()

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
  filter(is_not_empty_na(cell_type)) %>%
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
  filter(is_not_empty_na(cell_type)) %>%  
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
# update based on new entries in cell_info 

# get list of cell_types from cell_info that are in franken_meta
# filter out annotations from us
cell_info_comm <- recent_cell_info %>%
  filter(user_id != 355) %>% # uploads under Alex's ID
  filter(user_id != 52) %>% # uploads under Helen's ID
  filter((id < 366715) | (id > 412188)) # rows that we contributed

cell_info_ct_anno <- cell_info_comm %>%
  filter(tag %in% franken_meta_ct$cell_type)

# other cell name annotations
cell_info_otherNames_anno <- cell_info_comm %>%
  filter((tag2 == "neuron identity") & !(tag %in% franken_meta_ct$cell_type))


# join cell_info with franken_meta based on cell_type 
cell_info_fm_join <- cell_info_ct_anno %>%
  mutate(cell_type = tag) %>%
  left_join(franken_meta_both_joinCT, by = "cell_type") %>%
  distinct(id, .keep_all = TRUE)

# mod for joining with banc_meta
# align root_id with banc_meta updated ID
cell_info_fm_join_mod <- cell_info_fm_join %>%
  mutate(root_id = pt_root_id) %>%
  mutate(root_id = as.character(root_id)) %>%
  mutate(cell_info_id = as.character(id)) %>%
  mutate(supervoxel_id = as.character(pt_supervoxel_id)) %>%
  mutate(position = sapply(pt_position, function(x) paste(x, collapse = ", "))) %>%
  mutate(cell_info_user_id = as.character(user_id))

cell_info_fm_join_mod <- banc_updateids(cell_info_fm_join_mod)
  

cell_info_fm_join_mod <- cell_info_fm_join_mod %>%
  select(root_id, cell_info_id, cell_info_user_id, cell_type, flow, 
         super_class, cell_class, cell_sub_class, cell_function, 
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter_verified, 
         neuropeptide_verified, sexually_dimorphic)

# filter for rows of banc_meta that having matching root_id, for join
banc_meta_filt <- banc_meta_updatedID %>%
  filter(root_id %in% cell_info_fm_join_mod$root_id)

# root_ids that are missing from banc_meta
# stage for upload to banc_meta
missing_root_id <- cell_info_fm_join %>%
  mutate(root_id = as.character(pt_root_id)) %>%
  # filter for missing ids
  filter(!(root_id %in% banc_meta_updatedID$root_id)) %>%
  # modify columns to match banc_meta
  mutate(cell_info_id = as.character(id)) %>%
  mutate(cell_info_user_id = as.character(user_id)) %>%
  mutate(supervoxel_id = as.character(pt_supervoxel_id)) %>%
  mutate(position = pt_position) %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ", "))) %>%
  select(root_id, supervoxel_id, position, cell_info_id, cell_info_user_id, 
         cell_type, flow, super_class, cell_class, cell_sub_class, cell_function, 
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter_verified, 
         neuropeptide_verified, sexually_dimorphic)

# upload to banc_meta  
# WARNING: script not live, uncomment to run
# banctable_append_rows(df = missing_root_id, table='banc_meta',base='banc_meta')


# deal with root_ids already present in banc_meta
# Get shared columns (excluding the join key)
shared_cols <- intersect(names(banc_meta_filt), names(cell_info_fm_join_mod))
shared_cols <- shared_cols[shared_cols != "root_id"]

# Get the actual matching root_id values between the two dataframes
matching_root_ids <- intersect(banc_meta_filt$root_id, 
                               cell_info_fm_join_mod$root_id)


banc_meta_filt_join <- banc_meta_filt %>%
  left_join(cell_info_fm_join_mod, by = "root_id", suffix = c("_banc", "_cell_info")) %>%
  # Create conflict detection - only for rows with matching root_id values
  rowwise() %>%
  mutate(
    bm_cell_info_conflict = {
      if(!(root_id %in% matching_root_ids)) {
        FALSE
      } else {
        conflicts <- c()
        for(col in shared_cols) {
          banc_col <- paste0(col, "_banc")
          cell_info_col <- paste0(col, "_cell_info")
          
          if(banc_col %in% names(cur_data()) && cell_info_col %in% names(cur_data())) {
            banc_val <- cur_data()[[banc_col]]        # Use cur_data()[[]] instead of get()
            cell_info_val <- cur_data()[[cell_info_col]]  # Use cur_data()[[]] instead of get()
            
            # Conflict if both exist and are different
            conflict <- !(is_na_or_empty(banc_val)) & !(is_na_or_empty(cell_info_val)) & banc_val != cell_info_val
            conflicts <- c(conflicts, conflict)
          }
        }
        any(conflicts, na.rm = TRUE)
      }
    }
  ) %>%
  ungroup()

banc_meta_filt_join <- banc_meta_filt_join %>%
  distinct(`_id`, .keep_all = TRUE)

# For checking conflicts
banc_meta_filt_join_conflicts <- banc_meta_filt_join %>%
  filter(bm_cell_info_conflict) %>%
  select(root_id, `_id`, super_class_banc,
         super_class_cell_info, cell_class_banc, cell_class_cell_info,
         cell_type_banc, cell_type_cell_info, cell_type_banc, 
         cell_info_user_id_cell_info, cell_type_source)

overwrite_id <- c()


# Now create the final columns - conditional overwriting
for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  cell_info_col <- paste0(col, "_cell_info")
  
  if(banc_col %in% names(banc_meta_filt_join)) {
    if(cell_info_col %in% names(banc_meta_filt_join)) {
      # Determine source once per row: use cell_info if:
      # 1. root_id is in the matching set AND
      # 2. (No conflicts between banc and cell_info 
      #     OR annotation comes from trusted optic lobe annotator (by cell_info_user_id)
      #     OR _id is in overwrite_id (manually determined))
      use_cell_info_source <- ((banc_meta_filt_join$root_id %in% matching_root_ids) &
                                 ((!banc_meta_filt_join$bm_cell_info_conflict) |
                                    (banc_meta_filt_join$cell_info_user_id_cell_info %in% user_id_overwrite_optic) |
                                    (banc_meta_filt_join$`_id` %in% overwrite_id)))
      
      # Apply the same source decision to this column
      banc_meta_filt_join[[col]] <- ifelse(use_cell_info_source,
                                           banc_meta_filt_join[[cell_info_col]],  # Always cell_info for this row
                                           banc_meta_filt_join[[banc_col]])       # Always banc for this row
    } else {
      # Only banc column exists
      banc_meta_filt_join[[col]] <- banc_meta_filt_join[[banc_col]]
    }
  }
}


# Remove suffixed columns
banc_meta_filt_join <- banc_meta_filt_join %>%
  select(-ends_with("_banc"), -ends_with("_cell_info"))

banc_meta_filt_join_select <- banc_meta_filt_join %>%
  select(`_id`, cell_info_id, cell_info_user_id, flow, super_class, cell_class,
         cell_sub_class, cell_type, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         hemilineage, neurotransmitter_verified, neuropeptide_verified, 
         sexually_dimorphic)

# update banc_meta
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_filt_join_select, table='banc_meta',base='banc_meta')


# write last time cell_info checked into banc_meta_cell_info_update_file
cat(cell_info_read_datetime, "\n", file = banc_meta_cell_info_update_file, append = TRUE)



###

# update based on what's not already in banc_meta, even if entry was put in
#  cell_info before the last time cell_info was checked

cell_info <- banc_cave_query("cell_info", live = FALSE, version = as.integer(882))

# get list of cell_types from cell_info that are in franken_meta
# filter out annotations from us
cell_info_comm <- cell_info %>%
  filter(user_id != 355) %>% # uploads under Alex's ID
  filter(user_id != 52) %>% # uploads under Helen's ID
  filter((id < 366715) | (id > 412188)) # rows that we contributed

# don't reconsider those entries of cell_info already in banc_meta
# Extract all IDs from the comma-separated cell_info_id column in banc_meta
bm_cell_info_ids <- banc_meta$cell_info_id %>%
  na.omit() %>%
  strsplit(",") %>%
  unlist() %>%
  trimws() %>%
  as.integer()

# Filter cell_info_comm to exclude rows with IDs in banc_meta
cell_info_comm_notBM <- cell_info_comm %>%
  filter(!id %in% bm_cell_info_ids)


cell_info_ct_anno <- cell_info_comm_notBM %>%
  filter(tag %in% franken_meta_ct$cell_type)

# other cell name annotations
cell_info_otherNames_anno <- cell_info_comm_notBM %>%
  filter((tag2 == "neuron identity") & !(tag %in% franken_meta_ct$cell_type))



# join cell_info with franken_meta based on cell_type 
cell_info_fm_join <- cell_info_ct_anno %>%
  mutate(cell_type = tag) %>%
  left_join(franken_meta_both_joinCT, by = "cell_type") %>%
  distinct(id, .keep_all = TRUE)

# mod for joining with banc_meta
cell_info_fm_join_mod <- cell_info_fm_join %>%
  mutate(root_id = pt_root_id) %>%
  mutate(root_id = as.character(root_id)) %>%
  mutate(cell_info_id = as.character(id)) %>%
  mutate(cell_info_user_id = as.character(user_id)) %>%
  select(root_id, cell_info_id, cell_info_user_id, cell_type, flow, 
         super_class, cell_class, cell_sub_class, cell_function, 
         cell_function_detailed, peripheral_target_type, body_part_sensory,
         body_part_effector, hemilineage, neurotransmitter_verified, 
         neuropeptide_verified, sexually_dimorphic)

# filter for rows of banc_meta that having matching root_id, for join
banc_meta_filt <- banc_meta_updatedID %>%
  filter(root_id %in% cell_info_fm_join_mod$root_id)

# those root_ids that have multiple annotations, exclude for now and deal with
#  manually
dup_cell_info <- cell_info_fm_join_mod %>%
  group_by(root_id) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(count) %>%
  filter(count > 1)

banc_meta_filt <- banc_meta_filt %>%
  filter(!(root_id %in% dup_cell_info$root_id))

# Get shared columns (excluding the join key)
shared_cols <- intersect(names(banc_meta_filt), names(cell_info_fm_join_mod))
shared_cols <- shared_cols[shared_cols != "root_id"]

# Get the actual matching root_id values between the two dataframes
matching_root_ids <- intersect(banc_meta_filt$root_id, 
                               cell_info_fm_join_mod$root_id)


# Columns to exclude from conflict detection
conflict_exclude_cols <- c("cell_info_id", "cell_info_user_id")

banc_meta_filt_join <- banc_meta_filt %>%
  left_join(cell_info_fm_join_mod, by = "root_id", suffix = c("_banc", "_cell_info")) %>%
  # Create conflict detection - only for rows with matching root_id values
  rowwise() %>%
  mutate(
    bm_cell_info_conflict = {
      if(!(root_id %in% matching_root_ids)) {
        FALSE
      } else {
        conflicts <- c()
        for(col in shared_cols[!(shared_cols %in% conflict_exclude_cols)]) {
          banc_col <- paste0(col, "_banc")
          cell_info_col <- paste0(col, "_cell_info")
          
          if(banc_col %in% names(cur_data()) && cell_info_col %in% names(cur_data())) {
            banc_val <- cur_data()[[banc_col]]
            cell_info_val <- cur_data()[[cell_info_col]]
            
            # Conflict if both exist and are different
            conflict <- !(is_na_or_empty(banc_val)) & !(is_na_or_empty(cell_info_val)) & banc_val != cell_info_val
            conflicts <- c(conflicts, conflict)
          }
        }
        any(conflicts, na.rm = TRUE)
      }
    }
  ) %>%
  ungroup()


banc_meta_filt_join <- banc_meta_filt_join %>%
  distinct(`_id`, .keep_all = TRUE)

# For checking conflicts
banc_meta_filt_join_conflicts <- banc_meta_filt_join %>%
  filter(bm_cell_info_conflict) %>%
  select(root_id, `_id`, super_class_banc,
         super_class_cell_info, cell_class_banc, cell_class_cell_info,
         cell_type_banc, cell_type_cell_info, cell_type_banc, 
         cell_info_id_banc, cell_info_id_cell_info, cell_info_user_id_banc, cell_info_user_id_cell_info)


banc_meta_filt_join_noConflicts <- banc_meta_filt_join %>%
  filter(!bm_cell_info_conflict)  %>%
  filter(is_na_or_empty(cell_type_banc)) %>%
  select(root_id, `_id`, super_class_banc,
         super_class_cell_info, cell_class_banc, cell_class_cell_info,
         cell_type_banc, cell_type_cell_info, cell_type_banc, 
         cell_info_user_id_cell_info)

banc_meta_filt_join_optic_conflicts <- banc_meta_filt_join %>%
  filter(grepl("optic", super_class_cell_info) | grepl("photoreceptor", cell_class_cell_info) | 
           grepl("optic", super_class_banc) | grepl("photoreceptor", cell_class_banc)) %>%
  select(root_id, `_id`, super_class_banc,
         super_class_cell_info, cell_class_banc, cell_class_cell_info,
         cell_type_banc, cell_type_cell_info, cell_type_banc, 
         cell_info_id_banc, cell_info_id_cell_info, cell_info_user_id_banc, cell_info_user_id_cell_info)

overwrite_id <- c()

# ids to exclude - manual check
wrong_id <- c()

banc_meta_filt_join <- banc_meta_filt_join %>%
  filter(!(`_id` %in% wrong_id))


# Now create the final columns - conditional overwriting
# concatenate cell_info_id and cell_info_user_id when no conflict
for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  cell_info_col <- paste0(col, "_cell_info")
  
  if(banc_col %in% names(banc_meta_filt_join)) {
    if(cell_info_col %in% names(banc_meta_filt_join)) {
      
      # Determine source once per row: use cell_info if:
      # 1. root_id is in the matching set AND
      # 2. (No conflicts between banc and cell_info 
      #     OR annotation comes from trusted optic lobe annotator (by cell_info_user_id)
      #     OR _id is in overwrite_id (manually determined))
      use_cell_info_source <- ((banc_meta_filt_join$root_id %in% matching_root_ids) &
                                 ((!banc_meta_filt_join$bm_cell_info_conflict) |
                                    (banc_meta_filt_join$cell_info_user_id_cell_info %in% user_id_overwrite_optic) |
                                    (banc_meta_filt_join$`_id` %in% overwrite_id)))
      
      if(col %in% c("cell_info_id", "cell_info_user_id")) {
        # For these columns, concatenate when there is no conflict
        no_conflict_match <- (banc_meta_filt_join$root_id %in% matching_root_ids) &
          (!banc_meta_filt_join$bm_cell_info_conflict)
        
        banc_meta_filt_join[[col]] <- mapply(
          function(banc_val, cell_info_val, concat, use_ci) {
            if(isTRUE(concat)) {
              # Concatenate both values, dropping NAs/empty strings
              vals <- c(as.character(banc_val), as.character(cell_info_val))
              vals <- vals[!is.na(vals) & vals != "NA" & vals != ""]
              if(length(vals) == 0) NA_character_ else paste(vals, collapse = ", ")
            } else if(isTRUE(use_ci)) {
              # Conflict but trusted user or overwrite_id - use cell_info
              cell_info_val
            } else {
              # Not in matching set or unresolved conflict - use banc
              banc_val
            }
          },
          banc_meta_filt_join[[banc_col]],
          banc_meta_filt_join[[cell_info_col]],
          no_conflict_match,
          use_cell_info_source,
          SIMPLIFY = TRUE
        )
      } else {
        # Apply the same source decision to this column
        banc_meta_filt_join[[col]] <- ifelse(use_cell_info_source,
                                             banc_meta_filt_join[[cell_info_col]], # Always cell_info for this row
                                             banc_meta_filt_join[[banc_col]]) # Always banc for this row
      }
    } else {
      # Only banc column exists
      banc_meta_filt_join[[col]] <- banc_meta_filt_join[[banc_col]]
    }
  }
}


# Remove suffixed columns
banc_meta_filt_join_select <- banc_meta_filt_join %>%
  select(-ends_with("_banc"), -ends_with("_cell_info"))

banc_meta_filt_join_select <- banc_meta_filt_join_select %>%
  filter((`_id` %in% overwrite_id) | (!bm_cell_info_conflict)) %>%
  select(`_id`, cell_info_id, cell_info_user_id, flow, super_class, cell_class,
         cell_sub_class, cell_type, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         hemilineage, neurotransmitter_verified, neuropeptide_verified, 
         sexually_dimorphic)

# update banc_meta
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_filt_join_select, table='banc_meta',base='banc_meta')


# resolve conflicts

banc_meta_filt_join_conflict_optic <- banc_meta_filt_join %>%
  filter(bm_cell_info_conflict) %>%
  filter(grepl("optic_lobe", super_class_cell_info)) %>%
  filter(cell_info_user_id_cell_info %in% user_id_overwrite_optic)
  


# resolve optic lobe conflicts
for(col in shared_cols) {
  banc_col <- paste0(col, "_banc")
  cell_info_col <- paste0(col, "_cell_info")
  
  if(banc_col %in% names(banc_meta_filt_join_conflict_optic)) {
    if(cell_info_col %in% names(banc_meta_filt_join_conflict_optic)) {
      
      # Extract ID columns for comparison
      banc_id <- banc_meta_filt_join_conflict_optic$cell_info_id_banc
      cell_info_id_val <- banc_meta_filt_join_conflict_optic$cell_info_id_cell_info
      
      # Trusted annotator condition: user is trusted AND
      # (banc cell_info_id is empty OR cell_info_id_cell_info > cell_info_id_banc)
      trusted_annotator_condition <-
        (banc_meta_filt_join_conflict_optic$cell_info_user_id_cell_info %in% user_id_overwrite_optic) &
        (is.na(banc_id) | banc_id == "" |
           (as.integer(cell_info_id_val) > as.integer(banc_id)))
      
      use_cell_info_source <- ((banc_meta_filt_join_conflict_optic$root_id %in% matching_root_ids) &
                                 ((!banc_meta_filt_join_conflict_optic$bm_cell_info_conflict) |
                                    trusted_annotator_condition |
                                    (banc_meta_filt_join_conflict_optic$`_id` %in% overwrite_id)))
      
      banc_meta_filt_join_conflict_optic[[col]] <- ifelse(use_cell_info_source,
                                                          banc_meta_filt_join_conflict_optic[[cell_info_col]],
                                                          banc_meta_filt_join_conflict_optic[[banc_col]])
    } else {
      banc_meta_filt_join_conflict_optic[[col]] <- banc_meta_filt_join_conflict_optic[[banc_col]]
    }
  }
}

# Remove suffixed columns
banc_meta_filt_join_conflict_optic_select <- banc_meta_filt_join_conflict_optic %>%
  select(-ends_with("_banc"), -ends_with("_cell_info"))

banc_meta_filt_join_conflict_optic_select <- banc_meta_filt_join_conflict_optic_select %>%
  select(root_id, `_id`, cell_info_id, cell_info_user_id, flow, super_class, cell_class,
         cell_sub_class, cell_type, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector,
         hemilineage, neurotransmitter_verified, neuropeptide_verified, 
         sexually_dimorphic)

banc_meta_filt_join_conflict_optic_test <- banc_meta_filt_join_conflict_optic %>%
  select(root_id, `_id`, cell_type_banc, cell_type_cell_info, cell_info_id_banc,
         cell_info_id_cell_info)

test_cell_info <- cell_info %>%
  filter(id == 345240 | id == 151919)

# update banc_meta
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_filt_join_conflict_optic_select, table='banc_meta',base='banc_meta')
  
