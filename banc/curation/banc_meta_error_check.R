# script to check for inconsistencies in banc_meta

source("banc/curation/banc_curation_init.R")

# get banc_meta
banc_meta <- banctable_query()

# get franken meta
franken_meta <- franken_meta()

# select only the relevant columns of banc_meta
# make sure columns with multiple entries are sorted
banc_meta_sel <- banc_meta %>%
  select(region, hemilineage, nerve, flow, super_class,
         cell_class, cell_sub_class, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector, 
         cell_type, neuropeptide_verified, neurotransmitter_verified, 
         sexually_dimorphic,
         root_id, `_id`, status, cell_type_source) %>%
  mutate(nerve = str_replace(nerve, "^(left_|right_)", "")) %>%
  rowwise() %>%
  mutate(cell_function = sort_status(cell_function)) %>%
  mutate(cell_function_detailed = sort_status(cell_function_detailed)) %>%
  mutate(peripheral_target_type = sort_status(peripheral_target_type)) %>%
  mutate(body_part_sensory = sort_status(body_part_sensory)) %>%
  mutate(body_part_effector = sort_status(body_part_effector)) %>%
  mutate(neuropeptide_verified = sort_status(neuropeptide_verified)) %>%
  mutate(neurotransmitter_verified = sort_status(neurotransmitter_verified)) %>%
  mutate(nerve = sort_status(nerve)) %>%
  ungroup()


# select only relevant columns of franken_meta
# make sure columns with multiple entries are sorted
franken_meta_sel <- franken_meta %>%
  select(region, hemilineage, nerve, flow, super_class,
         cell_class, cell_sub_class, cell_function, cell_function_detailed,
         peripheral_target_type, body_part_sensory, body_part_effector, 
         cell_type, neuropeptide, neurotransmitter,
         sexually_dimorphic) %>%
  mutate(nerve = str_replace(nerve, "^(left_|right_)", "")) %>%
  rowwise() %>%
  mutate(cell_function = sort_status(cell_function)) %>%
  mutate(cell_function_detailed = sort_status(cell_function_detailed)) %>%
  mutate(peripheral_target_type = sort_status(peripheral_target_type)) %>%
  mutate(body_part_sensory = sort_status(body_part_sensory)) %>%
  mutate(body_part_effector = sort_status(body_part_effector)) %>%
  mutate(neuropeptide_verified = sort_status(neuropeptide)) %>%
  mutate(neurotransmitter_verified = sort_status(neurotransmitter)) %>%
  mutate(nerve = sort_status(nerve)) %>%
  ungroup() %>%
  select(-neurotransmitter, -neuropeptide)

# make a version of franken_meta_sel that removes duplicate rows
# this is the list of all possible valid combinations
franken_meta_opts <- franken_meta_sel %>%
  mutate(across(everything(), ~ if_else(.x == "", NA_character_, .x))) %>%
  distinct()


# identify all rows of banc_meta_sel where the combination of entries is not
#  found in franken_meta_opts
# operate only on rows that have a cell type

# Get shared columns (excluding _id, status, root_id, cell_type_source)
shared_cols <- setdiff(
  intersect(names(banc_meta_sel), names(franken_meta_opts)),
  c("_id", "root_id", "status", "cell_type_source")
)

banc_meta_err_rows <- banc_meta_sel %>%
  filter(is_not_empty_na(cell_type)) %>%
  filter((cell_type %in% franken_meta_opts$cell_type)) %>%
  mutate(across(all_of(shared_cols), ~ if_else(.x == "", NA_character_, .x))) %>%
  anti_join(
    franken_meta_opts %>%
      mutate(across(all_of(shared_cols), ~ if_else(.x == "", NA_character_, .x))),
    by = shared_cols
  )

# Get shared columns (excluding _id, status, root_id, cell_type_source)
# Exclude join key of cell_type
shared_cols <- setdiff(
  intersect(names(banc_meta_sel), names(franken_meta_opts)),
  c("_id", "root_id", "status", "cell_type_source", "cell_type")
)

# returns details of which columns have errors for each row, assuming cell_type
#  is correct
banc_meta_err_rows_detailed <- banc_meta_err_rows %>%
  left_join(franken_meta_opts, by = "cell_type", suffix = c("_banc", "_franken")) %>%
  # Create conflict detection
  rowwise() %>%
  mutate(
    conflict_info = list({
      if(!(cell_type %in% franken_meta_opts$cell_type)) {
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
            
            # Conflict they are different
            if(coalesce(banc_val,"") != coalesce(franken_val,"")) {
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
    bm_fm_conflict = conflict_info$has_conflict,
    conflict_cols = conflict_info$conflict_cols,
    bm_empty = conflict_info$has_bm_empty
  ) %>%
  ungroup() %>%
  select(-conflict_info)

banc_meta_err_rows_detailed <- banc_meta_err_rows_detailed %>%
  distinct(`_id`, .keep_all = TRUE)


# arrange columns
interleaved_cols <- as.vector(rbind(
  paste0(shared_cols, "_banc"),
  paste0(shared_cols, "_franken")
))

# Reorder columns
banc_meta_err_rows_detailed <- banc_meta_err_rows_detailed %>%
  select(
    `_id`, root_id, cell_type, bm_fm_conflict, conflict_cols, bm_empty, 
    status, cell_type_source, 
    all_of(interleaved_cols)
  )

  


# write to CSV
csv_file_path <- ""

write.csv(banc_meta_err_rows_detailed, 
          file = file.path(csv_file_path, "banc_meta_error_rows_detailed.csv"),
          row.names = FALSE)

write.csv(franken_meta_opts, 
          file = file.path(csv_file_path, "franken_meta_cell_types.csv"),
          row.names = FALSE)


