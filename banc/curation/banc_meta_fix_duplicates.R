# deal with duplicated rows in banc_meta

source("banc/curation/banc_curation_init.R")

banc_meta <- banctable_query()
banc_meta <- banc_updateids(banc_meta) %>%
  distinct(`_id`, .keep_all = TRUE)

# get franken meta
franken_meta <- franken_meta()

# replace semi-colon separated list with commas
franken_meta <- franken_meta %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, ";\\s*", ",")))

# get list of all cell_types 
franken_meta_ct <- franken_meta %>%
  distinct(cell_type)

# deal with those with DELETE_THIS_DUPLICATE tag
banc_meta_delete_dup_status <- banc_meta %>%
  filter((grepl("DELETE_THIS_DUPLICATE",status)) |
           grepl("DELETE_THIS_DUPLICATE",status_text))

banc_meta_delete_dup_wDup <- banc_meta %>%
  filter(root_id %in% banc_meta_delete_dup_status$root_id)

# Count how many times each root_id appears in banc_meta_delete_dup_wDup
root_id_dup_counts <- banc_meta_delete_dup_wDup %>%
  count(root_id)

# Find root_ids that appear only once (i.e. no duplicate found)
no_duplicate <- root_id_dup_counts %>%
  filter(n < 2) %>%
  pull(root_id)

# remove DELETE_THIS_DUPLICATE status tag for those that have no duplicate
banc_meta_rmv_del_status <- banc_meta %>%
  filter(root_id %in% no_duplicate) %>%
  rowwise() %>%
  mutate(status = subtract_status(status, "DELETE_THIS_DUPLICATE")) %>%
  mutate(status = subtract_status(status, "DELETE_THIS_DUPLICATE")) %>%
  ungroup() %>%
  select(`_id`, status, status_text)

# update banc_meta
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_rmv_del_status, 
                      # table='banc_meta',base='banc_meta')

# manually delete those rows using web SeaTable interface

# re-run banctable_query() after update
banc_meta <- banctable_query()

banc_meta_delete_dup_status <- banc_meta %>%
  filter((grepl("DELETE_THIS_DUPLICATE",status)) |
           grepl("DELETE_THIS_DUPLICATE",status_text))


###

# remove _X columns except for _id from consideration
banc_meta <- banc_meta %>%
  select(-`_locked`, -`_locked_by`, -`_archived`, -`_creator`, -`_ctime`,
         -`_last_modifier`, -`_mtime`)

# standardize position columns
banc_meta <- banc_meta %>%
  mutate(
    position = standardize_position(position),
    nucleus_position = standardize_position(nucleus_position),
    nucleus_position_nm = standardize_position(nucleus_position_nm),
    root_position = standardize_position(root_position),
    root_position_nm = standardize_position(root_position_nm)
  )

# find all duplicated root_ids
dup_root_id <- banc_meta$root_id[duplicated(banc_meta$root_id) | duplicated(banc_meta$root_id, fromLast = TRUE)]
dup_root_id <- unique(dup_root_id)

test_dup_root_id <- banc_meta %>%
  group_by(root_id) %>%
  summarize(count = n(), .groups = "drop") %>%
  filter(count > 1)

# return banc_meta with only rows that have duplicated root_ids
banc_meta_dup <- banc_meta %>%
  filter(root_id %in% dup_root_id)


###
# new way of dealing with duplicates now that publication version cell_ids_id
#  exist and are stable
# if 1 row has cell_ids_tag and the others don't, then use the row with cell_ids_tag
#  if multiple rows have cell_ids_tag and the value is different, keep and flag
#   for additional eval
#  if multiple rows have cell_ids_tag and the value is the same or none of the 
#   rows have cell_ids_tag, feed that root_id into the standard pipeline

# Initialize delete_this_row column in banc_meta_dup
banc_meta_dup$delete_this_row <- FALSE

# Initialize dup_root_id_sln
dup_root_id_sln <- data.frame(root_id = character(), sln = character(), stringsAsFactors = FALSE)

# Process each duplicated root_id
for (rid in dup_root_id) {
  
  # Get indices of rows belonging to this root_id
  row_idx <- which(banc_meta_dup$root_id == rid)
  
  # Get cell_ids_tag values for those rows
  tags <- banc_meta_dup$cell_ids_tag[row_idx]
  
  # Logical vector: which rows have a non-NA, non-empty tag
  has_tag <- !is.na(tags) & nzchar(tags)
  
  if (sum(has_tag) == 1) {
    # Case 1: exactly one row has a tag — delete the rest
    banc_meta_dup$delete_this_row[row_idx[!has_tag]] <- TRUE
    banc_meta_dup$delete_this_row[row_idx[has_tag]]  <- FALSE
    sln_val <- "cell_ids"
    
  } else if (sum(has_tag) > 1) {
    # Multiple rows have tags
    tag_values <- tags[has_tag]
    banc_meta_dup$delete_this_row[row_idx] <- FALSE
    
    if (length(unique(tag_values)) > 1) {
      # Case 2: tags present but different
      sln_val <- "multi_cell_ids"
    } else {
      # Case 3: tags present and identical
      sln_val <- "reg_pipeline_w_id"
    }
    
  } else {
    # Case 4: no rows have a tag
    banc_meta_dup$delete_this_row[row_idx] <- FALSE
    sln_val <- "reg_pipeline_no_id"
  }
  
  dup_root_id_sln <- rbind(
    dup_root_id_sln,
    data.frame(root_id = rid, sln = sln_val, stringsAsFactors = FALSE)
  )
}

# deal with each category of duplicates
# Case 1: exactly one row has a tag — delete the rest
cell_ids_root_id <- dup_root_id_sln %>%
  filter(sln == "cell_ids")

banc_meta_dup_tag_upload <- banc_meta_dup %>%
  filter(root_id %in% cell_ids_root_id$root_id) %>%
  filter(delete_this_row == TRUE) %>%
  rowwise() %>%
  mutate(status = append_status(status, "DELETE_THIS_DUPLICATE")) %>%
  mutate(status_text = append_status(status_text, "DELETE_THIS_DUPLICATE")) %>%
  ungroup() %>%
  select(root_id, status, status_text, `_id`)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_dup_tag_upload, table='banc_meta',base='banc_meta')

# Case 2: tags present but different 
multi_cell_ids_root_id <- dup_root_id_sln %>%
  filter(sln == "multi_cell_ids")

banc_meta_dup_multi_cell_ids <- banc_meta_dup %>%
  filter(root_id %in% multi_cell_ids_root_id$root_id)
# deal with manually

# Case 3: tags present and identical
reg_pipeline_w_id_root_id <- dup_root_id_sln %>%
  filter(sln == "reg_pipeline_w_id")
# none

# Case 4: no rows have a tag
reg_pipeline_no_id_root_id <- dup_root_id_sln %>%
  filter(sln == "reg_pipeline_no_id")

banc_meta_dup_no_id <- banc_meta_dup %>%
  filter(root_id %in% reg_pipeline_no_id_root_id$root_id) %>%
  select(-delete_this_row)
dup_data <- banc_meta_dup_no_id



###
dup_data <- banc_meta_dup



# columns to exclude from conflict checking
exclude_cols <- c("status", "cell_type_source", "status_text", 
                  "cell_type_source_text", "_id")
numeric_cols <- c("output_connections", "input_connections", "l2_nodes", 
                  "segregation_index", "pd_width", "input_side_index", 
                  "output_side_index", "mitochondria_volume", "l2_cable_length_um",
                  "volume_nm3", "soma_dcv_count", "soma_dcv_density", 
                  "mitochondria", "neurotransmitter_score", "banc_nblast",
                  "neurotransmitter_score_v1", "neurotransmitter_score_v2",
                  "neurotransmitter_score_v3", "malecns_nblast", 
                  "fafb_nblast", "manc_nblast", "fanc_nblast", "hemibrain_nblast",
                  "cell_ids_id", "cell_ids_id_626")
special_merge_cols <- c("status", "cell_type_source", "other_names", "manc_match",
                        "fafb_match", "status_text", "cell_type_source_text")
notes_merge_cols <- c("notes")


# Helper functions
is_na_or_empty_numeric <- function(x, col_name) {
  if(col_name %in% numeric_cols) {
    return(is.na(x) | x == 0)
  } else {
    return(is.na(x) | x == "" | (is.character(x) & trimws(x) == ""))
  }
}

get_empty_value <- function(col_name) {
  if(col_name %in% numeric_cols) {
    # return(0)
    return(NA)
  } else {
    return(NA)
  }
}

# Initialize outputs
merge_conflict_root_id <- data.frame(root_id = character(), conflict_cols = character(), stringsAsFactors = FALSE)
banc_meta_merged_dup <- data.frame()

# Process each group
dup_groups <- dup_data %>%
  group_by(root_id) %>%
  group_split()

for(group in dup_groups) {
  # Check for conflicts (only in non-excluded columns)
  regular_cols <- setdiff(names(group), exclude_cols)
  conflict_cols <- c()
  
  for(col in regular_cols) {
    valid_values <- group[[col]][!sapply(group[[col]], is_na_or_empty_numeric, col_name = col)]
    if(length(unique(valid_values)) > 1) {
      conflict_cols <- c(conflict_cols, col)
    }
  }
  
  if(length(conflict_cols) > 0) {
    # Add to conflict dataframe
    new_conflict <- data.frame(
      root_id = group$root_id[1],
      conflict_cols = paste(conflict_cols, collapse = ", "),
      stringsAsFactors = FALSE
    )
    merge_conflict_root_id <- bind_rows(merge_conflict_root_id, new_conflict)
  } else {
    # Merge the group
    merged_row <- group[1, ]  # Template
    
    # Get column categories
    all_cols <- names(group)
    excluded_cols <- intersect(all_cols, exclude_cols)
    cols_for_regular_merge <- c(regular_cols, setdiff(excluded_cols, special_merge_cols))
    
    # Merge regular columns (including excluded columns that don't need special merge)
    for(col in cols_for_regular_merge) {
      valid_values <- group[[col]][!sapply(group[[col]], is_na_or_empty_numeric, col_name = col)]
      if(length(valid_values) > 0) {
        merged_row[[col]] <- valid_values[1]
      } else {
        merged_row[[col]] <- get_empty_value(col)
      }
    }
    
    # Handle special merge columns with append_status
    for(col in special_merge_cols) {
      if(col %in% names(group)) {
        merged_val <- group[[col]][1]
        for(i in 2:nrow(group)) {
          merged_val <- append_status(merged_val, group[[col]][i])
        }
        merged_row[[col]] <- merged_val
      }
    }
    
    banc_meta_merged_dup <- bind_rows(banc_meta_merged_dup, merged_row)
  }
}

# summarize merge conflicts
merge_conflict_summary <- merge_conflict_root_id %>%
  group_by(conflict_cols) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(desc(count))


# generate rows for duplicates with conflicts
# Helper functions
count_non_empty <- function(row_data, all_cols) {
  count <- 0
  for(col in all_cols) {
    if(!is_na_or_empty_numeric(row_data[[col]], col)) {
      count <- count + 1
    }
  }
  return(count)
}

count_status_elements <- function(status_val) {
  if(is.na(status_val) || status_val == "") {
    return(0)
  }
  return(length(str_split(status_val, ",")[[1]]))
}

merge_notes <- function(notes_vector) {
  valid_notes <- notes_vector[!is_na_or_empty_numeric(notes_vector, "notes")]
  
  if(length(valid_notes) == 0) {
    return(NA)
  } else if(length(valid_notes) == 1) {
    return(valid_notes[1])
  } else {
    unique_notes <- unique(valid_notes)
    return(paste(unique_notes, collapse = "\n"))
  }
}


# Process conflict root_ids
conflict_root_ids <- merge_conflict_root_id$root_id

# pre-process conflict rows of banc_meta for resolution
# conflict_data <- banc_meta_dup_join %>%
#   rowwise() %>%
#   mutate(status = subtract_status(status, "DUPLICATED")) %>%
#   mutate(status_text = subtract_status(status_text, "DUPLICATED")) %>%
#   ungroup()
# 
# 
# conflict_data <- banc_meta_dup %>% 
#   filter(root_id %in% conflict_root_ids) %>%
#   rowwise() %>%
#   mutate(status = subtract_status(status, "DUPLICATED")) %>%
#   mutate(status_text = subtract_status(status_text, "DUPLICATED")) %>%
#   ungroup()

conflict_data <- dup_data %>% 
  filter(root_id %in% conflict_root_ids) %>%
  rowwise() %>%
  mutate(status = subtract_status(status, "DUPLICATED")) %>%
  mutate(status_text = subtract_status(status_text, "DUPLICATED")) %>%
  ungroup()

banc_meta_conflict_sel <- data.frame()

# Process each group
conflict_groups <- conflict_data %>%
  group_by(root_id) %>%
  group_split()

for(group in conflict_groups) {
  # Get column categories
  all_cols <- names(group)
  excluded_cols <- intersect(all_cols, exclude_cols)
  cols_for_regular_merge <- c(
    setdiff(all_cols, exclude_cols),
    setdiff(excluded_cols, c(special_merge_cols, notes_merge_cols))
  )
  
  # Identify conflict vs non-conflict columns
  conflict_cols <- c()
  non_conflict_cols <- c()
  
  for(col in cols_for_regular_merge) {
    valid_values <- group[[col]][!sapply(group[[col]], is_na_or_empty_numeric, col_name = col)]
    if(length(unique(valid_values)) > 1) {
      conflict_cols <- c(conflict_cols, col)
    } else {
      non_conflict_cols <- c(non_conflict_cols, col)
    }
  }
  
  # Find best row for conflict resolution
  best_row_idx <- 1
  if(length(conflict_cols) > 0) {
    row_scores <- data.frame(
      row_idx = 1:nrow(group),
      non_empty_count = sapply(1:nrow(group), function(i) count_non_empty(group[i, ], all_cols)),
      status_count = sapply(1:nrow(group), function(i) {
        if("status" %in% names(group)) {
          count_status_elements(group$status[i])
        } else {
          0
        }
      })
    )
    
    best_row <- row_scores %>%
      arrange(desc(non_empty_count), desc(status_count), row_idx) %>%
      slice(1)
    
    best_row_idx <- best_row$row_idx
  }
  
  # Create merged row
  merged_row <- group[1, ]
  
  # Non-conflict columns: regular merge
  for(col in non_conflict_cols) {
    valid_values <- group[[col]][!sapply(group[[col]], is_na_or_empty_numeric, col_name = col)]
    if(length(valid_values) > 0) {
      merged_row[[col]] <- valid_values[1]
    } else {
      merged_row[[col]] <- get_empty_value(col)
    }
  }
  
  # Conflict columns: use best row
  for(col in conflict_cols) {
    merged_row[[col]] <- group[[col]][best_row_idx]
  }
  
  # Special merge columns (status, cell_type_source)
  for(col in special_merge_cols) {
    if(col %in% names(group)) {
      merged_val <- group[[col]][1]
      for(i in 2:nrow(group)) {
        merged_val <- append_status(merged_val, group[[col]][i])
      }
      merged_row[[col]] <- merged_val
    }
  }
  
  # Notes merge columns (concatenate with line breaks)
  for(col in notes_merge_cols) {
    if(col %in% names(group)) {
      merged_row[[col]] <- merge_notes(group[[col]])
    }
  }
  
  banc_meta_conflict_sel <- bind_rows(banc_meta_conflict_sel, merged_row)
}


# formerly duplicated rows, for upload
banc_meta_merged_upload <- bind_rows(banc_meta_merged_dup, banc_meta_conflict_sel)

# bug with SeaTable upload with numeric columns
# banc_meta_merged_upload_num <- banc_meta_merged_upload %>%
#   mutate(across(all_of(numeric_cols), ~ ifelse(is.na(.x), get_empty_value(cur_column()), .x)))

# duplicated root_id rows for deletion - flag with status tag DELETE_THIS_DUPLICATE
banc_meta_duplicate_delete_stage <- dup_data %>%
  filter(!(`_id` %in% banc_meta_merged_upload$`_id`)) %>%
  rowwise() %>%
  mutate(status = append_status(status, "DELETE_THIS_DUPLICATE")) %>%
  mutate(status_text = append_status(status_text, "DELETE_THIS_DUPLICATE")) %>%
  ungroup() %>%
  select(status, status_text, `_id`)

# modify banc_meta with both of these changes
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_merged_upload, table='banc_meta',base='banc_meta')
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_duplicate_delete_stage, table='banc_meta',base='banc_meta')

# because of bug in append_status, NA pushed to cell_type_source and other_names
#  - correct this
banc_meta_cts_wNA <- banc_meta %>%
  filter(grepl("NA", cell_type_source)) %>%
  rowwise() %>%
  mutate(cell_type_source = subtract_status(cell_type_source, "NA")) %>%
  ungroup() %>%
  mutate(cell_type_source = case_when(
    cell_type_source == "" ~ NA,
    TRUE ~ cell_type_source
  )) %>%
  select(cell_type_source, malecns_match, `_id`)

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_cts_wNA, table='banc_meta',base='banc_meta')

banc_meta_otherNames_wNA <- banc_meta %>%
  filter(grepl("NA", other_names)) %>%
  rowwise() %>%
  mutate(other_names = subtract_status(other_names, "NA")) %>%
  ungroup() %>%
  mutate(other_names = case_when(
    other_names == "" ~ NA,
    TRUE ~ other_names
  )) %>%
  select(other_names, malecns_match, `_id`)  

# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_otherNames_wNA, table='banc_meta',base='banc_meta')
           
