# script to generate cell_representative_point CAVE table
#

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


banc_version <- 882

# get CAVE tables
backbone_proofread <- banc_cave_query("backbone_proofread", live = FALSE, 
                                      version = as.integer(banc_version))
cell_info <- banc_cave_query("cell_info", live = FALSE, 
                             version = as.integer(banc_version))
proofreading_notes <- banc_cave_query("proofreading_notes", live = FALSE, 
                                      version = as.integer(banc_version))
peripheral_nerves <- banc_cave_query("peripheral_nerves", live = FALSE, 
                                     version = as.integer(banc_version))
neck_connective_y92500 <- banc_cave_query("neck_connective_y92500", live = FALSE, 
                                          version = as.integer(banc_version))
neck_connective_y121000 <- banc_cave_query("neck_connective_y121000", live = FALSE, 
                                           version = as.integer(banc_version))
# somas CAVE tables
somas_v1a <- banc_cave_query(table = "somas_v1a", live = FALSE, 
                             version = as.integer(banc_version))
somas_v1b <- banc_cave_query(table = "somas_v1b", live = FALSE, 
                             version = as.integer(banc_version))

# correct points in somas_v1a with somas_v1b
somas_v1b_mod <- somas_v1b %>%
  select(id, pt_supervoxel_id, pt_root_id, pt_position, created, valid)

somas_v1a_mod <- somas_v1a %>%
  select(id, pt_supervoxel_id, pt_root_id, pt_position, created, valid) %>%
  filter(!(id %in% somas_v1b_mod$id))

somas_all <- bind_rows(somas_v1a_mod, somas_v1b_mod)


# check if there are any entries since the last materialization
backbone_proofread_live <- banc_cave_query("backbone_proofread", live = 2)
cell_info_live <- banc_cave_query("cell_info", live = 2)
proofreading_notes_live <- banc_cave_query("proofreading_notes", live = 2)

# append in extra rows
backbone_proofread_live_miss <- backbone_proofread_live %>%
  filter(created > max(backbone_proofread$created, na.rm = TRUE))
backbone_proofread_all <- bind_rows(backbone_proofread, backbone_proofread_live_miss)

cell_info_live_miss <- cell_info_live %>%
  filter(created > max(cell_info$created, na.rm = TRUE))
cell_info_all <- bind_rows(cell_info, cell_info_live_miss)

proofreading_notes_live_miss <- proofreading_notes_live %>%
  filter(created > max(proofreading_notes$created, na.rm = TRUE))
proofreading_notes_all <- bind_rows(proofreading_notes, proofreading_notes_live_miss)


# only valid IDs from each table
backbone_proofread_all <- backbone_proofread_all %>%
  filter((valid == "TRUE"))
cell_info_ct_all <- cell_info_all %>%
  filter((valid == "TRUE")) %>%
  filter(tag %in% franken_meta_ct$cell_type)
roughly_proofread <- proofreading_notes_all %>%
  filter((valid == "TRUE") & (tag == "roughly proofread"))
peripheral_nerves <- peripheral_nerves %>%
  filter((valid == "TRUE"))
somas_all <- somas_all %>%
  filter((valid == "TRUE"))
neck_connective_y92500 <- neck_connective_y92500 %>%
  filter((valid == "TRUE"))
neck_connective_y121000 <- neck_connective_y121000 %>%
  filter((valid == "TRUE"))


# only unique root_ids; if duplicated, pick the most recent one; 
# also, modify types
backbone_proofread_uni <- backbone_proofread_all %>%
  group_by(pt_root_id) %>%
  slice_max(created, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  rename(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_backbone = pt_position) %>%
  mutate(supervoxel_backbone = pt_supervoxel_id) %>%
  mutate(position_backbone = standardize_position(position_backbone)) %>%
  select(root_id, position_backbone, supervoxel_backbone)

# only unique root_ids; if duplicated, pick the one that matches the cell type 
#  in banc_meta; otherwise, pick the most recent one
cell_info_ct_uni <- cell_info_ct_all %>%
  mutate(pt_root_id_char = as.character(pt_root_id)) %>%
  left_join(banc_meta %>% select(root_id, cell_type),
            by = c("pt_root_id_char" = "root_id")) %>%
  mutate(
    type_match_priority = case_when(
      !is.na(cell_type) & cell_type != "" &
        !is.na(tag) & tag != "" &
        tag == cell_type ~ 1,
      TRUE ~ 2
    )
  ) %>%
  group_by(pt_root_id) %>%
  arrange(type_match_priority, desc(created), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-pt_root_id_char, -cell_type, -type_match_priority) %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_cell_info = pt_position) %>%
  mutate(supervoxel_cell_info = pt_supervoxel_id) %>%
  mutate(position_cell_info = standardize_position(position_cell_info)) %>%
  select(root_id, position_cell_info, supervoxel_cell_info)

roughly_proofread_uni <- roughly_proofread %>%
  group_by(pt_root_id) %>%
  slice_max(created, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_rough_proof = pt_position) %>%
  mutate(supervoxel_rough_proof = pt_supervoxel_id) %>%
  mutate(position_rough_proof = standardize_position(position_rough_proof)) %>%
  select(root_id, position_rough_proof, supervoxel_rough_proof)

peripheral_nerves_uni <- peripheral_nerves %>%
  group_by(pt_root_id) %>%
  slice_max(created, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_nerves = pt_position) %>%
  mutate(supervoxel_nerves = pt_supervoxel_id) %>%
  mutate(position_nerves = standardize_position(position_nerves)) %>%
  select(root_id, position_nerves, supervoxel_nerves)

neck_connective_y92500_id <- neck_connective_y92500 %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_neck_y92500 = pt_position) %>%
  mutate(supervoxel_neck_y92500 = pt_supervoxel_id) %>%
  mutate(position_neck_y92500 = standardize_position(position_neck_y92500)) %>%
  select(root_id, position_neck_y92500, supervoxel_neck_y92500)

neck_connective_y121000_id <- neck_connective_y121000 %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_neck_y12100 = pt_position) %>%
  mutate(supervoxel_neck_y12100 = pt_supervoxel_id) %>%
  mutate(position_neck_y12100 = standardize_position(position_neck_y12100)) %>%  
  select(root_id, position_neck_y12100, supervoxel_neck_y12100)

# temp: remove not real nuclei (will be removed from CAVE table later)
not_real_nuclei_id <- c("73044200840495136",
                        "73887387813413131",
                        "72833575912735349",
                        "72831925571551381",
                        "72761557834007000",
                        "73113607042238803",
                        "73042964091241586",
                        "73254139110360940",
                        "73114706822299761",
                        "73113264048832769",
                        "72904012839518677",
                        "72972663798104282",
                        "72832063480267234",
                        "72902432224444677",
                        "72972594541756546",
                        "72973831559446903", 
                        "72973831559446932",
                        "73324576506906241",
                        "72903668973699481",
                        "72763344607511085",
                        "73114706486755501",
                        "72973763578167741")

somas_all_id <- somas_all %>%
  mutate(across(where(~ inherits(.x, "list")), ~ sapply(.x, paste, collapse = ","))) %>%
  mutate(root_id = pt_root_id) %>%
  mutate(across(where(~ inherits(.x, "integer64")), as.character)) %>%
  mutate(position_soma = pt_position) %>%
  mutate(supervoxel_soma = pt_supervoxel_id) %>%
  mutate(position_soma = standardize_position(position_soma)) %>%  
  filter(!(id %in% not_real_nuclei_id)) %>%
  select(root_id, position_soma, supervoxel_soma)





# update banc_meta with latest in what cells are proofread
banc_meta_proof <- banc_meta %>%
  mutate(proofread = case_when(
    root_id %in% backbone_proofread_uni$root_id ~ "TRUE",
    TRUE ~ proofread
  )) %>%
  mutate(roughly_proofread = case_when(
    root_id %in% roughly_proofread_uni$root_id ~ "TRUE",
    TRUE ~ roughly_proofread
  ))


# generate a version of banc_meta that is only the cells to be included:
# backbone_proofread, roughly_proofread, has cell_type or match to another 
#  dataset, status tag "INCLUDE_THIS_DUPLICATE", cell_class == "astrocyte"
#  status tag "HAS_MANUAL_ANNOTATION", cell_type_source is not empty and is not 
#  just wilson_lab, or has cell_type from cell_info
# NOT glia (unless astrocyte)
banc_meta_valid_cells <- banc_meta_proof %>%
  filter(((proofread == "TRUE") & !(grepl("glia", super_class))) | 
           ((roughly_proofread == "TRUE") & !(grepl("glia", super_class))) |
           grepl("astrocyte", cell_class) |
           grepl("INCLUDE_THIS_DUPLICATE", status) |
           grepl("HAS_MANUAL_ANNOTATION", status) |
           is_not_empty_na(cell_type) | 
           is_not_empty_na(fafb_match) |
           is_not_empty_na(manc_match) |
           is_not_empty_na(fanc_match) |
           is_not_empty_na(malecns_match) |
           is_not_empty_na(hemibrain_match) |
           is_not_empty_na(cell_info_id) |
           (is_not_empty_na(cell_type_source) & (cell_type_source != "wilson_lab")))


# select only necessary columns
banc_meta_mod <- banc_meta_valid_cells %>%
  rename(position_banc_meta = position) %>%
  rename(supervoxel_banc_meta = supervoxel_id) %>%
  mutate(cell_ids_tag = NA) %>%
  select(`_id`, root_id, position_banc_meta, supervoxel_banc_meta, 
         proofread, roughly_proofread,
         status,
         nucleus_position, nucleus_supervoxel_id,
         cell_ids_position, cell_ids_position_nm, cell_ids_tag)

# join with other sources of position info
banc_meta_join <- banc_meta_mod %>%
  left_join(backbone_proofread_uni, by = "root_id") %>%
  left_join(roughly_proofread_uni, by = "root_id") %>%
  left_join(cell_info_ct_uni, by = "root_id") %>%
  left_join(peripheral_nerves_uni, by = "root_id")


# pick the appropriate position and cell_ids_supervoxel_id for cell_ids_position 
#  based on the following rules, in order:
# set which_source as:
# if status contains INCLUDE_THIS_DUPLICATE, use banc_meta
# if proofread == TRUE AND 
#  position_backbone isn't nucleus_position AND
#  AND supervoxel_backbone isn't nucleus_supervoxel_id AND
#  position_backbone isn't in somas_all_id$position_soma AND
#  supervoxel_backbone isn't in somas_all_id$supervoxel_soma
#   which_source == "backbone"
# if position_cell_info is not empty AND 
#  position_cell_info isn't nucleus_position AND
#  supervoxel_cell_info isn't nucleus_supervoxel_id AND
#  position_cell_info isn't in somas_all_id$position_soma AND
#  supervoxel_cell_info isn't in somas_all_id$supervoxel_soma
#   which_source == "cell_info"
# if position_banc_meta is not empty AND 
#  position_banc_meta isn't nucleus_position AND
#  supervoxel_banc_meta isn't nucleus_supervoxel_id AND
#  position_banc_meta isn't in somas_all_id$position_soma AND
#  supervoxel_banc_meta isn't in somas_all_id$supervoxel_soma
#   which_source == "banc_meta"
# if position_nerves is not empty
#   which_source == "nerves"
# if proofread != TRUE AND roughly_proofread == TRUE AND
#  position_rough_proof isn't nucleus_position AND
#  supervoxel_rough_proof isn't nucleus_supervoxel_id AND
#  position_rough_proof isn't in somas_all_id$position_soma AND
#  supervoxel_rough_proof isn't in somas_all_id$supervoxel_soma
#   which_source == "rough_proof"
# if proofread == TRUE,
#   which_source == "backbone"

banc_meta_join_merge <- banc_meta_join %>%
  mutate(which_source = case_when(
    (grepl("INCLUDE_THIS_DUPLICATE", status) ~ "banc_meta"),
    ((proofread == "TRUE") & 
       is_not_empty_na(position_backbone) &
       (is.na(nucleus_position) | position_backbone != nucleus_position) &
       (is.na(nucleus_supervoxel_id) | supervoxel_backbone != nucleus_supervoxel_id) &
       (!(position_backbone %in% somas_all_id$position_soma)) &
       (!(supervoxel_backbone %in% somas_all_id$supervoxel_soma))
    ) ~ "backbone",
    ((is_not_empty_na(position_cell_info)) & 
       (is.na(nucleus_position) | position_cell_info != nucleus_position) &
       (is.na(nucleus_supervoxel_id) | supervoxel_cell_info != nucleus_supervoxel_id) &
       (!(position_cell_info %in% somas_all_id$position_soma)) &
       (!(supervoxel_cell_info %in% somas_all_id$supervoxel_soma))
    ) ~ "cell_info",
    ((is_not_empty_na(position_banc_meta)) & 
       (is.na(nucleus_position) | position_banc_meta != nucleus_position) &
       (is.na(nucleus_supervoxel_id) | supervoxel_banc_meta != nucleus_supervoxel_id) &
       (!(position_banc_meta %in% somas_all_id$position_soma)) &
       (!(supervoxel_banc_meta %in% somas_all_id$supervoxel_soma))
    ) ~ "banc_meta", 
    (is_not_empty_na(position_nerves)) ~ "nerves",
    ((proofread != "TRUE") & (roughly_proofread == "TRUE") & 
       is_not_empty_na(position_rough_proof) &
       (is.na(nucleus_position) | position_rough_proof != nucleus_position) &
       (is.na(nucleus_supervoxel_id) | supervoxel_rough_proof != nucleus_supervoxel_id) &
       (!(position_rough_proof %in% somas_all_id$position_soma)) &
       (!(supervoxel_rough_proof %in% somas_all_id$supervoxel_soma))
    ) ~ "rough_proof",  
    ((proofread == "TRUE") & is_not_empty_na(position_backbone)) ~ "backbone",
    ((roughly_proofread == "TRUE") & is_not_empty_na(position_rough_proof)) ~ "rough_proof",
    (is_not_empty_na(position_banc_meta)) ~ "banc_meta",
    TRUE ~ NA
  ))

# check for missing positions
missing_source <- banc_meta_join_merge %>%
  filter(is_na_or_empty(which_source))


# get the appropriate position and supervoxel
banc_meta_join_merge_pos <- banc_meta_join_merge %>%
  mutate(
    cell_ids_position = case_when(
      which_source == "backbone"   ~ position_backbone,
      which_source == "cell_info"  ~ position_cell_info,
      which_source == "banc_meta"  ~ position_banc_meta,
      which_source == "nerves"     ~ position_nerves,
      which_source == "rough_proof" ~ position_rough_proof,
      TRUE ~ NA_character_
    ),
    cell_ids_supervoxel_id = case_when(
      which_source == "backbone"   ~ supervoxel_backbone,
      which_source == "cell_info"  ~ supervoxel_cell_info,
      which_source == "banc_meta"  ~ supervoxel_banc_meta,
      which_source == "nerves"     ~ supervoxel_nerves,
      which_source == "rough_proof" ~ supervoxel_rough_proof,
      TRUE ~ NA_character_
    )
  )
  
# check that row got a position and supervoxel_id
missing_val <- banc_meta_join_merge_pos %>%
  filter(is_na_or_empty(cell_ids_position) | is_na_or_empty(cell_ids_supervoxel_id))



# get cell_ids_position_nm
non_na_mask <- !is.na(banc_meta_join_merge_pos$cell_ids_position) & 
  banc_meta_join_merge_pos$cell_ids_position != ""

positions_to_convert <- banc_meta_join_merge_pos$cell_ids_position[non_na_mask]

# Pass character vector directly - xyzmatrix handles the parsing
coords_nm <- banc_raw2nm(positions_to_convert)

position_nm <- paste(
  format(coords_nm[, "X"], scientific = FALSE, trim = TRUE),
  format(coords_nm[, "Y"], scientific = FALSE, trim = TRUE),
  format(coords_nm[, "Z"], scientific = FALSE, trim = TRUE),
  sep = ", "
)

banc_meta_join_merge_pos$cell_ids_position_nm <- NA_character_
banc_meta_join_merge_pos$cell_ids_position_nm[non_na_mask] <- position_nm

# get cell_ids_tag and cell_ids_id
banc_meta_join_merge_pos <- banc_meta_join_merge_pos %>%
  mutate(
    # Generate unique 6-digit integer IDs starting from 0
    cell_ids_tag = sprintf("%06d", seq(0, 0 + n() - 1))
  ) %>%
  # convert tag to id
  mutate(cell_ids_id = as.integer(cell_ids_tag))


# update position and supervoxel_id in main banc_meta columns to align with
#  cell_ids_position and cell_ids_supervoxel_id
# also, select columns
banc_meta_final_update <- banc_meta_join_merge_pos %>%
  mutate(position = cell_ids_position) %>%
  mutate(supervoxel_id = cell_ids_supervoxel_id) %>%
  select(`_id`, supervoxel_id, position, cell_ids_position, cell_ids_position_nm,
         cell_ids_supervoxel_id, cell_ids_tag, cell_ids_id)

# update banc_meta
# WARNING: script not live, uncomment to run
# banctable_update_rows(banc_meta_final_update, table='banc_meta',base='banc_meta')




# generate cell_representative_point table

banc_meta_final_update_sv <- banc_meta_final_update %>%
  mutate(supervoxel_id_int = as.integer64.character(cell_ids_supervoxel_id))

# no scientific notation
options(scipen = 999)

# Create cell_ids dataframe - in nm
cell_ids <- banc_meta_final_update_sv %>%
  transmute(
    id = format(cell_ids_id, scientific = FALSE, trim = TRUE),
    valid = TRUE,
    pt = cell_ids_position_nm,
    supervoxel_id = supervoxel_id_int
  )

# Save to CSV files
refTable_csv_filepath <- ""

# Save cell_ids separately
cell_ids_path <- file.path(refTable_csv_filepath, "cell_representative_point_nm_wSV.csv")
write.csv(cell_ids, cell_ids_path, row.names = FALSE)
cat("Saved cell_ids to", cell_ids_path, "\n")

