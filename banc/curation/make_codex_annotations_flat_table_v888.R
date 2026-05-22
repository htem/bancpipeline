# generate a flat table for codex

library(arrow)
source("banc/curation/banc_curation_init.R")

# get banc_meta
banc_meta <- banctable_query()

banc_meta_newIDs_filt <- banc_meta %>%
  filter(is_not_empty_na(cell_ids_tag))

# load in from CSVs
crp_filepath <- "cell_representative_point_nm_wSV.csv"

cell_rep_pt <- read_csv(crp_filepath, col_types = cols(
  .default       = col_character(),
  valid          = col_logical(),
  supervoxel_id  = col_character()   # read as string first to avoid scientific notation
)) %>%
  mutate(supervoxel_id = as.integer64(supervoxel_id))

# get position in voxels
positions_to_convert <- cell_rep_pt$pt

coords_raw <- banc_nm2raw(positions_to_convert)

position_raw <- paste(
  format(coords_raw[, "X"], scientific = FALSE, trim = TRUE),
  format(coords_raw[, "Y"], scientific = FALSE, trim = TRUE),
  format(coords_raw[, "Z"], scientific = FALSE, trim = TRUE),
  sep = ", "
)

cell_rep_pt$position <- position_raw

cell_rep_pt <- cell_rep_pt %>%
  select(-valid)

# pull v888 root_id from cell rep pt parquet
crp_filepath <- "neuron_annotations_v888_cell_representative_point.parquet"
crp <- read_parquet(crp_filepath)

root_888 <- crp %>%
  select(pt_root_id, id)

cell_rep_pt_id <- cell_rep_pt %>%
  mutate(id = as.integer(id)) %>%
  left_join(root_888, by = "id")

cell_rep_pt_id_mod <- cell_rep_pt_id %>%
  mutate(cell_ids_id = as.integer(id)) %>%
  rename(position_nm = pt)

# join with banc_meta
banc_meta_newIDs_filt <- banc_meta_newIDs_filt %>%
  select(-supervoxel_id, -position, -root_id) %>%
  left_join(cell_rep_pt_id_mod, by = "cell_ids_id")


test_dup_root_id <- crp %>%
  group_by(pt_root_id) %>%
  summarize(count = n(), .groups = "drop") %>%
  filter(count > 1)

###

# create a user_id column of banc_meta from cell_type_source and cell_info_user_id
# Define the mapping
source_to_user_id <- c(
  "ache_lab" = 1776, "aelysia" = 4462, "agrawal_lab" = 3953, "azevedo" = 188,
  "bidaye_lab" = 1238, "dacks_lab" = 1163, "daly_lab" = 1234, "dickerson_lab" = 1044,
  "eichler_lab" = 92, "jefferis_lab" = 60, "lee_lab" = 158, "pankratz_lab" = 1152,
  "seeds_hampel_lab" = 125, "seung_lab" = 1081, "suver_lab" = 3910, "wilson_lab" = 355,
  "zandawala_lab" = 3153, "tomke" = 847, "yao_lab" = 1205, "huston_lab" = 52, 
  "bates" = 355, "yang" = 52, "jasper" = 2660
)

# Vectorized function to create user_id column
create_user_id_vectorized <- function(df) {
  df %>%
    mutate(
      user_id = pmap_chr(list(cell_type_source, cell_info_user_id), function(cts, ciu) {
        all_ids <- c()
        
        # Process cell_type_source
        if(!is.na(cts) && cts != "") {
          sources <- trimws(unlist(strsplit(cts, ",")))
          sources <- sources[sources != ""]
          source_ids <- as.character(source_to_user_id[sources[sources %in% names(source_to_user_id)]])
          all_ids <- c(all_ids, source_ids)
        }
        
        # Add cell_info_user_id
        if(!is.na(ciu) && ciu != "") {
          all_ids <- c(all_ids, as.character(ciu))
        }
        
        # Return unique IDs
        if(length(all_ids) > 0) {
          paste(unique(all_ids), collapse = ",")
        } else {
          NA_character_
        }
      })
    )
}

# Create user_id column
banc_meta_newIDs_filt <- create_user_id_vectorized(banc_meta_newIDs_filt)

test<-banc_meta_newIDs_filt %>%
  arrange(cell_ids_id) %>%
  select(pt_root_id, position, position_nm)

# get banc_meta to join to cell_rep_pt
cell_rep_pt_bm <- banc_meta_newIDs_filt %>%
  select(pt_root_id, cell_ids_id, position, position_nm, supervoxel_id, flow, 
         super_class, cell_class, cell_sub_class, cell_type, 
         region, side, cell_function, cell_function_detailed, peripheral_target_type, 
         body_part_sensory, body_part_effector, nerve, hemilineage, sexually_dimorphic,
         neurotransmitter_verified, neuropeptide_verified, neurotransmitter_predicted, 
         neurotransmitter_score, fafb_match, fafb_cell_type, manc_match, manc_cell_type, 
         hemibrain_match, hemibrain_cell_type, fanc_match, fanc_cell_type, 
         malecns_match, malecns_cell_type, fafb_alignment_cell_type,
         other_names,  user_id) 

# generate CSV
codex_table_flat <- cell_rep_pt_bm %>%
  mutate(across(
    c(fafb_cell_type, manc_cell_type, hemibrain_cell_type, fanc_cell_type, malecns_cell_type),
    ~ if_else(str_detect(.x, "auto:"), NA_character_, .x)
  )) %>%
  rename(target_id = cell_ids_id) %>%
  rename(fafb_783_match_id = fafb_match) %>%
  rename(fafb_783_cell_type = fafb_cell_type) %>%
  rename(manc_121_match_id = manc_match) %>%
  rename(manc_121_cell_type = manc_cell_type) %>%
  rename(fanc_1116_match_id = fanc_match) %>%
  rename(fanc_1116_cell_type = fanc_cell_type) %>%
  rename(hemibrain_121_match_id = hemibrain_match) %>%
  rename(hemibrain_121_cell_type = hemibrain_cell_type) %>%
  rename(malecns_09_match_id = malecns_match) %>%
  rename(malecns_09_cell_type = malecns_cell_type) %>%
  rename(alignment_cell_type = fafb_alignment_cell_type) %>%
  arrange(target_id)


# write to CSV
codex_table_filepath <- "codex_annotations_flat_table.csv"
# Temporarily disable scientific notation
old_scipen <- getOption("scipen")
options(scipen = 999)


# WARNING: script not live, uncomment to run
# Write the CSV
# write.csv(codex_table_flat, codex_table_filepath,
#           row.names = FALSE,    # Exclude row names
#           quote = TRUE)

# Restore original scipen setting
options(scipen = old_scipen)




###
# make the CAVE version too (in old format)

# Define columns
comma_split_cols <- c("peripheral_target_type", "body_part_sensory", "body_part_effector", 
                      "nerve", "neurotransmitter_verified", "neuropeptide_verified", 
                      "fafb_783_match_id", "manc_121_match_id", "other_names", 
                      "cell_function", "cell_function_detailed", "user_id")

all_cols <- c("region", 
              "side", 
              "flow", 
              "super_class", 
              "cell_class", 
              "cell_sub_class", 
              "cell_type", 
              "cell_function", 
              "cell_function_detailed", 
              "peripheral_target_type", 
              "body_part_sensory", 
              "body_part_effector", 
              "nerve", 
              "hemilineage", 
              "sexually_dimorphic",
              "neurotransmitter_verified", 
              "neuropeptide_verified", 
              "neurotransmitter_predicted",
              "neurotransmitter_score",
              "other_names",
              "fafb_783_cell_type", 
              "manc_121_cell_type", 
              "hemibrain_121_cell_type", 
              "fanc_1116_cell_type",
              "malecns_09_cell_type",
              "alignment_cell_type",
              "fafb_783_match_id", 
              "manc_121_match_id", 
              "hemibrain_121_match_id",
              "fanc_1116_match_id",
              "malecns_09_match_id",
              "user_id")

# Create codex_annotations using map_dfr
codex_table_flat <- codex_table_flat %>%
  rename(id = target_id)
codex_annotations <- map_dfr(all_cols, function(col) {
  if(col %in% names(codex_table_flat)) {
    base_df <- codex_table_flat %>%
      select(id, !!sym(col)) %>%
      filter(!is.na(!!sym(col)) & !!sym(col) != "")
    
    if(nrow(base_df) > 0) {
      if(col %in% comma_split_cols) {
        base_df %>%
          mutate(!!sym(col) := as.character(!!sym(col))) %>%
          separate_rows(!!sym(col), sep = ",\\s*") %>%
          mutate(!!sym(col) := trimws(!!sym(col))) %>%
          filter(!!sym(col) != "" & !is.na(!!sym(col))) %>%
          transmute(
            valid = TRUE,
            target_id = id,
            classification_system = col,
            cell_type = as.character(!!sym(col))
          )
      } else {
        base_df %>%
          transmute(
            valid = TRUE,
            target_id = id,
            classification_system = col,
            cell_type = as.character(!!sym(col))
          )
      }
    } else {
      tibble()
    }
  } else {
    tibble()
  }
})

# add id column to codex_annotations
codex_annotations_wId <- codex_annotations %>%
  mutate(id = seq(0, nrow(codex_annotations) - 1)) %>%
  select(id, valid, target_id, classification_system, cell_type)

# Save to CSV files
refTable_csv_filepath <- ""

# WARNING: script not live, uncomment to run
# Save codex_annotations
codex_path <- file.path(refTable_csv_filepath, "codex_annotations.csv")
# write.csv(codex_annotations_wId, codex_path, row.names = FALSE)
# cat("Saved codex_annotations to", codex_path, "\n")
