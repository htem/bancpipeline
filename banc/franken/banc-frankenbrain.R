#' banc-frankenbrain — Build the franken-brain meta CSV (FAFB + MANC merge).
#'
#' Merges FlyWire (FAFB) and MANC metadata into a single cross-dataset
#' "franken" metadata table, using BANC-mediated cell_type bridges where
#' available. The output is the canonical reference loaded everywhere via
#' `franken_meta()`.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/flywire_meta.csv` (FAFB)
#'   - `<banc.meta.save.path>/manc_meta.csv` (MANC)
#'   - SeaTable `banc_meta` (cell_type bridges)
#'
#' @section Writes:
#'   - `<banc.meta.save.path>/frankenbrain_v<X.Y>_meta.csv`
#'
#' @section Notes:
#'   - The MANC mapping treats sensory_neuron / intrinsic_neuron / VNC
#'     endocrine as `efferent`. FAFB sensory_descending is mapped the same.

###########################
### BUILD FRANKEN-BRAIN ###
###########################
source("banc/banc-startup.R")

# # Run locally
# banc.meta.save.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
# banc.connectivity.save.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity"
# banc.dropbox.connectivity.save.path <- "/Users/abates/HMS Dropbox/Alexander Bates/neuroanat/connectomes"
# banc.path <- "/Users/papers/BANC-project/data/meta"

# manc: sensory_neuron super_class, intrinsic_neuron, vnc endorcine need to be in efferent
# fafb: sensory_descending

############
### Main ###
############

# Read meta
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))

# Get BANC-mediated cell type level mapping
bancr:::banctable_updateids()
bc <- banctable_query() 
bc.meta.neck <- bc %>%
  dplyr::filter(!grepl("DELETE",status),!grepl("NOT_A_NEURON",status),!grepl("GLIA",status),
                !grepl("MERGE_ERROR",super_class), !grepl("efferent|motor",super_class)) %>%
  dplyr::mutate(cell_class = gsub("auto:","",cell_class),
                super_class = gsub("auto:","",super_class),
                cell_type = ifelse(grepl("auto:",cell_type),NA,cell_type),
                fafb_cell_type = ifelse(grepl("auto:",fafb_cell_type),NA,fafb_cell_type),
                manc_cell_type = ifelse(grepl("auto:",manc_cell_type),NA,manc_cell_type)) %>%
  dplyr::filter(cell_class %in% c("^AN$","^DN$","^sensory ascending$", "^descending$", "^ascending$")|
                super_class %in% c("^AN$","^DN$","^sensory ascending$", "^descending$", "^ascending$", "sensroy_ascending", "sensory_descending")|
                grepl("neck",region)) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    !is.na(an_type) ~ paste0("AN_",an_type),
    !is.na(dn_type) ~ paste0("DN_",dn_type),
    TRUE ~ cell_sub_class
  )) %>% dplyr::mutate(
    super_class = dplyr::case_when(
      super_class=="sensory_descending" ~ "descending",
      TRUE ~ super_class
    )) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
bc.meta.neck$fafb_cell_type <- fw.meta$cell_type[match(bc.meta.neck$fafb_match,fw.meta$root_783)]
bc.meta.neck$manc_cell_type <- mc.meta$cell_type[match(bc.meta.neck$manc_match,mc.meta$bodyid)]
bc.fafb.matches.df.valid <- bc.meta.neck %>%
  dplyr::distinct(pt_root_id = root_id,
                  pt_supervoxel_id = supervoxel_id,
                  pt_position = position,
                  query_root_id =  root_id,
                  match_id = fafb_match,
                  match_cell_type = fafb_cell_type) %>%
  dplyr::filter(!is.na(pt_root_id),!is.na(match_id)) %>%
  dplyr::mutate(validation = TRUE, score = 1)
bc.manc.matches.df.valid <- bc.meta.neck %>%
  dplyr::distinct(pt_root_id = root_id,
                  pt_supervoxel_id = supervoxel_id,
                  pt_position = position,
                  query_root_id =  root_id,
                  match_id = manc_match,
                  match_cell_type = manc_cell_type) %>%
  dplyr::filter(!is.na(pt_root_id),!is.na(match_id)) %>%
  dplyr::mutate(validation = TRUE, score = 1)

# Matches adjusted for side
fafb.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE) %>%
  dplyr::mutate(score=0.9) %>%
  dplyr::filter(valid=='t') %>%
  dplyr::mutate(validation = TRUE) %>%
  dplyr::select(-valid) %>%
  dplyr::rename_with(~sub("^query_id$", "query_root_id", .x))
fafb.matches.df.valid <- rbind(bc.fafb.matches.df.valid,fafb.matches.df.valid) %>%
  dplyr::distinct(pt_root_id,
                  match_id,
                  .keep_all=TRUE)
banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather")) %>%
  rbind(fafb.matches.df.valid) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(!grepl("\\_",match_id), score > 0) %>%
  dplyr::rename(fafb_match=match_id)
banc.meta.fafb.nb <- banc_updateids(banc.meta.fafb.nb, 
                                    root.column = "pt_root_id", 
                                    supervoxel.column = "pt_supervoxel_id", 
                                    position.column = "pt_position")
manc.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE) %>%
  dplyr::mutate(score=0.9)  %>%
  dplyr::filter(valid=='t') %>%
  dplyr::mutate(validation = TRUE) %>%
  dplyr::select(-valid) %>%
  dplyr::rename_with(~sub("^query_id$", "query_root_id", .x))
manc.matches.df.valid <- rbind(bc.manc.matches.df.valid,manc.matches.df.valid) %>%
  dplyr::distinct(pt_root_id,
                  match_id,
                  .keep_all=TRUE)
banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather")) %>%
  rbind(manc.matches.df.valid) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(!grepl("\\_",match_id), score > 0) %>%
  dplyr::rename(manc_match=match_id)
banc.meta.manc.nb <- banc_updateids(banc.meta.manc.nb, 
                                    root.column = "pt_root_id", 
                                    supervoxel.column = "pt_supervoxel_id", 
                                    position.column = "pt_position")

# Synergise meta data
fw.meta <- fw.meta %>%
  dplyr::mutate(super_class = dplyr::case_when(
    cell_class=="DN" ~ "descending_neuron",
    grepl("motor", cell_class) ~ "brain_motor_neuron",
    grepl("motor", super_class) ~ "brain_motor_neuron",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    cell_class=="AN" & grepl("sensory", super_class) ~ "sensory_ascending",
    cell_class=="AN" ~ "ascending_neuron",
    cell_class=="DN" ~ "descending_neuron",
    grepl("motor", cell_class) ~ "head_motor_neuron",
    grepl("motor", super_class) ~ "head_motor_neuron",
    TRUE ~ cell_class
  )) %>%
  dplyr::mutate(flow = dplyr::case_when(
    cell_class%in%c("descending_neuron") ~ "brain_vnc_transfer",
    cell_class%in%c("ascending_neuron") ~ "vnc_brain_transfer",
    flow=="intrinsic" ~ "brain_intrinsic",
    flow=="afferent" ~ "brain_afferent",
    flow=="efferent" ~ "brain_efferent",
    TRUE ~ flow
  )) %>%
  dplyr::distinct(root_783, .keep_all = TRUE)
mc.meta <- mc.meta %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    cell_class=="intrinsic_neuron" ~ "vnc_intrinsic_neuron",
    cell_class=="motor_neuron" ~ "vnc_motor_neuron",
    cell_class=="efferent_neuron" ~ "vnc_efferent_neuron",
    TRUE ~ cell_class
  ))  %>%
  dplyr::mutate(flow = dplyr::case_when(
    cell_class%in%c("descending_neuron") ~ "brain_vnc_transfer",
    cell_class%in%c("ascending_neuron") ~ "vnc_brain_transfer",
    flow=="intrinsic" ~ "vnc_intrinsic",
    flow=="afferent" ~ "vnc_afferent",
    flow=="transfer" ~ "vnc_afferent",
    flow=="efferent" ~ "vnc_efferent",
    TRUE ~ flow
  )) %>%
  dplyr::distinct(bodyid, .keep_all = TRUE)

# Get neck neurons in other data sets
fw.meta.neck <- fw.meta %>% 
  dplyr::filter(super_class %in% c("descending","ascending","sensory_ascending"))
mc.meta.neck <- mc.meta %>% 
  dplyr::filter(super_class %in% c("descending","ascending","sensory_ascending"))

# Function to perform iterative matching
iterative_matching <- function(data, id_col1, id_col2, score_col, cell_type_col, pre_matched = NULL) {
  library(dplyr)
  
  # Validation checks
  if(!all(c(id_col1, id_col2, score_col, cell_type_col) %in% base::names(data))) {
    missing_cols <- base::setdiff(c(id_col1, id_col2, score_col, cell_type_col), base::names(data))
    base::stop("Missing columns in main data: ", base::paste(missing_cols, collapse = ", "))
  }
  
  if(!is.null(pre_matched) && !all(c(id_col1, id_col2, cell_type_col) %in% base::names(pre_matched))) {
    missing_cols <- base::setdiff(c(id_col1, id_col2, cell_type_col), base::names(pre_matched))
    base::stop("Missing columns in pre_matched data: ", base::paste(missing_cols, collapse = ", "))
  }
  
  # Initial statistics
  base::message("Initial data:")
  base::message("Unique ", id_col1, ": ", base::length(base::unique(data[[id_col1]])))
  base::message("Unique ", id_col2, ": ", base::length(base::unique(data[[id_col2]])))
  base::message("Total rows: ", base::nrow(data))
  
  # Initialize matched with pre-existing matches if provided
  if(!is.null(pre_matched)) {
    matched <- pre_matched
    base::message("Pre-existing matches: ", base::nrow(pre_matched))
    # Remove pre-matched entries from the data
    remaining <- data %>% 
      dplyr::anti_join(matched, by = c(id_col1, id_col2, cell_type_col)) %>%
      dplyr::arrange(dplyr::desc(.data[[score_col]]))
  } else {
    matched <- base::data.frame()
    # Use all data if no pre-matched entries
    remaining <- data %>% 
      dplyr::arrange(dplyr::desc(.data[[score_col]]))
  }
  
  # Track matched IDs
  matched_id1 <- base::unique(c(matched[[id_col1]], base::character(0)))
  matched_id2 <- base::unique(c(matched[[id_col2]], base::character(0)))
  
  # Phase 1: Find 1:1 matches
  iteration <- 1
  while(base::nrow(remaining) > 0) {
    base::message("\nIteration ", iteration, " (1:1 matching phase)")
    
    # Remove rows with already matched IDs
    remaining <- remaining %>%
      dplyr::filter(!(.data[[id_col1]] %in% matched_id1)) %>%
      dplyr::filter(!(.data[[id_col2]] %in% matched_id2))
    
    if(base::nrow(remaining) == 0) break
    
    # Find best 1:1 match
    best_match <- remaining %>%
      dplyr::slice_head(n = 1)
    
    # Check if the match respects cell type boundaries
    if(base::nrow(matched) > 0) {
      existing_match_1 <- matched %>% 
        dplyr::filter(.data[[id_col1]] == best_match[[id_col1]]) %>%
        dplyr::pull(.data[[cell_type_col]])
      
      existing_match_2 <- matched %>% 
        dplyr::filter(.data[[id_col2]] == best_match[[id_col2]]) %>%
        dplyr::pull(.data[[cell_type_col]])
      
      if(base::length(existing_match_1) > 0 && existing_match_1 != best_match[[cell_type_col]]) {
        base::message("Skipping match due to cell type mismatch for ", id_col1)
        next
      }
      
      if(base::length(existing_match_2) > 0 && existing_match_2 != best_match[[cell_type_col]]) {
        base::message("Skipping match due to cell type mismatch for ", id_col2)
        next
      }
    }
    
    # Update matched IDs and results
    matched_id1 <- base::c(matched_id1, best_match[[id_col1]])
    matched_id2 <- base::c(matched_id2, best_match[[id_col2]])
    matched <- dplyr::bind_rows(matched, best_match)
    
    base::message("1:1 match made: ", base::nrow(best_match))
    iteration <- iteration + 1
  }
  
  # Track matched IDs
  matched_id1 <- base::unique(c(matched[[id_col1]], base::character(0)))
  matched_id2 <- base::unique(c(matched[[id_col2]], base::character(0)))
  
  # Reset remaining to include all original unmatched entries
  remaining <- data %>%
    dplyr::filter(!(.data[[id_col1]] %in% matched_id1)) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]))
  
  # Phase 2: Allow 1:many matches for remaining unmatched entries within cell types
  if(base::nrow(remaining) > 0) {
    base::message("\nStarting 1:many matching phase for remaining unmatched entries")
    
    many_matches <- remaining %>%
      dplyr::group_by(.data[[cell_type_col]], .data[[id_col1]]) %>%
      dplyr::slice_max(.data[[score_col]], n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    # Check if the many matches respect cell type boundaries
    valid_many_matches <- many_matches %>%
      dplyr::left_join(
        matched %>% 
          dplyr::select(!!id_col2, matched_cell_type = !!cell_type_col) %>%
          dplyr::distinct(),
        by = id_col2
      ) %>%
      dplyr::filter(is.na(.data$matched_cell_type) | .data$matched_cell_type == .data[[cell_type_col]]) %>%
      dplyr::select(-matched_cell_type)
    
    matched <- dplyr::bind_rows(matched, valid_many_matches)
    base::message("Additional 1:many matches made: ", base::nrow(valid_many_matches))
  }
  
  # Track matched IDs
  matched_id1 <- base::unique(c(matched[[id_col1]], base::character(0)))
  matched_id2 <- base::unique(c(matched[[id_col2]], base::character(0)))
  
  # Reset remaining to include all original unmatched entries
  remaining <- data %>%
    dplyr::filter(!(.data[[id_col2]] %in% matched_id2)) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]))
  
  # Phase 3: Allow 1:many matches for remaining unmatched entries within cell types
  if(base::nrow(remaining) > 0) {
    base::message("\nStarting 1:many matching phase for remaining unmatched entries")
    
    many_matches <- remaining %>%
      dplyr::group_by(.data[[cell_type_col]], .data[[id_col2]]) %>%
      dplyr::slice_max(.data[[score_col]], n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    # Check if the many matches respect cell type boundaries
    valid_many_matches <- many_matches %>%
      dplyr::left_join(
        matched %>% 
          dplyr::select(!!id_col1, matched_cell_type = !!cell_type_col) %>%
          dplyr::distinct(),
        by = id_col1
      ) %>%
      dplyr::filter(is.na(.data$matched_cell_type) | .data$matched_cell_type == .data[[cell_type_col]]) %>%
      dplyr::select(-matched_cell_type)
    
    matched <- dplyr::bind_rows(matched, valid_many_matches)
    base::message("Additional 1:many matches made: ", base::nrow(valid_many_matches))
  }
  
  # Sort and summarize results
  matched <- matched %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]))
  
  # Final statistics
  base::message("\nFinal results:")
  base::message("Total matches found: ", base::nrow(matched))
  base::message("Unique ", id_col1, " matched: ", base::length(base::unique(matched[[id_col1]])))
  base::message("Unique ", id_col2, " matched: ", base::length(base::unique(matched[[id_col2]])))
  base::message("Unique ", id_col1, " unmatched: ", base::length(base::unique(setdiff(data[[id_col1]],matched[[id_col1]]))))
  base::message("Unique ", id_col2, " unmatched: ", base::length(base::unique(setdiff(data[[id_col2]],matched[[id_col2]]))))
  
  # Count and report 1:many matches
  many_matches <- matched %>%
    dplyr::count(.data[[id_col2]]) %>%
    dplyr::filter(n > 1)
  if(base::nrow(many_matches) > 0) {
    base::message("\n1:many matches:")
    base::message("Number of ", id_col2, " with multiple matches: ", base::nrow(many_matches))
    base::message("Average matches per reused ", id_col2, ": ", base::mean(many_matches$n))
  }
  return(matched)
}

# Get BANC cell type matches
bc.meta.neck.cts <- bc.meta.neck %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type) & grepl("^AN|^S", manc_cell_type) ~ manc_cell_type,
    is.na(cell_type) & !grepl("^AN|^S", manc_cell_type) ~ fafb_cell_type,
    TRUE ~ cell_type
  )) %>%
  dplyr::mutate(fafb_cell_type = fw.meta$cell_type[match(fafb_match,fw.meta$root_783)]) %>%
  dplyr::mutate(manc_cell_type = mc.meta$cell_type[match(manc_match,mc.meta$bodyid)]) %>%
  dplyr::distinct(root_id, region, side,
                  super_class, cell_type, 
                  fafb_cell_type, manc_cell_type, cell_sub_class) %>%
  dplyr::arrange(super_class,
                 cell_type, fafb_cell_type, manc_cell_type, 
                 cell_sub_class, root_id)

# Step 1: Join FAFB metadata with matching scores
fafb_joined_data <- banc.meta.fafb.nb %>%
  dplyr::distinct(pt_root_id, fafb_match, .keep_all = TRUE) %>%
  dplyr::full_join(bc.meta.neck.cts  %>% dplyr::select(root_id, region, side, super_class, fafb_cell_type, cell_sub_class),
                   by = c("pt_root_id" = "root_id")) %>%
  dplyr::full_join(fw.meta.neck[, c("root_783")], by = c("fafb_match" = "root_783")) %>%
  dplyr::mutate(fafb_super_class = fw.meta$super_class[match(fafb_match,fw.meta$root_783)]) %>%
  dplyr::filter(match_cell_type == fafb_cell_type,
                gsub("sensory_|efferent_","",super_class) == gsub("sensory_|efferent_","",fafb_super_class),
                (pt_root_id %in% bc.meta.neck$root_id | fafb_match %in% fw.meta.neck$root_783))  %>%
  dplyr::mutate(fafb_cell_type=match_cell_type,
                fafb_nblast_score=score) %>%
  dplyr::select(-score)

# Perform iterative matching for FAFB to BANC
fafb_banc_matches <- iterative_matching(fafb_joined_data, 
                                        id_col1 = "pt_root_id", 
                                        id_col2 = "fafb_match", 
                                        score_col = "fafb_nblast_score", 
                                        cell_type_col = "match_cell_type")

# Fill in missing with NBLAST matches for FAFB
fw.id.missing <- setdiff(fw.meta.neck$root_783, fafb_banc_matches$fafb_match)
banc.fafb.id.missing <- setdiff(bc.meta.neck$root_id, fafb_banc_matches$pt_root_id)
fafb.missing <- banc.meta.fafb.nb %>%
  dplyr::left_join(bc.meta.neck.cts %>% dplyr::select(root_id, region, side, super_class, fafb_cell_type, cell_sub_class), 
                   by = c("pt_root_id" = "root_id")) %>%
  dplyr::mutate(fafb_super_class = fw.meta$super_class[match(fafb_match,fw.meta$root_783)]) %>%
  dplyr::filter(fafb_match %in% fw.id.missing | pt_root_id %in% banc.fafb.id.missing,
                match_cell_type %in% fw.meta.neck$cell_type,
                gsub("sensory_|efferent_","",super_class) == gsub("sensory_|efferent_","",fafb_super_class)) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::distinct(fafb_match, .keep_all = TRUE) %>%
  dplyr::mutate(fafb_cell_type=match_cell_type,
                fafb_nblast_score=score) %>%
  dplyr::select(-score)
fafb.missing.df <- iterative_matching(fafb.missing,  
                                      id_col1 = "fafb_match", 
                                      id_col2 = "pt_root_id", 
                                      score_col = "fafb_nblast_score", 
                                      cell_type_col = "match_cell_type", 
                                      pre_matched = fafb_banc_matches)
fafb_banc_matches <- unique(plyr::rbind.fill(fafb_banc_matches,fafb.missing.df))

# Join MANC metadata with matching scores
manc_joined_data <- banc.meta.manc.nb %>%
  dplyr::distinct(pt_root_id, manc_match, score, .keep_all = TRUE) %>%
  dplyr::full_join(bc.meta.neck.cts %>% dplyr::select(root_id, region, side, super_class, manc_cell_type, cell_sub_class), 
                   by = c("pt_root_id" = "root_id")) %>%
  dplyr::full_join(mc.meta.neck[, c("bodyid")], by = c("manc_match" = "bodyid")) %>%
  dplyr::mutate(manc_super_class = mc.meta$super_class[match(manc_match,mc.meta$bodyid)]) %>%
  dplyr::filter(match_cell_type == manc_cell_type,
                gsub("sensory_|efferent_","",super_class) == gsub("sensory_|efferent_","",manc_super_class),
                pt_root_id %in% bc.meta.neck$root_id | manc_match %in% mc.meta.neck$bodyid) %>%
  dplyr::mutate(manc_cell_type=match_cell_type,
                manc_nblast_score=score) %>%
  dplyr::select(-score)

# Step 4: Perform iterative matching for BANC to MANC
banc_manc_matches <- iterative_matching(manc_joined_data, 
                                        id_col1 = "pt_root_id", 
                                        id_col2 =  "manc_match", 
                                        score_col = "manc_nblast_score", 
                                        cell_type_col = "match_cell_type")

# Step 5: Fill in missing with NBLAST matches for MANC
manc.id.missing <- setdiff(mc.meta.neck$bodyid, banc_manc_matches$manc_match)
banc.manc.id.missing <- setdiff(bc.meta.neck$root_id, banc_manc_matches$pt_root_id)
manc.missing <- banc.meta.manc.nb %>%
  dplyr::left_join(bc.meta.neck.cts %>% dplyr::select(root_id, region, side, super_class, manc_cell_type, cell_sub_class), 
                   by = c("pt_root_id" = "root_id")) %>%
  dplyr::mutate(manc_super_class = mc.meta$super_class[match(manc_match,mc.meta$bodyid)]) %>%
  dplyr::filter(manc_match %in% manc.id.missing | pt_root_id %in% banc.manc.id.missing,
                match_cell_type %in% mc.meta.neck$cell_type,
                gsub("sensory_|efferent_","",super_class) == gsub("sensory_|efferent_","",manc_super_class)) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::distinct(manc_match, .keep_all = TRUE) %>%
  dplyr::mutate(manc_cell_type=match_cell_type,
                manc_nblast_score=score) %>%
  dplyr::select(-score)
manc.missing.df <- iterative_matching(manc.missing, 
                                      id_col1 = "manc_match", 
                                      id_col2 = "pt_root_id", 
                                      score_col = "manc_nblast_score", 
                                      cell_type_col = "match_cell_type",
                                      pre_matched = banc_manc_matches)
banc_manc_matches <- unique(plyr::rbind.fill(banc_manc_matches,manc.missing.df))

# Step 6: Create the final mapping of FAFB to BANC to MANC
final_mapping <- 
  fafb_banc_matches[,c("pt_root_id","side","fafb_match","fafb_cell_type","fafb_nblast_score")] %>%
  dplyr::full_join(banc_manc_matches[,c("pt_root_id","side","manc_match","manc_cell_type","manc_nblast_score")], 
                   relationship = "many-to-many",
                   by = "pt_root_id") %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    is.na(side.x) ~ side.y,
    TRUE ~ side.x
  )) %>%
  dplyr::select(pt_root_id, 
                fafb_match, manc_match, 
                fafb_cell_type, manc_cell_type, 
                fafb_nblast_score, manc_nblast_score) %>%
  dplyr::left_join(fw.meta[,c("side","super_class","cell_class","flow","root_id","sez_class","hemilineage")], by = c("fafb_match"="root_id")) %>%
  dplyr::left_join(mc.meta[,c("side","super_class","cell_class","flow","bodyid","hemilineage","serial_motif")], by = c("manc_match"="bodyid")) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    is.na(super_class.x) ~ super_class.y,
    TRUE ~ super_class.x
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    is.na(cell_class.x) ~ cell_class.y,
    TRUE ~ cell_class.x
  )) %>%
  dplyr::mutate(flow = dplyr::case_when(
    is.na(flow.x) ~ flow.y,
    TRUE ~ flow.x
  )) %>%
  dplyr::mutate(hemilineage = dplyr::case_when(
    is.na(hemilineage.x) ~ hemilineage.y,
    TRUE ~ hemilineage.x
  )) %>%
  dplyr::mutate(side = dplyr::case_when(
    !is.na(side.x) ~ side.x,
    is.na(side.x) ~ side.y,
    is.na(side.y) ~ side.x,
    super_class=="ascending" ~ side.y,
    super_class=="sensory_ascending" ~ side.y,
    super_class=="descending" ~ side.x,
    TRUE ~ side.x
  )) %>%
  # region="neck_connective" disabled — ascending/descending neurons are now
  # identified via super_class; cervical_connective lives as a tract entry.
  # dplyr::mutate(region="neck_connective") %>%
  dplyr::left_join(bc[,c("root_id","an_type","dn_type","side")], by = c("pt_root_id"="root_id")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    is.na(an_type) ~ paste0("DN_",dn_type),
    is.na(dn_type) ~ paste0("AN_",an_type),
    TRUE ~ NA
  )) %>%
  dplyr::mutate(side = dplyr::case_when(
    is.na(side.x) ~ side.y,
    is.na(side.y) ~ side.x,
    TRUE ~ side.x
  )) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pt_root_id, fafb_match, manc_match,
                  fafb_cell_type, manc_cell_type, 
                  fafb_nblast_score, manc_nblast_score,
                  side, region, flow, super_class, cell_class, cell_sub_class, 
                  sez_class, serial_motif,
                  hemilineage) %>%
  dplyr::arrange(super_class, pt_root_id, fafb_match, manc_match, fafb_cell_type, manc_cell_type)

# Remove any rows where any of the matches is NA
final_mapping_tw <- final_mapping %>%
  dplyr::filter(!is.na(manc_match), !is.na(fafb_match), !is.na(pt_root_id)) %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(fafb_cell_type) | grepl("^AN|^S", manc_cell_type) ~ manc_cell_type,
    TRUE ~ fafb_cell_type
  ))

# Print summary of matches
cat("FAFB to BANC matches:", nrow(fafb_banc_matches), "\n")
cat("MANC to BANC matches:", nrow(banc_manc_matches), "\n")
cat("Final three-way matches:", nrow(final_mapping_tw), "\n")

# Check for distinctness
cat("Distinct FAFB IDs in FAFB-BANC matches:", fafb_banc_matches %>% dplyr::distinct(fafb_match) %>% nrow(), "\n")
cat("Distinct MANC IDs in MANC-BANC matches:", banc_manc_matches %>% dplyr::distinct(manc_match) %>% nrow(), "\n")
cat("Distinct BANC IDs in FAFB-BANC matches:", fafb_banc_matches %>% dplyr::distinct(pt_root_id) %>% nrow(), "\n")
cat("Distinct BANC IDs in MANC-BANC matches:", banc_manc_matches %>% dplyr::distinct(pt_root_id) %>% nrow(), "\n")

# Check for non-optimal missing
cat("FAFB neurons lacking human-verified BANC match:", length(fw.id.missing), "\n")
cat("BANC neurons lacking human-verified FAFB match:", length(banc.fafb.id.missing), "\n")
cat("MANC neurons lacking human-verified BANC match:", length(manc.id.missing), "\n")
cat("BANC neurons lacking human-verified MANC match:", length(banc.manc.id.missing), "\n")

# Check for missing
cat("FAFB neurons missing from neck bridge:", sum(!fw.meta.neck$root_783 %in% final_mapping$fafb_match), "\n")
cat("MANC neurons missing from neck bridge:", sum(!mc.meta.neck$bodyid %in% final_mapping$manc_match), "\n")
cat("BANC neurons missing from neck bridge:", sum(!bc.meta.neck$root_id %in% final_mapping$pt_root_id), "\n")

# Matches per BANC neuron
banc.load <- final_mapping %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(fafb_ids = length(unique(fafb_match)),
                manc_ids = length(unique(manc_match))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pt_root_id, fafb_ids, manc_ids) %>%
  dplyr::arrange(dplyr::desc(manc_ids), dplyr::desc(fafb_ids))
head(as.data.frame(banc.load), 100)

banc.load <- final_mapping %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(fafb_ids = length(unique(fafb_match)),
                manc_ids = length(unique(manc_match))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pt_root_id, fafb_ids, manc_ids) %>%
  dplyr::arrange(dplyr::desc(fafb_ids), dplyr::desc(manc_ids))
head(as.data.frame(banc.load), 100)

#####################################
### Load edgelists and manipulate ###
#####################################

# Combine meta data
fw.mapping <- final_mapping %>%
  dplyr::distinct(pt_root_id, fafb_match, .keep_all = TRUE) %>%
  dplyr::select(banc_id = pt_root_id, fafb_match, 
                fafb_cell_type, manc_cell_type, 
                fafb_nblast_score, manc_nblast_score,
                side, region, flow, 
                super_class, cell_class, cell_sub_class, 
                sez_class, serial_motif,
                hemilineage) %>%
  dplyr::filter(!is.na(fafb_match), !is.na(banc_id))
mc.mapping <- final_mapping %>%
  dplyr::distinct(pt_root_id, manc_match, .keep_all = TRUE) %>%
  dplyr::select(banc_id = pt_root_id, manc_match, 
                fafb_cell_type, manc_cell_type, 
                fafb_nblast_score, manc_nblast_score,
                side, region, flow, 
                super_class, cell_class, cell_sub_class, 
                sez_class, serial_motif,
                hemilineage) %>%
  dplyr::filter(!is.na(manc_match), !is.na(banc_id))
meta <- fw.meta %>%
  dplyr::mutate(soma = is.na(soma_x)) %>%
  dplyr::left_join(fw.mapping, by = c("root_783" = "fafb_match")) %>%
  dplyr::mutate(region = dplyr::case_when(
    !is.na(region.y)&region.y!="" ~ region.y,
    TRUE ~ region.x
  )) %>%
  dplyr::mutate(side = dplyr::case_when(
    !is.na(side.y)&side.y!="" ~ side.y,
    TRUE ~ side.x
  )) %>%
  dplyr::mutate(flow = dplyr::case_when(
    !is.na(flow.y)&flow.y!="" ~ flow.y,
    TRUE ~ flow.x
  )) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    !is.na(super_class.y)&super_class.y!="" ~ super_class.y,
    TRUE ~ super_class.x
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    !is.na(cell_class.y)&cell_class.y!="" ~ cell_class.y,
    TRUE ~ cell_class.x
  )) %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    !is.na(cell_sub_class.y)&cell_sub_class.y!="" ~ cell_sub_class.y,
    TRUE ~ cell_sub_class.x
  )) %>%
  dplyr::mutate(hemilineage = dplyr::case_when(
    !is.na(hemilineage.y)&hemilineage.y!="" ~ hemilineage.y,
    TRUE ~ hemilineage.x
  )) %>%
  dplyr::mutate(sez_class = dplyr::case_when(
    !is.na(sez_class.y)&sez_class.y!="" ~ sez_class.y,
    TRUE ~ sez_class.x
  )) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y")) %>%
  plyr::rbind.fill(mc.meta %>% 
                     dplyr::mutate(region = "vnc") %>%
                     dplyr::left_join(mc.mapping, by = c("bodyid" = "manc_match")) %>%
                     dplyr::mutate(region = dplyr::case_when(
                       !is.na(region.y)&region.y!="" ~ region.y,
                       TRUE ~ region.x
                     )) %>%
                     dplyr::mutate(side = dplyr::case_when(
                       !is.na(side.y)&side.y!="" ~ side.y,
                       TRUE ~ side.x
                     )) %>%
                     dplyr::mutate(flow = dplyr::case_when(
                       !is.na(flow.y)&flow.y!="" ~ flow.y,
                       TRUE ~ flow.x
                     )) %>%
                     dplyr::mutate(super_class = dplyr::case_when(
                       !is.na(super_class.y)&super_class.y!="" ~ super_class.y,
                       TRUE ~ super_class.x
                     )) %>%
                     dplyr::mutate(cell_class = dplyr::case_when(
                       !is.na(cell_class.y)&cell_class.y!="" ~ cell_class.y,
                       TRUE ~ cell_class.x
                     )) %>%
                     dplyr::mutate(cell_sub_class = dplyr::case_when(
                       !is.na(cell_sub_class.y)&cell_sub_class.y!="" ~ cell_sub_class.y,
                       TRUE ~ cell_sub_class.x
                     )) %>%
                     dplyr::mutate(hemilineage = dplyr::case_when(
                       !is.na(hemilineage.y)&hemilineage.y!="" ~ hemilineage.y,
                       TRUE ~ hemilineage.x
                     )) 
  ) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(id = dplyr::coalesce(banc_id, root_783, bodyid)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-supervoxel_id, -proofread, -status, 
                -pos_x, -pos_y, -pos_z, -nucleus_id, 
                -soma_x, -soma_y, -soma_z, -tosoma_location, -soma_location, -root_location,
                -hemibrain_match, -notes, -type) %>%
  dplyr::mutate(across(where(is.character), ~replace_na(., "")),
                across(where(is.factor), ~as.factor(replace_na(as.character(.), ""))),
                across(where(is.numeric), ~replace_na(., 0)),
                across(where(is.logical), ~.)) %>%
  dplyr::select(id, fafb_id = root_783, manc_id = bodyid, banc_id, dataset,
                side, region, cell_function, seed, 
                hemilineage, nerve, origin, target, serial_motif, neuromere, morphology_group,
                flow, super_class, cell_class, cell_sub_class, cell_type, sez_class,
                top_nt, top_nt_p, known_nt, known_nt_source,
                soma_dcv_count, soma_dcv_density,
                cell_dcv_count, cell_dcv_density,
                #total_outputs, axon_outputs, dend_outputs,
                #total_inputs, dend_inputs, dend_outputs,
                #total_outputs_density, axon_outputs_density, dend_outputs_density,
                #total_inputs_density, dend_inputs_density, dend_outputs_density,
                #total_length, axon_length, dend_length, pn_length, pd_length,
                soma_volume, cable_length
  )  %>%
  dplyr::mutate(dataset = dplyr::case_when(
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-FAFB-MANC",
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!="" ~ "BANC-FAFB",
    !is.na(banc_id)&banc_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-MANC",
    !is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "FAFB-MANC",
    !is.na(fafb_id)&fafb_id!="" ~ "FAFB",
    !is.na(manc_id)&manc_id!="" ~ "MANC",
    !is.na(banc_id)&banc_id!="" ~ "BANC",
    TRUE ~ dataset
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    cell_class=="Kenyon_cell" ~ "kenyon_cell",
    cell_class=="sensory_ascending" & cell_function=="gustatory" ~ "chemosensory",
    cell_class=="ascending_neuron" ~ "ascending",
    cell_class=="descending_neuron" ~ "descending",
    TRUE ~ cell_function
  )) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    super_class=="sensory_neuron" ~ "sensory",
    super_class=="intrinsic_neuron" & region=="vnc" ~ "vnc_central_other",
    super_class=="central" & region!="vnc" ~ "brain_central_other",
    super_class=="endocrine" & region=="vnc" ~ "vnc_endocrine",
    super_class=="endocrine" & region!="vnc" ~ "brain_endocrine",
    cell_class=="ascending_neuron" ~ "ascending",
    cell_class=="descending_neuron" ~ "descending",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(seed = dplyr::case_when(
    !is.na(cell_sub_class) ~ seed,
    cell_sub_class=="" ~ seed,
    cell_function=="ascending" ~ paste0(cell_sub_class,"_",side),
    cell_function=="descending" ~ paste0(cell_sub_class,"_",side),
    TRUE ~ seed
  ))  %>%
  dplyr::mutate(region = dplyr::case_when(
    super_class == "visual_projection" ~ "optic",
    super_class == "visual_centrifugal" ~ "optic",
    # "neck_connective" disabled as a region value — ascending/descending
    # neurons are identified via super_class; cervical_connective is a tract.
    TRUE ~ region
  )) %>%
  dplyr::mutate(seed = gsub("__","_",seed)) %>%
  dplyr::mutate(seed = ifelse(grepl("unknown",seed),'',seed)) %>%
  dplyr::mutate(cell_function = ifelse(grepl("unknown",cell_function),'',cell_function))

## Copy
# file.copy(from = "/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783/flywire_783_data.sqlite", 
#           to = "/n/data1/hms/neurobio/wilson/banc/connectivity/fafb_783_data.sqlite",
#           overwrite = TRUE)

# Get FAFB data, not split
fw.elist.simp <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "flywire_783_edgelist.feather")) %>%
  dplyr::left_join(fw.mapping %>% dplyr::distinct(fafb_match, .keep_all = TRUE),
                   by = c("pre" = "fafb_match")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pre = dplyr::coalesce(banc_id, pre)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-banc_id) %>%
  dplyr::left_join(fw.mapping %>% dplyr::distinct(fafb_match, .keep_all = TRUE),
                   by = c("post" = "fafb_match")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(post = dplyr::coalesce(banc_id, post)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-banc_id)

# Get MANC data, not split
mc.elist.simp <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "manc_1.2.1_edgelist.feather")) %>%
  dplyr::left_join(mc.mapping %>% dplyr::distinct(manc_match, .keep_all = TRUE),
                   by = c("pre" = "manc_match")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pre = dplyr::coalesce(banc_id, pre)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-banc_id) %>%
  dplyr::left_join(mc.mapping %>% dplyr::distinct(manc_match, .keep_all = TRUE),
                   by = c("post" = "manc_match")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(post = dplyr::coalesce(banc_id, post)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-banc_id)

# Simplify
elist.simp <- fw.elist.simp %>%
  dplyr::select(post, pre, count) %>%
  rbind(mc.elist.simp %>% dplyr::select(post, pre, count)) %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre, post, .keep_all = TRUE) %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(norm = round(count/post_count,6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(pre, post, count, norm, post_count, pre_count)

# Assess neck bridge
library(ggsankey)

# Step 1: Count the total number of neurons in each pool
step1 <- bind_rows(
  bc.meta.neck %>% distinct(id = root_id) %>% mutate(node = "BANC"),
  mc.meta.neck %>% distinct(id = bodyid) %>% mutate(node = "MANC"),
  fw.meta.neck %>% distinct(id = root_783) %>% mutate(node = "FAFB")) %>% 
  dplyr::mutate(x = 1,
                next_x = NA,
                next_node = NA)

# Step 2: Count the matches and no_matches
step2 <- final_mapping %>%
  # BANC perspective
  dplyr::distinct(id = pt_root_id, .keep_all = TRUE) %>%
  dplyr::full_join(bc.meta.neck %>% 
                     dplyr::distinct(id = root_id),
                   by = "id") %>%
  dplyr::mutate(node=dplyr::case_when(
    !is.na(fafb_match)&!is.na(manc_match) ~ "BANC-FAFB-MANC",
    !is.na(fafb_match) ~ "BANC-FAFB",
    !is.na(manc_match) ~ "BANC-MANC",
    TRUE ~ "BANC-none"
  ),
  next_node="BANC") %>%
  dplyr::select(id, node, next_node) %>%
  # FAFB perspecive
  rbind(
    final_mapping %>%
      # BANC perspective
      dplyr::distinct(id = fafb_match, .keep_all = TRUE) %>%
      dplyr::full_join(fw.meta.neck %>% 
                         dplyr::distinct(id = root_783),
                       by = "id") %>%
      dplyr::mutate(node=dplyr::case_when(
        !is.na(pt_root_id)&!is.na(manc_match) ~ "BANC-FAFB-MANC",
        !is.na(pt_root_id) ~ "BANC-FAFB",
        TRUE ~ "FAFB-none"
      ),
      next_node="FAFB") %>%
      dplyr::select(id, node, next_node)
  ) %>%
  # MANC perspecive
  rbind(
    final_mapping %>%
      # BANC perspective
      dplyr::distinct(id = manc_match, .keep_all = TRUE) %>%
      dplyr::full_join(mc.meta.neck %>% 
                         dplyr::distinct(id = bodyid),
                       by = "id") %>%
      dplyr::mutate(node=dplyr::case_when(
        !is.na(pt_root_id)&!is.na(fafb_match) ~ "BANC-FAFB-MANC",
        !is.na(pt_root_id) ~ "BANC-MANC",
        TRUE ~ "MANC-none"
      ),
      next_node="MANC") %>%
      dplyr::select(id, node, next_node)
  ) %>%
  dplyr::mutate(x = 2,
                next_x = 1)

# 1:1, 1:many or many-to-many
step3 <- final_mapping %>%
  # BANC-FAFB-MANC perspective for FAFB
  dplyr::filter(!is.na(fafb_match), !is.na(manc_match)) %>%
  dplyr::group_by(id=pt_root_id) %>%
  dplyr::mutate(
    rank = length(unique(na.omit(fafb_match))),
    node =  paste0("FAFB_", cut(rank, 
                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                include.lowest = TRUE,
                                right = FALSE,
                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
    next_node = "BANC-FAFB-MANC") %>%
  dplyr::select(id, rank, node, next_node) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  # BANC-FAFB-MANC perspective for MANC
  rbind(
    final_mapping %>%
      dplyr::filter(!is.na(fafb_match), !is.na(manc_match)) %>%
      dplyr::group_by(id=pt_root_id) %>%
      dplyr::mutate(rank = length(unique(na.omit(manc_match))),
                    node =  paste0("MANC_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-FAFB-MANC") %>%
      dplyr::select(id, rank, node, next_node) %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # BANC-FAFB-MANC perspective for BANC
  rbind(
    final_mapping %>%
      dplyr::filter(!is.na(fafb_match), !is.na(manc_match)) %>%
      dplyr::group_by(id=fafb_match) %>%
      dplyr::mutate(rank = length(unique(na.omit(pt_root_id))),
                    node =  paste0("BANC_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-FAFB-MANC") %>%
      dplyr::select(id, rank, node, next_node) %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # BANC-FAFB perspective
  rbind(
    final_mapping %>%
      dplyr::filter(!is.na(fafb_match), is.na(manc_match)) %>%
      dplyr::group_by(id=pt_root_id) %>%
      dplyr::mutate(rank = length(unique(na.omit(fafb_match))),
                    node =  paste0("FAFB_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-FAFB") %>%
      dplyr::select(id, node, next_node)  %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # BANC-FAFB perspective for BANC
  rbind(
    final_mapping %>%
      dplyr::filter(!is.na(fafb_match), is.na(manc_match)) %>%
      dplyr::group_by(id=fafb_match) %>%
      dplyr::mutate(rank = length(unique(na.omit(pt_root_id))),
                    node =  paste0("BANC_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-FAFB") %>%
      dplyr::select(id, rank, node, next_node) %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # BANC-MANC perspective
  rbind(
    final_mapping %>%
      dplyr::filter(is.na(fafb_match), !is.na(manc_match)) %>%
      dplyr::group_by(id=pt_root_id) %>%
      dplyr::mutate(rank = length(unique(na.omit(manc_match))),
                    node =  paste0("MANC_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-MANC") %>%
      dplyr::select(id, node, next_node)  %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # BANC-MANC perspective for BANC
  rbind(
    final_mapping %>%
      dplyr::filter(is.na(fafb_match), !is.na(manc_match)) %>%
      dplyr::group_by(id=manc_match) %>%
      dplyr::mutate(rank = length(unique(na.omit(pt_root_id))),
                    node =  paste0("BANC_", cut(rank, 
                                                breaks = c(1, 2, seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5)),
                                                include.lowest = TRUE,
                                                right = FALSE,
                                                labels = c("1", paste0("[", seq(1, ceiling(max(100, na.rm = TRUE)/5)*5 - 5, by=5), 
                                                                       ",", seq(5, ceiling(max(100, na.rm = TRUE)/5)*5, by=5), ")")))),
                    next_node = "BANC-MANC") %>%
      dplyr::select(id, rank, node, next_node) %>%
      dplyr::distinct(id, .keep_all = TRUE)
  ) %>%
  # Steps
  dplyr::arrange(rank) %>%
  dplyr::mutate(x = 3,
                next_x = 2,
                node = factor(node, levels = unique(node))) %>%
  dplyr::arrange(node) %>%
  dplyr::select(-rank)

# Combine all data
sankey_data <- dplyr::bind_rows(step1,step2,step3)

# Create plot
g.sank <- ggplot(sankey_data, aes(x = x, 
                                  next_x = next_x, 
                                  node = node, 
                                  next_node = next_node,
                                  fill = node, 
                                  label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black") +
  geom_sankey_label(size = 3) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Save
ggsave(plot = g.sank,
       filename = "inst/images/frankenbrain_neck_bridge_v.1.7.png", 
       width = 12, height = 8, dpi = 300)

# ## Sanity check
# fw.elist.raw <- dplyr::tbl(con, "edgelist_simple") %>%
#   dplyr::collect() 
# pre.missing <- setdiff(fw.meta$root_783,fw.elist.raw$pre)
# post.missing <- setdiff(fw.meta$root_783,fw.elist.raw$post)
# pre.missing.meta <- subset(fw.meta, root_783 %in% pre.missing)
# post.missing.meta <- subset(fw.meta, root_783 %in% post.missing)

# ## Sanity check 2
# dnao1 <- subset(meta, cell_type=="DNa01")
# elist.simp.dna01 <- subset(elist.simp, post %in% dnao1$id) %>%
#   dplyr::filter(count > 25)
# table(meta$cell_type[match(elist.simp.dna01$pre, meta$id)])

# # Get FAFB data, split
# con <- DBI::dbConnect(RSQLite::SQLite(),
#                       file.path(banc.connectivity.save.path,"fafb_783_data.sqlite"))
# fw.dends <- subset(fw.meta, cell_class == "descending_neuron" | grepl("efferent",flow))$root_783
# fw.axons <- subset(fw.meta, cell_class == "ascending_neuron" | grepl("afferent",flow))$root_783
# fw.elist <- dplyr::tbl(con, "edgelist") %>%
#   dplyr::collect() %>%
#   dplyr::mutate(pre_label = dplyr::case_when(
#     pre %in% fw.dends ~ "dendrite",
#     pre %in% fw.axons ~ "axon",
#     TRUE ~ pre_label
#   )) %>%
#   dplyr::mutate(post_label = dplyr::case_when(
#     post %in% fw.dends ~ "dendrite",
#     post %in% fw.axons ~ "axon",
#     TRUE ~ post_label
#   )) %>%
#   dplyr::left_join(fw.mapping, by = c("pre" = "fafb_match")) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(pre = dplyr::coalesce(banc_id, pre)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-banc_id) %>%
#   dplyr::left_join(fw.mapping, by = c("post" = "fafb_match")) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(post = dplyr::coalesce(banc_id, post)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-banc_id, -index) %>%
#   dplyr::rename(pre_top_nt=pre_conf_nt, 
#                 pre_top_nt_p=pre_conf_nt_p,
#                 post_top_nt=post_conf_nt,
#                 post_top_nt_p=post_conf_nt_p) 
# dbDisconnect(con)

# #  Get MANC data, split
# con <- DBI::dbConnect(RSQLite::SQLite(),
#                       file.path(banc.connectivity.save.path,"manc_1.2.1_data.sqlite"))
# mc.dends <- subset(mc.meta, grepl("ascending_neuron",cell_class) | grepl("efferent",flow))$bodyid
# mc.axons <- subset(mc.meta, grepl("descending_neuron",cell_class) | grepl("afferent",flow))$bodyid
# mc.elist <- dplyr::tbl(con, "edgelist") %>%
#   dplyr::collect() %>%
#   dplyr::mutate(pre_label = dplyr::case_when(
#     pre %in% mc.dends ~ "dendrite",
#     pre %in% mc.axons ~ "axon",
#     TRUE ~ pre_label
#   )) %>%
#   dplyr::mutate(post_label = dplyr::case_when(
#     post %in% mc.dends ~ "dendrite",
#     post %in% mc.axons ~ "axon",
#     TRUE ~ post_label
#   )) %>%
#   dplyr::left_join(mc.mapping, by = c("pre" = "manc_match")) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(pre = dplyr::coalesce(banc_id, pre)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-banc_id) %>%
#   dplyr::left_join(mc.mapping, by = c("post" = "manc_match")) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(post = dplyr::coalesce(banc_id, post)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-banc_id)
# dbDisconnect(con)

# # Combine, split
# elist <- fw.elist %>%
#   rbind(mc.elist) %>%
#   dplyr::group_by(pre, post, pre_label, post_label) %>%
#   dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(pre, post, pre_label, post_label, .keep_all = TRUE) %>%
#   dplyr::group_by(post) %>%
#   dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(pre) %>%
#   dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(post, post_label) %>%
#   dplyr::mutate(post_label_count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(pre, post, pre_label, post_label) %>%
#   dplyr::mutate(norm = count/post_count) %>%
#   dplyr::ungroup()

# # Simplify
# elist.simp <- elist %>%
#   dplyr::group_by(pre, post) %>%
#   dplyr::mutate(count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-pre_label, -post_label, -post_label_count) %>%
#   dplyr::distinct(pre, post, .keep_all = TRUE) %>%
#   dplyr::group_by(post) %>%
#   dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(pre) %>%
#   dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(norm = count/post_count)

# Calculate number of connections
input.connections<- elist.simp %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(input_connections = sum(count, na.rm = TRUE)) %>%
  dplyr::mutate(input_degree = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(post, input_connections, input_degree) %>%
  dplyr::arrange(dplyr::desc(input_connections))
output.connections<- elist.simp %>%
  dplyr::group_by(pre) %>%
  dplyr::mutate(output_connections = sum(count, na.rm = TRUE)) %>%
  dplyr::mutate(output_degree = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre, output_connections, output_degree) %>%
  dplyr::arrange(dplyr::desc(output_connections))
franken.pre <- unique(elist.simp$pre)
franken.post <- unique(elist.simp$post)
franken.meta <- meta %>%
  dplyr::left_join(input.connections, by = c("id"="post")) %>%
  dplyr::left_join(output.connections, by = c("id"="pre")) %>%
  dplyr::filter(id %in% franken.pre | id %in% franken.post) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    # Special handling for id columns - concatenate unique non-NA values
    dplyr::across(c(fafb_id, manc_id), 
                  ~paste(unique(na.omit(.)), collapse = ",")),
    # Standard handling for all other columns - first non-NA value
    dplyr::across(!c(fafb_id, manc_id), 
                  ~if(all(is.na(.))) {
                    NA
                  } else {
                    first(na.omit(.))
                  })) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(bc.meta.neck[,c("root_id","an_type","dn_type","side")], by = c("banc_id"="root_id")) %>%
  dplyr::mutate(side = dplyr::case_when(
    is.na(side.y)  ~ side.x,
    is.na(side.x)  ~ side.y,
    TRUE ~ side.x
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    super_class == "ascending" ~ paste0("AN_",an_type),
    super_class == "descending" ~ paste0("DN_",dn_type),
    super_class %in% c("ascending","descending") ~ NA,
    TRUE ~ cell_sub_class
  )) %>%
  dplyr::mutate(seed = dplyr::case_when(
    super_class == "ascending" ~ paste0("AN_",an_type,"_",side),
    super_class == "descending" ~ paste0("DN_",dn_type,"_",side),
    TRUE ~ seed
  )) %>%
  dplyr::select(-side.x, -side.y, -an_type, -dn_type) %>%
  dplyr::mutate(dataset = dplyr::case_when(
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-FAFB-MANC",
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!="" ~ "BANC-FAFB",
    !is.na(banc_id)&banc_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-MANC",
    !is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "FAFB-MANC",
    !is.na(fafb_id)&fafb_id!="" ~ "FAFB",
    !is.na(manc_id)&manc_id!="" ~ "MANC",
    !is.na(banc_id)&banc_id!="" ~ "BANC",
    TRUE ~ dataset
  ))

# Save the frankenbrain
franken.meta[franken.meta==""] <- NA
arrow::write_feather(elist.simp,
                     file.path(banc.connectivity.save.path, "frankenbrain_v.1.7_edgelist_simple.feather"))
arrow::write_feather(final_mapping,
                     file.path(banc.connectivity.save.path, "frankenbrain_v.1.7_neck_bridge.feather"))
arrow::write_feather(franken.meta,
                     file.path(banc.connectivity.save.path, "frankenbrain_v.1.7_meta.feather"))

# Update in BANC data
arrow::write_feather(final_mapping,
                     file.path(banc.connectivity.save.path, "banc_neck_bridge.feather"))

# Send
readr::write_csv(meta, file = "~/BANC-project/data/meta/frankenbrain_v.1.7_meta.csv")
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/frankenbrain_v.1.7_meta.feather gs://lee-lab_brain-and-nerve-cord-fly-connectome/frankenbrain")
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/frankenbrain_v.1.7_edgelist_simple.feather gs://lee-lab_brain-and-nerve-cord-fly-connectome/frankenbrain")
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/frankenbrain_v.1.7_neck_bridge.feather gs://lee-lab_brain-and-nerve-cord-fly-connectome/frankenbrain")

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'
sync_files(path(A, "connectivity/"), path(B, "connectivity/"), extensions = c("csv", "feather", "parquet"), move.old = FALSE)
remove_empty_dirs(path(A, "connectivity/"))

# Define the remote name
remote_name <- "hms"
local_files <- "/n/data1/hms/neurobio/wilson/banc/connectivity/"
local_files <- list.files(local_files, pattern = "feather$|parquet$", full.names = TRUE)
for(local_file in local_files){
  message("Working on: ", local_file)
  remote_path <- "neuroanat/connectomes/"  # Replace with your desired path on Dropbox
  system(paste("rclone copy", local_file, paste0(remote_name, ":", remote_path)))
  cat("File transfer complete.\n")  
}

# # Check
franken.meta <- arrow::read_feather(
  file.path(banc.dropbox.connectivity.save.path, "frankenbrain_v.1.7_meta.feather"))


#################
#### Seatable ###
################# 

# NEW BANC IDS
meta.bc <- meta %>%
  dplyr::filter(!is.na(banc_id)&banc_id!="") %>%
  dplyr::rename(neuron_id = id)
meta.bc.fafb.ids <- unique(meta.bc$fafb_id)
meta.bc.manc.ids <- unique(meta.bc$manc_id)

# GET OLD
franken.meta.old <- franken_meta("SELECT _id, neuron_id, banc_id, fafb_id, manc_id, dataset from banc_meta")
franken.meta.old.bc <- franken.meta.old %>%
  dplyr::filter(fafb_id %in% meta.bc.fafb.ids|manc_id %in% meta.bc.manc.ids|(!is.na(banc_id)&banc_id!="")) %>%
  dplyr::select(`_id`, neuron_id, fafb_id, manc_id, dataset)

# CORRECT IDS
franken.meta.update <- franken.meta.old.bc %>%
  dplyr::left_join(meta %>%
                     dplyr::filter(!is.na(fafb_id)) %>%
                     dplyr::select(banc_id, fafb_id) %>%
                     dplyr::distinct(fafb_id, .keep_all = TRUE),
                   by = "fafb_id") %>%
  dplyr::left_join(meta %>%
                     dplyr::filter(!is.na(manc_id)) %>%
                     dplyr::select(banc_id, manc_id) %>%
                     dplyr::distinct(manc_id, .keep_all = TRUE),
                   by = "manc_id") %>%
  dplyr::mutate(banc_id = dplyr::case_when(
    is.na(banc_id.x) ~ banc_id.y,
    TRUE ~ banc_id.x
  )) %>%
  dplyr::mutate(dataset = dplyr::case_when(
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-FAFB-MANC",
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!="" ~ "BANC-FAFB",
    !is.na(banc_id)&banc_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-MANC",
    !is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "FAFB-MANC",
    !is.na(fafb_id)&fafb_id!="" ~ "FAFB",
    !is.na(manc_id)&manc_id!="" ~ "MANC",
    !is.na(banc_id)&banc_id!="" ~ "BANC",
    TRUE ~ dataset
  )) %>%
  dplyr::mutate(neuron_id = dplyr::case_when(
    !is.na(banc_id) ~ banc_id,
    !is.na(fafb_id) ~ fafb_id,
    !is.na(manc_id) ~ manc_id,
    TRUE ~ neuron_id
  )) %>%
  dplyr::select(`_id`, neuron_id, banc_id, fafb_id, manc_id, dataset)

# UPDATE
franken.meta.update[is.na(franken.meta.update)] <- ""
banctable_update_rows(base='cns_meta', 
                      table = 'franken_meta', 
                      df = franken.meta.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

# ################
# ### Seatable ###
# ################ add banc supervoxel_id
# 
# # FAFB orig
# fw.meta.orig <- readr::read_tsv("/Users/papers/BANC-project/data/meta/fafb_schlegel_et_al_2024_meta.tsv", 
#                                  col_types = banc.col.types)
# fw.meta.orig <- fw.meta.orig %>%
#   dplyr::mutate(soma = !is.na(soma_x)&soma_x!="") %>%
#   dplyr::select(-soma_x, -soma_y, -soma_z, -pos_x, -pos_y, -pos_z, -known_nt_source, -known_nt)
# colnames(fw.meta.orig) <- paste0("FAFB_",colnames(fw.meta.orig))
# 
# # MANC orig
# malevnc:::choose_malevnc_dataset('MANC')
# mc.find <- neuprintr::neuprint_search("Traced",field="status",dataset="manc:v1.2.1")
# mc.ids <- unique(mc.find$bodyid)
# mc.meta.orig <- manc_neuprint_meta(mc.ids)
# readr::write_csv(mc.meta.orig,"/Users/papers/BANC-project/data/meta/manc_v.1.2_meta.csv")
# mc.meta.orig <- mc.meta.orig %>%
#   dplyr::select(-avgLocation, -somaLocation,
#                 -post, -pre, -downstream, -upstream, -synweight)
# colnames(mc.meta.orig) <- paste0("MANC_",colnames(mc.meta.orig))
# 
# # Save for basis of CNS seatable
# cns.meta <- meta %>%
#   dplyr::filter(super_class!="glia", dataset!="") %>%
#   dplyr::distinct(fafb_id, manc_id, .keep_all = TRUE) %>%
#   dplyr::select(-soma_volume, -cable_length, -top_nt_p, -soma_dcv_count, -soma_dcv_density, -cell_dcv_count, -cell_dcv_density) %>%
#   dplyr::left_join(fw.meta.orig, by = c("fafb_id"="FAFB_root_id")) %>%
#   dplyr::left_join(mc.meta.orig, by = c("manc_id"="MANC_bodyid")) %>%
#   dplyr::rename(neuron_id = id) %>%
#   as.data.frame()
# cns.meta[is.na(cns.meta)] <- ""
# readr::write_csv(cns.meta, file.path(banc.path,"frankenbrain_v.1.7_meta.csv"))
# banctable_append_rows(base = 'cns_meta',
#                       table = "franken_meta",
#                       df = subset(cns.meta, region=="optic_lobes"),
#                       chunksize = 300,
#                       bigdata = TRUE)
# banctable_move_to_bigdata(base='cns_meta',
#                           table = "franken_meta",
#                           where = "`region` = 'optic_lobes'",
#                           invert = FALSE)
# # fbc <- banctable_query("SELECT * FROM franken_meta", 
# #                        base = "cns_meta") 

