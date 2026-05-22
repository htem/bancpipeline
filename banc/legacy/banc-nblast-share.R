###########################################
### Compile NBLAST from various sources ###
###########################################
source("banc/banc-startup.R")
message("##### Building BANC NBLAST data files #####")
version <- "elastix_tpsreg_240721"

# Register cores
cl <- setup_parallel()

##########################
### DELETE WRONG MATCH ###
##########################

# Read latest seatable data
bancr:::banctable_updateids()
bc <- banctable_query()
conflicts <- data.frame()

# ### Conflicts files ###
# conflicts <- suppressWarnings(readr::read_csv(file = file.path("tracing","conflicts","2024-12-04_matching_conflicts.csv"), 
#                                                       col_types = banc.col.types,
#                                                       show_col_types = FALSE)) %>%
#   dplyr::select(root_id, current, suggested, fafb_match, manc_match) %>%
#   dplyr::mutate(manc_png_match = dplyr::case_when(
#     !is.na(manc_match) & suggested==manc_match ~ current,
#     !is.na(manc_match) & current==manc_match ~ suggested,
#     TRUE ~ NA
#   )) %>%
#   dplyr::mutate(fafb_png_match = dplyr::case_when(
#     !is.na(fafb_match) & suggested==fafb_match ~ current,
#     !is.na(fafb_match) & current==fafb_match ~ suggested,
#     TRUE ~ NA
#   )) %>%
#   dplyr::distinct(root_id, fafb_png_match, manc_png_match) %>%
#   dplyr::mutate(valid = 'f')

### MANC ###

# Read wrong matches
bc.manc.png.wrong <- bc %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(grepl("MANC_PNG_MATCH_WRONG", status), 
                !is.na(manc_png_match),
                !is.na(root_id)) %>%
  dplyr::mutate(valid = 'f') %>%
  dplyr::select(`_id`,status,root_id,
                manc_match,manc_png_match,
                cell_type, manc_cell_type,
                valid) %>%
  plyr::rbind.fill(conflicts) %>%
  dplyr::filter(!is.na(manc_png_match),
                !is.na(root_id)) %>%
  as.data.frame()

# Set to false in our matching file
if(nrow(bc.manc.png.wrong)){
  banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.feather")) %>%
    dplyr::anti_join(bc.manc.png.wrong, by = c("match_id"="manc_png_match", "pt_root_id"="root_id"))
  manc.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"),
                                           col_types = banc.col.types, 
                                           show_col_types = FALSE) %>%
    dplyr::left_join(bc.manc.png.wrong %>%
                       dplyr::select(root_id,manc_png_match,valid), 
                     by = c("match_id"="manc_png_match", "pt_root_id"="root_id")) %>%
    dplyr::mutate(valid = dplyr::case_when(
      !is.na(valid.y) ~ valid.y,
      TRUE ~ valid.x
    )) %>%
    dplyr::select(-valid.x,-valid.y)
  banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                               col_types = banc.col.types, 
                               show_col_types = FALSE) %>%
    dplyr::mutate(manc_auto_cell_type = ifelse(root_id %in% bc.manc.png.wrong$root_id, NA, manc_auto_cell_type)) %>%
    dplyr::mutate(manc_auto_cell_type = ifelse(supervoxel_id %in% bc.manc.png.wrong$supervoxel_id, NA, manc_auto_cell_type))
  arrow::write_feather(banc.meta.manc.nb, file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.feather"))
  readr::write_csv(manc.matches.df.valid,file = file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"))  
  readr::write_csv(banc.meta,file = file.path(banc.meta.save.path,"banc_meta.csv"))  
  
  # Remove saved png file
  A <- '/n/data1/hms/neurobio/wilson/banc/'
  B <- '/n/files/Neurobio/wilsonlab/banc/'
  sync_files(path(B, "matching/manc/correct/"), path(A, "matching/manc/correct/"), extensions = c("png"), move.old = FALSE)
  inventory <- list.files(banc.manc.correct.match.path, recursive = TRUE, full.names = TRUE)
  for(i in 1:nrow(bc.manc.png.wrong)){
    query <- bc.manc.png.wrong[i,"root_id"]
    match <- bc.manc.png.wrong[i,"manc_png_match"]
    poss <- inventory[grepl(match,inventory)]
    ids <- regmatches(basename(poss), regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE))
    ids <- banc_latestid(ids)
    del <- which(ids==query)
    if(sum(del)){
      file.remove(poss[del])
      message("deleting: ", paste(basename(poss[del]),collapse=", "),"\n") 
    }
  }
  sync_files(path(A, "matching/manc/correct/"), path(B, "matching/manc/correct/"), extensions = c("png"), move.old = FALSE)
  
  # Update seatable
  bc.new <- bc.manc.png.wrong %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      cell_type = dplyr::case_when(
        is.na(manc_cell_type)|manc_cell_type=="" ~ cell_type,
        is.na(cell_type)|cell_type=="" ~ cell_type,
        is.na(manc_png_match)|manc_png_match=="" ~ cell_type,
        is.na(manc_match)|manc_match=="" ~ cell_type,
        manc_match!=manc_png_match  ~ cell_type,
        (cell_type==manc_cell_type) & (manc_match==manc_png_match) ~ "",
        TRUE ~ cell_type
      ),
      manc_cell_type = dplyr::case_when(
        is.na(manc_cell_type)|manc_cell_type=="" ~ manc_cell_type,
        is.na(manc_png_match)|manc_png_match=="" ~ manc_cell_type,
        is.na(manc_match)|manc_match=="" ~ "",
        manc_match!=manc_png_match  ~ manc_cell_type,
        (manc_match==manc_png_match) ~ "",
        TRUE ~ manc_cell_type
      ),
      manc_match = dplyr::case_when(
        is.na(manc_png_match)|manc_png_match=="" ~ manc_match,
        is.na(manc_match)|manc_match=="" ~ manc_match,
        (manc_match==manc_png_match) ~ "",
        TRUE ~ manc_match
      )
    ) %>%
    dplyr::mutate(manc_png_match = '' ) %>%
    dplyr::mutate(status = subtract_status(status,"MANC_PNG_MATCH_WRONG")) %>%
    dplyr::distinct(`_id`, manc_match, manc_png_match, cell_type, manc_cell_type, status) %>%
    as.data.frame()
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.new, 
                        append_allowed = FALSE, 
                        chunksize = 1000)
}

### maleCNS ###

# Read wrong matches
bc.malecns.png.wrong <- bc %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(grepl("MALECNS_PNG_MATCH_WRONG", status),
                !is.na(malecns_png_match),
                !is.na(root_id)) %>%
  dplyr::mutate(valid = 'f') %>%
  dplyr::select(`_id`,status,root_id,
                malecns_match,malecns_png_match,
                cell_type, malecns_cell_type,
                valid) %>%
  plyr::rbind.fill(conflicts) %>%
  dplyr::filter(!is.na(malecns_png_match),
                !is.na(root_id)) %>%
  as.data.frame()

# Set to false in our matching file
if(nrow(bc.malecns.png.wrong)){
  banc.meta.malecns.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.feather")) %>%
    dplyr::anti_join(bc.malecns.png.wrong, by = c("match_id"="malecns_png_match", "pt_root_id"="root_id"))
  malecns.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_malecns_reviewed_matches.csv"),
                                              col_types = banc.col.types,
                                              show_col_types = FALSE) %>%
    dplyr::left_join(bc.malecns.png.wrong %>%
                       dplyr::select(root_id,malecns_png_match,valid),
                     by = c("match_id"="malecns_png_match", "pt_root_id"="root_id")) %>%
    dplyr::mutate(valid = dplyr::case_when(
      !is.na(valid.y) ~ valid.y,
      TRUE ~ valid.x
    )) %>%
    dplyr::select(-valid.x,-valid.y)
  banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                               col_types = banc.col.types,
                               show_col_types = FALSE) %>%
    dplyr::mutate(malecns_auto_cell_type = ifelse(root_id %in% bc.malecns.png.wrong$root_id, NA, malecns_auto_cell_type)) %>%
    dplyr::mutate(malecns_auto_cell_type = ifelse(supervoxel_id %in% bc.malecns.png.wrong$supervoxel_id, NA, malecns_auto_cell_type))
  arrow::write_feather(banc.meta.malecns.nb, file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.feather"))
  readr::write_csv(malecns.matches.df.valid,file = file.path(banc.meta.save.path,"banc_malecns_reviewed_matches.csv"))
  readr::write_csv(banc.meta,file = file.path(banc.meta.save.path,"banc_meta.csv"))

  # Remove saved png file
  A <- '/n/data1/hms/neurobio/wilson/banc/'
  B <- '/n/files/Neurobio/wilsonlab/banc/'
  sync_files(path(B, "matching/malecns/correct/"), path(A, "matching/malecns/correct/"), extensions = c("png"), move.old = FALSE)
  inventory <- list.files(banc.malecns.correct.match.path, recursive = TRUE, full.names = TRUE)
  for(i in 1:nrow(bc.malecns.png.wrong)){
    query <- bc.malecns.png.wrong[i,"root_id"]
    match <- bc.malecns.png.wrong[i,"malecns_png_match"]
    poss <- inventory[grepl(match,inventory)]
    ids <- regmatches(basename(poss), regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE))
    ids <- banc_latestid(ids)
    del <- which(ids==query)
    if(sum(del)){
      file.remove(poss[del])
      message("deleting: ", paste(basename(poss[del]),collapse=", "),"\n")
    }
  }
  sync_files(path(A, "matching/malecns/correct/"), path(B, "matching/malecns/correct/"), extensions = c("png"), move.old = FALSE)

  # Update seatable
  bc.new <- bc.malecns.png.wrong %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      cell_type = dplyr::case_when(
        is.na(malecns_cell_type)|malecns_cell_type=="" ~ cell_type,
        is.na(cell_type)|cell_type=="" ~ cell_type,
        is.na(malecns_png_match)|malecns_png_match=="" ~ cell_type,
        is.na(malecns_match)|malecns_match=="" ~ cell_type,
        malecns_match!=malecns_png_match  ~ cell_type,
        (cell_type==malecns_cell_type) & (malecns_match==malecns_png_match) ~ "",
        TRUE ~ cell_type
      ),
      malecns_cell_type = dplyr::case_when(
        is.na(malecns_cell_type)|malecns_cell_type=="" ~ malecns_cell_type,
        is.na(malecns_png_match)|malecns_png_match=="" ~ malecns_cell_type,
        is.na(malecns_match)|malecns_match=="" ~ "",
        malecns_match!=malecns_png_match  ~ malecns_cell_type,
        (malecns_match==malecns_png_match) ~ "",
        TRUE ~ malecns_cell_type
      ),
      malecns_match = dplyr::case_when(
        is.na(malecns_png_match)|malecns_png_match=="" ~ malecns_match,
        is.na(malecns_match)|malecns_match=="" ~ malecns_match,
        (malecns_match==malecns_png_match) ~ "",
        TRUE ~ malecns_match
      )
    ) %>%
    dplyr::mutate(malecns_png_match = '' ) %>%
    dplyr::mutate(status = subtract_status(status,"MALECNS_PNG_MATCH_WRONG")) %>%
    dplyr::distinct(`_id`, malecns_match, malecns_png_match, cell_type, malecns_cell_type, status) %>%
    as.data.frame()
  banctable_update_rows(base='banc_meta',
                        table = "banc_meta",
                        df = bc.new,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

### FAFB ###

# Read wrong matches
bc.fafb.png.wrong <- bc %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(grepl("FAFB_PNG_MATCH_WRONG", status),
                !is.na(fafb_png_match),
                !is.na(root_id)) %>%
  dplyr::mutate(valid = 'f') %>%
  dplyr::select(`_id`,status,root_id,
                fafb_match,fafb_png_match,
                cell_type, fafb_cell_type,
                valid) %>%
  plyr::rbind.fill(conflicts) %>%
  dplyr::filter(!is.na(fafb_png_match),
                !is.na(root_id)) %>%
  as.data.frame()

# Set to false in our matching file
if(nrow(bc.fafb.png.wrong)){
  banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_fafb_783_nblast.feather")) %>%
    dplyr::anti_join(bc.fafb.png.wrong, by = c("match_id"="fafb_png_match", "pt_root_id"="root_id"))
  fafb.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"),
                                           col_types = banc.col.types, 
                                           show_col_types = FALSE) %>%
    dplyr::left_join(bc.fafb.png.wrong %>%
                       dplyr::select(root_id,fafb_png_match,valid), 
                     by = c("match_id"="fafb_png_match", "pt_root_id"="root_id")) %>%
    dplyr::mutate(valid = dplyr::case_when(
      !is.na(valid.y) ~ valid.y,
      TRUE ~ valid.x
    )) %>%
    dplyr::select(-valid.x,-valid.y)
  banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                               col_types = banc.col.types, 
                               show_col_types = FALSE) %>%
    dplyr::mutate(fafb_auto_cell_type = ifelse(root_id %in% bc.fafb.png.wrong$root_id, NA, fafb_auto_cell_type)) %>%
    dplyr::mutate(fafb_auto_cell_type = ifelse(supervoxel_id %in% bc.fafb.png.wrong$supervoxel_id, NA, fafb_auto_cell_type))
  arrow::write_feather(banc.meta.fafb.nb, file.path(banc.meta.save.path,"banc_fafb_783_nblast.feather"))
  readr::write_csv(fafb.matches.df.valid,file = file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"))  
  readr::write_csv(banc.meta,file = file.path(banc.meta.save.path,"banc_meta.csv"))  
  
  # Remove saved png file
  A <- '/n/data1/hms/neurobio/wilson/banc/'
  B <- '/n/files/Neurobio/wilsonlab/banc/'
  sync_files(path(B, "matching/fafb/correct/"), path(A, "matching/fafb/correct/"), extensions = c("png"), move.old = FALSE)
  inventory <- list.files(banc.fafb.correct.match.path, recursive = TRUE, full.names = TRUE)
  for(i in 1:nrow(bc.fafb.png.wrong)){
    query <- bc.fafb.png.wrong[i,"root_id"]
    match <- bc.fafb.png.wrong[i,"fafb_png_match"]
    poss <- inventory[grepl(match,inventory)]
    ids <- regmatches(basename(poss), regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE))
    ids <- banc_latestid(ids)
    del <- which(ids==query)
    if(sum(del)){
      file.remove(poss[del])
      message("deleting: ", paste(basename(poss[del]),collapse=", "),"\n") 
    }
  }
  sync_files(path(A, "matching/fafb/correct/"), path(B, "matching/fafb/correct/"), extensions = c("png"), move.old = FALSE)
  
  # Update seatable
  bc.new <- bc.fafb.png.wrong %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      cell_type = dplyr::case_when(
        is.na(fafb_cell_type)|fafb_cell_type=="" ~ cell_type,
        is.na(cell_type)|cell_type=="" ~ cell_type,
        is.na(fafb_png_match)|fafb_png_match=="" ~ cell_type,
        is.na(fafb_match)|fafb_match=="" ~ cell_type,
        fafb_match!=fafb_png_match  ~ cell_type,
        (cell_type==fafb_cell_type) & (fafb_match==fafb_png_match) ~ "",
        TRUE ~ cell_type
      ),
      fafb_cell_type = dplyr::case_when(
        is.na(fafb_cell_type)|fafb_cell_type=="" ~ fafb_cell_type,
        is.na(fafb_png_match)|fafb_png_match=="" ~ fafb_cell_type,
        is.na(fafb_match)|fafb_match=="" ~ "",
        fafb_match!=fafb_png_match  ~ fafb_cell_type,
        (fafb_match==fafb_png_match) ~ "",
        TRUE ~ fafb_cell_type
      ),
      fafb_match = dplyr::case_when(
        is.na(fafb_png_match)|fafb_png_match=="" ~ fafb_match,
        is.na(fafb_match)|fafb_match=="" ~ fafb_match,
        (fafb_match==fafb_png_match) ~ "",
        TRUE ~ fafb_match
      )
    ) %>%
    dplyr::mutate(fafb_png_match = '' ) %>%
    dplyr::mutate(status = subtract_status(status,"FAFB_PNG_MATCH_WRONG")) %>%
    dplyr::distinct(`_id`, fafb_match, fafb_png_match,cell_type, fafb_cell_type, status) %>%
    as.data.frame()
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.new, 
                        append_allowed = FALSE, 
                        chunksize = 1000)
}


### hemibrain ###

# Read wrong matches
bc.hemibrain.png.wrong <- bc %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(grepl("HEMIBRAIN_PNG_MATCH_WRONG", status), 
                !is.na(hemibrain_png_match),
                !is.na(root_id)) %>%
  dplyr::mutate(valid = 'f') %>%
  dplyr::select(`_id`,status,root_id,
                hemibrain_match,hemibrain_png_match,
                cell_type, hemibrain_cell_type,
                valid) %>%
  plyr::rbind.fill(conflicts) %>%
  dplyr::filter(!is.na(hemibrain_png_match),
                !is.na(root_id)) %>%
  as.data.frame()

# Set to false in our matching file
if(nrow(bc.hemibrain.png.wrong)){
  banc.meta.hemibrain.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.feather")) %>%
    dplyr::anti_join(bc.hemibrain.png.wrong, by = c("match_id"="hemibrain_png_match", "pt_root_id"="root_id"))
  hemibrain.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"),
                                                col_types = banc.col.types, 
                                                show_col_types = FALSE) %>%
    dplyr::left_join(bc.hemibrain.png.wrong %>%
                       dplyr::select(root_id,hemibrain_png_match,valid), 
                     by = c("match_id"="hemibrain_png_match", "pt_root_id"="root_id")) %>%
    dplyr::mutate(valid = dplyr::case_when(
      !is.na(valid.y) ~ valid.y,
      TRUE ~ valid.x
    )) %>%
    dplyr::select(-valid.x,-valid.y)
  banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                               col_types = banc.col.types, 
                               show_col_types = FALSE) %>%
    dplyr::mutate(hemibrain_auto_cell_type = ifelse(root_id %in% bc.hemibrain.png.wrong$root_id, NA, hemibrain_auto_cell_type)) %>%
    dplyr::mutate(hemibrain_auto_cell_type = ifelse(supervoxel_id %in% bc.hemibrain.png.wrong$supervoxel_id, NA, hemibrain_auto_cell_type))
  arrow::write_feather(banc.meta.hemibrain.nb, file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.feather"))
  readr::write_csv(hemibrain.matches.df.valid,file = file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"))  
  readr::write_csv(banc.meta,file = file.path(banc.meta.save.path,"banc_meta.csv"))  
  
  # Remove saved png file
  A <- '/n/data1/hms/neurobio/wilson/banc/'
  B <- '/n/files/Neurobio/wilsonlab/banc/'
  sync_files(path(B, "matching/hemibrain/correct/"), path(A, "matching/hemibrain/correct/"), extensions = c("png"), move.old = FALSE)
  inventory <- list.files(banc.hemibrain.correct.match.path, recursive = TRUE, full.names = TRUE)
  for(i in 1:nrow(bc.hemibrain.png.wrong)){
    query <- bc.hemibrain.png.wrong[i,"root_id"]
    match <- bc.hemibrain.png.wrong[i,"hemibrain_png_match"]
    poss <- inventory[grepl(match,inventory)]
    ids <- regmatches(basename(poss), regexpr("(?<=root_id_)\\d+", basename(poss), perl = TRUE))
    ids <- banc_latestid(ids)
    del <- which(ids==query)
    if(sum(del)){
      file.remove(poss[del])
      message("deleting: ", paste(basename(poss[del]),collapse=", "),"\n") 
    }
  }
  sync_files(path(A, "matching/hemibrain/correct/"), path(B, "matching/hemibrain/correct/"), extensions = c("png"), move.old = FALSE)
  
  # Update seatable
  bc.new <- bc.hemibrain.png.wrong %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      cell_type = dplyr::case_when(
        is.na(hemibrain_cell_type)|hemibrain_cell_type=="" ~ cell_type,
        is.na(cell_type)|cell_type=="" ~ cell_type,
        is.na(hemibrain_png_match)|hemibrain_png_match=="" ~ cell_type,
        is.na(hemibrain_match)|hemibrain_match=="" ~ cell_type,
        hemibrain_match!=hemibrain_png_match  ~ cell_type,
        (cell_type==hemibrain_cell_type) & (hemibrain_match==hemibrain_png_match) ~ "",
        TRUE ~ cell_type
      ),
      hemibrain_cell_type = dplyr::case_when(
        is.na(hemibrain_cell_type)|hemibrain_cell_type=="" ~ hemibrain_cell_type,
        is.na(hemibrain_png_match)|hemibrain_png_match=="" ~ hemibrain_cell_type,
        is.na(hemibrain_match)|hemibrain_match=="" ~ "",
        hemibrain_match!=hemibrain_png_match  ~ hemibrain_cell_type,
        (hemibrain_match==hemibrain_png_match) ~ "",
        TRUE ~ hemibrain_cell_type
      ),
      hemibrain_match = dplyr::case_when(
        is.na(hemibrain_png_match)|hemibrain_png_match=="" ~ hemibrain_match,
        is.na(hemibrain_match)|hemibrain_match=="" ~ hemibrain_match,
        (hemibrain_match==hemibrain_png_match) ~ "",
        TRUE ~ hemibrain_match
      )
    ) %>%
    dplyr::mutate(hemibrain_png_match = '' ) %>%
    dplyr::mutate(status = subtract_status(status,"hemibrain_PNG_MATCH_WRONG")) %>%
    dplyr::distinct(`_id`, hemibrain_match, hemibrain_png_match, cell_type, hemibrain_cell_type, status) %>%
    as.data.frame()
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.new, 
                        append_allowed = FALSE, 
                        chunksize = 1000)
}

##########################
### COMPILE MATCH DATA ###
##########################

# Older results
banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_fafb_783_nblast.feather"))
banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.feather"))
banc.meta.fanc.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_fanc_1116_nblast.feather"))
banc.meta.hemibrain.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.feather"))
banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_mirror_nblast.feather"))
mirror.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_mirror_reviewed_matches.csv"),
                                           col_types = banc.col.types, 
                                           show_col_types = FALSE)
fafb.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"),
                                         col_types = banc.col.types, 
                                         show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id))
manc.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"),
                                         col_types = banc.col.types, 
                                         show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))
hemibrain.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"),
                                              col_types = banc.col.types,
                                              show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))
fanc.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fanc_reviewed_matches.csv"),
                                         col_types = banc.col.types,
                                         show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))
banc.meta.malecns.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.feather"))
malecns.matches.df.valid <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_malecns_reviewed_matches.csv"),
                                            col_types = banc.col.types,
                                            show_col_types = FALSE)  %>%
  dplyr::filter(!grepl("\\_",match_id))

##############################
### Update verified matches ##
##############################

banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(fafb.matches.df.valid[,c("pt_root_id","match_id","valid")], 
                   by=c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))

banc.meta.manc.nb <- banc.meta.manc.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(manc.matches.df.valid[,c("pt_root_id","match_id","valid")], 
                   by=c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))

banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(hemibrain.matches.df.valid[,c("pt_root_id","match_id","valid")], 
                   by=c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))

banc.meta.fanc.nb <- banc.meta.fanc.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(fanc.matches.df.valid[,c("pt_root_id","match_id","valid")],
                   by=c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))

banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
  dplyr::select(-valid) %>%
  dplyr::left_join(malecns.matches.df.valid[,c("pt_root_id","match_id","valid")],
                   by=c("pt_root_id", "match_id"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))

#####################
### Get meta data ###
#####################

# Read BANC meta local
banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                             col_types = banc.col.types, 
                             show_col_types = FALSE)
banc.meta <- banc_updateids(banc.meta)
readr::write_csv(x = banc.meta, file=file.path(banc.meta.save.path,"banc_ids.csv"))
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
hb.ids <- unique(hb.meta$bodyid)
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(fc.meta$cell_id)
mcns.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"malecns_09_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
mcns.ids <- unique(mcns.meta$malecns_09_id)

# Get NBLAST versions
fafb.nblast.folder <- file.path(banc.nblast.fafb.save.path,"results",version)
mirror.nblast.folder <- file.path(banc.nblast.mirror.save.path,"results")
manc.nblast.folder <- file.path(banc.nblast.manc.save.path,"results",version)
fanc.nblast.folder <- file.path(banc.nblast.fanc.save.path,"results",version)
hemibrain.nblast.folder <- file.path(banc.nblast.hemibrain.save.path,"results_with_mirrored",version)
malecns.nblast.folder <- file.path(banc.nblast.malecns.save.path,"results","navis_tpsreg_250206")

#######################################
### Collect manual matching results ###
#######################################

# Mirror
message("Working on BANC-mirror reviewed matches ...")
matches1 <- list.files(file.path(banc.mirror.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.mirror.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.mirror.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.mirror.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.mirror.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.mirror.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.mirror.correct.match.path,"blue_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.mirror.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.mirror.correct.match.path,"investigate"), pattern = "png")
mirror.matches <- unique(c(matches1,
                           matches2,
                           matches3,
                           matches4,
                           matches5))
mirror.matches.df <- data.frame()
for(file in mirror.matches){
  mdf <- data.frame(root_id = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*","",file), 
                    supervoxel_id = ifelse(grepl("supervoxel_id",file),gsub(".*_supervoxel_id_|_hit.*|_query.*","",file),NA),
                    banc_match = gsub(".*_hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*","",file))
  mirror.matches.df <- rbind(mirror.matches.df,mdf)
}
mirror.matches.df <- mirror.matches.df %>%
  dplyr::anti_join(mirror.matches.df.valid %>% dplyr::filter(valid=='t'), by = c("banc_match"="match_id", "root_id"="query_id"))
if(nrow(mirror.matches.df)){
  
  # Add reverse matches
  mirror.matches.df.missing <- mirror.matches.df %>%
    dplyr::filter(!banc_match%in%root_id) %>%
    dplyr::rename(root_id = banc_match,
                  banc_match = root_id)
  mirror.matches.df <- rbind(mirror.matches.df,mirror.matches.df.missing) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  
  # Update IDs
  mirror.matches.df$query <- mirror.matches.df$root_id
  mirror.matches.df$match_id <- mirror.matches.df$banc_match
  mirror.matches.df$root_id <- banc_updateids(mirror.matches.df$root_id)
  mirror.matches.df$banc_match <- banc_updateids(mirror.matches.df$banc_match)
  
  # Update ids
  mirror.matches.df$position <- banc.meta$position[match(mirror.matches.df$root_id,banc.meta$root_id)]
  mirror.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(mirror.matches.df$root_id,banc.meta$root_id)]
  mirror.matches.df$banc_match_position <- banc.meta$position[match(mirror.matches.df$banc_match,banc.meta$root_id)]
  mirror.matches.df$banc_match_supervoxel_id <- banc.meta$supervoxel_id[match(mirror.matches.df$banc_match,banc.meta$root_id)]
  
  # Organise
  mirror.matches.df <- mirror.matches.df %>%
    dplyr::rename(pt_root_id = root_id, 
                  pt_supervoxel_id = supervoxel_id, 
                  pt_position = position, 
                  match_root_id = banc_match,
                  match_supervoxel_id = banc_match_supervoxel_id,
                  match_position = banc_match_position,
                  query_id = query) %>%
    dplyr::select(pt_root_id, 
                  pt_supervoxel_id, 
                  pt_position, 
                  match_root_id,
                  match_supervoxel_id,
                  match_position,
                  query_id,
                  match_id) %>%
    plyr::rbind.fill(mirror.matches.df.valid) %>%
    dplyr::mutate(valid = 't') %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::filter(pt_root_id!=match_root_id) %>%
    dplyr::distinct(pt_root_id, match_root_id, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(match_supervoxel_id))
  
  # Save
  mirror.matches.df <- banc_updateids(mirror.matches.df,
                                      root.column = "pt_root_id",
                                      position.column = "pt_position",
                                      supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(mirror.matches.df, file.path(banc.meta.save.path,"banc_mirror_reviewed_matches.csv"))
}  

# FAFB
message("Working on BANC-FAFB reviewed matches ...")
matches1 <- list.files(file.path(banc.fafb.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.fafb.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.fafb.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.fafb.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.fafb.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.fafb.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.fafb.correct.match.path,"blue_wrong"), pattern = "png")
green_wrong <- list.files(file.path(banc.fafb.correct.match.path,"green_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.fafb.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.fafb.correct.match.path,"investigate"), pattern = "png")
fafb.matches <- unique(c(matches1,
                         matches2,
                         matches3,
                         matches4,
                         matches5))
fafb.matches.df <- data.frame()
for(file in fafb.matches){
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*","",file), 
                    supervoxel_ids = ifelse(grepl("supervoxel_id",file),gsub(".*_supervoxel_id_|_hit.*|_query.*","",file),NA),
                    fafb_match = gsub(".*root_783_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png","",file))
  fafb.matches.df <- rbind(fafb.matches.df,mdf)
}
fafb.matches.df <- fafb.matches.df %>%
  dplyr::filter(!grepl("\\_",fafb_match)) %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids)  %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
fafb.matches.df <- fafb.matches.df %>%
  dplyr::anti_join(fafb.matches.df.valid %>% dplyr::filter(valid=='t'), by = c("fafb_match"="match_id", "root_id"="query_id"))
if(nrow(fafb.matches.df)){
  
  # Update IDs
  fafb.matches.df$query <- fafb.matches.df$root_id
  fafb.matches.df <- banc_updateids(fafb.matches.df, 
                                    root.column = "root_id",
                                    supervoxel.column = "supervoxel_id")
  
  # Assign needed data
  fafb.matches.df$position <- banc.meta$position[match(fafb.matches.df$root_id,banc.meta$root_id)]
  fafb.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(fafb.matches.df$root_id,banc.meta$root_id)]
  
  # Organise
  fafb.matches.df <- fafb.matches.df %>%
    dplyr::left_join(fw.meta[,c("root_783","cell_type")],by=c("fafb_match"="root_783")) %>%
    dplyr::rename(pt_root_id = root_id, 
                  pt_supervoxel_id = supervoxel_id, 
                  pt_position = position, 
                  query_id = query,
                  match_id = fafb_match,
                  match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, 
                  pt_supervoxel_id, 
                  pt_position, 
                  query_id,
                  match_id,
                  match_cell_type) %>%
    dplyr::mutate(valid = 't') %>%
    rbind(fafb.matches.df.valid) %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)
  
  # Save
  fafb.matches.df <- banc_updateids(fafb.matches.df, 
                                    root.column = "pt_root_id", 
                                    position.column = "pt_position", 
                                    supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(fafb.matches.df, file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"))
}

# MANC
message("Working on BANC-MANC reviewed matches ...")
matches1 <- list.files(file.path(banc.manc.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.manc.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.manc.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.manc.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.manc.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.manc.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.manc.correct.match.path,"blue_wrong"), pattern = "png")
green_wrong <- list.files(file.path(banc.manc.correct.match.path,"green_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.manc.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.manc.correct.match.path,"investigate"), pattern = "png")
manc.matches <- unique(c(matches1,
                         matches2,
                         matches3,
                         matches4,
                         matches5))
manc.matches.df <- data.frame()
for(file in manc.matches){
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*","",file), 
                    supervoxel_ids = ifelse(grepl("supervoxel_id",file),gsub(".*_supervoxel_id_|_hit.*|_query.*","",file),NA),
                    manc_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png","",file))
  manc.matches.df <- rbind(manc.matches.df,mdf)
}
manc.matches.df <- manc.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids)  %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
manc.matches.df <- manc.matches.df %>%
  dplyr::anti_join(manc.matches.df.valid %>% dplyr::filter(valid=='t'), 
                   by = c("manc_match"="match_id", "root_id"="query_id"))
if(nrow(manc.matches.df)){
  
  # Update IDs
  manc.matches.df$query <- manc.matches.df$root_id
  manc.matches.df <- banc_updateids(manc.matches.df, 
                                    root.column = "root_id",
                                    supervoxel.column = "supervoxel_id")
  
  # Assign needed data
  manc.matches.df$position <- banc.meta$position[match(manc.matches.df$root_id,banc.meta$root_id)]
  manc.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(manc.matches.df$root_id,banc.meta$root_id)]
  
  # Organise
  manc.matches.df <- manc.matches.df %>%
    dplyr::left_join(mc.meta[,c("bodyid","cell_type")],by=c("manc_match"="bodyid")) %>%
    dplyr::rename(pt_root_id = root_id, 
                  pt_supervoxel_id = supervoxel_id, 
                  pt_position = position, 
                  query_id = query,
                  match_id = manc_match,
                  match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, 
                  pt_supervoxel_id, 
                  pt_position, 
                  query_id,
                  match_id,
                  match_cell_type) %>%
    dplyr::mutate(valid = 't') %>%
    rbind(manc.matches.df.valid) %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)
  
  # Save
  manc.matches.df <- banc_updateids(manc.matches.df, 
                                    root.column = "pt_root_id", 
                                    position.column = "pt_position", 
                                    supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(manc.matches.df, file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"))
}

# hemibrain
message("Working on BANC-hemibrain reviewed matches ...")
matches1 <- list.files(file.path(banc.hemibrain.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.hemibrain.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.hemibrain.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.hemibrain.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.hemibrain.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.hemibrain.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"blue_wrong"), pattern = "png")
green_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"green_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.hemibrain.correct.match.path,"investigate"), pattern = "png")
hemibrain.matches <- unique(c(matches1,
                              matches2,
                              matches3,
                              matches4))
hemibrain.matches.df <- data.frame()
for(file in hemibrain.matches){
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id_.*","",file), 
                    supervoxel_ids = ifelse(grepl("supervoxel_id",file),gsub(".*_supervoxel_id_|_hit.*|_query.*","",file),NA),
                    hemibrain_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*|\\.png","",file))
  hemibrain.matches.df <- rbind(hemibrain.matches.df,mdf)
}
hemibrain.matches.df <- hemibrain.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
hemibrain.matches.df <- hemibrain.matches.df %>%
  dplyr::anti_join(hemibrain.matches.df.valid %>% dplyr::filter(valid=='t'), by = c("hemibrain_match"="match_id","root_id"="query_id"))
if(nrow(hemibrain.matches.df)){
  
  # Update IDs
  hemibrain.matches.df$query <- hemibrain.matches.df$root_id
  hemibrain.matches.df <- banc_updateids(hemibrain.matches.df, 
                                         root.column = "root_id",
                                         supervoxel.column = "supervoxel_id")
  
  # Assign needed data
  hemibrain.matches.df$position <- banc.meta$position[match(hemibrain.matches.df$root_id,banc.meta$root_id)]
  hemibrain.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(hemibrain.matches.df$root_id,banc.meta$root_id)]
  
  # Organise
  hemibrain.matches.df <- hemibrain.matches.df %>%
    dplyr::left_join(hb.meta[,c("bodyid","cell_type")],by=c("hemibrain_match"="bodyid")) %>%
    dplyr::rename(pt_root_id = root_id, 
                  pt_supervoxel_id = supervoxel_id, 
                  pt_position = position, 
                  query_id = query,
                  match_id = hemibrain_match,
                  match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, 
                  pt_supervoxel_id, 
                  pt_position, 
                  query_id,
                  match_id,
                  match_cell_type) %>%
    dplyr::mutate(valid = 't') %>%
    rbind(hemibrain.matches.df.valid) %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)
  
  # Save
  hemibrain.matches.df <- banc_updateids(hemibrain.matches.df, 
                                         root.column = "pt_root_id", 
                                         position.column = "pt_position", 
                                         supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(hemibrain.matches.df, file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"))
}

# FANC
message("Working on BANC-fanc reviewed matches ...")
matches1 <- list.files(file.path(banc.fanc.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.fanc.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.fanc.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.fanc.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.fanc.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.fanc.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.fanc.correct.match.path,"blue_wrong"), pattern = "png")
green_wrong <- list.files(file.path(banc.fanc.correct.match.path,"green_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.fanc.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.fanc.correct.match.path,"investigate"), pattern = "png")
fanc.matches <- unique(c(matches1,
                         matches2,
                         matches3,
                         matches4))
fanc.matches.df <- data.frame()
for(file in fanc.matches){
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id_.*","",file), 
                    supervoxel_ids = gsub(".*supervoxel_id_([0-9]+).*", "\\1", file),
                    fanc_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel_id_.*|\\.png","",file))
  fanc.matches.df <- rbind(fanc.matches.df,mdf)
}
fanc.matches.df <- fanc.matches.df %>%
  tidyr::separate_rows(root_ids, sep = "_") %>%
  dplyr::rename(root_id = root_ids) %>%
  tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
  dplyr::rename(supervoxel_id = supervoxel_ids)
fanc.matches.df <- fanc.matches.df %>%
  dplyr::anti_join(fanc.matches.df.valid %>% dplyr::filter(valid=='t'), 
                   by = c("fanc_match"="match_id","root_id"="query_id"))
if(nrow(fanc.matches.df)){
  
  # Update IDs
  fanc.matches.df$query <- fanc.matches.df$root_id
  fanc.matches.df <- banc_updateids(fanc.matches.df, 
                                    root.column = "root_id",
                                    supervoxel.column = "supervoxel_id")
  
  # Assign needed data
  fanc.matches.df$position <- banc.meta$position[match(fanc.matches.df$root_id,banc.meta$root_id)]
  fanc.matches.df$fanc_match <- fc.meta$cell_id[match(fanc.matches.df$fanc_match,fc.meta$root_id)]
  
  # Organise
  fanc.matches.df <- fanc.matches.df %>%
    dplyr::left_join(fc.meta[,c("cell_id","cell_type")],by=c("fanc_match"="cell_id")) %>%
    dplyr::rename(pt_root_id = root_id, 
                  pt_supervoxel_id = supervoxel_id, 
                  pt_position = position, 
                  query_id = query,
                  match_id = fanc_match,
                  match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id, 
                  pt_supervoxel_id, 
                  pt_position, 
                  query_id,
                  match_id,
                  match_cell_type) %>%
    dplyr::mutate(valid = 't') %>%
    rbind(fanc.matches.df.valid) %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)
  
  # Save
  fanc.matches.df <- banc_updateids(fanc.matches.df, 
                                    root.column = "pt_root_id", 
                                    position.column = "pt_position", 
                                    supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(fanc.matches.df, file.path(banc.meta.save.path,"banc_fanc_reviewed_matches.csv"))
}

# maleCNS
message("Working on BANC-maleCNS reviewed matches ...")
matches1 <- list.files(file.path(banc.malecns.correct.match.path,"1_perfect"), pattern = "png")
matches2 <- list.files(file.path(banc.malecns.correct.match.path,"2_confident"), pattern = "png")
matches3 <- list.files(file.path(banc.malecns.correct.match.path,"3_good"), pattern = "png")
matches4 <- list.files(file.path(banc.malecns.correct.match.path,"4_likely"), pattern = "png")
matches5 <- list.files(file.path(banc.malecns.correct.match.path,"5_possible"), pattern = "png")
matches6 <- list.files(file.path(banc.malecns.correct.match.path,"6_no_match"), pattern = "png")
blue_wrong <- list.files(file.path(banc.malecns.correct.match.path,"blue_wrong"), pattern = "png")
green_wrong <- list.files(file.path(banc.malecns.correct.match.path,"green_wrong"), pattern = "png")
red_wrong <- list.files(file.path(banc.malecns.correct.match.path,"red_wrong"), pattern = "png")
investigate <- list.files(file.path(banc.malecns.correct.match.path,"investigate"), pattern = "png")
malecns.matches <- unique(c(matches1,
                            matches2,
                            matches3,
                            matches4,
                            matches5))
malecns.matches.df <- data.frame()
for(file in malecns.matches){
  mdf <- data.frame(root_ids = gsub(".*_root_id_|_nucleus_id_.*|_supervoxel_id.*","",file),
                    supervoxel_ids = ifelse(grepl("supervoxel_id",file),stringr::str_match(file, "supervoxel_id_([0-9]+)_")[, 2],NA),
                    malecns_match = gsub(".*bodyid_|.*hit_id_|_hit_nucleus_id_.*|_hit_supervoxel.*|_query.*|\\.png","",file))
  malecns.matches.df <- rbind(malecns.matches.df,mdf)
}
if(nrow(malecns.matches.df)){
  malecns.matches.df <- malecns.matches.df %>%
    tidyr::separate_rows(root_ids, sep = "_") %>%
    dplyr::rename(root_id = root_ids) %>%
    tidyr::separate_rows(supervoxel_ids, sep = "_") %>%
    dplyr::rename(supervoxel_id = supervoxel_ids)
  malecns.matches.df <- malecns.matches.df %>%
    dplyr::anti_join(malecns.matches.df.valid %>% dplyr::filter(valid=='t'),
                     by = c("malecns_match"="match_id", "root_id"="query_id"))
}
if(nrow(malecns.matches.df)){
 
  # Update IDs
  malecns.matches.df$query <- malecns.matches.df$root_id
  malecns.matches.df <- banc_updateids(malecns.matches.df,
                                       root.column = "root_id",
                                       supervoxel.column = "supervoxel_id")

  # Assign needed data
  malecns.matches.df$position <- banc.meta$position[match(malecns.matches.df$root_id,banc.meta$root_id)]
  malecns.matches.df$supervoxel_id <- banc.meta$supervoxel_id[match(malecns.matches.df$root_id,banc.meta$root_id)]

  # Organise
  malecns.matches.df <- malecns.matches.df %>%
    dplyr::left_join(mcns.meta[,c("malecns_09_id","cell_type")],by=c("malecns_match"="malecns_09_id")) %>%
    dplyr::rename(pt_root_id = root_id,
                  pt_supervoxel_id = supervoxel_id,
                  pt_position = position,
                  query_id = query,
                  match_id = malecns_match,
                  match_cell_type = cell_type) %>%
    dplyr::select(pt_root_id,
                  pt_supervoxel_id,
                  pt_position,
                  query_id,
                  match_id,
                  match_cell_type) %>%
    dplyr::mutate(valid = 't') %>%
    rbind(malecns.matches.df.valid) %>%
    dplyr::arrange(dplyr::desc(pt_position)) %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)

  # Save
  malecns.matches.df <- banc_updateids(malecns.matches.df,
                                       root.column = "pt_root_id",
                                       position.column = "pt_position",
                                       supervoxel.column = "pt_supervoxel_id")
  readr::write_csv(malecns.matches.df, file.path(banc.meta.save.path,"banc_malecns_reviewed_matches.csv"))
}

##########################
### Mirror NBLAST list ###
##########################

# Get new NBLAST files
message("Working on BANC mirror NBLAST matches ...")
mirror.nblast.files <- list.files(mirror.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
mirror.nblast.files <- mirror.nblast.files[grepl("supervoxel_id",mirror.nblast.files)]
done <- unique(banc.meta.mirror.nb$query_id)
mirror.nblast.files <- mirror.nblast.files[!gsub(".*_root_id_|.csv","",basename(mirror.nblast.files))%in%done]

# Read nblast mirror files and combine
if(length(mirror.nblast.files)){
  by.query.mirror <- foreach::foreach(mfile = mirror.nblast.files) %do% {
    try({
      id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
        dplyr::group_by(query) %>%
        dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
        dplyr::ungroup()
      mdf
    })
  }
  by.query.mirror <- by.query.mirror[unlist(lapply(by.query.mirror,is.data.frame))]
  
  # Combine with meta data
  banc.meta.mirror.scores.all <- do.call(rbind,by.query.mirror) %>%
    dplyr::mutate(banc_nblast_match = root_id,
                  root_id = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(nb)) %>%
    dplyr::filter(nb<=1)
  
  # Only do full update for new data
  banc.meta.mirror.scores.todo <- banc.meta.mirror.scores.all %>%
    dplyr::anti_join(banc.meta.mirror.nb, by = c("banc_nblast_match"="match_id", "query"="query_id"))
  if(nrow(banc.meta.mirror.scores.todo)){
    
    # Update ids
    queries.all <- banc.meta.mirror.scores.todo$root_id
    queries <- unique(queries.all)
    sps <- banc.meta$supervoxel_id[match(queries,banc.meta$root_id)]
    root_ids <- queries
    root_ids[is.na(sps)] <- banc_updateids(queries[is.na(sps)])
    root_ids.all <- root_ids[match(queries.all,queries)]
    banc.meta.mirror.scores.todo$root_id <- root_ids.all
    banc.meta.mirror.scores.todo$position <- banc.meta$position[match(banc.meta.mirror.scores.todo$root_id,banc.meta$root_id)]
    banc.meta.mirror.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.mirror.scores.todo$root_id,banc.meta$root_id)]
    
    # Update match ids
    targets.all <- banc.meta.mirror.scores.todo$banc_nblast_match
    targets <- unique(targets.all)
    sps <- banc.meta$supervoxel_id[match(targets,banc.meta$root_id)]
    banc_nblast_matches <- targets
    banc_nblast_matches[is.na(sps)] <- banc_updateids(targets[is.na(sps)])
    banc_nblast_matches.all <- banc_nblast_matches[match(targets.all,targets)]
    banc.meta.mirror.scores.todo$banc_nblast_match_latest <- banc_nblast_matches.all
    banc.meta.mirror.scores.todo$banc_match_position <- banc.meta$position[match(banc.meta.mirror.scores.todo$banc_nblast_match_latest,banc.meta$root_id)]
    banc.meta.mirror.scores.todo$banc_match_supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.mirror.scores.todo$banc_nblast_match_latest,banc.meta$root_id)]
    
    # Organise
    banc.meta.mirror.scores <- banc.meta.mirror.scores.todo %>%
      dplyr::rename(pt_root_id = root_id, 
                    pt_supervoxel_id = supervoxel_id, 
                    pt_position = position, 
                    match_root_id = banc_nblast_match_latest,
                    match_supervoxel_id = banc_match_supervoxel_id,
                    match_position = banc_match_position,
                    query_id = query,
                    match_id = banc_nblast_match,
                    score = nb) %>%
      dplyr::select(pt_root_id, 
                    pt_supervoxel_id, 
                    pt_position, 
                    match_root_id,
                    match_supervoxel_id,
                    match_position,
                    query_id,
                    match_id,
                    score) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      plyr::rbind.fill(banc.meta.mirror.nb) %>%
      dplyr::select(-valid) %>%
      dplyr::filter(pt_root_id!=match_root_id) %>%
      dplyr::distinct(pt_root_id, match_root_id, .keep_all = TRUE) %>%
      dplyr::left_join(mirror.matches.df[,c("pt_root_id","match_id","valid")], 
                       by=c("pt_root_id","match_id")) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
    
    # Save
    banc.meta.mirror.scores <- banc_updateids(banc.meta.mirror.scores)
    banc.meta.mirror.scores$root_626 <- banc_rootid(banc.meta.mirror.scores$pt_supervoxel_id, version="626")
    banc.meta.mirror.scores$match_root_626 <- banc_rootid(banc.meta.mirror.scores$match_supervoxel_id, version="626")
    arrow::write_feather(banc.meta.mirror.scores, file.path(banc.meta.save.path,"banc_mirror_nblast.feather"))
  }
}

########################
### FANC NBLAST list ###
########################

# Read NBLAST FANC files and combine
message("Working on BANC-FANC NBLAST matches ...")
fanc.nblast.files <- list.files(fanc.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.fanc.nb$query_id)
fanc.nblast.files <- fanc.nblast.files[!gsub(".*_|.csv","",basename(fanc.nblast.files))%in%done]
if(length(fanc.nblast.files)){
  by.query.fanc <- foreach::foreach(mfile = fanc.nblast.files) %do% {
    id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
    mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
    mdf$query <- id
    mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
      dplyr::group_by(query) %>%
      dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
      dplyr::ungroup()
    mdf
  }
  if(!length(by.query.fanc)){
    break
  }
  by.query.fanc <- by.query.fanc[unlist(lapply(by.query.fanc,is.data.frame))]
  
  # Combine with meta data
  fanc.nblast.scores.all <- do.call(plyr::rbind.fill,by.query.fanc) %>%
    dplyr::mutate(root_id = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(nb)) %>%
    dplyr::distinct(nb,fanc_id,cell_id,cell_type,query,root_id,supervoxel_id,cell_type) %>%
    dplyr::filter(nb<=1)
  
  # # update validated matches
  # banc.meta.fanc.nb <- banc.meta.fanc.nb %>%
  #   dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)  %>%
  #   dplyr::select(-valid) %>%
  #   dplyr::left_join(fanc.matches.df[,c("query_id","match_cell_type","valid")],
  #                    by=c("query_id","match_cell_type"),
  #                    relationship = "many-to-many") %>%
  #   dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
  
  # Only do full update for new data
  banc.meta.fanc.scores.todo <- fanc.nblast.scores.all %>%
    dplyr::anti_join(banc.meta.fanc.nb, by = c("cell_id"="match_id", "query"="query_id"))
  if(nrow(banc.meta.fanc.scores.todo)){
    
    # Update match ids
    queries.all <- banc.meta.fanc.scores.todo$query
    queries <- unique(queries.all)
    root_ids <- banc_updateids(queries)
    root_ids.all <- root_ids[match(queries.all,queries)]
    banc.meta.fanc.scores.todo$root_id <- root_ids.all
    banc.meta.fanc.scores.todo$position <- banc.meta$position[match(banc.meta.fanc.scores.todo$supervoxel_id,banc.meta$supervoxel_id)]
    #banc.meta.fanc.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.fanc.scores.todo$root_id,banc.meta$root_id)]
    
    # Organise
    banc.meta.fanc.nb <- banc.meta.fanc.scores.todo %>%
      dplyr::rename(pt_root_id = root_id,
                    pt_supervoxel_id = supervoxel_id,
                    pt_position = position,
                    query_id = query,
                    match_id = cell_id,
                    match_cell_type = cell_type,
                    score = nb) %>%
      dplyr::select(pt_root_id,
                    pt_supervoxel_id,
                    pt_position,
                    query_id,
                    match_id,
                    match_cell_type,
                    score) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      plyr::rbind.fill(banc.meta.fanc.nb) %>%
      dplyr::select(-valid) %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)  %>%
      dplyr::left_join(fanc.matches.df[,c("pt_root_id","match_id","valid")],
                       by=c("pt_root_id","match_id")) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
    
    # Save
    fanc.nblast.scores <- banc_updateids(banc.meta.fanc.nb, root.column = 'pt_root_id', supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
    fanc.nblast.scores$root_626 <- banc_rootid(fanc.nblast.scores$pt_supervoxel_id, version="626")
    arrow::write_feather(fanc.nblast.scores, file.path(banc.meta.save.path,"banc_fanc_1116_nblast.feather"))
  }
}

########################
### FAFB NBLAST list ###
########################

# Read nblast fafb files and combine
message("Working on BANC-FAFB NBLAST matches ...")
fafb.nblast.files <- list.files(fafb.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.fafb.nb$query_id)
fafb.nblast.files <- fafb.nblast.files[!gsub(".*_|.csv","",basename(fafb.nblast.files))%in%done]
batch_size <- 250
total_files <- length(fafb.nblast.files)
num_batches <- ceiling(total_files / batch_size)
fafb.nblast.files <- sample(fafb.nblast.files)
for (batch_idx in 1:num_batches) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, total_files)
  current_batch <- fafb.nblast.files[start_idx:end_idx]
  message(sprintf("Processing batch %d/%d (files %d to %d)", 
                  batch_idx, num_batches, start_idx, end_idx))
  if(length(current_batch)){
    by.query.fafb <- NULL
    by.query.fafb <- foreach::foreach(mfile = current_batch, 
                                      .errorhandling = 'pass') %do% {
                                        id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
                                        mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
                                        if(length(mdf)){
                                          numeric_columns <- sapply(mdf, is.numeric)
                                          mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
                                          mdf$query <- id
                                          mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
                                            dplyr::group_by(query) %>%
                                            dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
                                            dplyr::ungroup()
                                          mdf 
                                        }else{
                                          NULL
                                        }
                                      }
    if(!length(by.query.fafb)){
      break
    }
    by.query.fafb <- by.query.fafb[unlist(lapply(by.query.fafb,is.data.frame))]
    
    # Combine with meta data
    fafb.nblast.scores.all <- do.call(plyr::rbind.fill,by.query.fafb) %>%
      dplyr::mutate(root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::distinct(nb,root_783,nucleus_id,cell_type,query,root_id) %>%
      dplyr::filter(nb<=1)
    
    # update validated matches
    banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)  %>%
      dplyr::select(-valid) %>%
      dplyr::left_join(fafb.matches.df[,c("query_id","match_cell_type","valid")], 
                       by=c("query_id","match_cell_type"),
                       relationship = "many-to-many") %>%
      dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
    
    # Only do full update for new data
    banc.meta.fafb.scores.todo <- fafb.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.fafb.nb, by = c("root_783"="match_id", "query"="query_id"))
    if(nrow(banc.meta.fafb.scores.todo)){
      
      # Update match ids
      queries.all <- banc.meta.fafb.scores.todo$root_id
      queries <- unique(queries.all)
      root_ids <- banc_updateids(queries)
      root_ids.all <- root_ids[match(queries.all,queries)]
      banc.meta.fafb.scores.todo$root_id <- root_ids.all
      banc.meta.fafb.scores.todo$position <- banc.meta$position[match(banc.meta.fafb.scores.todo$root_id,banc.meta$root_id)]
      banc.meta.fafb.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.fafb.scores.todo$root_id,banc.meta$root_id)]
      
      # Organise
      banc.meta.fafb.nb <- banc.meta.fafb.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, 
                      pt_supervoxel_id = supervoxel_id, 
                      pt_position = position, 
                      query_id = query,
                      match_id = root_783,
                      match_cell_type = cell_type,
                      score = nb) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::select(pt_root_id, 
                      pt_supervoxel_id, 
                      pt_position, 
                      query_id,
                      match_id,
                      match_cell_type,
                      score) %>%
        plyr::rbind.fill(banc.meta.fafb.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE)  %>%
        dplyr::left_join(fafb.matches.df[,c("pt_root_id","match_id","valid")], 
                         by=c("pt_root_id","match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid),'f',valid)) %>%
        dplyr::filter(!is.na(pt_supervoxel_id), 
                      !is.na(score), 
                      score<=1)
      
      # Save
      fafb.nblast.scores <- banc_updateids(banc.meta.fafb.nb, root.column = 'pt_root_id', supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      fafb.nblast.scores$root_626 <- banc_rootid(fafb.nblast.scores$pt_supervoxel_id, version="626")
      arrow::write_feather(fafb.nblast.scores, file.path(banc.meta.save.path,"banc_fafb_783_nblast.feather"))
    }
  }
}

########################
### MANC NBLAST list ###
########################

# Read nblast manc files and combine
message("Working on BANC-MANC NBLAST matches ...")
manc.nblast.files <- list.files(manc.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- banc.meta.manc.nb$query_id
manc.nblast.files <- manc.nblast.files[!gsub(".*_|.csv","",basename(manc.nblast.files))%in%done]
batch_size <- 250 
total_files <- length(manc.nblast.files)
num_batches <- ceiling(total_files / batch_size)
for (batch_idx in 1:num_batches) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, total_files)
  current_batch <- manc.nblast.files[start_idx:end_idx]
  message(sprintf("Processing batch %d/%d (files %d to %d)", 
                  batch_idx, num_batches, start_idx, end_idx))
  if(length(current_batch)){
    by.query.manc <- foreach::foreach(mfile = current_batch) %do% {
      id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
        dplyr::group_by(query) %>%
        dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
        dplyr::ungroup()
      mdf
    }
    if(!length(by.query.manc)){
      break
    }
    by.query.manc <- by.query.manc[unlist(lapply(by.query.manc,is.data.frame))]
    
    # Combine with meta data
    manc.nblast.scores.all <- do.call(plyr::rbind.fill,by.query.manc)
    manc.nblast.scores.all$root_location <- NULL
    manc.nblast.scores.all <- manc.nblast.scores.all %>%
      dplyr::mutate(root_id = query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(nb)) %>%
      dplyr::filter(nb<=1)
    
    # update validated matches
    banc.meta.manc.nb <- banc.meta.manc.nb %>%
      dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)  %>%
      dplyr::select(-valid) %>%
      dplyr::left_join(manc.matches.df[,c("query_id","match_cell_type","valid")], 
                       by=c("query_id","match_cell_type"),
                       relationship = "many-to-many") %>%
      dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
    
    # Only do full update for new data
    banc.meta.manc.scores.todo <- manc.nblast.scores.all %>%
      dplyr::anti_join(banc.meta.manc.nb, by = c("bodyid"="match_id", "query"="query_id"))
    if(nrow(banc.meta.manc.scores.todo)){
      
      # Update ids
      queries.all <- banc.meta.manc.scores.todo$root_id
      queries <- unique(queries.all)
      root_ids <- banc_updateids(queries)
      root_ids.all <- root_ids[match(queries.all,queries)]
      banc.meta.manc.scores.todo$root_id <- root_ids.all
      banc.meta.manc.scores.todo$position <- banc.meta$position[match(banc.meta.manc.scores.todo$root_id,banc.meta$root_id)]
      banc.meta.manc.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.manc.scores.todo$root_id,banc.meta$root_id)]
      
      # Organise
      banc.meta.manc.nb <- banc.meta.manc.scores.todo %>%
        dplyr::rename(pt_root_id = root_id, 
                      pt_supervoxel_id = supervoxel_id, 
                      pt_position = position, 
                      query_id = query,
                      match_id = bodyid,
                      match_cell_type = cell_type,
                      score = nb) %>%
        dplyr::select(pt_root_id, 
                      pt_supervoxel_id, 
                      pt_position, 
                      query_id,
                      match_id,
                      match_cell_type,
                      score) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        plyr::rbind.fill(banc.meta.manc.nb) %>%
        dplyr::select(-valid) %>%
        dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
        dplyr::left_join(manc.matches.df[,c("pt_root_id","match_id","valid")], 
                         by=c("pt_root_id","match_id"),
                         relationship = "many-to-many") %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::mutate(valid = ifelse(is.na(valid),'f',valid)) %>%
        dplyr::filter(!is.na(pt_supervoxel_id),!is.na(pt_root_id))
      
      # Save
      manc.nblast.scores <- banc_updateids(banc.meta.manc.nb, root.column = 'pt_root_id', supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
      manc.nblast.scores$root_626 <- banc_rootid(manc.nblast.scores$pt_supervoxel_id, version="626")
      arrow::write_feather(manc.nblast.scores, file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.feather")) 
    }
  }
}

#############################
### Hemibrain NBLAST list ###
#############################

# Read nblast hemibrain files and combine
message("Working on BANC-hemibrain NBLAST matches ...")
hemibrain.nblast.files <- list.files(hemibrain.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- unique(banc.meta.hemibrain.nb$query_id)
hemibrain.nblast.files <- hemibrain.nblast.files[!gsub(".*_|.csv","",basename(hemibrain.nblast.files))%in%done]
if(length(hemibrain.nblast.files)){
  by.query.hemibrain <- foreach::foreach(mfile = hemibrain.nblast.files) %do% {
    id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
    svid <- gsub(".*supervoxel_id_|_root_id_.*|\\.csv","",basename(mfile))
    mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
    mdf$query <- id
    mdf$query_supervoxel_id <- svid
    mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
      dplyr::group_by(query) %>%
      dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
      dplyr::ungroup()
    mdf
  }
  if(!length(by.query.hemibrain)){
    break
  }
  by.query.hemibrain <- by.query.hemibrain[unlist(lapply(by.query.hemibrain,is.data.frame))]
  
  # Combine with meta data
  hemibrain.nblast.scores.all <- do.call(plyr::rbind.fill,by.query.hemibrain) %>%
    dplyr::mutate(root_id = query,
                  supervoxel_id=query_supervoxel_id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(nb))  %>%
    dplyr::filter(nb<=1)
  
  # # update validated matches
  banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
    dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE)  %>%
    dplyr::select(-valid) %>%
    dplyr::left_join(hemibrain.matches.df[,c("query_id","match_cell_type","valid")],
                     by=c("query_id","match_cell_type"),
                     relationship = "many-to-many") %>%
    dplyr::mutate(valid = ifelse(is.na(valid),'f',valid))
  
  # Only do full update for new data
  banc.meta.hemibrain.scores.todo <- hemibrain.nblast.scores.all %>%
    dplyr::anti_join(banc.meta.hemibrain.nb, by = c("bodyid"="match_id", "query"="query_id"))
  if(nrow(banc.meta.hemibrain.scores.todo)){
    
    # Update ids
    queries.all <- banc.meta.hemibrain.scores.todo$supervoxel_id
    queries <- unique(queries.all)
    root_ids <- banc_rootid(queries)
    root_ids.all <- root_ids[match(queries.all,queries)]
    banc.meta.hemibrain.scores.todo$root_id <- root_ids.all
    banc.meta.hemibrain.scores.todo$position <- banc.meta$position[match(banc.meta.hemibrain.scores.todo$root_id,banc.meta$root_id)]
    
    # Organise
    banc.meta.hemibrain.nb <- banc.meta.hemibrain.scores.todo %>%
      dplyr::rename(pt_root_id = root_id, 
                    pt_supervoxel_id = supervoxel_id, 
                    pt_position = position, 
                    query_id = query,
                    match_id = bodyid,
                    match_cell_type = cell_type,
                    score = nb) %>%
      dplyr::select(pt_root_id, 
                    pt_supervoxel_id, 
                    pt_position, 
                    query_id,
                    match_id,
                    match_cell_type,
                    score) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      plyr::rbind.fill(banc.meta.hemibrain.nb) %>%
      dplyr::select(-valid) %>%
      dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
      dplyr::left_join(hemibrain.matches.df[,c("pt_root_id","match_id","valid")],
                       by=c("pt_root_id","match_id"),
                       relationship = "many-to-many") %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(valid = ifelse(is.na(valid),'f',valid)) %>%
      dplyr::filter(!is.na(pt_supervoxel_id))
    
    # Save
    hemibrain.nblast.scores <- banc_updateids(banc.meta.hemibrain.nb, root.column = 'pt_root_id', supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
    hemibrain.nblast.scores$root_626 <- banc_rootid(hemibrain.nblast.scores$pt_supervoxel_id, version="626")
    arrow::write_feather(hemibrain.nblast.scores, file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.feather"))
  }
}

############################
### maleCNS NBLAST list ###
############################

# Read nblast maleCNS files and combine
message("Working on BANC-maleCNS NBLAST matches ...")
malecns.nblast.files <- list.files(malecns.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
done <- banc.meta.malecns.nb$query_id
malecns.nblast.files <- malecns.nblast.files[!gsub(".*_|.csv","",basename(malecns.nblast.files))%in%done]
batch_size <- 250
total_files <- length(malecns.nblast.files)
if(total_files){
  num_batches <- ceiling(total_files / batch_size)
  for (batch_idx in 1:num_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, total_files)
    current_batch <- malecns.nblast.files[start_idx:end_idx]
    message(sprintf("Processing batch %d/%d (files %d to %d)",
                    batch_idx, num_batches, start_idx, end_idx))
    if(length(current_batch)){
      by.query.malecns <- foreach::foreach(mfile = current_batch) %do% {
        id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
        sv_id <- str_match(basename(mfile), "supervoxel_id_([0-9]+)_root_id")[, 2]
        mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
        numeric_columns <- sapply(mdf, is.numeric)
        mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
        mdf$query <- id
        mdf$supervoxel_id <- sv_id
        mdf <- dplyr::arrange(mdf, dplyr::desc(nb)) %>%
          dplyr::group_by(query) %>%
          dplyr::filter(nb >= 0.3 | nb >= nth(sort(nb, decreasing = TRUE), 5)) %>%
          dplyr::ungroup()
        mdf
      }
      if(!length(by.query.malecns)){
        break
      }
      by.query.malecns <- by.query.malecns[unlist(lapply(by.query.malecns,is.data.frame))]
      
      # Combine with meta data
      malecns.nblast.scores.all <- do.call(plyr::rbind.fill,by.query.malecns)
      malecns.nblast.scores.all$root_location <- NULL
      malecns.nblast.scores.all <- malecns.nblast.scores.all %>%
        dplyr::mutate(root_id = query) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(nb)) %>%
        dplyr::filter(nb<=1)
      
      # update validated matches
      if(nrow(malecns.matches.df)){
        banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
          dplyr::distinct(pt_root_id, match_id, .keep_all = TRUE) %>%
          dplyr::select(-valid) %>%
          dplyr::left_join(malecns.matches.df[,c("query_id","match_cell_type","valid")],
                           by=c("query_id","match_cell_type"),
                           relationship = "many-to-many") %>%
          dplyr::mutate(valid = ifelse(is.na(valid),'f',valid)) 
      }
      
      # Only do full update for new data
      banc.meta.malecns.scores.todo <- malecns.nblast.scores.all %>%
        dplyr::anti_join(banc.meta.malecns.nb, by = c("bodyid"="match_id", "query"="query_id"))
      if(nrow(banc.meta.malecns.scores.todo)){
        
        # Update ids
        queries.all <- banc.meta.malecns.scores.todo$root_id
        queries <- unique(queries.all)
        #root_ids <- banc_updateids(queries)
        root_ids.all <- root_ids[match(queries.all,queries)]
        banc.meta.malecns.scores.todo$root_id <- root_ids.all
        banc.meta.malecns.scores.todo$position <- banc.meta$position[match(banc.meta.malecns.scores.todo$root_id,banc.meta$root_id)]
        banc.meta.malecns.scores.todo$supervoxel_id <- banc.meta$supervoxel_id[match(banc.meta.malecns.scores.todo$root_id,banc.meta$root_id)]
        banc.meta.malecns.scores.todo <- subset(banc.meta.malecns.scores.todo, !is.na(position))
        
        # Organise
        banc.meta.malecns.nb <- banc.meta.malecns.scores.todo %>%
          dplyr::rename(pt_root_id = root_id,
                        pt_supervoxel_id = supervoxel_id,
                        pt_position = position,
                        query_id = query,
                        match_id = bodyid,
                        match_cell_type = cell_type,
                        score = nb) %>%
          dplyr::select(pt_root_id,
                        pt_supervoxel_id,
                        pt_position,
                        query_id,
                        match_id,
                        match_cell_type,
                        score) %>%
          dplyr::arrange(dplyr::desc(score)) %>%
          plyr::rbind.fill(banc.meta.malecns.nb) %>%
          dplyr::select(-valid) %>%
          dplyr::distinct(pt_root_id, query_id, match_id, .keep_all = TRUE) %>%
          dplyr::left_join(malecns.matches.df[,c("pt_root_id","match_id","valid")],
                           by=c("pt_root_id","match_id"),
                           relationship = "many-to-many") %>%
          dplyr::arrange(dplyr::desc(score)) %>%
          dplyr::mutate(valid = ifelse(is.na(valid),'f',valid)) %>%
          dplyr::filter(!is.na(pt_supervoxel_id),!is.na(pt_root_id))
        
        # Save
        malecns.nblast.scores <- banc_updateids(banc.meta.malecns.nb, root.column = 'pt_root_id', supervoxel.column = 'pt_supervoxel_id', position.column = 'pt_position')
        malecns.nblast.scores$root_626 <- banc_rootid(malecns.nblast.scores$pt_supervoxel_id, version="626")
        arrow::write_feather(malecns.nblast.scores, file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.feather"))
      }
    }
  }
}

# Stop parallel backend
stopImplicitCluster()

# Push to google storage
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_fafb_783_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_fanc_1116_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_mirror_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_malecns_v0.9_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_manc_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_fafb_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_hemibrain_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_mirror_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_fanc_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/nblast", banc.meta.save.path, "banc_malecns_reviewed_matches.csv"))

#################
### Send v626 ###
#################

# Push to google storage
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_fafb_783_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_fanc_1116_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_mirror_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_malecns_v0.9_nblast.feather"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_manc_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_fafb_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_hemibrain_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_mirror_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_fanc_reviewed_matches.csv"))
system(sprintf("gsutil cp %s/%s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/nblast", banc.meta.save.path, "banc_malecns_reviewed_matches.csv"))

##########################
### Plot our inventory ###
##########################

# Step 1: Process the data
bc <- bc %>% dplyr::distinct(root_id, .keep_all = TRUE)
bc.cts <- bc %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON",status)) %>%
  dplyr::filter(!is.na(root_id), 
                !is.na(supervoxel_id), 
                !is.na(position), 
                root_id!="", 
                root_id!="0") %>%
  dplyr::arrange(
    fafb_cell_type,
    manc_cell_type,
    hemibrain_cell_type,
    cell_type
  ) %>%
  dplyr::filter(!duplicated(root_id)) %>%
  dplyr::mutate(
    fafb_cell_type = dplyr::case_when(
      grepl("auto",fafb_cell_type) ~ NA,
      is.na(fafb_match) ~ NA,
      TRUE ~ fafb_cell_type
    ),
    manc_cell_type = dplyr::case_when(
      grepl("auto",manc_cell_type) ~ NA,
      is.na(manc_match) ~ NA,
      TRUE ~ manc_cell_type
    ),
    hemibrain_cell_type = dplyr::case_when(
      grepl("auto",hemibrain_cell_type) ~ NA,
      is.na(hemibrain_match) ~ NA,
      TRUE ~ hemibrain_cell_type
    )
  ) %>%
  dplyr::filter(!is.na(fafb_match)|!is.na(manc_match)|!is.na(hemibrain_match)) %>%
  dplyr::select(root_id, supervoxel_id, position, fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::arrange(fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  dplyr::distinct(root_id, fafb_cell_type, manc_cell_type, hemibrain_cell_type, .keep_all = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, fafb_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, manc_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(root_id, hemibrain_cell_type) %>%
  dplyr::mutate(duplicate_flag = dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!duplicate_flag) %>%
  dplyr::rename(pt_root_id = root_id, 
                pt_supervoxel_id = supervoxel_id, 
                pt_position = position) %>%
  dplyr::select(-duplicate_flag) %>%
  reshape2::melt(id = c("pt_root_id",  "pt_supervoxel_id", "pt_position"),
                 value.name = "tag",
                 variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA","unknown","Uknown","fragment","Fragment","None","none","no_match"),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity",
                user_id = 355,
                valid = "t") %>%
  dplyr::mutate(tag2 = gsub("_"," ",tag2)) %>%
  dplyr::filter(!grepl("auto",tag), !grepl("no_match",tag), 
                !grepl("auto",tag2), !grepl("no_match",tag2), 
                !is.na(tag), !is.na(tag2)) %>%
  dplyr::distinct(pt_root_id, pt_supervoxel_id, pt_position, tag2, tag, user_id, valid)
processed_data <- bc.cts %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(both = any(grepl("fafb",tag2))&any(grepl("manc",tag2)),
                both2 = any(grepl("fafb",tag2))&any(grepl("hemibrain",tag2))) %>%
  dplyr::mutate(tag2 = dplyr::case_when(
    both ~ "manual fafb+manc cell type",
    both2 ~ "manual fafb+hemibrain cell type",
    TRUE ~ tag2
  )) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(bc[,c("root_id","region")], by = c("pt_root_id"='root_id')) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(tag2 = ifelse(is.na(tag2),"unmatched",tag2)) %>%
  dplyr::group_by(tag2, region) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(tag2, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

# Step 2: Define custom colors
tag2_colors <- scales::hue_pal()(length(unique(processed_data$tag2)))
names(tag2_colors) <- unique(processed_data$tag2)

# Calculate the positions for the labels
processed_data <- processed_data %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Step 3: Create the plot
pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = tag2)) +
  facet_wrap(vars(region)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 2) +
  geom_segment(aes(x = 4, 
                   y = position, 
                   xend = 4.05, 
                   yend = position), 
               color = "black", 
               size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tag2_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "match status", title = "distribution of entries by match status")

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/cell_types/bc_verified_match_status_distribution.png", width = 12, height = 12, dpi = 300)

# Step 1: Process the data
# Sort out data to share, NBLAST high:
nb.thresh <- 0.5
fw.nblast.scores <- banc.meta.fafb.nb %>% 
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(fafb_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
mc.nblast.scores <- banc.meta.manc.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(manc_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
hb.nblast.scores <- banc.meta.hemibrain.nb %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::filter(score>=nb.thresh) %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::rename(hemibrain_cell_type = match_cell_type) %>%
  dplyr::mutate(valid == "f")
nblast.scores <- fw.nblast.scores %>%
  dplyr::full_join(mc.nblast.scores, by = c("pt_root_id","pt_supervoxel_id","pt_position")) %>%
  dplyr::full_join(hb.nblast.scores, by = c("pt_root_id","pt_supervoxel_id","pt_position")) %>%
  dplyr::filter(!is.na(pt_root_id), !is.na(pt_supervoxel_id), !is.na(pt_position), pt_root_id!="", pt_root_id!="0") %>%
  dplyr::select(pt_root_id, pt_supervoxel_id, pt_position, 
                fafb_cell_type, manc_cell_type, hemibrain_cell_type) %>%
  reshape2::melt(id = c("pt_root_id",  "pt_supervoxel_id", "pt_position"),
                 value.name = "tag",
                 variable.name = "tag2") %>%
  dplyr::filter(!tag %in% c("NA","unknown","Uknown","fragment","Fragment","None","none","no_match"),
                !is.na(tag)) %>%
  dplyr::mutate(tag2 = "neuron identity",
                tag = paste0(tag,"?"),
                user_id = 355,
                valid = "t") %>% 
  dplyr::filter(!is.na(tag), !is.na(tag2))
nblast.scores <- bancr::banc_updateids(nblast.scores, root.column = "pt_root_id",
                                       supervoxel.column = "pt_supervoxel_id",
                                       position.column = "pt_position")
nblast.scores <- nblast.scores %>% 
  dplyr::group_by(pt_root_id, tag2) %>%
  dplyr::filter(!duplicated(pt_position), 
                !duplicated(pt_root_id),
                pt_root_id!="0") %>%
  dplyr::ungroup()
processed_data <- nblast.scores %>%
  dplyr::group_by(pt_root_id) %>%
  dplyr::mutate(both = any(grepl("fafb",tag2))&any(grepl("manc",tag2)),
                both2 = any(grepl("fafb",tag2))&any(grepl("hemibrain",tag2))) %>%
  dplyr::mutate(tag2 = dplyr::case_when(
    both ~ "manual fafb+manc cell type",
    both2 ~ "manual fafb+hemibrain cell type",
    TRUE ~ tag2
  )) %>%
  dplyr::full_join(bc[,c("root_id","region")], by = c("pt_root_id"='root_id')) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(tag2 = ifelse(is.na(tag2),"unmatched",tag2)) %>%
  dplyr::group_by(tag2, region) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(
    proportion = count / sum(count),
    label = paste0(tag2, "\n(", count, ", ", scales::percent(proportion, accuracy = 0.1), ")")
  )

# Step 2: Define custom colors
# You can replace this with your own color palette if desired
tag2_colors <- scales::hue_pal()(length(unique(processed_data$tag2)))
names(tag2_colors) <- unique(processed_data$tag2)

# Calculate the positions for the labels
processed_data <- processed_data %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    position = (ymax + ymin) / 2,
    position_label = 1.05
  )

# Step 3: Create the plot
pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = tag2)) +
  facet_wrap(vars(region)) +
  geom_rect() +
  ggrepel::geom_text_repel(aes(label = label, x = 4, y = position), 
                           nudge_x = 0.5,
                           size = 2) +
  geom_segment(aes(x = 4, 
                   y = position, 
                   xend = 4.05, 
                   yend = position), 
               color = "black", 
               size = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tag2_colors) +
  xlim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "match status", title = "distribution of entries by match status")

# Display the plot
print(pie_chart)

# Save
ggsave(plot = pie_chart,
       filename = "inst/images/cell_types/bc_high_nblast_match_status_distribution.png", width = 12, height = 12, dpi = 300)


