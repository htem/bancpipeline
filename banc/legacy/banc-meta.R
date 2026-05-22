##############################################
### Compile meta data from various sources ###
##############################################
source("banc/banc-startup.R")
message("##### Building BANC meta data #####")
version <- "elastix_tpsreg_0714"
version <- "elastix_tpsreg_240721"

# Register cores
cl <- setup_parallel()

# Decide on thresholds
fafb.nblast.threshold <- 0.3
manc.nblast.threshold <- 0.3
malecns.nblast.threshold <- 0.3
hemibrain.nblast.threshold <- 0.3
fanc.nblast.threshold <- 0.3

# Get NBLAST versions
fafb.nblast.folder <- file.path(banc.nblast.fafb.save.path,"results",version)
mirror.nblast.folder <- file.path(banc.nblast.mirror.save.path,"results")
manc.nblast.folder <- file.path(banc.nblast.manc.save.path,"results",version)
malecns.nblast.folder <- file.path(banc.nblast.malecns.save.path,"results","navis_tpsreg_250206")

# # Get meta data
# fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
# fw.ids <- unique(fw.meta$root_783)
# mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
# mc.ids <- unique(mc.meta$bodyid)
franken.meta <- franken_meta()
fw.meta <- franken.meta %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::rename(root_783=fafb_id)
fw.ids <- unique(fw.meta$root_783)
mc.meta <- franken.meta %>%
  dplyr::filter(!is.na(manc_id)) %>%
  dplyr::rename(bodyid=manc_id)
mc.ids <- unique(mc.meta$manc_id)
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
hb.ids <- unique(hb.meta$bodyid)
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(fc.meta$cell_id)
mcns.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"malecns_09_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
mcns.ids <- unique(mcns.meta$malecns_09_id)

# Direct us to the BANC dataset
choose_banc()

# Read IDs
banc.ids <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                            col_types = banc.col.types, 
                            show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id),!is.na(position))

# Read root nodes
soma.positions <- readr::read_csv(file=file.path(banc.save.path,"banc_root_positions.csv"), 
                                  col_types = banc.col.types, 
                                  show_col_types = FALSE)
bc.regions <- soma.positions %>%
  dplyr::select(root_id, region) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
banc.ids <- dplyr::left_join(banc.ids, bc.regions, by = "root_id")
if("region.y"%in%colnames(banc.ids)){
  banc.ids <- banc.ids %>%
    dplyr::mutate(region = dplyr::case_when(
      is.na(region.x) ~ region.y,
      is.na(region.y) ~ region.x,
      TRUE ~ region.x
    )) 
}

# Update with L2 skeleton stats
banc.metrics.files <- list.files(banc.metrics.save.path, pattern = "\\.csv", full.names = TRUE)

# Read metrics files and combine
proofread <- banc_backbone_proofread()
by.query.metrics <- foreach::foreach(mfile = banc.metrics.files) %do% {
  mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
  numeric_columns <- sapply(mdf, is.numeric)
  mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 3)
  mdf[,c("l2_n_trees","l2_endpoints","l2_branchpoints","l2_segments","l2_root")] <- NULL
  mdf
}
by.query.metrics <- by.query.metrics[unlist(lapply(by.query.metrics,is.data.frame))]
banc.meta <- do.call(plyr::rbind.fill, by.query.metrics)
banc.meta <- banc.meta %>%
  dplyr::mutate(l2_cable_length_um = round(l2_cable_length/1000,2)) %>%
  dplyr::distinct(root_id, l2_nodes, l2_cable_length_um) %>%
  dplyr::arrange(dplyr::desc(l2_nodes))
banc.meta <- banc.meta %>%
  dplyr::left_join(banc.ids[,c("root_id","supervoxel_id")] 
                   %>% dplyr::distinct(root_id,.keep_all=TRUE), 
                   by = "root_id") %>%
  dplyr::left_join(banc.ids %>% 
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::mutate(root_id=as.character(root_id),
                                   supervoxel_id = as.character(supervoxel_id)) %>%
                     dplyr::select(root_id, supervoxel_id),
                   by=c("root_id")) %>%
  dplyr::mutate(supervoxel_id = dplyr::case_when(
    is.na(supervoxel_id.x) ~ supervoxel_id.y,
    is.na(supervoxel_id.y) ~ supervoxel_id.x,
    TRUE ~ supervoxel_id.x
  )) %>%
  dplyr::select(-supervoxel_id.x, -supervoxel_id.y) %>%
  dplyr::left_join(banc.ids %>% dplyr::distinct(root_id,.keep_all=TRUE), 
                   by = "root_id") %>%
  dplyr::mutate(supervoxel_id = dplyr::case_when(
    is.na(supervoxel_id.x) ~ supervoxel_id.y,
    is.na(supervoxel_id.y) ~ supervoxel_id.x,
    TRUE ~ supervoxel_id.x
  )) %>%
  dplyr::select(-supervoxel_id.x, -supervoxel_id.y) %>%
  dplyr::filter(!is.na(supervoxel_id))
banc.meta <- banc_updateids(banc.meta, 
                            root.column = "root_id", 
                            supervoxel.column = "supervoxel_id",
                            position.column = "position")

# Read nblast mirror files and combine
banc.meta.mirror.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_mirror_nblast.csv"), 
                                                        col_types = hemibrainr:::sql_col_types))
banc.meta.mirror.nb <- banc.meta.mirror.nb %>%
  dplyr::rename(banc_nblast_match=match_root_id,
                banc_nblast_match_supervoxel_id=match_supervoxel_id,
                root_id=pt_root_id,
                banc_nblast = score) %>%
  dplyr::arrange(dplyr::desc(banc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, 
                banc_nblast_match, 
                banc_nblast,
                banc_nblast_match_supervoxel_id) %>%
  # dplyr::mutate(banc_match = ifelse(valid=='t'),banc_nblast_match,'') %>%
  # dplyr::mutate(banc_match_supervoxel_id = ifelse(valid=='t'),banc_nblast_match_supervoxel_id,'') %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
banc.meta$banc_nblast_match <- NULL
banc.meta$banc_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, banc.meta.mirror.nb, by = c("root_id"))

# Update BANC NBLAST match column
banc_nblast_match_missing <- is.na(banc.meta$banc_nblast_match)|banc.meta$banc_nblast_match==""|banc.meta$banc_nblast_match=="0"
banc.meta$banc_nblast_match[!banc_nblast_match_missing] <- banc_rootid(banc.meta$banc_nblast_match_supervoxel_id[!banc_nblast_match_missing])
banc.meta$banc_nblast_match[banc.meta$banc_nblast_match=="0"] <- ""

# # Update BANC match column 
# banc_match_missing <- is.na(banc.meta$banc_match)|banc.meta$banc_match==""|banc.meta$banc_match=="0"
# banc.meta$banc_match[!banc_match_missing] <- banc_rootid(banc.meta$banc_match_supervoxel_id[!banc_match_missing])
# banc.meta$banc_match[banc.meta$banc_match=="0"] <- ""

# Read nblast fafb files and combine
banc.meta.fafb.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_fafb_783_nblast.csv"), 
                                                      col_types = hemibrainr:::sql_col_types))
banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
  dplyr::rename(fafb_nblast_match=match_id, 
                root_id=pt_root_id,
                fafb_nblast=score) %>%
  dplyr::arrange(dplyr::desc(fafb_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, 
                fafb_nblast_match, 
                fafb_nblast) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(fw.meta[,c("root_783","hemilineage", "nerve", "cell_function",
                              "top_nt", "flow", "super_class", "cell_class", 
                              "cell_type", "hemibrain_type")] %>%
                     dplyr::distinct(root_783, .keep_all = TRUE), 
                   by = c("fafb_nblast_match"="root_783")) %>%
  dplyr::mutate(
    fafb_ito_lee_hemilineage = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(hemilineage) ~ paste0("auto:",hemilineage),
      TRUE ~ NA
    ),
    fafb_nerve = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(nerve) ~ paste0("auto:",nerve),
      TRUE ~ NA
    ),
    fafb_top_nt = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(top_nt) ~ paste0("auto:",top_nt),
      TRUE ~ NA
    ),
    fafb_flow = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(flow) ~ paste0("auto:",flow),
      TRUE ~ NA
    ),
    fafb_cell_function= dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(cell_function) ~ paste0("auto:",cell_function),
      TRUE ~ NA
    ),
    fafb_super_class = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(super_class) ~ paste0("auto:",super_class),
      TRUE ~ NA
    ),
    fafb_cell_class = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(cell_class) ~ paste0("auto:",cell_class),
      TRUE ~ NA
    ),
    fafb_auto_cell_type = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    ),
    # fafb_hemibrain_match = dplyr::case_when(
    #   fafb_nblast>=fafb.nblast.threshold&!is.na(hemibrain_match) ~ paste0("auto:",hemibrain_match),
    #   TRUE ~ NA
    # ),
    fafb_hemibrain_type = dplyr::case_when(
      fafb_nblast>=fafb.nblast.threshold&!is.na(hemibrain_type) ~ paste0("auto:",hemibrain_type),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-c("hemilineage", 
                   "nerve",
                   "top_nt", 
                   "flow", 
                   "cell_function",
                   "super_class", 
                   "cell_class", 
                   "cell_type", 
                   "hemibrain_type")) %>%
  dplyr::arrange(dplyr::desc(fafb_nblast))

# Combine with previous data
banc.meta$fafb_nblast_match <- NULL
banc.meta$fafb_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, 
                              banc.meta.fafb.nb, 
                              by = c("root_id"))

# Read nblast fanc files and combine
banc.meta.fanc.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_fanc_1116_nblast.csv"), 
                                                      col_types = hemibrainr:::sql_col_types))
banc.meta.fanc.nb <- banc.meta.fanc.nb %>%
  dplyr::rename(fanc_nblast_match=match_id, 
                root_id=pt_root_id,
                fanc_nblast=score) %>%
  dplyr::arrange(dplyr::desc(fanc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, 
                fanc_nblast_match, 
                fanc_nblast) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(fc.meta[,c("cell_id", "super_class", "cell_class", 
                              "cell_sub_class", "nerve", "hemilineage",
                              "cell_type")] , 
                   by = c("fanc_nblast_match"="cell_id")) %>%
  dplyr::mutate(
    fanc_ito_lee_hemilineage = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(hemilineage) ~ paste0("auto:",hemilineage),
      TRUE ~ NA
    ),
    fanc_nerve = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(nerve) ~ paste0("auto:",nerve),
      TRUE ~ NA
    ),
    fanc_super_class = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(super_class) ~ paste0("auto:",super_class),
      TRUE ~ NA
    ),
    fanc_cell_class = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(cell_class) ~ paste0("auto:",cell_class),
      TRUE ~ NA
    ),
    fanc_cell_sub_class = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(cell_sub_class) ~ paste0("auto:",cell_sub_class),
      TRUE ~ NA
    ),
    fanc_auto_cell_type = dplyr::case_when(
      fanc_nblast>=fanc.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-c("hemilineage", 
                   "nerve",
                   "super_class", 
                   "cell_class", 
                   "cell_type", 
                   "cell_sub_class")) %>%
  dplyr::arrange(dplyr::desc(fanc_nblast))

# Combine with previous data
banc.meta$fanc_nblast_match <- NULL
banc.meta$fanc_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, 
                              banc.meta.fanc.nb, 
                              by = c("root_id"))

# Read nblast manc files and combine
banc.meta.manc.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.csv"), 
                                                      col_types = hemibrainr:::sql_col_types))
banc.meta.manc.nb <- banc.meta.manc.nb %>%
  dplyr::rename(manc_nblast_match=match_id, 
                root_id=pt_root_id,
                manc_nblast=score) %>%
  dplyr::arrange(dplyr::desc(manc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, 
                manc_nblast_match, 
                manc_nblast) %>%
  dplyr::left_join(mc.meta[,c("bodyid",
                              "hemilineage",
                              "top_nt",
                              "nerve", 
                              "cell_function",
                              "cell_class", 
                              "cell_type")] %>%
                     dplyr::distinct(bodyid, .keep_all = TRUE), 
                   by = c("manc_nblast_match"="bodyid")) %>%
  dplyr::mutate(
    manc_hemilineage = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(hemilineage) ~ paste0("auto:",hemilineage),
      TRUE ~ NA
    ),
    manc_top_nt = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(top_nt) ~ paste0("auto:",top_nt),
      TRUE ~ NA
    ),
    manc_cell_class = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(cell_class) ~ paste0("auto:",cell_class),
      TRUE ~ NA
    ),
    manc_auto_cell_type = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    ),
    cell_function_cell_function= dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(cell_function) ~ paste0("auto:",cell_function),
      TRUE ~ NA
    ),
    manc_nerve = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(nerve) ~ paste0("auto:",nerve),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-c("hemilineage",
                   "top_nt", 
                   "cell_class", 
                   "cell_type", 
                   "cell_function",
                   "nerve")) %>%
  dplyr::arrange(dplyr::desc(manc_nblast))

# Combine with previous data
banc.meta$manc_nblast_match <- NULL
banc.meta$manc_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta,
                              banc.meta.manc.nb,
                              by = c("root_id"))

# Read nblast maleCNS files and combine
banc.meta.malecns.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.csv"),
                                                         col_types = hemibrainr:::sql_col_types))
banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
  dplyr::rename(malecns_nblast_match=match_id,
                root_id=pt_root_id,
                malecns_nblast=score) %>%
  dplyr::arrange(dplyr::desc(malecns_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id,
                malecns_nblast_match,
                malecns_nblast) %>%
  dplyr::left_join(mcns.meta[,c("malecns_09_id",
                                "hemilineage",
                                "neurotransmitter_predicted",
                                "nerve",
                                "cell_function",
                                "cell_class",
                                "cell_type")] %>%
                     dplyr::distinct(malecns_09_id, .keep_all = TRUE),
                   by = c("malecns_nblast_match"="malecns_09_id")) %>%
  dplyr::mutate(
    malecns_hemilineage = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(hemilineage) ~ paste0("auto:",hemilineage),
      TRUE ~ NA
    ),
    malecns_top_nt = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(neurotransmitter_predicted) ~ paste0("auto:",neurotransmitter_predicted),
      TRUE ~ NA
    ),
    malecns_cell_class = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(cell_class) ~ paste0("auto:",cell_class),
      TRUE ~ NA
    ),
    malecns_auto_cell_type = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    ),
    malecns_cell_function = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(cell_function) ~ paste0("auto:",cell_function),
      TRUE ~ NA
    ),
    malecns_nerve = dplyr::case_when(
      malecns_nblast>=malecns.nblast.threshold&!is.na(nerve) ~ paste0("auto:",nerve),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-c("hemilineage",
                   "neurotransmitter_predicted",
                   "cell_class",
                   "cell_type",
                   "cell_function",
                   "nerve")) %>%
  dplyr::arrange(dplyr::desc(malecns_nblast))

# Combine with previous data
banc.meta$malecns_nblast_match <- NULL
banc.meta$malecns_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta,
                              banc.meta.malecns.nb,
                              by = c("root_id"))

# Read nblast hemibrain files and combine
banc.meta.hemibrain.nb <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.csv"), 
                                                           col_types = hemibrainr:::sql_col_types))
banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
  dplyr::rename(hemibrain_nblast_match=match_id, 
                root_id=pt_root_id,
                hemibrain_nblast=score) %>%
  dplyr::arrange(dplyr::desc(hemibrain_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, 
                hemibrain_nblast_match, 
                hemibrain_nblast) %>%
  dplyr::mutate(hemibrain_nblast_match = gsub("m","",hemibrain_nblast_match)) %>%
  dplyr::left_join(hb.meta[,c("bodyid",
                              "top_nt",
                              "cell_type")], 
                   by = c("hemibrain_nblast_match"="bodyid")) %>%
  dplyr::mutate(
    hemibrain_auto_cell_type = dplyr::case_when(
      hemibrain_nblast>=hemibrain.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    ),
    hemibrain_top_nt= dplyr::case_when(
      hemibrain_nblast>=hemibrain.nblast.threshold&!is.na(top_nt) ~ paste0("auto:",top_nt),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-c("cell_type","top_nt")) %>%
  dplyr::arrange(dplyr::desc(hemibrain_nblast))

# Combine with previous data
banc.meta$hemibrain_nblast_match <- NULL
banc.meta$hemibrain_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, banc.meta.hemibrain.nb, by = c("root_id"))

# Save
banc.meta <- banc.meta %>%
  dplyr::mutate(
    fafb_nblast = dplyr::case_when(
      is.na(fafb_nblast) ~ NA,
      fafb_nblast>1 ~ NA,
      TRUE ~ fafb_nblast
    ),
    fanc_nblast = dplyr::case_when(
      is.na(fanc_nblast) ~ NA,
      fanc_nblast>1 ~ NA,
      TRUE ~ fanc_nblast
    ),
    manc_nblast = dplyr::case_when(
      is.na(manc_nblast) ~ NA,
      manc_nblast>1 ~ NA,
      TRUE ~ manc_nblast
    ),
    hemibrain_nblast = dplyr::case_when(
      is.na(hemibrain_nblast) ~ NA,
      hemibrain_nblast>1 ~ NA,
      TRUE ~ hemibrain_nblast
    ),
    malecns_nblast = dplyr::case_when(
      is.na(malecns_nblast) ~ NA,
      malecns_nblast>1 ~ NA,
      TRUE ~ malecns_nblast
    ),
  )
readr::write_csv(banc.meta, file.path(banc.meta.save.path,"banc_meta.csv"))

# Announce
message("##### BANCpipeline: banc metrics updated #####")
message(sprintf("##### we have metrics for : %s neurons", nrow(banc.meta)))

# Stop cores
stop_parallel(cl)

##########################
#### MANC NBLAST PLOT ####
##########################

# Get our calls
manc.nblast.threshold.borderline <- 0.2
manc.nblast.threshold <- 0.3
manc.nblast.threshold.very_good <- 0.5

# Remove 'auto' prefix and empty entries
bc_cleaned <- banc.meta %>%
  dplyr::filter(grepl("neck|vnc",region)) %>%
  dplyr::mutate(
    manc_nblast = as.numeric(manc_nblast),
    assignments_match = dplyr::case_when(
      is.na(manc_nblast) ~ "no_nblast",
      manc_nblast < manc.nblast.threshold.borderline ~ "bad",
      manc_nblast >= manc.nblast.threshold.very_good ~ "confident",
      manc_nblast >= manc.nblast.threshold ~ "good",
      manc_nblast >= manc.nblast.threshold.borderline ~ "borderline",
      TRUE ~ "bad"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = c("no_nblast","bad","borderline","good","confident")))

# Define colours
col.values = c("no_nblast" = "lightgrey",
               "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
               "borderline"= hemibrainr:::hemibrain_bright_colors[["orange"]],
               "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
               "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col.values) +
  theme_minimal() +
  labs(
    title = "BANC-MANC top NBLAST scores per BANC neuron",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of manc_top NBLAST scores
plot2 <- ggplot(bc_cleaned, aes(x = manc_nblast, 
                                y= (..count..)/nrow(bc_cleaned))) +
  geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-MANC top NBLAST scores",
    x = "BANC-MANC NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$manc_nblast, na.rm = TRUE), by = 0.1)) +
  ggplot2::geom_vline(xintercept = 0, color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = manc.nblast.threshold.borderline, color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = manc.nblast.threshold, color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = manc.nblast.threshold.very_good, color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = manc.nblast.threshold.borderline, y = Inf, label = sprintf("borderline: >%s",manc.nblast.threshold.borderline), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
  ggplot2::annotate("text", x = manc.nblast.threshold, y = Inf, label = sprintf("good: >%s",manc.nblast.threshold), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
  ggplot2::annotate("text", x = manc.nblast.threshold.very_good, y = Inf, label = sprintf("confident: >%s",manc.nblast.threshold.very_good), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/nblast/manc_nblast_match_banc.png", width = 18, height = 8, dpi = 300)

#############################
#### maleCNS NBLAST PLOT ####
#############################

# Get our calls
malecns.nblast.threshold.borderline <- 0.2
malecns.nblast.threshold <- 0.3
malecns.nblast.threshold.very_good <- 0.5

# Remove 'auto' prefix and empty entries
bc_cleaned <- banc.meta %>%
  dplyr::mutate(
    malecns_nblast = as.numeric(malecns_nblast),
    assignments_match = dplyr::case_when(
      is.na(malecns_nblast) ~ "no_nblast",
      malecns_nblast < malecns.nblast.threshold.borderline ~ "bad",
      malecns_nblast >= malecns.nblast.threshold.very_good ~ "confident",
      malecns_nblast >= malecns.nblast.threshold ~ "good",
      malecns_nblast >= malecns.nblast.threshold.borderline ~ "borderline",
      TRUE ~ "bad"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match,
                                         levels = c("no_nblast","bad","borderline","good","confident")))

# Define colours
col.values = c("no_nblast" = "lightgrey",
               "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
               "borderline"= hemibrainr:::hemibrain_bright_colors[["orange"]],
               "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
               "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col.values) +
  theme_minimal() +
  labs(
    title = "BANC-maleCNS top NBLAST scores per BANC neuron",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))),
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of malecns_top NBLAST scores
plot2 <- ggplot(bc_cleaned, aes(x = malecns_nblast,
                                y= (..count..)/nrow(bc_cleaned))) +
  geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-maleCNS top NBLAST scores",
    x = "BANC-maleCNS NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$malecns_nblast, na.rm = TRUE), by = 0.1)) +
  ggplot2::geom_vline(xintercept = 0, color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = malecns.nblast.threshold.borderline, color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = malecns.nblast.threshold, color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = malecns.nblast.threshold.very_good, color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = malecns.nblast.threshold.borderline, y = Inf, label = sprintf("borderline: >%s",malecns.nblast.threshold.borderline),
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
  ggplot2::annotate("text", x = malecns.nblast.threshold, y = Inf, label = sprintf("good: >%s",malecns.nblast.threshold),
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
  ggplot2::annotate("text", x = malecns.nblast.threshold.very_good, y = Inf, label = sprintf("confident: >%s",malecns.nblast.threshold.very_good),
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1,
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/nblast/malecns_nblast_match_banc.png", width = 18, height = 8, dpi = 300)

##########################
#### FAFB NBLAST PLOT ####
##########################

# Get our calls
fafb.nblast.threshold.borderline <- 0.2
fafb.nblast.threshold <- 0.3
fafb.nblast.threshold.very_good <- 0.5

# Remove 'auto' prefix and empty entries
bc_cleaned <- banc.meta %>%
  dplyr::filter(grepl("brain|neck|optic",region)) %>%
  dplyr::mutate(
    fafb_nblast = as.numeric(fafb_nblast),
    assignments_match = dplyr::case_when(
      is.na(fafb_nblast) ~ "no_nblast",
      fafb_nblast < fafb.nblast.threshold.borderline ~ "bad",
      fafb_nblast >= fafb.nblast.threshold.very_good ~ "confident",
      fafb_nblast >= fafb.nblast.threshold ~ "good",
      fafb_nblast >= fafb.nblast.threshold.borderline ~ "borderline",
      TRUE ~ "bad"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = c("no_nblast","bad","borderline","good","confident")))

# Define colours
col.values = c("no_nblast" = "lightgrey",
               "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
               "borderline"= hemibrainr:::hemibrain_bright_colors[["orange"]],
               "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
               "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col.values) +
  theme_minimal() +
  labs(
    title = "BANC-FAFB top NBLAST scores per BANC neuron",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of fafb_top NBLAST scores
plot2 <- ggplot(bc_cleaned, aes(x = fafb_nblast)) +
  geom_density(aes(y = after_stat(density))) +
  geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-FAFB top NBLAST scores",
    x = "BANC-FAFB NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$fafb_nblast, na.rm = TRUE), by = 0.1)) +
  ggplot2::geom_vline(xintercept = 0, color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = fafb.nblast.threshold.borderline, color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = fafb.nblast.threshold, color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = fafb.nblast.threshold.very_good, color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = fafb.nblast.threshold.borderline, y = Inf, label = sprintf("borderline: >%s",fafb.nblast.threshold.borderline), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
  ggplot2::annotate("text", x = fafb.nblast.threshold, y = Inf, label = sprintf("good: >%s",fafb.nblast.threshold), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
  ggplot2::annotate("text", x = fafb.nblast.threshold.very_good, y = Inf, label = sprintf("confident: >%s",fafb.nblast.threshold.very_good), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/nblast/fafb_nblast_match_banc.png", width = 18, height = 8, dpi = 300)


###############################
#### Hemibrain NBLAST PLOT ####
###############################

# Get our calls
hemibrain.nblast.threshold.borderline <- 0.2
hemibrain.nblast.threshold <- 0.3
hemibrain.nblast.threshold.very_good <- 0.5

# Remove 'auto' prefix and empty entries
bc_cleaned <- banc.meta %>%
  dplyr::filter(grepl("brain|neck ",region)) %>%
  dplyr::mutate(
    hemibrain_nblast = as.numeric(hemibrain_nblast),
    assignments_match = dplyr::case_when(
      is.na(hemibrain_nblast) ~ "no_nblast",
      hemibrain_nblast < hemibrain.nblast.threshold.borderline ~ "bad",
      hemibrain_nblast >= hemibrain.nblast.threshold.very_good ~ "confident",
      hemibrain_nblast >= hemibrain.nblast.threshold ~ "good",
      hemibrain_nblast >= hemibrain.nblast.threshold.borderline ~ "borderline",
      TRUE ~ "bad"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = c("no_nblast","bad","borderline","good","confident")))

# Define colours
col.values = c("no_nblast" = "lightgrey",
               "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
               "borderline"= hemibrainr:::hemibrain_bright_colors[["orange"]],
               "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
               "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col.values) +
  theme_minimal() +
  labs(
    title = "BANC-Hemibrain top NBLAST scores per BANC neuron",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of hemibrain_top NBLAST scores
plot2 <- ggplot(bc_cleaned, aes(x = hemibrain_nblast, 
                                y= (..count..)/nrow(bc_cleaned))) +
  geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-Hemibrain top NBLAST scores",
    x = "BANC-Hemibrain NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$hemibrain_nblast, na.rm = TRUE), by = 0.1)) +
  ggplot2::geom_vline(xintercept = 0, color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = hemibrain.nblast.threshold.borderline, color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = hemibrain.nblast.threshold, color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = hemibrain.nblast.threshold.very_good, color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = hemibrain.nblast.threshold.borderline, y = Inf, label = sprintf("borderline: >%s",hemibrain.nblast.threshold.borderline), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
  ggplot2::annotate("text", x = hemibrain.nblast.threshold, y = Inf, label = sprintf("good: >%s",hemibrain.nblast.threshold), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
  ggplot2::annotate("text", x = hemibrain.nblast.threshold.very_good, y = Inf, label = sprintf("confident: >%s",hemibrain.nblast.threshold.very_good), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# save
ggsave(plot = combined_plot,
       filename = "inst/images/nblast/hemibrain_nblast_match_banc.png", width = 18, height = 8, dpi = 300)


############################
#### Mirror NBLAST PLOT ####
############################

# Get our calls
banc.nblast.threshold.borderline <- 0.2
banc.nblast.threshold <- 0.3
banc.nblast.threshold.very_good <- 0.5

# Remove 'auto' prefix and empty entries
bc_cleaned <- banc.meta %>%
  dplyr::mutate(
    banc_nblast = as.numeric(banc_nblast),
    assignments_match = dplyr::case_when(
      is.na(banc_nblast) ~ "no_nblast",
      banc_nblast < banc.nblast.threshold.borderline ~ "bad",
      banc_nblast >= banc.nblast.threshold.very_good ~ "confident",
      banc_nblast >= banc.nblast.threshold ~ "good",
      banc_nblast >= banc.nblast.threshold.borderline ~ "borderline",
      TRUE ~ "bad"
    ),
  ) %>%
  dplyr::mutate(assignments_match=factor(assignments_match, 
                                         levels = c("no_nblast","bad","borderline","good","confident")))

# Define colours
col.values = c("no_nblast" = "lightgrey",
               "bad" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
               "borderline"= hemibrainr:::hemibrain_bright_colors[["orange"]],
               "good" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
               "confident" = hemibrainr:::hemibrain_bright_colors[["green"]])

# Calculate the percentage of matching and non-matching assignments
match_summary <- bc_cleaned %>%
  dplyr::group_by(assignments_match) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

# Plot 1: Stacked bar plot of matching and non-matching assignments
plot1 <- ggplot(match_summary, aes(x = "", y = percentage, fill = assignments_match)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col.values) +
  theme_minimal() +
  labs(
    title = "BANC-Mirror top NBLAST scores per BANC neuron",
    x = "",
    y = "percentage",
    fill = "assignments match"
  ) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", percentage, scales::comma(count))), 
            position = position_stack(vjust = 0.5), size = 3.5)

# Plot 2: Normalized density plot of banc_top NBLAST scores
plot2 <- ggplot(bc_cleaned, aes(x = banc_nblast, 
                                y= (..count..)/nrow(bc_cleaned))) +
  geom_density(alpha = 0.5, fill = hemibrainr:::hemibrain_bright_colors[["marine"]]) +
  theme_minimal() +
  labs(
    title = "normalised density of BANC-Mirror top NBLAST scores",
    x = "BANC-Mirror NBLAST score",
    y = "density",
    fill = "assignments match"
  ) +
  scale_x_continuous(breaks = seq(0, max(bc_cleaned$banc_nblast, na.rm = TRUE), by = 0.1)) +
  ggplot2::geom_vline(xintercept = 0, color = hemibrainr:::hemibrain_bright_colors[["cerise"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = banc.nblast.threshold.borderline, color = hemibrainr:::hemibrain_bright_colors[["orange"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = banc.nblast.threshold, color = hemibrainr:::hemibrain_bright_colors[["yellow"]], linetype = "dashed", size = 1) +
  ggplot2::geom_vline(xintercept = banc.nblast.threshold.very_good, color = hemibrainr:::hemibrain_bright_colors[["green"]], linetype = "dashed", size = 1) +
  ggplot2::annotate("text", x = banc.nblast.threshold.borderline, y = Inf, label = sprintf("borderline: >%s",banc.nblast.threshold.borderline), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["orange"]]) +
  ggplot2::annotate("text", x = banc.nblast.threshold, y = Inf, label = sprintf("good: >%s",banc.nblast.threshold), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["yellow"]]) +
  ggplot2::annotate("text", x = banc.nblast.threshold.very_good, y = Inf, label = sprintf("confident: >%s",banc.nblast.threshold.very_good), 
                    vjust = 2.2, hjust = 1.1, color = hemibrainr:::hemibrain_bright_colors[["green"]])

# Combine the plots into a single image
combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, 
                                         top = "")

# Display the combined plot
print(combined_plot)

# Save
ggsave(plot = combined_plot,
       filename = "inst/images/nblast/banc_nblast_match_banc.png", width = 18, height = 8, dpi = 300)

