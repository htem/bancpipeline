#' banc-meta — Compile `banc_meta.csv` from L2 metrics + NBLAST CSVs + IDs.
#'
#' Builds `banc_meta.csv` by combining: (1) L2 skeleton metrics from
#' per-neuron CSVs, (2) `banc_ids.csv` + root positions, (3) compiled
#' NBLAST CSVs to derive `auto:*` cell types above score thresholds.
#' Plots have moved to `banc/nblast/banc-nblast-plot.R`.
#'
#' @section Reads:
#'   - `<banc.metrics.save.path>/*.csv` (per-neuron L2 metrics)
#'   - `<banc.meta.save.path>/banc_ids.csv`
#'   - compiled NBLAST CSVs under `<banc.nblast.save.path>/`
#'
#' @section Writes:
#'   - `<banc.meta.save.path>/banc_meta.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_update.sh`

##############################################
### Compile meta data from various sources ###
###
### Builds banc_meta.csv by combining:
###   1. L2 skeleton metrics from per-neuron CSVs
###   2. banc_ids.csv + root positions
###   3. Compiled NBLAST CSVs → derive auto:*
###      cell types above score thresholds
###
### Plots have been moved to banc-nblast-plot.R
##############################################
source("banc/banc-startup.R")

local({

message("##### Building BANC meta data #####")
t_start <- Sys.time()
version <- banc.nblast.version

# Helper: ensure columns exist in a data frame, filling missing with NA
ensure_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  df
}

# Register cores
cl <- setup_parallel()
on.exit(stop_parallel(cl), add = TRUE)

# Decide on thresholds
fafb.nblast.threshold <- 0.2
manc.nblast.threshold <- 0.2
malecns.nblast.threshold <- 0.2
hemibrain.nblast.threshold <- 0.2
fanc.nblast.threshold <- 0.2

# Get NBLAST versions
fafb.nblast.folder <- file.path(banc.nblast.fafb.save.path,"results",version)
mirror.nblast.folder <- file.path(banc.nblast.mirror.save.path,"results")
manc.nblast.folder <- file.path(banc.nblast.manc.save.path,"results",version)
malecns.nblast.folder <- file.path(banc.nblast.malecns.save.path,"results",banc.nblast.malecns.version)

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

# Direct us to the BANC dataset
choose_banc()

# Read IDs
banc.ids <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                            col_types = banc.col.types, 
                            show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id),!is.na(position))

# Read regions
bc.regions <- arrow::read_feather(file.path(banc.save.path, "banc_regions.feather")) %>%
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

# Update with L2 skeleton stats — only read files for current root_ids
banc.metrics.files <- file.path(banc.metrics.save.path, paste0(banc.ids$root_id, ".csv"))
banc.metrics.files <- banc.metrics.files[file.exists(banc.metrics.files)]

# Read metrics files and combine
message(sprintf("Reading %d L2 metrics files (of %d in SeaTable)...",
                length(banc.metrics.files), nrow(banc.ids)))
proofread <- banc_backbone_proofread()
n_mfiles <- length(banc.metrics.files)
by.query.metrics <- foreach::foreach(mfile = banc.metrics.files, i = seq_along(banc.metrics.files)) %do% {
  if (i == 1 || i %% 500 == 0 || i == n_mfiles)
    message(sprintf("  Metrics: %d/%d (%.0f%%)", i, n_mfiles, 100 * i / n_mfiles))
  tryCatch({
    mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 3)
    mdf[,c("l2_n_trees","l2_endpoints","l2_branchpoints","l2_segments","l2_root")] <- NULL
    mdf
  }, error = function(e) {
    message(sprintf("  Error reading %s: %s", basename(mfile), e$message))
    NULL
  })
}
by.query.metrics <- by.query.metrics[unlist(lapply(by.query.metrics,is.data.frame))]
banc.meta <- dplyr::bind_rows(by.query.metrics)
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
message("Reading mirror NBLAST data...")
banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_mirror_nblast.feather"))
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
message("Reading FAFB NBLAST data...")
banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_fafb_783_nblast.feather"))
banc.meta.fafb.nb <- banc.meta.fafb.nb %>%
  dplyr::rename(fafb_nblast_match=match_id,
                root_id=pt_root_id,
                fafb_nblast=score) %>%
  # Prefer fresh entries (NBLAST was computed on the CURRENT root_id) over
  # stale entries that fast_updateids_df migrated forward after a post-NBLAST
  # proofread edit. Within each tier, take the highest score.
  dplyr::mutate(.is_fresh = root_id == query_root_id) %>%
  dplyr::arrange(dplyr::desc(.is_fresh), dplyr::desc(fafb_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(-.is_fresh) %>%
  dplyr::select(root_id, 
                fafb_nblast_match, 
                fafb_nblast) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(fw.meta[, intersect(c("root_783","hemilineage", "nerve", "cell_function",
                              "top_nt", "flow", "super_class", "cell_class",
                              "cell_type", "hemibrain_type"), colnames(fw.meta))] %>%
                     dplyr::distinct(root_783, .keep_all = TRUE),
                   by = c("fafb_nblast_match"="root_783")) %>%
  ensure_cols(c("hemilineage", "nerve", "cell_function", "top_nt", "flow",
                "super_class", "cell_class", "cell_type", "hemibrain_type")) %>%
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
  dplyr::select(-dplyr::any_of(c("hemilineage",
                   "nerve",
                   "top_nt",
                   "flow",
                   "cell_function",
                   "super_class",
                   "cell_class",
                   "cell_type",
                   "hemibrain_type"))) %>%
  dplyr::arrange(dplyr::desc(fafb_nblast))

# Combine with previous data
banc.meta$fafb_nblast_match <- NULL
banc.meta$fafb_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, 
                              banc.meta.fafb.nb, 
                              by = c("root_id"))

# Read nblast fanc files and combine
message("Reading FANC NBLAST data...")
banc.meta.fanc.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_fanc_1116_nblast.feather"))
banc.meta.fanc.nb <- banc.meta.fanc.nb %>%
  dplyr::rename(fanc_nblast_match=match_id,
                root_id=pt_root_id,
                fanc_nblast=score) %>%
  # Prefer fresh entries (NBLAST was computed on the CURRENT root_id) over
  # stale entries that fast_updateids_df migrated forward after a post-NBLAST
  # proofread edit. Within each tier, take the highest score.
  dplyr::mutate(.is_fresh = root_id == query_root_id) %>%
  dplyr::arrange(dplyr::desc(.is_fresh), dplyr::desc(fanc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(-.is_fresh) %>%
  dplyr::select(root_id, 
                fanc_nblast_match, 
                fanc_nblast) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(fc.meta[, intersect(c("cell_id", "super_class", "cell_class",
                              "cell_sub_class", "nerve", "hemilineage",
                              "cell_type"), colnames(fc.meta))] ,
                   by = c("fanc_nblast_match"="cell_id")) %>%
  ensure_cols(c("hemilineage", "nerve", "super_class", "cell_class", "cell_sub_class", "cell_type")) %>%
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
  dplyr::select(-dplyr::any_of(c("hemilineage",
                   "nerve",
                   "super_class",
                   "cell_class",
                   "cell_type",
                   "cell_sub_class"))) %>%
  dplyr::arrange(dplyr::desc(fanc_nblast))

# Combine with previous data
banc.meta$fanc_nblast_match <- NULL
banc.meta$fanc_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta, 
                              banc.meta.fanc.nb, 
                              by = c("root_id"))

# Read nblast manc files and combine
message("Reading MANC NBLAST data...")
banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather"))
banc.meta.manc.nb <- banc.meta.manc.nb %>%
  dplyr::rename(manc_nblast_match=match_id,
                root_id=pt_root_id,
                manc_nblast=score) %>%
  # Prefer fresh entries (NBLAST was computed on the CURRENT root_id) over
  # stale entries that fast_updateids_df migrated forward after a post-NBLAST
  # proofread edit. Within each tier, take the highest score.
  dplyr::mutate(.is_fresh = root_id == query_root_id) %>%
  dplyr::arrange(dplyr::desc(.is_fresh), dplyr::desc(manc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(-.is_fresh) %>%
  dplyr::select(root_id, 
                manc_nblast_match, 
                manc_nblast) %>%
  dplyr::left_join(mc.meta[, intersect(c("bodyid",
                              "hemilineage",
                              "top_nt",
                              "nerve",
                              "cell_function",
                              "cell_class",
                              "cell_type"), colnames(mc.meta))] %>%
                     dplyr::distinct(bodyid, .keep_all = TRUE),
                   by = c("manc_nblast_match"="bodyid")) %>%
  ensure_cols(c("hemilineage", "top_nt", "nerve", "cell_function", "cell_class", "cell_type")) %>%
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
    manc_cell_function = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(cell_function) ~ paste0("auto:",cell_function),
      TRUE ~ NA
    ),
    manc_nerve = dplyr::case_when(
      manc_nblast>=manc.nblast.threshold&!is.na(nerve) ~ paste0("auto:",nerve),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-dplyr::any_of(c("hemilineage",
                   "top_nt",
                   "cell_class",
                   "cell_type",
                   "cell_function",
                   "nerve"))) %>%
  dplyr::arrange(dplyr::desc(manc_nblast))

# Combine with previous data
banc.meta$manc_nblast_match <- NULL
banc.meta$manc_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta,
                              banc.meta.manc.nb,
                              by = c("root_id"))

message("Reading maleCNS NBLAST data...")
# Read nblast maleCNS files and combine
banc.meta.malecns.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_malecns_v0.9_nblast.feather"))
banc.meta.malecns.nb <- banc.meta.malecns.nb %>%
  dplyr::rename(malecns_nblast_match=match_id,
                root_id=pt_root_id,
                malecns_nblast=score) %>%
  # Prefer fresh entries (NBLAST was computed on the CURRENT root_id) over
  # stale entries that fast_updateids_df migrated forward after a post-NBLAST
  # proofread edit. Within each tier, take the highest score.
  dplyr::mutate(.is_fresh = root_id == query_root_id) %>%
  dplyr::arrange(dplyr::desc(.is_fresh), dplyr::desc(malecns_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(-.is_fresh) %>%
  dplyr::select(root_id,
                malecns_nblast_match,
                malecns_nblast) %>%
  dplyr::left_join(mcns.meta[, intersect(c("malecns_09_id",
                                "hemilineage",
                                "neurotransmitter_predicted",
                                "nerve",
                                "cell_function",
                                "cell_class",
                                "cell_type"), colnames(mcns.meta))] %>%
                     dplyr::distinct(malecns_09_id, .keep_all = TRUE),
                   by = c("malecns_nblast_match"="malecns_09_id")) %>%
  ensure_cols(c("hemilineage", "neurotransmitter_predicted", "nerve", "cell_function", "cell_class", "cell_type")) %>%
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
  dplyr::select(-dplyr::any_of(c("hemilineage",
                   "neurotransmitter_predicted",
                   "cell_class",
                   "cell_type",
                   "cell_function",
                   "nerve"))) %>%
  dplyr::arrange(dplyr::desc(malecns_nblast))

# Combine with previous data
banc.meta$malecns_nblast_match <- NULL
banc.meta$malecns_nblast <- NULL
banc.meta <- dplyr::left_join(banc.meta,
                              banc.meta.malecns.nb,
                              by = c("root_id"))

message("Reading hemibrain NBLAST data...")
# Read nblast hemibrain files and combine
banc.meta.hemibrain.nb <- arrow::read_feather(file.path(banc.meta.save.path, "banc_hemibrain_v1.2.1_nblast.feather"))
banc.meta.hemibrain.nb <- banc.meta.hemibrain.nb %>%
  dplyr::rename(hemibrain_nblast_match=match_id,
                root_id=pt_root_id,
                hemibrain_nblast=score) %>%
  # Prefer fresh entries (NBLAST was computed on the CURRENT root_id) over
  # stale entries that fast_updateids_df migrated forward after a post-NBLAST
  # proofread edit. Within each tier, take the highest score.
  dplyr::mutate(.is_fresh = root_id == query_root_id) %>%
  dplyr::arrange(dplyr::desc(.is_fresh), dplyr::desc(hemibrain_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(-.is_fresh) %>%
  dplyr::select(root_id, 
                hemibrain_nblast_match, 
                hemibrain_nblast) %>%
  dplyr::mutate(hemibrain_nblast_match = gsub("^m|m$","",hemibrain_nblast_match)) %>%
  dplyr::left_join(hb.meta[, intersect(c("bodyid", "hemilineage", "top_nt", "cell_type"),
                                        colnames(hb.meta))] %>%
                     dplyr::distinct(bodyid, .keep_all = TRUE),
                   by = c("hemibrain_nblast_match"="bodyid")) %>%
  ensure_cols(c("hemilineage", "top_nt", "cell_type")) %>%
  dplyr::mutate(
    hemibrain_hemilineage = dplyr::case_when(
      hemibrain_nblast>=hemibrain.nblast.threshold&!is.na(hemilineage) ~ paste0("auto:",hemilineage),
      TRUE ~ NA
    ),
    hemibrain_auto_cell_type = dplyr::case_when(
      hemibrain_nblast>=hemibrain.nblast.threshold&!is.na(cell_type) ~ paste0("auto:",cell_type),
      TRUE ~ NA
    ),
    hemibrain_top_nt= dplyr::case_when(
      hemibrain_nblast>=hemibrain.nblast.threshold&!is.na(top_nt) ~ paste0("auto:",top_nt),
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-dplyr::any_of(c("hemilineage","cell_type","top_nt"))) %>%
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
message(sprintf("##### BANCpipeline: banc metrics updated [%s] #####",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
message(sprintf("##### we have metrics for : %s neurons", nrow(banc.meta)))

})
