####################################
### COLLATE BANC BRAIN INFLUENCE ###
####################################
source("banc/banc-startup.R")
overwrite <- TRUE
n_cores <- max(1, parallel::detectCores() - 1)

####################
### WRANGLE DATA ###
####################

# Get current meta
banc.meta <- banctable_query() %>%
  dplyr::arrange(dplyr::desc(output_connections)) %>%
  dplyr::filter(!is.na(root_626)) %>%
  dplyr::distinct(root_id,.keep_all = TRUE)

# Get orig meta
meta <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "banc_610_meta.feather"))

# Master folder
banc.influence.save.path <- "/n/data1/hms/neurobio/wilson/banc/influence/"
influence.master <- file.path(banc.influence.save.path,"banc_610/")
seed.groups <- list.files(influence.master, recursive = FALSE, full.names = TRUE)
seed.sets <- list.files(seed.groups, recursive = FALSE, full.names = TRUE)

# Output parquet file
parquet_file <- file.path(banc.influence.save.path, "influence_banc_626.parquet")
if(overwrite && file.exists(parquet_file)){
  file.remove(parquet_file)
}

# CSV column types for individual influence result files
inf.col.types <- readr::cols(
  id = readr::col_character(),
  is_seed = readr::col_logical(),
  .default = readr::col_number())

# Collect all results across seed groups
all_results <- list()

# Iterate over seed groups
for(sg in seed.groups){
  message("Working on: ", basename(sg))
  seed.level <- tolower(gsub(".*from_|\\.csv*|_influence.*","",basename(sg)))
  do.not.process <- grepl("_13$|_13_",sg)
  if(do.not.process){
    next
  }

  # Overwrite?
  if(!overwrite && file.exists(parquet_file)){
    influence.db <- arrow::read_parquet(parquet_file, col_select = c("level", "seed")) %>%
      dplyr::filter(level == seed.level)
    if(nrow(influence.db)){
      seeds.done <- unique(influence.db$seed)
    }else{
      seeds.done <- NULL
    }
    rm('influence.db')
  }else{
    seeds.done <- NULL
  }

  # Get seed files to process
  seed.files <- list.files(sg, recursive = FALSE, full.names = TRUE, pattern = "\\.csv")
  seed.seeds <- gsub(".*from_|\\.csv*|_influence.*","",basename(seed.files))
  seed.files <- seed.files[!seed.seeds%in%seeds.done]
  if(!length(seed.files)){
    next
  }

  message(sprintf("  Processing %d files in parallel (%d cores)", length(seed.files), n_cores))

  # Process CSVs in parallel
  results <- parallel::mclapply(seed.files, function(csv){
    seed.seed <- gsub(".*from_|\\.csv*|_influence.*","",basename(csv))
    data <- readr::read_csv(csv,
                            progress = FALSE,
                            col_types = inf.col.types)
    data <- data[,-1]
    colnames(data) <- c("id","is_seed","influence")
    n.seeds <- sum(data$is_seed, na.rm = TRUE)
    id.seeds <- subset(data, is_seed)$id
    syns.seeds <- subset(banc.meta, root_id %in% id.seeds |
                           supervoxel_id %in% subset(meta, root_id %in% id.seeds)$supervoxel_id)$output_connections
    n.syns.seed <- sum(syns.seeds, na.rm = TRUE)
    if(n.syns.seed==0){
      n.syns.seed <- 1
    }

    # Skip if no influence
    if(all(data$influence=="0")){
      warning("no values for: ", csv)
      return(NULL)
    }

    # supervoxel id
    data <- dplyr::left_join(data,
                             meta %>%
                               dplyr::distinct(root_id, .keep_all = TRUE) %>%
                               dplyr::select(root_id, supervoxel_id),
                             by = c("id"="root_id")) %>%
      dplyr::left_join(banc.meta %>%
                         dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
                         dplyr::select(root_626, supervoxel_id),
                       by = c("supervoxel_id")) %>%
      dplyr::ungroup()

    # Apply seed label
    data$seed <- seed.seed
    data$level <- seed.level

    # Perform scaling operations
    data$influence_original <- as.numeric(data$influence)
    data$influence_norm_original <- as.numeric(data$influence)/n.seeds
    data$influence_syn_norm_original <- as.numeric(data$influence)/n.syns.seed
    data <- data %>%
      calculate_influence_norms()
    data
  }, mc.cores = n_cores)

  # Filter out NULLs and collect
  results <- results[!sapply(results, is.null)]
  if(length(results)){
    all_results <- c(all_results, results)
  }
  message(sprintf("  Done: %d results collected", length(results)))
}

# Combine and write as parquet
if(length(all_results)){
  influence_data <- dplyr::bind_rows(all_results)
  message(sprintf("Writing %d rows to parquet", nrow(influence_data)))
  write_connectome_data(influence_data, parquet_file, format = "parquet")
}


# Get current meta
banc.meta <- banctable_query("SELECT root_id, root_626, supervoxel_id, position from banc_meta") %>%
  dplyr::filter(!is.na(root_626)) %>%
  dplyr::distinct(root_id,.keep_all = TRUE)

# Get orig meta
meta <- arrow::read_feather(
  file.path(banc.connectivity.save.path, "banc_610_meta.feather"))

# Get ID mapping
mappings <- banc.meta %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
  dplyr::select(root_626, supervoxel_id) %>%
  dplyr::left_join(meta %>%
                     dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
                     dplyr::select(root_610 = root_id, supervoxel_id),
                   by = "supervoxel_id") %>%
  dplyr::distinct(root_610, root_626)
meta <- meta %>%
  dplyr::left_join(mappings, by=c("root_id"="root_610"))

# Write meta as feather
arrow::write_feather(meta,
                     file.path(banc.influence.save.path, "influence_banc_626_meta.feather"))

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'
sync_files(path(A, "influence/"), path(B, "influence/"), extensions = c("feather", "parquet"), move.old = FALSE)
remove_empty_dirs(path(A, "influence/"))

# Define the remote name
remote_name <- "hms"
local_files <- "/n/data1/hms/neurobio/wilson/banc/influence/"
local_files <- list.files(local_files, pattern = "feather$|parquet$", full.names = TRUE)
local_files
for(local_file in local_files){
  message("Working on: ", local_file)
  remote_path <- "neuroanat/influence/"
  system(paste("rclone copy", local_file, paste0(remote_name, ":", remote_path)))
  cat("File transfer complete.\n")
}

# Send
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/influence",
               parquet_file))
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/influence",
               file.path(banc.influence.save.path, "influence_banc_626_meta.feather")))
message("influence scores sent to google bucket")

# Connect to old influence data for export
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.influence.save.path,"influence_banc_610.sqlite"))
influence.db <- dplyr::tbl(con, "banc_610_influence") %>%
  dplyr::filter(!is_seed,
                level %in% c("seed_02")) %>%
  dplyr::collect()
dbDisconnect(con)

# Connect to old influence meta
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.influence.save.path,"influence_banc_610.sqlite"))
influence.meta <- dplyr::tbl(con, "meta") %>%
  dplyr::collect()
dbDisconnect(con)

# Write files to google bucket
seed03.df <- influence.db %>%
  dplyr::select(source = seed,
                root_610 = id,
                adjusted_influence = influence_norm_log) %>%
  dplyr::left_join(mappings,
                   by = c("root_610")) %>%
  dplyr::distinct() %>%
  dplyr::select(-root_610)

# Send
readr::write_csv(x = seed03.df,
                 file= file.path(banc.influence.save.path,"influence_downstream_of_sensory_cell_sub_classes.csv"))
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/influence",
               file.path(banc.influence.save.path,"influence_downstream_of_sensory_cell_sub_classes.csv")))
message("influence scores sent to google bucket")
