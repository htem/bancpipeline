#######################################
### COLLATE FRANKEN BRAIN INFLUENCE ###
#######################################
source("banc/banc-startup.R")
library(progress)


######################################
### INFLUENCE CALCULATION FUNCTION ###
######################################

# function
calculate_influence_norms <- function(influence.df,
                                      const = -24){
  inf.threshold <- exp(const)
  if(!"target"%in%colnames(influence.df)){
    influence.df$target <- influence.df$id
    orig.target = FALSE
  }else{
    orig.target = TRUE
  }
  if(!"influence_syn_original"%in%colnames(influence.df)){
    influence.df$influence_syn_norm <- 1
  }
  if(!"influence_original"%in%colnames(influence.df)){
    influence.df$influence_original <- influence.df$influence
  }
  if(!"influence_norm_original"%in%colnames(influence.df)){
    influence.df$influence_norm_original <- influence.df$influence_norm
  }
  influence.df <- influence.df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(no_seeds = influence_original/influence_norm_original,
                  no_seeds = ifelse(is.na(no_seeds),1,no_seeds),
                  no_synapses = influence_original/influence_syn_norm) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(no_targets = length(unique(id))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      no_seeds = as.numeric(no_seeds),
      no_targets = as.numeric(no_targets)
    ) %>%
    dplyr::group_by(seed) %>%
    dplyr::mutate(influence_per_seed = influence_original*no_seeds,
                  influence_per_synapse = influence_original*no_synapses) %>%
    dplyr::group_by(target, seed) %>%
    dplyr::mutate(influence = sum(influence_original,na.rm = TRUE),
                  total_seeds = sum(no_seeds,na.rm=TRUE),
                  total_synapses = sum(no_synapses,na.rm=TRUE),
                  influence_norm = sum(influence_per_seed,na.rm = TRUE)/(total_seeds*no_targets),
                  influence = ifelse(influence<inf.threshold,inf.threshold,influence),
                  influence_norm = ifelse(influence_norm<inf.threshold,inf.threshold,influence_norm),
                  influence_norm = sum(influence_norm,na.rm = TRUE),
                  influence_syn_norm =  sum(influence_per_synapse,na.rm = TRUE)/(no_targets*total_synapses),
                  influence_syn_norm = ifelse(influence_syn_norm<inf.threshold,inf.threshold,influence_syn_norm),
                  influence_syn_norm = sum(influence_syn_norm,na.rm = TRUE),
                  influence_norm_log = log(influence_norm),
                  influence_log = log((influence/no_targets)),
                  influence_syn_norm_log = log(influence_syn_norm)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(influence_norm_log = influence_norm_log-const,
                  influence_log = influence_log-const,
                  influence_syn_norm_log = influence_syn_norm_log-const) %>%
    dplyr::group_by(seed) %>%
    dplyr::mutate(influence_log = ifelse(is.na(influence),0,influence_log)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(target, seed, .keep_all = TRUE) %>%
    dplyr::mutate(influence = signif(influence,4),
                  influence_log = signif(influence_log,4),
                  influence_norm = signif(influence_norm,4),
                  influence_norm_log = signif(influence_norm_log,4)
    ) %>%
    dplyr::distinct(target,
                    seed, 
                    .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(influence = signif(influence,4),
                  influence_log = signif(influence_log,4),
                  influence_norm = signif(influence_norm,4),
                  influence_syn_norm = signif(influence_syn_norm,4),
                  influence_norm_log = signif(influence_norm_log,4),
                  influence_syn_norm_log = signif(influence_syn_norm_log,4)
    ) 
  if(!orig.target){
    influence.df$target <- NULL
  }
  influence.df
}

# SQL column types
inf.col.types <- readr::cols(
  id = readr::col_character(),
  is_seed = readr::col_logical(),
  .default = readr::col_number())

####################
### WRANGLE DATA ###
####################

# Get current meta
banc.meta <- banctable_query() %>%
  dplyr::filter(!is.na(root_626)) %>%
  dplyr::distinct(root_id,.keep_all = TRUE)

# Get orig meta
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.connectivity.save.path,"banc_610_data.sqlite"))
meta <- dplyr::tbl(con, "meta") %>%
  dplyr::collect()
dbDisconnect(con)

# Get ID mapping
mappings <- banc.meta %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%  
  dplyr::select(root_626, supervoxel_id) %>%
  dplyr::left_join(meta %>%
                     dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
                     dplyr::select(root_610 = root_id, supervoxel_id),
                   by = "supervoxel_id") %>%
  dplyr::mutate(change = root_626!=root_610)

# Master folder
seed13.save.path <- "/n/data1/hms/neurobio/wilson/banc/influence/banc_610/seed_13"

# Rename files
changed.mappings <- mappings %>%
  dplyr::filter(change)
changed.ids <- changed.mappings$root_610
for(changed.id in changed.ids){
  file.to.change <- file.path(seed13.save.path,paste0(changed.id,"_influence.csv"))
  new.name <- file.path(seed13.save.path,paste0(changed.mappings$root_626[match(changed.id,changed.mappings$root_610)],"_influence.csv"))
  if(file.exists(file.to.change)){
    file.rename(file.to.change, new.name)
  }
}

# Send into google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/banc/influence/banc_610/seed_13/ gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/all_influence_downstream")
message("seed influence sent to google bucket")

# Send skeletons into google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/banc/l1/ gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/l2")
message("seed influence sent to google bucket")

# Wrangle output influence
banc.eff.meta <- banc.meta %>%
  dplyr::filter(flow=="efferent") %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
  is.na(cell_sub_class) ~ cell_class,
  TRUE ~ cell_sub_class
  )) %>%
  dplyr::left_join(mappings %>%
                     dplyr::distinct(root_610,supervoxel_id),
                   by = "supervoxel_id")
post.ids <- na.omit(unique(banc.eff.meta$root_610))
pre.ids <- na.omit(unique(mappings$root_610))
progressr::with_progress({
  p <- progressr::progressor(steps = length(pre.ids))
  inf.data <- foreach(pre.id = pre.ids, 
                      .packages = c("dplyr", "readr"), 
                      .errorhandling = 'pass') %do% {
                        
                        # Read all data
                        csv <- file.path(seed13.save.path,paste0(pre.id,"_influence.csv"))
                        if(!file.exists(csv)){
                          message("no results for: ", pre.id)
                          return(NULL)
                        }
                        data <- readr::read_csv(csv,
                                                show_col_types = FALSE, 
                                                progress = FALSE,
                                                col_types = inf.col.types, )
                        data <- data[,-1]
                        colnames(data) <- c("id","is_seed","influence")
                        n.seeds <- sum(data$is_seed, na.rm = TRUE)
                        id.seeds <- subset(data, is_seed)$id
                        syns.seeds <- subset(banc.meta, root_626 %in% id.seeds | 
                                               supervoxel_id %in% subset(banc.meta, root_id %in% id.seeds)$supervoxel_id)$output_connections
                        n.syns.seed <- sum(syns.seeds, na.rm = TRUE)
                        if(n.syns.seed==0){
                          n.syns.seed <- 1
                        }
                        
                        # Skip if no influence
                        if(all(data$influence=="0")){
                          message("No values for: ", csv)
                          return(NULL)
                        }
                        
                        # Apply seed label
                        data <- data %>%
                          dplyr::filter(id %in% post.ids) %>%
                          dplyr::left_join(banc.eff.meta %>%
                                             dplyr::distinct(root_610, .keep_all = TRUE) %>%
                                             dplyr::select(root_610, root_626, target = cell_sub_class),
                                           by = c("id"="root_610")) %>%
                          dplyr::mutate(seed = gsub(".*from_|\\.csv*|_influence.*","",basename(csv)),
                                        level = "seed_13") 
                        
                        # Perform scaling operations
                        data$influence_original <- as.numeric(data$influence)
                        data$influence_norm_original <- as.numeric(data$influence)/n.seeds
                        data$influence_syn_norm_original <- as.numeric(data$influence)/n.syns.seed
                        data <- calculate_influence_norms(data) %>%
                          dplyr::distinct(root_626 = seed, 
                                        target,
                                        adjusted_influence = influence_norm_log)
                        p()
                        data
                      } 
})

# combine
eff.inf.data <- do.call(rbind,inf.data)

# Send
readr::write_csv(x = eff.inf.data, 
                 file= file.path(banc.influence.save.path,"influence_upstream_of_effector_cell_sub_classes.csv"))
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v626/influence",
               file.path(banc.influence.save.path,"influence_upstream_of_effector_cell_sub_classes.csv")))
message("influence scores sent to google bucket")

####################
### WRANGLE DATA ###
####################
overwrite <- TRUE

# Get current meta
banc.meta <- banctable_query() %>%
  dplyr::arrange(dplyr::desc(output_connections)) %>%
  dplyr::filter(!is.na(root_626)) %>%
  dplyr::distinct(root_id,.keep_all = TRUE)

# Get orig meta
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.connectivity.save.path,"banc_610_data.sqlite"))
meta <- dplyr::tbl(con, "meta") %>%
  dplyr::collect()
dbDisconnect(con)

# Master folder
banc.influence.save.path <- "/n/data1/hms/neurobio/wilson/banc/influence/"
influence.master <- file.path(banc.influence.save.path,"banc_610/")
seed.groups <- list.files(influence.master, recursive = FALSE, full.names = TRUE)
seed.groups <- seed.groups[grepl("_13$",seed.groups)]
seed.sets <- list.files(seed.groups, recursive = FALSE, full.names = TRUE)
seed.sets <- seed.sets[grepl("_13$",seed.sets)]

# Create .sql file
table_name <- "influence_all_by_all"
dir.create(banc.influence.save.path, showWarnings = FALSE)
if(overwrite){
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        file.path(banc.influence.save.path,"banc_626_all_by_all_data.sqlite"))
  DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", table_name))
  dbDisconnect(con)  
}

# SQL column types
inf.col.types <- readr::cols(
  id = readr::col_character(),
  is_seed = readr::col_logical(),
  .default = readr::col_number())

# Iterate over seed groups
for(sg in seed.groups){
  message("Working on: ", basename(sg))
  influence.list <- list()
  seed.level <- tolower(gsub(".*from_|\\.csv*|_influence.*","",basename(sg)))
  
  # Overwrite?
  if(!overwrite){
    con <- DBI::dbConnect(RSQLite::SQLite(),
                          file.path(banc.influence.save.path,"banc_626_all_by_all_data.sqlite"))
    influence.db <- dplyr::tbl(con, table_name) %>%
      dplyr::filter(level == seed.level) %>%
      dplyr::collect()
    dbDisconnect(con)
    if(nrow(influence.db)){
      seeds.done <- unique(influence.db$seed)
    }else{
      seeds.done <- NULL
    }
    rm('influence.db') 
  }else{
    seeds.done <- NULL
  }
  
  # Set up error bar
  seed.files <- list.files(sg, recursive = FALSE, full.names = TRUE, pattern = "\\.csv")
  seed.seeds <- gsub(".*from_|\\.csv*|_influence.*","",basename(seed.files))
  seed.files <- seed.files[!seed.seeds%in%seeds.done]
  if(!length(seed.files)){
    next
  }
  
  # Batch for parallel processing
  multiplier <- 10000
  numCores <- 1
  seed.files <- sample(seed.files)
  upper <- ifelse((numCores*multiplier)<length(seed.files),numCores*multiplier,length(seed.files))
  batches <- split(seed.files, round(seq(from = 1, to = upper, length.out = length(seed.files))))
  
  # Register cores
  cl <- setup_parallel()

  # Iterate
  with_progress({
    p <- progressor(steps = length(batches))
    for(batch in 1:length(batches)){
      batch.seed.files <- batches[[batch]]
      p()
      all_data <- foreach(csv = batch.seed.files, 
                          .combine = dplyr::bind_rows, 
                          .packages = c("readr", "dplyr", "DBI", "RSQLite")) %do% {
                            
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
                          }
      # Write to sql data base
      con <- DBI::dbConnect(RSQLite::SQLite(),
                            file.path(banc.influence.save.path,"banc_626_all_by_all_data.sqlite"))
      DBI::dbWriteTable(con,
                        name = table_name,
                        value = all_data,
                        overwrite = FALSE,
                        append = TRUE)
      dbDisconnect(con)
    }
  })
}

# Write to sql data base
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.influence.save.path,"banc_626_all_by_all_data.sqlite"))
DBI::dbWriteTable(con,
                  name = "meta",
                  value = meta,
                  overwrite = TRUE,
                  append = FALSE)
dbDisconnect(con) 

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'
sync_files(path(A, "influence/"), path(B, "influence/"), extensions = c("sqlite"), move.old = FALSE)
remove_empty_dirs(path(A, "influence/"))

# Define the remote name
remote_name <- "hms"
local_files <- "/n/data1/hms/neurobio/wilson/banc/influence/"
local_files <- list.files(local_files, pattern = "sqlite$", full.names = TRUE)
for(local_file in local_files){
  message("Working on: ", local_file)
  remote_path <- "neuroanat/influence/"
  system(paste("rclone copy", local_file, paste0(remote_name, ":", remote_path)))
  cat("File transfer complete.\n")  
}

# Send
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/influence",
               file.path(banc.influence.save.path,"banc_626_all_by_all_data.sqlite")))
message("influence scores sent to google bucket")

