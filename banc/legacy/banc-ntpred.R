################################
### Make BANC NT predictions ###
################################
source("banc/banc-startup.R")
poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin","tyramine")
bancsynapses.nts <- "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v2/banc_nt_prediction_w_sizethresh_5_11052025.parquet"
bancsynapses.nts <- "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v2/banc_nt_prediction_w_sizethresh_5_11102025.parquet"

# Direct us to the BANC dataset
bancr::choose_banc()

# Read IDs
banc.meta <- banctable_query("SELECT _id, root_id, root_626, supervoxel_id, position, cell_type, neurotransmitter_predicted, neurotransmitter_score from banc_meta")

# # Read IDs
# banc.ntpred <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ntpred.csv"),
#                                col_types = banc.col.types, 
#                                show_col_types = FALSE)
# banc.ntpred <- banc.ntpred %>% 
#   dplyr::filter(!is.na(neurotransmitter_predicted),
#                 root_id %in% banc.root.ids) %>%
#   dplyr::arrange(dplyr::desc(count)) %>%
#   dplyr::distinct(root_id, .keep_all = TRUE)
# banc.nt.done <- unique(banc.ntpred$root_id)

# # Which have been done already?
# banc.nt.todo <- setdiff(banc.root.ids, banc.nt.done)
# message(length(banc.nt.todo), " nt predictions to fetch, fetched: ", length(banc.nt.done))

# Get version
version.path <- file.path(banc.save.path,"v746")
banc.meta$root_746 <- banc_rootid(banc.meta$supervoxel_id,version = "746")
banc.meta <- banc.meta %>%
  dplyr::mutate(root_746 = banc_rootid(banc.meta$supervoxel_id,version = "746"),
                root_746 = ifelse(is.na(root_746)|root_746=="0",root_id,root_746))
proof.ids <- na.omit(unique(banc.meta$root_746))
proof.ids <- proof.ids[proof.ids!="0"]

# Batch for parallel processing
banc.root.ids <- na.omit(unique(banc.meta$root_746))
banc.nt.todo <- banc.root.ids
multiplier <- 1000
numCores <- 1
banc.nt.todo <- sample(banc.nt.todo)
upper <- ifelse((numCores*multiplier)<length(banc.nt.todo),numCores*multiplier,length(banc.nt.todo))
batches <- split(banc.nt.todo, round(seq(from = 1, to = upper, length.out = length(banc.nt.todo))))

# Register cores
cl <- setup_parallel()

# Set up progress bar
iterations <- length(batches)
pb <- utils::txtProgressBar(max = iterations, style = 3)
progress <- function(n) utils::setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Connect to .nt database
synapses.nt <- arrow::read_parquet(bancsynapses.nts)
max.prob <- max(synapses.nt$probability, na.rm = TRUE)

#  Get synapses
desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
col_types <- cols(
  id = col_character(),
  size = col_double(),
  pre_root_id = col_character(),
  post_root_id = col_character(),
  ctr_x = col_double(),
  ctr_y = col_double(),
  ctr_z = col_double(),
  .default = col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

# dir.create(version.path)
# system(sprintf("gsutil cp gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v746/synapses_v2_human_readable.csv.gz %s ",
#                file.path(version.path,"synapses_v2_human_readable.csv.gz")))
# system(sprintf("gunzip %s",file.path(version.path,"synapses_v2_human_readable.csv.gz")))
banc.syns <- vroom::vroom(file.path(version.path,"synapses_v2_human_readable.csv"),
                          col_names = column_names,
                          col_select = dplyr::all_of(desired_columns),
                          col_types = col_types,
                          skip = 1) %>%
  dplyr::rename(X=ctr_x,
                Y=ctr_y,
                Z=ctr_z) %>%
  dplyr::filter(pre_root_id %in% proof.ids,
                pre_root_id!=post_root_id,
  ) %>%
  tibble::as_tibble()

# Join
banc.syns.nt <- banc.syns %>%
  dplyr::left_join(synapses.nt %>%
                     dplyr::mutate(id = as.character(id)), 
                   by = "id") %>%
  dplyr::filter(!is.na(predicted_nt))

# Calculate the weighted sum by pre_pt_root_id and predicted_nt, plus count totals
transmitters <- c("acetylcholine", "dopamine", "gaba", "glutamate", 
                  "histamine", "octopamine", "serotonin", "tyramine")
total_counts <- banc.syns.nt %>%
  dplyr::filter(!is.na(predicted_nt)) %>%
  dplyr::group_by(pre_root_id) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop")

# Summed scores (and also counts, though not used below)
per_nt_scores <- banc.syns.nt %>%
  dplyr::filter(!is.na(predicted_nt)) %>%
  dplyr::group_by(pre_root_id, predicted_nt) %>%
  dplyr::summarise(
    neurotransmitter_score = sum(probability, na.rm = TRUE),
    neurotransmitter_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    pre_root_id,
    predicted_nt = transmitters,
    fill = list(neurotransmitter_score = 0, neurotransmitter_count = 0)
  )

# Find the predicted nt and the confidence percentage for each neuron
nt_pred_and_conf <- per_nt_scores %>%
  dplyr::group_by(pre_root_id) %>%
  dplyr::mutate(
    total_confidence = sum(neurotransmitter_score),
    max_score = max(neurotransmitter_score),
    most_confident_nt = predicted_nt[which.max(neurotransmitter_score)],
    confidence_pct = dplyr::if_else(total_confidence > 0,
                                    (max_score / total_confidence) * 100,
                                    0)
  ) %>%
  dplyr::slice_head(n = 1) %>%  # get just one row per neuron; all mutant cols same across NTs
  dplyr::ungroup() %>%
  dplyr::select(
    pre_root_id,
    neurotransmitter_predicted = most_confident_nt,
    neurotransmitter_score = confidence_pct,
  )

# Now make your full summary table (single row per neuron)
dat <- per_nt_scores %>%
  dplyr::select(pre_root_id, predicted_nt, neurotransmitter_count) %>%
  tidyr::pivot_wider(
    id_cols = pre_root_id,
    names_from = predicted_nt,
    values_from = neurotransmitter_count,
    names_prefix = ""
  ) %>%
  dplyr::left_join(nt_pred_and_conf, by = "pre_root_id") %>%
  dplyr::left_join(total_counts, by = "pre_root_id")

# Save
dat <- dat %>%
  dplyr::filter(pre_root_id %in% proof.ids,
                pre_root_id!="0") %>%
  dplyr::rename(root_id = pre_root_id) %>%
  dplyr::left_join(banc.meta %>%
                     dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
                     dplyr::select(root_id, supervoxel_id, position),
                   by = "root_id") %>%
  dplyr::mutate(neurotransmitter_score = round(neurotransmitter_score,4)/100)

# Get cell type average
dat <- dat %>%
  dplyr::left_join(banc.meta %>%
                     dplyr::distinct(root_id=root_746,cell_type),
                   by = "root_id")

# Calculate
cell_type_nt <- dat |>
  dplyr::filter(
    !is.na(cell_type),
    !is.na(neurotransmitter_predicted),
    !is.na(neurotransmitter_score)
  ) |>
  dplyr::group_by(cell_type, neurotransmitter_predicted) |>
  dplyr::summarise(
    total_count     = sum(count, na.rm = TRUE),
    cell_type_neurotransmitter_score  = sum(neurotransmitter_score * count, na.rm = TRUE) /
      sum(count, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  dplyr::slice_max(cell_type_neurotransmitter_score, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::rename(
    cell_type_neurotransmitter_predicted = neurotransmitter_predicted
  ) %>%
  dplyr::mutate(cell_type_neurotransmitter_score=round(cell_type_neurotransmitter_score,4)) %>%
  dplyr::select(cell_type,cell_type_neurotransmitter_predicted,cell_type_neurotransmitter_score)
dat <- dat %>%
  dplyr::left_join(cell_type_nt,by="cell_type")
  
# Send
readr::write_csv(dat, 
                   file = file.path(banc.meta.save.path,"banc_neurotransmitter_746_prediction.csv"), 
                   append = FALSE) 
system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v746/banc_neurotransmitter_prediction.csv",
       file.path(banc.meta.save.path,"banc_neurotransmitter_746_prediction.csv")))
message("NT predictions sent to google bucket")

#######################
### UPDATE SEATABLE ###
#######################

bc <- banctable_query("SELECT _id, root_id, supervoxel_id, cell_type, neurotransmitter_predicted_v1, neurotransmitter_predicted_v2, neurotransmitter_predicted_v3 from banc_meta")
bc$root_746 <- banc_rootid(bc$supervoxel_id,version = "746")
bc$neurotransmitter_predicted <- NA
bc$neurotransmitter_predicted <- NA
bc.update <- bc %>%
  dplyr::mutate(root_746 = ifelse(is.na(root_746)|root_746=="0",root_id,root_746)) %>%
  dplyr::select(`_id`, root_id, root_746, supervoxel_id, cell_type, neurotransmitter_predicted_v1, neurotransmitter_predicted_v2, neurotransmitter_predicted_v3) %>%
  dplyr::left_join(dat %>%
                     dplyr::filter(count>=10) %>%
                     dplyr::distinct(root_746=root_id,.keep_all=TRUE) %>%
                     dplyr::select(root_746=root_id,
                                     neurotransmitter_score,
                                     neurotransmitter_predicted),
                   by = "root_746") %>%
  dplyr::left_join(cell_type_nt %>%
                     dplyr::distinct(cell_type,
                                     cell_type_neurotransmitter_predicted),
                   by = "cell_type") %>%
  dplyr::mutate(neurotransmitter_score = round(neurotransmitter_score,4)) %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(`_id`)) %>%
  dplyr::select(`_id`, cell_type, cell_type_neurotransmitter_predicted, neurotransmitter_predicted, neurotransmitter_score, neurotransmitter_predicted_v1, neurotransmitter_predicted_v2, neurotransmitter_predicted_v3) %>%
  as.data.frame()
bc.update$neurotransmitter_score <- as.numeric(bc.update$neurotransmitter_score)
bc.update$neurotransmitter_score[is.na(bc.update$neurotransmitter_score)] <- 0
# bc.update$neurotransmitter_score_v2 <- as.numeric(bc.update$neurotransmitter_score_v2)
# bc.update$neurotransmitter_score_v2[is.na(bc.update$neurotransmitter_score_v2)] <- 0
# bc.update$neurotransmitter_score_v3 <- as.numeric(bc.update$neurotransmitter_score_v3)
# bc.update$neurotransmitter_score_v3[is.na(bc.update$neurotransmitter_score_v3)] <- 0
# bc.update$neurotransmitter_predicted <- bc.update$neurotransmitter_predicted_v3
# bc.update$neurotransmitter_score <- bc.update$neurotransmitter_score_v3
bc.update$cell_type <- NULL
if(nrow(bc.update)){
  banctable_update_rows(base='banc_meta', 
                        table = 'banc_meta', 
                        df = bc.update[,],
                        append_allowed = FALSE, 
                        chunksize = 1000)  
}

