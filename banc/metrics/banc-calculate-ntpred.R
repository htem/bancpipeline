#' banc-calculate-ntpred — Per-neuron neurotransmitter predictions from per-synapse classifier.
#'
#' Joins external `synister_banc` per-synapse NT parquet (K. Dasari) to the
#' versioned v2/v3 synapse table; computes per-neuron NT + confidence plus
#' cell-type consensus.
#'
#' @section Reads:
#'   - `synister_banc/banc_nt_prediction_w_sizethresh_{5,10}_*.parquet`
#'   - `banc_<ver>_synapses_<src>.parquet`, SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `banc_ntpred_<src>.feather`
#'   - `banc_neurotransmitter_<ver>_prediction_<src>.csv` (+ GCS mirrors)
#'
#' @section CLI:
#'   --source {v2,v3}
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Schema:
#'   banc_888_neurotransmitter_prediction_v2.md.
#'
#' @section Paper:
#'   Methods §"Neurotransmitter prediction".

###########################################################
### Calculate neurotransmitter predictions for BANC neurons
###
### Reads synapse-level NT predictions from parquet,
### joins with the current-version v2 or v3 synapse table, computes per-neuron
### predicted NT and confidence scores, plus cell-type consensus predictions.
###
### Source selection: --source v2|v3 CLI arg OR BANC_SYN_SOURCE env var
### (default: banc.synapse.source.default from banc-startup.R; usually v3).
###
### Outputs (source-suffixed):
###   banc.save.path/banc_ntpred_<src>.feather
###   banc.meta.save.path/banc_neurotransmitter_<version>_prediction_<src>.csv
###   GCS: v<ver>/banc_neurotransmitter_prediction_<src>.csv
###   GCS: compiled_data/banc_<ver>/banc_<ver>_neurotransmitter_prediction_<src>.csv
###########################################################
source("banc/banc-startup.R")

local({

###########################
### Source selection    ###
###########################
.parse_source <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == "--source")
  if (length(i) == 1 && length(args) >= i + 1) return(tolower(args[i + 1]))
  env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
  if (!is.na(env) && nzchar(env)) return(tolower(env))
  if (exists("banc.synapse.source.default")) return(tolower(banc.synapse.source.default))
  "v2"
}
syn_source <- .parse_source()
if (!syn_source %in% c("v2", "v3")) {
  stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", syn_source))
}

message(sprintf("### banc: calculating NT predictions [source=%s] ###", syn_source))
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Configuration      ###
###########################

poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin","tyramine")

# Per-synapse NT prediction parquet — paths match the kd193 inference results
# and the v2.0 / v3.0 layout under gs://lee-lab.../synapses/.
bancsynapses.nts <- switch(syn_source,
  v2 = "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v2/banc_nt_prediction_w_sizethresh_5_11102025.parquet",
  v3 = "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v3/banc_nt_prediction_v3_w_sizethresh_10_05042026.parquet"
)
stopifnot(file.exists(bancsynapses.nts))

###########################
### Read current state  ###
###########################

banc.meta <- banctable_query("SELECT _id, root_id, root_626, root_888, root_888, supervoxel_id, position, cell_type, neurotransmitter_predicted, neurotransmitter_score from banc_meta")

# Get current-version root IDs
version.path <- file.path(banc.save.path, paste0("v", banc.version))
banc.meta <- banc.meta %>%
  dplyr::mutate(root_888 = banc_rootid(banc.meta$supervoxel_id, version = banc.version),
                root_888 = ifelse(is.na(root_888) | root_888 == "0", root_id, root_888))
proof.ids <- na.omit(unique(banc.meta$root_888))
proof.ids <- proof.ids[proof.ids != "0"]

banc.root.ids <- na.omit(unique(banc.meta$root_888))
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.root.ids <- intersect(banc.root.ids, banc.test.ids)

###########################
### Read synapse data   ###
###########################

# Read NT prediction parquet (source-specific path resolved above)
synapses.nt <- arrow::read_parquet(bancsynapses.nts)
max.prob <- max(synapses.nt$probability, na.rm = TRUE)

# Read current-version synapse CSV — same 15-col schema for v2 and v3,
# only difference is the source file name. v2 has a header row (skip=1);
# v3 has no header (skip=0). Both are produced by banc-calculate-connectivity.R
# (downloaded + decompressed from GCS).
desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
col_types <- vroom::cols(
  id = vroom::col_character(),
  size = vroom::col_double(),
  pre_root_id = vroom::col_character(),
  post_root_id = vroom::col_character(),
  ctr_x = vroom::col_double(),
  ctr_y = vroom::col_double(),
  ctr_z = vroom::col_double(),
  .default = vroom::col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

syn_csv <- file.path(version.path,
                     sprintf("synapses_%s_human_readable.csv", syn_source))
skip_n  <- if (syn_source == "v2") 1L else 0L      # v2 has header, v3 doesn't
stopifnot(file.exists(syn_csv))

banc.syns <- vroom::vroom(syn_csv,
                          col_names = column_names,
                          col_select = dplyr::all_of(desired_columns),
                          col_types = col_types,
                          skip = skip_n) %>%
  dplyr::rename(X = ctr_x, Y = ctr_y, Z = ctr_z) %>%
  dplyr::filter(pre_root_id %in% proof.ids,
                pre_root_id != post_root_id) %>%
  tibble::as_tibble()

###########################
### Compute predictions ###
###########################

# Join synapses with NT predictions
banc.syns.nt <- banc.syns %>%
  dplyr::left_join(synapses.nt %>%
                     dplyr::mutate(id = as.character(id)),
                   by = "id") %>%
  dplyr::filter(!is.na(predicted_nt))

# Total synapse counts per neuron
transmitters <- c("acetylcholine", "dopamine", "gaba", "glutamate",
                  "histamine", "octopamine", "serotonin", "tyramine")
total_counts <- banc.syns.nt %>%
  dplyr::filter(!is.na(predicted_nt)) %>%
  dplyr::group_by(pre_root_id) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop")

# Per-NT scores and counts
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

# Find predicted NT and confidence percentage for each neuron
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
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    pre_root_id,
    neurotransmitter_predicted = most_confident_nt,
    neurotransmitter_score = confidence_pct,
  )

# Full summary table (single row per neuron)
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

dat <- dat %>%
  dplyr::filter(pre_root_id %in% proof.ids,
                pre_root_id != "0") %>%
  dplyr::rename(root_id = pre_root_id) %>%
  dplyr::left_join(banc.meta %>%
                     dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
                     dplyr::select(root_id = root_888, supervoxel_id, position),
                   by = "root_id") %>%
  dplyr::mutate(neurotransmitter_score = round(neurotransmitter_score, 4) / 100)

# Cell type consensus
dat <- dat %>%
  dplyr::left_join(banc.meta %>%
                     dplyr::distinct(root_id = root_888, cell_type),
                   by = "root_id")

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
  dplyr::mutate(cell_type_neurotransmitter_score = round(cell_type_neurotransmitter_score, 4)) %>%
  dplyr::select(cell_type, cell_type_neurotransmitter_predicted, cell_type_neurotransmitter_score)

dat <- dat %>%
  dplyr::left_join(cell_type_nt, by = "cell_type")

###########################
### Save results        ###
###########################

# Save feather (source-suffixed). Legacy `banc_ntpred.feather` (un-suffixed)
# is the v2-source artifact; keep writing it for one cycle as a back-compat
# symlink target so existing consumers (banc-update-ntpred.R, banc-data.R
# Section 5 NT join in v2 branch) don't break before we wire their source-
# aware variants. Drop after the next consumer-side refactor.
feather_file <- file.path(banc.save.path,
                          sprintf("banc_ntpred_%s.feather", syn_source))
arrow::write_feather(dat, feather_file)
message(sprintf("Wrote %s (%d neurons)", basename(feather_file), nrow(dat)))
if (syn_source == "v2") {
  legacy_feather <- file.path(banc.save.path, "banc_ntpred.feather")
  arrow::write_feather(dat, legacy_feather)
  message("Also wrote legacy banc_ntpred.feather for back-compat")
}

# Save CSV + GCS upload (source-suffixed in filename)
ntpred_csv <- file.path(banc.meta.save.path,
                         sprintf("banc_neurotransmitter_%s_prediction_%s.csv",
                                 banc.version, syn_source))
readr::write_csv(dat, file = ntpred_csv, append = FALSE)
# Source-suffixed (_v2_ / _v3_) — these predictions are computed from the
# named synapse table; include the source in the filename so downstream
# consumers know which synapse table backs the NT calls. Ship only to the
# canonical compiled_data/banc_<ver>/ location; the legacy flat
# lee-lab/v<ver>/ path was decommissioned 2026-05-13.
system(sprintf("gsutil cp %s gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/banc_%s_neurotransmitter_prediction_%s.csv",
               ntpred_csv, banc.version, banc.version, syn_source))
message(sprintf("NT predictions (%s) sent to google bucket", syn_source))

# Save cell type consensus separately (source-suffixed)
cell_type_nt_file <- file.path(banc.save.path,
                                sprintf("banc_ntpred_cell_type_%s.feather",
                                        syn_source))
arrow::write_feather(cell_type_nt, cell_type_nt_file)
if (syn_source == "v2") {
  legacy_ct <- file.path(banc.save.path, "banc_ntpred_cell_type.feather")
  arrow::write_feather(cell_type_nt, legacy_ct)
}

message(sprintf("### banc: NT prediction calculation complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
