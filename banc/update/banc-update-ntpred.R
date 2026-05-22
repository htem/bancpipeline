#' banc-update-ntpred — Push neurotransmitter predictions to SeaTable.
#'
#' Reads `banc_ntpred_<src>.feather` (from `banc-calculate-ntpred.R`) and
#' pushes per-neuron predictions to SeaTable's `neurotransmitter_predicted_v<N>`
#' column (v1 reserved for older detector; v2/v3 populated from
#' `--source v2/v3` runs). The `cell_type_neurotransmitter_predicted`
#' column is shared (cell-type consensus is source-agnostic).
#'
#' @section Reads:
#'   - `<banc.save.path>/banc_ntpred_<src>.feather`
#'   - SeaTable `banc_meta`
#'   - env var `BANC_SYN_SOURCE`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `neurotransmitter_predicted_v<2|3>`,
#'     `neurotransmitter_score_v<2|3>`,
#'     `cell_type_neurotransmitter_predicted`
#'
#' @section CLI:
#'   --source {v2,v3}   default banc.synapse.source.default
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_v888_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v850_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v890_rebuild.sh`

###########################################################
### Push neurotransmitter predictions to BANC seatable
###
### Reads banc_ntpred_<src>.feather (from banc-calculate-ntpred.R) and pushes
### per-neuron predictions to SeaTable's neurotransmitter_predicted_v<N>
### column (v1 is reserved for an older detector; we now populate v2 and v3
### from --source v2/v3 runs). The cell_type_neurotransmitter_predicted
### column is shared (cell-type consensus is source-agnostic).
###
### Source selection: --source v2|v3 CLI arg OR BANC_SYN_SOURCE env var
### (default: banc.synapse.source.default from banc-startup.R).
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
# SeaTable column to write per-neuron prediction into. v1 is reserved.
sea_col <- paste0("neurotransmitter_predicted_", syn_source)

message(sprintf("### banc: pushing NT predictions to seatable [source=%s → %s] ###",
                syn_source, sea_col))
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read predictions    ###
###########################

feather_file <- file.path(banc.save.path,
                          sprintf("banc_ntpred_%s.feather", syn_source))
if (!file.exists(feather_file)) {
  message(sprintf("%s not found — run banc-calculate-ntpred.R --source %s first.",
                  feather_file, syn_source))
  return(invisible())
}
dat <- arrow::read_feather(feather_file) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

cell_type_nt_file <- file.path(banc.save.path,
                                sprintf("banc_ntpred_cell_type_%s.feather",
                                        syn_source))
if (file.exists(cell_type_nt_file)) {
  cell_type_nt <- arrow::read_feather(cell_type_nt_file)
} else {
  cell_type_nt <- data.frame(cell_type = character(0),
                             cell_type_neurotransmitter_predicted = character(0))
}

###########################
### Build update        ###
###########################

bc <- banctable_query("SELECT _id, root_id, supervoxel_id, cell_type, neurotransmitter_predicted_v1, neurotransmitter_predicted_v2, neurotransmitter_predicted_v3 from banc_meta")
bc$root_888 <- banc_rootid(bc$supervoxel_id, version = banc.version)
bc.update <- bc %>%
  dplyr::mutate(root_888 = ifelse(is.na(root_888) | root_888 == "0", root_id, root_888)) %>%
  dplyr::select(`_id`, root_id, root_888, supervoxel_id, cell_type,
                neurotransmitter_predicted_v1, neurotransmitter_predicted_v2,
                neurotransmitter_predicted_v3) %>%
  dplyr::left_join(dat %>%
                     dplyr::filter(count >= 10) %>%
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::select(root_888 = root_id,
                                   neurotransmitter_score,
                                   neurotransmitter_predicted),
                   by = "root_888") %>%
  dplyr::left_join(cell_type_nt %>%
                     dplyr::distinct(cell_type,
                                     cell_type_neurotransmitter_predicted),
                   by = "cell_type") %>%
  dplyr::mutate(neurotransmitter_score = round(neurotransmitter_score, 4)) %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(`_id`)) %>%
  as.data.frame()

# Route per-neuron prediction into the right source-versioned SeaTable column,
# leaving the other (v1 / sibling source) columns untouched.
bc.update[[sea_col]] <- bc.update$neurotransmitter_predicted
bc.update$neurotransmitter_predicted <- NULL

bc.update$neurotransmitter_score <- as.numeric(bc.update$neurotransmitter_score)
bc.update$neurotransmitter_score[is.na(bc.update$neurotransmitter_score)] <- 0
bc.update$cell_type <- NULL

# Final select preserves the two sibling-source columns we don't write so the
# update doesn't clobber them, plus the cell-type consensus + our target col.
keep_cols <- c("_id", "cell_type_neurotransmitter_predicted",
               "neurotransmitter_score",
               "neurotransmitter_predicted_v1",
               "neurotransmitter_predicted_v2",
               "neurotransmitter_predicted_v3")
# Ensure the target col holds the freshly-written values (override the SeaTable
# echo from the earlier SELECT).
bc.update[[sea_col]] <- bc.update[[sea_col]]
bc.update <- bc.update[, intersect(keep_cols, colnames(bc.update)), drop = FALSE]

###########################
### Push to seatable    ###
###########################

message(sprintf("Pushing NT predictions for %d neurons to seatable", nrow(bc.update)))

if (nrow(bc.update)) {
  banctable_update_rows(base = 'banc_meta',
                        table = 'banc_meta',
                        df = bc.update,
                        append_allowed = FALSE,
                        chunksize = 1000)
}

message(sprintf("### banc: NT prediction push complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
