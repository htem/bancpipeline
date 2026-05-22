#' banc-dcv-density — Push per-neuron soma DCV count + density to SeaTable.
#'
#' Reads the per-neuron DCV-in-soma parquet from GCS, joins to SeaTable on
#' `nucleus_id`, and pushes `soma_dcv_count` (3D detection count) and
#' `soma_dcv_density` (dcv_voxels / soma_voxels). No `root_id` needed.
#'
#' @section Reads:
#'   - GCS `lee-lab_..._/dcv/banc_dcv_predictions_in_soma_per_neuron_*.parquet`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `soma_dcv_count`, `soma_dcv_density`

###########################################################
### DCV density: read soma DCV predictions from GCS,
### join to seatable on nucleus_id, and push
### soma_dcv_count and soma_dcv_density
###
### Data source: GCS parquet with per-neuron DCV
### predictions in soma. Columns:
###   nucleus_id, dcv_3d_count, dcv_density_by_dcv_voxels
###
### SeaTable columns updated (numeric):
###   soma_dcv_count   — 3D DCV detection count
###   soma_dcv_density — dcv_voxels / soma_voxels
###
### Join key: nucleus_id (no root_id needed)
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: DCV density update ###")
t_start <- Sys.time()

###########################
### Read DCV parquet    ###
###########################

gcs_parquet <- "lee-lab_brain-and-nerve-cord-fly-connectome/dcv/banc_dcv_predictions_in_soma_per_neuron_02262026.parquet"
message("Reading DCV predictions from GCS: ", gcs_parquet)

gcs <- arrow::GcsFileSystem$create()
dcv <- arrow::read_parquet(gcs$path(gcs_parquet))
message(sprintf("  Loaded %d rows from parquet", nrow(dcv)))

# Rename to seatable column names
dcv <- dcv %>%
  dplyr::mutate(nucleus_id = as.character(nucleus_id)) %>%
  dplyr::rename(soma_dcv_count = dcv_3d_count,
                soma_dcv_density = dcv_density_by_dcv_voxels) %>%
  dplyr::filter(!is.na(nucleus_id), nucleus_id != "")

# De-duplicate: keep highest DCV density per nucleus_id
n.dup <- sum(duplicated(dcv$nucleus_id))
if (n.dup > 0) {
  dcv <- dcv %>%
    dplyr::group_by(nucleus_id) %>%
    dplyr::slice_max(soma_dcv_density, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  message(sprintf("  De-duplicated %d rows with duplicate nucleus_ids", n.dup))
}

# Keep only the columns we need
dcv <- dcv %>%
  dplyr::select(nucleus_id, soma_dcv_count, soma_dcv_density)

message(sprintf("  %d neurons with DCV data", nrow(dcv)))

###########################
### Join to seatable    ###
###########################

bc <- banctable_query("SELECT _id, nucleus_id, soma_dcv_count, soma_dcv_density FROM banc_meta") %>%
  dplyr::filter(!is.na(nucleus_id), nucleus_id != "") %>%
  dplyr::mutate(nucleus_id = as.character(nucleus_id)) %>%
  dplyr::distinct(nucleus_id, .keep_all = TRUE)

message(sprintf("  %d seatable rows with nucleus_id", nrow(bc)))

# Join DCV data onto seatable rows
bc.merged <- bc %>%
  dplyr::select(`_id`, nucleus_id) %>%
  dplyr::left_join(dcv, by = "nucleus_id") %>%
  dplyr::filter(!is.na(soma_dcv_count) | !is.na(soma_dcv_density))

# Only push rows where values actually changed
bc.changed <- bc.merged %>%
  dplyr::anti_join(bc, by = c("_id", "soma_dcv_count", "soma_dcv_density"))

###########################
### Push to seatable    ###
###########################

if (nrow(bc.changed)) {
  push.df <- bc.changed %>%
    dplyr::select(`_id`, soma_dcv_count, soma_dcv_density) %>%
    as.data.frame()

  message(sprintf("Updating soma_dcv_count and soma_dcv_density for %d neurons", nrow(push.df)))

  # Push both columns together (2 data columns avoids single-column
  # serialization bug in banc_df2updatepayload where column names are lost).
  # bancr's jsonlite::toJSON(na = "null") handles any remaining NAs.
  banctable_update_rows(base = 'banc_meta',
                        table = "banc_meta",
                        df = push.df,
                        append_allowed = FALSE,
                        chunksize = 1000)
} else {
  message("No DCV density changes to push")
}

message(sprintf("### banc: DCV density update complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
