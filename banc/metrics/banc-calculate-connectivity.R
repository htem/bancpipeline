#' banc-calculate-connectivity — Build the versioned BANC synapse parquet + simple edgelist.
#'
#' v3 pulls CAVE-ingested synapses from GCS (not the local v3 parquet — see
#' Notes); v2 pulls the human-readable CSV.gz. Both branches dedupe autapses,
#' filter by size and SeaTable inclusion set.
#'
#' @section Reads:
#'   - GCS `synapses_<src>_human_readable.csv.gz`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `banc_<ver>_synapses_<src>.parquet`
#'   - `banc_<ver>_edgelist_simple_<src>.feather`
#'
#' @section CLI:
#'   --source {v2,v3}
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   BANC-project/R/startup/banc-edgelist.R (loads simple edgelist via `construct_path`).
#'
#' @section Schema:
#'   banc_888_edgelist_simple_v2.md.
#'
#' @section Paper:
#'   Methods §"Synapse detection".
#'
#' @section Notes:
#'   - Pre-2026-04-20 v3 read a locally-produced parquet; ~30% had spurious "0"
#'     root_ids from svid_cache contamination, so we pivoted to the GCS export.

###############################################
### Calculate versioned BANC connectivity   ###
###############################################
### Builds a versioned simple edgelist from either:
###   - v2 (synapses_v2_human_readable_*.csv.gz from GCS), or
###   - v3 (synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet
###         from the v<ver> GCS export — CAVE-ingested, trusted root_ids).
### Source selected via CLI arg `--source v2|v3` or BANC_SYN_SOURCE env var;
### defaults to banc.synapse.source.default.
###
### Output files (in banc.connectivity.save.path):
###   banc_{version}_synapses_v2.parquet        (v2 branch — local raw cache)
###   banc_{version}_edgelist_simple_{v2|v3}.feather
###   banc_{version}_meta.feather
###
### Requires: gsutil (both branches now — v3 also pulls from GCS).
###
### History: pre-2026-04-20, v3 branch read from a locally-produced
### banc_<ver>_synapses_v3.parquet (output of banc-synapses-v3-optimised.R).
### That file was found to have ~30% spurious "0" pre/post root_ids from
### svid_cache contamination — pivoted to GCS export. Local v3 parquet
### remains useful for spatial neuropil annotation only.
###############################################
source("banc/banc-startup.R")

local({

#####################################
### CONFIG — source selection   ###
#####################################

.parse_source <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == "--source")
  if (length(i) == 1 && length(args) >= i + 1) return(tolower(args[i + 1]))
  env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
  if (!is.na(env) && nzchar(env)) return(tolower(env))
  if (exists("banc.synapse.source.default")) return(tolower(banc.synapse.source.default))
  "v3"
}
syn_source <- .parse_source()
if (!syn_source %in% c("v2", "v3")) {
  stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", syn_source))
}

message(sprintf("### banc: calculating versioned connectivity [v%s | source=%s] ###",
                banc.version, syn_source))
t_total_start <- Sys.time()

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R

# Synapse size threshold (minimum synapse size to include)
banc.size.threshold <- 5

# GCS bucket prefix for published synapse tables. Mirrors the preprint-era
# v626 layout — public release artifacts live under `neuron_connectivity/v<ver>/`.
# (Pre-2026-05-13 this pointed at a flat `gs://lee-lab.../v<ver>/` stub that
# was never populated for v888; we now publish v888 to neuron_connectivity/v888/.)
gcs_bucket <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_connectivity"

dir.create(banc.connectivity.save.path, recursive = TRUE, showWarnings = FALSE)

parquet_file <- file.path(banc.connectivity.save.path,
                          sprintf("banc_%s_synapses_%s.parquet",
                                  banc.version, syn_source))

if (syn_source == "v2") {
  #######################
  ### V2: DOWNLOAD CSV
  #######################

  version.path <- file.path(banc.save.path, paste0("v", banc.version))
  dir.create(version.path, recursive = TRUE, showWarnings = FALSE)

  csv_file <- file.path(version.path, "synapses_v2_human_readable.csv")
  gz_file <- paste0(csv_file, ".gz")

  if (!file.exists(csv_file)) {
    gcs_uri <- sprintf("%s/v%s/synapses_v2_human_readable.csv.gz",
                       gcs_bucket, banc.version)
    message("Downloading synapse CSV from GCS: ", gcs_uri)
    system(sprintf("gsutil cp %s %s", gcs_uri, gz_file))

    message("Decompressing...")
    system(sprintf("gunzip %s", gz_file))

    stopifnot(file.exists(csv_file))
    message("Download complete: ", csv_file)
  } else {
    message("Synapse CSV already exists: ", csv_file)
  }

  ###########################
  ### V2: READ & SAVE PARQUET
  ###########################

  if (!file.exists(parquet_file)) {
    message("Reading synapse CSV (this may take a few minutes)...")
    t_csv_start <- Sys.time()

    # Column definitions (15 columns in the CSV)
    column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                      'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                      'pre_root_id', 'post_supervoxel_id', 'post_root_id')
    # Keep presyn / postsyn endpoint coords (pre_x/y/z, post_x/y/z) and the
    # supervoxel IDs alongside the centroid (ctr_*) — these ride through
    # neuropil-inclusion and land in banc_<ver>_synapses_v2_enriched.parquet
    # so downstream consumers don't have to re-pull the GCS CSV.
    desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id',
                         'ctr_x', 'ctr_y', 'ctr_z',
                         'pre_x', 'pre_y', 'pre_z',
                         'post_x', 'post_y', 'post_z',
                         'pre_supervoxel_id', 'post_supervoxel_id')
    col_types <- vroom::cols(
      id = vroom::col_character(),
      size = vroom::col_double(),
      pre_root_id = vroom::col_character(),
      post_root_id = vroom::col_character(),
      pre_supervoxel_id  = vroom::col_character(),
      post_supervoxel_id = vroom::col_character(),
      ctr_x = vroom::col_double(),
      ctr_y = vroom::col_double(),
      ctr_z = vroom::col_double(),
      .default = vroom::col_double()
    )

    synapses <- vroom::vroom(csv_file,
                             col_names = column_names,
                             col_select = dplyr::all_of(desired_columns),
                             col_types = col_types,
                             skip = 1) %>%
      dplyr::rename(X = ctr_x, Y = ctr_y, Z = ctr_z) %>%
      dplyr::filter(pre_root_id != post_root_id) %>%    # remove autapses
      dplyr::distinct(id, .keep_all = TRUE) %>%          # deduplicate
      tibble::as_tibble()

    message(sprintf("Raw synapses: %s rows [read in %s]",
                    format(nrow(synapses), big.mark = ","),
                    format(round(difftime(Sys.time(), t_csv_start, units = "mins"), 1))))

    # Save raw synapse parquet (input for neuropil inclusion script)
    write_connectome_data(synapses, parquet_file, format = "parquet")
    message("Synapse parquet saved: ", parquet_file)
  } else {
    message("Synapse parquet already exists: ", parquet_file)
    synapses <- arrow::read_parquet(parquet_file)
  }
} else {
  #######################
  ### V3: DOWNLOAD CSV.GZ (CAVE-INGESTED EXPORT FROM GCS)
  #######################
  ### As of v890 we read the CSV.gz directly — same pattern as the v2 branch
  ### above. Reason: the curated parquet
  ### `synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet`
  ### is not always staged at GCS for new versions, but the canonical CSV.gz
  ### always is. The CSV is the CAVE-ingested export with trusted root_ids.
  ###
  ### Schema (15 cols, NO header): same as v2 — id, pre_x/y/z, post_x/y/z,
  ### ctr_x/y/z, size, pre_supervoxel_id, pre_root_id, post_supervoxel_id,
  ### post_root_id. We pull only id/size/pre_root_id/post_root_id.
  ###
  ### Local v3 parquet (banc_<ver>_synapses_v3.parquet, locally produced)
  ### remains useful for spatial neuropil annotations only — its root_ids
  ### were contaminated by svid-cache reuse and are NOT trusted here.

  version.path <- file.path(banc.save.path, paste0("v", banc.version))
  dir.create(version.path, recursive = TRUE, showWarnings = FALSE)

  csv_file <- file.path(version.path, "synapses_v3_human_readable.csv")
  gz_file <- paste0(csv_file, ".gz")

  if (!file.exists(csv_file)) {
    gcs_uri <- sprintf("%s/v%s/synapses_v3_human_readable.csv.gz",
                       gcs_bucket, banc.version)
    message("Downloading v3 synapse CSV from GCS: ", gcs_uri)
    st <- system(sprintf("gsutil cp %s %s", gcs_uri, gz_file))
    if (st != 0L) {
      stop(sprintf("gsutil cp failed (status %d) for %s", st, gcs_uri))
    }
    message("Decompressing...")
    system(sprintf("gunzip %s", gz_file))
    stopifnot(file.exists(csv_file))
    message("Download complete: ", csv_file)
  } else {
    message("v3 synapse CSV already exists: ", csv_file)
  }

  message("Reading v3 CSV (this may take ~15-20 min)...")
  t_csv_start <- Sys.time()

  # v3 CSV has 17 columns (vs v2's 15): same layout as v2 but with
  # mean_score + median_score inserted between `size` and `pre_supervoxel_id`.
  # Verified 2026-05-15 against v888 CSV. Reading with the v2 15-col schema
  # silently mis-binds pre/post_root_id to the score columns and produces an
  # empty edgelist.
  column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                    'ctr_x', 'ctr_y', 'ctr_z', 'size',
                    'mean_score', 'median_score',
                    'pre_supervoxel_id', 'pre_root_id',
                    'post_supervoxel_id', 'post_root_id')
  desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id')
  col_types <- vroom::cols(
    id = vroom::col_character(),
    size = vroom::col_double(),
    pre_root_id = vroom::col_character(),
    post_root_id = vroom::col_character(),
    .default = vroom::col_double()
  )

  synapses <- vroom::vroom(csv_file,
                           col_names = column_names,
                           col_select = dplyr::all_of(desired_columns),
                           col_types = col_types,
                           skip = 0) %>%             # no header row in v3 CSV
    dplyr::filter(pre_root_id != post_root_id) %>%    # remove autapses
    dplyr::distinct(id, .keep_all = TRUE) %>%         # deduplicate
    tibble::as_tibble()

  message(sprintf("Raw v3 synapses (CSV, post-autapse-removal): %s rows [read in %s]",
                  format(nrow(synapses), big.mark = ","),
                  format(round(difftime(Sys.time(), t_csv_start, units = "mins"), 1))))
}

##############################
### GET VERSION-SPECIFIC IDs #
##############################

message("Fetching neuron metadata from SeaTable...")

banc.meta <- banctable_query()

# Broad inclusion: every SeaTable row, version-specific root_id with fallback
# to the existing root_id column when banc_rootid returns NA/"0". Mirrors
# banc-data.R Section 1 so the edgelist's inclusion universe matches the meta
# distributed alongside it — no super_class / status pre-filter, so glia,
# trachea, and ambiguous classes all stay in the inclusion set. Both pre AND
# post must land in this set (& filter below).
message(sprintf("Converting %s supervoxel IDs to v%s root IDs...",
                format(nrow(banc.meta), big.mark = ","), banc.version))
t_rootid_start <- Sys.time()
banc.meta$version_root_id <- banc_rootid(banc.meta$supervoxel_id,
                                          version = banc.version)
banc.meta <- banc.meta %>%
  dplyr::mutate(version_root_id = ifelse(
    is.na(version_root_id) | version_root_id == "0",
    root_id, version_root_id)) %>%
  dplyr::filter(!is.na(version_root_id), version_root_id != "0")

inclusion.ids <- as.character(banc.meta$version_root_id)
message(sprintf("  Inclusion set: %s root IDs [%s]",
                format(length(inclusion.ids), big.mark = ","),
                format(round(difftime(Sys.time(), t_rootid_start, units = "mins"), 1))))

############################
### FILTER & BUILD ELIST ###
############################

# Filter synapses to inclusion set and apply size threshold
message("Filtering synapses...")
synapses.filtered <- synapses %>%
  dplyr::filter(pre_root_id %in% inclusion.ids,
                post_root_id %in% inclusion.ids)

if (banc.size.threshold > 1) {
  synapses.filtered <- synapses.filtered %>%
    dplyr::filter(size >= banc.size.threshold)
}

message(sprintf("Filtered synapses: %s rows (from %s)",
                format(nrow(synapses.filtered), big.mark = ","),
                format(nrow(synapses), big.mark = ",")))

# Build simple edgelist: pre, post, count, norm, post_count, pre_count
message("Building simple edgelist...")
t_elist_start <- Sys.time()
edgelist <- synapses.filtered %>%
  dplyr::group_by(pre_root_id, post_root_id) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(pre = pre_root_id, post = post_root_id) %>%
  dplyr::mutate(pre = as.character(pre),
                post = as.character(post)) %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(norm = round(count / post_count, 6)) %>%
  dplyr::select(pre, post, count, norm, post_count, pre_count)

message(sprintf("Edgelist: %s connections between %s neurons [built in %s]",
                format(nrow(edgelist), big.mark = ","),
                format(length(unique(c(edgelist$pre, edgelist$post))), big.mark = ","),
                format(round(difftime(Sys.time(), t_elist_start, units = "secs"), 0))))

############################
### BUILD META           ###
############################

# Minimal meta with version-specific root IDs (broad — every banc.meta row)
meta <- banc.meta %>%
  dplyr::mutate(root_id = as.character(version_root_id)) %>%
  dplyr::select(-version_root_id)
# meta[meta == ""] <- NA

############################
### SAVE                 ###
############################

edgelist_file <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_edgelist_simple_%s.feather",
                                   banc.version, syn_source))
meta_file <- file.path(banc.connectivity.save.path,
                       sprintf("banc_%s_meta.feather", banc.version))

arrow::write_feather(edgelist, edgelist_file)
message("Edgelist saved: ", edgelist_file)

arrow::write_feather(meta, meta_file)
message("Meta saved: ", meta_file)

message(sprintf("### banc: connectivity calculation complete for v%s source=%s [total: %s] ###",
                banc.version, syn_source,
                format(round(difftime(Sys.time(), t_total_start, units = "mins"), 1))))

})
