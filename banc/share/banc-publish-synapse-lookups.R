#' banc-publish-synapse-lookups — Publish version-stable syn_id → neuropil/region/side parquets.
#'
#' Synapse `id` is stable across BANC root_id versions; downstream consumers
#' join once and reuse.
#'
#' @section Reads:
#'   - per-version v1/v2/v3 neuropil parquets (inline paths)
#'
#' @section Writes:
#'   - GCS `synapses/v{1.1,2.0,3.0}/synapse_neuropil_lookup_v{1,2,3}.parquet`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_publish_lookups.sh`.
#'
#' @section Schema:
#'   synapse_neuropil_lookup_v2.md.
#'
#' @section Notes:
#'   - Deliberately does NOT source `banc/banc-startup.R`: the natverse + rgl
#'     stack plus the 12 GiB v1 CSV pushed peak RSS over 24 GiB.

###########################################################
### Publish per-version syn_id -> neuropil/region/side
### lookups to GCS.
###
### Synapse `id` is stable across segmentation versions; only
### root_id assignments change. So neuropil + region + side
### labels for a given synapse detection round (v1, v2, v3) can
### be published once and reused across root_id versions.
###
### Builds slim (id, neuropil, region, side) parquets for v1,
### v2, and v3 detection rounds and uploads to
###   gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses/v1.1/
###   gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses/v2.0/
###   gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses/v3.0/
### with filenames synapse_neuropil_lookup_v{1,2,3}.parquet.
###
### NOTE: deliberately does NOT source banc/banc-startup.R — bancr+nat+rgl
### startup adds ~2 GiB working set, and combined with the 12 GiB v1 CSV
### that pushed peak RSS over 24 GiB on the first attempt. Paths are
### inlined here.
###
### Usage: Rscript banc/share/banc-publish-synapse-lookups.R
###########################################################

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# v3 spatial parquet has embedded nul bytes in some string columns
# (region/side leftovers from the alpha-shape pipeline) — strip rather than
# crash on read.
options(arrow.skip_nul = TRUE)

banc.version <- "888"
banc.connectivity.save.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity"
banc.synapses.v3.path       <- "/n/data1/hms/neurobio/wilson/banc/synapses_v3"

local({

message("### banc: publish synapse neuropil lookups ###")
t_start <- Sys.time()

bucket <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses"
staging <- file.path(banc.connectivity.save.path, "publish_lookups")
dir.create(staging, recursive = TRUE, showWarnings = FALSE)

publish_one <- function(out_path, version_label, gcs_subdir, n_rows) {
  size_h <- format(structure(file.size(out_path), class = "object_size"),
                   units = "auto")
  message(sprintf("  Wrote %s (%s rows, %s)",
                  out_path, format(n_rows, big.mark = ","), size_h))
  target <- sprintf("%s/%s/synapse_neuropil_lookup_%s.parquet",
                    bucket, gcs_subdir, version_label)
  rc <- system(sprintf("gsutil cp %s %s", out_path, target))
  if (rc != 0) stop(sprintf("gsutil cp failed (rc=%d) for %s", rc, target))
  message(sprintf("  Uploaded -> %s", target))
}

#####################################################
### v1 — from CSV (id, neuropil, region, side)    ###
### Streamed via arrow::open_dataset to avoid      ###
### loading 12 GiB into R memory.                  ###
#####################################################
v1_csv <- file.path(banc.connectivity.save.path,
                    "banc_synapses_to_neuropils_v1.csv")
v1_out <- file.path(staging, "synapse_neuropil_lookup_v1.parquet")
if (file.exists(v1_out)) {
  v1_n <- arrow::open_dataset(v1_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  message(sprintf("v1: staged file already exists (%s rows) — skipping rebuild",
                  format(v1_n, big.mark = ",")))
  publish_one(v1_out, "v1", "v1.1", v1_n)
} else if (file.exists(v1_csv)) {
  message(sprintf("v1: streaming %s ...", v1_csv))
  v1_schema <- arrow::schema(id       = arrow::string(),
                             neuropil = arrow::string(),
                             region   = arrow::string(),
                             side     = arrow::string())
  v1_ds <- arrow::open_csv_dataset(v1_csv, schema = v1_schema, skip = 1L)
  arrow::write_parquet(v1_ds %>%
                         dplyr::select(id, neuropil, region, side) %>%
                         dplyr::collect(),
                       v1_out)
  v1_n <- arrow::open_dataset(v1_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  publish_one(v1_out, "v1", "v1.1", v1_n)
  rm(v1_ds); gc()
} else {
  message(sprintf("v1 CSV not found at %s — skipping", v1_csv))
}

#####################################################
### v2 — from banc_<ver>_synapses_v2_neuropils    ###
#####################################################
v2_pq <- file.path(banc.connectivity.save.path,
                   sprintf("banc_%s_synapses_v2_neuropils.parquet", banc.version))
v2_out <- file.path(staging, "synapse_neuropil_lookup_v2.parquet")
if (file.exists(v2_out)) {
  v2_n <- arrow::open_dataset(v2_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  message(sprintf("v2: staged file already exists (%s rows) — skipping rebuild",
                  format(v2_n, big.mark = ",")))
  publish_one(v2_out, "v2", "v2.0", v2_n)
} else if (file.exists(v2_pq)) {
  message(sprintf("v2: reading %s ...", v2_pq))
  v2_ds <- arrow::open_dataset(v2_pq)
  arrow::write_parquet(v2_ds %>%
                         dplyr::select(id, neuropil, region, side) %>%
                         dplyr::mutate(id = as.character(id)) %>%
                         dplyr::collect() %>%
                         dplyr::distinct(id, .keep_all = TRUE),
                       v2_out)
  v2_n <- arrow::open_dataset(v2_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  publish_one(v2_out, "v2", "v2.0", v2_n)
  rm(v2_ds); gc()
} else {
  message(sprintf("v2 parquet not found at %s — skipping", v2_pq))
}

#####################################################
### v3 — from banc_<ver>_synapses_v3.parquet      ###
### (rename syn_id -> id; spatial cols are clean) ###
#####################################################
v3_pq <- file.path(banc.synapses.v3.path,
                   sprintf("banc_%s_synapses_v3.parquet", banc.version))
v3_out <- file.path(staging, "synapse_neuropil_lookup_v3.parquet")
if (file.exists(v3_out)) {
  v3_n <- arrow::open_dataset(v3_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  message(sprintf("v3: staged file already exists (%s rows) — skipping rebuild",
                  format(v3_n, big.mark = ",")))
  publish_one(v3_out, "v3", "v3.0", v3_n)
} else if (file.exists(v3_pq)) {
  message(sprintf("v3: reading %s ...", v3_pq))
  v3_ds <- arrow::open_dataset(v3_pq)
  # syn_id is stored as float64 in the v3 source; as.character() on a double
  # uses scientific notation ("1.11e+13"), which won't join against integer
  # syn_ids in the synapse table. Cast double -> int64 -> string in arrow to
  # preserve the canonical integer representation.
  arrow::write_parquet(v3_ds %>%
                         dplyr::select(syn_id, neuropil, region, side) %>%
                         dplyr::rename(id = syn_id) %>%
                         dplyr::mutate(id = arrow::cast(
                           arrow::cast(id, arrow::int64()),
                           arrow::string())) %>%
                         dplyr::collect() %>%
                         dplyr::distinct(id, .keep_all = TRUE),
                       v3_out)
  v3_n <- arrow::open_dataset(v3_out) %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  publish_one(v3_out, "v3", "v3.0", v3_n)
  rm(v3_ds); gc()
} else {
  message(sprintf("v3 parquet not found at %s — skipping", v3_pq))
}

message(sprintf("### banc: lookup publish complete [%s mins] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
