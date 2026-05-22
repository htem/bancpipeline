#' banc-extract-synapse-lookups — Extract version-agnostic synapse → neuropil lookup.
#'
#' Synapse `id` is stable across BANC root_id versions; run ONCE before a
#' rebuild to preserve existing neuropil assignments.
#'
#' @section Reads:
#'   - `banc_<ver>_synapses_<src>_neuropils.parquet`
#'
#' @section Writes:
#'   - `synapse_neuropil_lookup.parquet`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.

###########################################################
### Extract version-agnostic synapse lookups
###
### Extracts synapse-level neuropil assignments (id -> neuropil,
### region, side) from existing versioned neuropil parquets
### and saves as a standalone lookup file that can be reused
### across version rebuilds.
###
### The synapse `id` column is stable across segmentation
### versions — only root_id assignments change. So neuropil
### and side assignments computed for one version are valid
### for all versions.
###
### Output:
###   synapse_neuropil_lookup.parquet  (in banc.connectivity.save.path)
###
### Usage: Run ONCE before a rebuild to preserve existing
###        neuropil assignments. Future runs of
###        banc-calculate-neuropil-inclusion.R will use this
###        lookup automatically.
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: extracting version-agnostic synapse lookups ###")
t_start <- Sys.time()

lookup_file <- file.path(banc.connectivity.save.path, "synapse_neuropil_lookup.parquet")

# Find existing neuropil parquet (check versions in reverse order). Tries the
# new explicit "_v2_" name first, falls back to the legacy unsuffixed name for
# any pre-rename files still on disk (e.g. older v746/v626 outputs).
existing <- NULL
for (v in c("888", "850", "821", "746", "626")) {
  for (suffix in c("v2_neuropils", "neuropils")) {
    candidate <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_synapses_%s.parquet", v, suffix))
    if (file.exists(candidate)) {
      existing <- candidate
      message(sprintf("Found existing neuropil data: %s (v%s)", candidate, v))
      break
    }
  }
  if (!is.null(existing)) break
}

# Fallback: CSV from banc-assess-synapses.R (banc.syns.nps.cleaned)
csv_fallback <- file.path(banc.connectivity.save.path, "banc_synapses_to_neuropils_v2.csv")

if (is.null(existing) && !file.exists(csv_fallback)) {
  message("No existing neuropil parquet or CSV found — lookup will be built from scratch by neuropil inclusion script")
  return(invisible())
}

# Read source data: parquet preferred, CSV fallback
if (!is.null(existing)) {
  message("Reading neuropil parquet...")
  nps <- arrow::read_parquet(existing, col_select = c("id", "neuropil", "region", "side"))
} else {
  message(sprintf("Reading neuropil CSV fallback: %s", csv_fallback))
  nps <- readr::read_csv(csv_fallback, col_types = readr::cols(.default = "c"),
                         show_col_types = FALSE) %>%
    dplyr::select(id, neuropil, region, side)
}
nps$id <- as.character(nps$id)
message(sprintf("  Read %s synapse rows", format(nrow(nps), big.mark = ",")))

# Deduplicate by id (shouldn't have dupes, but be safe)
nps <- nps[!duplicated(nps$id), ]

# If lookup already exists, merge (keep new data where overlapping)
if (file.exists(lookup_file)) {
  message("Merging with existing lookup...")
  old_lookup <- arrow::read_parquet(lookup_file)
  old_lookup$id <- as.character(old_lookup$id)

  # New data takes precedence over old
  new_ids <- nps$id
  old_only <- old_lookup[!old_lookup$id %in% new_ids, ]
  nps <- dplyr::bind_rows(nps, old_only)
  message(sprintf("  Merged: %s from new source + %s from old lookup",
                  format(length(new_ids), big.mark = ","),
                  format(nrow(old_only), big.mark = ",")))
}

arrow::write_parquet(nps, lookup_file)
message(sprintf("Saved lookup: %s (%s rows)",
                lookup_file, format(nrow(nps), big.mark = ",")))

message(sprintf("### banc: synapse lookup extraction complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "secs"), 0))))

})
