#' banc-nblast-search — Find BANC neurons matching a target cell_type via compiled NBLAST.
#'
#' Manual: set `search_cell_type` + `search_dataset` at the top; returns
#' the top-N BANC NBLAST hits.
#'
#' @section Reads:
#'   - per-dataset compiled NBLAST feathers

###########################################################
### NBLAST Search: find BANC neurons matching a cell type
### from another dataset
###
### Uses compiled NBLAST feather files to find the best
### BANC matches for a given cell type in a given dataset.
###
### Configure the search at the top, then run.
###########################################################
source("banc/banc-startup.R")

local({

#######################
### SEARCH CONFIG   ###
#######################

search_cell_type <- "SLP444"
search_dataset   <- "FAFB"       # FAFB, MANC, hemibrain, maleCNS, FANC
top_n            <- 10

#######################
### NBLAST FILES    ###
#######################

# Map dataset names to GCS feather filenames and key columns
nblast_files <- list(
  FAFB     = "banc_fafb_783_nblast.feather",
  MANC     = "banc_manc_v1.2.1_nblast.feather",
  hemibrain = "banc_hemibrain_v1.2.1_nblast.feather",
  maleCNS  = "banc_malecns_v0.9_nblast.feather",
  FANC     = "banc_fanc_1116_nblast.feather"
)

gcs_nblast_base <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast"
local_cache_dir <- file.path(banc.meta.save.path)

if (!search_dataset %in% names(nblast_files)) {
  stop("Unknown dataset: ", search_dataset,
       ". Choose from: ", paste(names(nblast_files), collapse = ", "))
}

#######################
### LOAD NBLAST     ###
#######################

nblast_filename <- nblast_files[[search_dataset]]
local_path <- file.path(local_cache_dir, nblast_filename)

# Try local first, then GCS
if (file.exists(local_path)) {
  message(sprintf("Reading NBLAST scores from local: %s", local_path))
  nblast <- arrow::read_feather(local_path)
} else {
  gcs_path <- file.path(gcs_nblast_base, nblast_filename)
  cache_path <- file.path("data", "cache", nblast_filename)
  if (file.exists(cache_path)) {
    message(sprintf("Reading NBLAST scores from cache: %s", cache_path))
    nblast <- arrow::read_feather(cache_path)
  } else {
    message(sprintf("Downloading NBLAST scores from GCS: %s", gcs_path))
    dir.create("data/cache", showWarnings = FALSE, recursive = TRUE)
    dl <- system2("gsutil", c("cp", gcs_path, cache_path),
                   stdout = FALSE, stderr = FALSE)
    if (dl != 0 || !file.exists(cache_path))
      stop("Failed to download ", gcs_path)
    nblast <- arrow::read_feather(cache_path)
    message(sprintf("  Cached to: %s", cache_path))
  }
}
message(sprintf("  %s NBLAST scores loaded (%s rows)",
                search_dataset, format(nrow(nblast), big.mark = ",")))

#######################
### FIND MATCHES    ###
#######################

# The feather has match_cell_type = cell type of the matched neuron (from the other dataset)
# We want BANC neurons whose best match is to a neuron of the search cell type

# Filter to rows matching the search cell type
matches <- nblast %>%
  dplyr::filter(match_cell_type == search_cell_type,
                score > 0)

if (nrow(matches) == 0) {
  message(sprintf("No NBLAST matches found for cell_type '%s' in %s.",
                  search_cell_type, search_dataset))
  message("Available cell types (sample):")
  cts <- sort(unique(nblast$match_cell_type))
  # Show types that partially match the search
  partial <- grep(search_cell_type, cts, ignore.case = TRUE, value = TRUE)
  if (length(partial) > 0) {
    message("  Partial matches: ", paste(head(partial, 20), collapse = ", "))
  } else {
    message("  ", paste(head(cts, 20), collapse = ", "), " ...")
  }
  return(invisible())
}

message(sprintf("  Found %d NBLAST hits for '%s'", nrow(matches), search_cell_type))

# Use root_888 as the BANC neuron ID (fall back to pt_root_id)
banc_id_col <- if ("root_888" %in% colnames(matches) &&
                    any(!is.na(matches$root_888))) "root_888" else "pt_root_id"

# For each BANC neuron, take the MAX score across all neurons of the search type
# (i.e., if there are 5 SLP444 neurons in FAFB, take the best match to any of them)
top_hits <- matches %>%
  dplyr::mutate(banc_id = as.character(.data[[banc_id_col]])) %>%
  dplyr::group_by(banc_id) %>%
  dplyr::summarise(
    max_score = max(score, na.rm = TRUE),
    best_match_id = match_id[which.max(score)],
    n_matches = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(max_score)) %>%
  dplyr::slice_head(n = top_n)

#######################
### ENRICH WITH META
#######################

# Get BANC meta for the hit neurons
banc.meta <- tryCatch(banctable_query_cached(), error = function(e) NULL)
if (!is.null(banc.meta)) {
  top_hits <- top_hits %>%
    dplyr::left_join(
      banc.meta %>%
        dplyr::mutate(root_id = as.character(root_id)) %>%
        dplyr::distinct(root_id, .keep_all = TRUE) %>%
        dplyr::select(root_id, dplyr::any_of(c(
          "cell_type", "super_class", "region", "side",
          "proofread", "roughly_proofread"))),
      by = c("banc_id" = "root_id")
    )
}

#######################
### REPORT          ###
#######################

message(sprintf("\n=== Top %d BANC matches for %s '%s' ===\n",
                min(top_n, nrow(top_hits)), search_dataset, search_cell_type))

for (i in seq_len(nrow(top_hits))) {
  h <- top_hits[i, ]
  ct_str <- if (!is.null(h$cell_type) && !is.na(h$cell_type) && h$cell_type != "")
    h$cell_type else "(untyped)"
  sc_str <- if (!is.null(h$super_class) && !is.na(h$super_class)) h$super_class else ""
  pr_str <- if (!is.null(h$proofread) && isTRUE(h$proofread)) "proofread" else "not proofread"

  message(sprintf("  %2d. %s  score=%.3f  type=%-20s  super_class=%-25s  %s  (matched %d %s neurons, best: %s)",
                  i, h$banc_id, h$max_score, ct_str, sc_str, pr_str,
                  h$n_matches, search_cell_type, h$best_match_id))
}

# Return the results invisibly for interactive use
assign("nblast_search_results", top_hits, envir = .GlobalEnv)
message(sprintf("\nResults saved to nblast_search_results (%d rows)", nrow(top_hits)))

})
