#' banc-alignment-discrepancies — Flag super_class discrepancies in optic-lobe alignment.
#'
#' Identifies BANC neurons where the alignment-assigned cell_type implies a
#' different `super_class` than the SeaTable curation, and emits per-row
#' neuroglancer links for review.
#'
#' @section Reads:
#'   - `data/optic_lobe/banc_optic_<side>_alignment_results.csv`
#'   - `data/optic_lobe/banc_optic_<side>_meta.csv`
#'   - SeaTable `banc_meta` (current curation)
#'
#' @section Writes:
#'   - `data/optic_lobe/banc_optic_discrepancies<suffix>.csv` (with NGL URLs)
#'
#' @section CLI:
#'   [side1] [side2]      sides to process (default `right left`)
#'   --suffix <name>      optional suffix for the output filename

###########################################################
### Optic Lobe Alignment: Super_class Discrepancy Detection
###
### Identifies BANC neurons where the connectivity-assigned
### cell type implies a different super_class than SeaTable.
### Generates neuroglancer links for review.
###
### Usage: Rscript alignment/banc-alignment-discrepancies.R [side]
###########################################################
source("banc/banc-startup.R")
source("alignment/alignment-data-sources.R")

local({

# Args: [side1] [side2] [--suffix name]
args <- commandArgs(trailingOnly = TRUE)
suffix <- ""
suffix_raw <- ""
if ("--suffix" %in% args) {
  suffix_idx <- which(args == "--suffix")
  suffix_raw <- args[suffix_idx + 1]
  suffix <- paste0("_", suffix_raw)
  args <- args[-c(suffix_idx, suffix_idx + 1)]
}
sides_to_process <- if (length(args) > 0) args else c("right", "left")
data_dir <- "data/optic_lobe"
target_name    <- "fafb"
region_name    <- "optic-lobe"
nblast_version <- "783"
syn_source     <- Sys.getenv("BANC_SYN_SOURCE",
                              unset = banc.synapse.source.default)

# Base neuroglancer scene (same as matching/ scripts)
ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5370752139264000"
ngl_base <- tryCatch({
  ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
  ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
  ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                      return = "text", cache = TRUE)
  fafbseg::ngl_decode_scene(
    fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))
}, error = function(e) {
  message("  Warning: Could not decode NGL base scene: ", e$message)
  NULL
})

# Find BANC segmentation layer and FAFB layers
banc_layer_idx <- NA_integer_
fafb_st_idx <- NA_integer_    # "seatable type" layer (red) — old FAFB match
fafb_conn_idx <- NA_integer_  # "connectivity type" layer (green) — new connectivity match
fafb_nb_idx <- NA_integer_    # "top nblast type" layer (yellow) — best NBLAST match
if (!is.null(ngl_base)) {
  ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))
  banc_layer_idx <- match("v626 neurons", ngl_ls$name)
  if (is.na(banc_layer_idx)) banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)
  if (is.na(banc_layer_idx)) banc_layer_idx <- grep("banc|segmentation", ngl_ls$name, ignore.case = TRUE)[1]
  fafb_st_idx <- grep("seatable.type", ngl_ls$name, ignore.case = TRUE)[1]
  fafb_conn_idx <- grep("connectivity.type", ngl_ls$name, ignore.case = TRUE)[1]
  fafb_nb_idx <- grep("nblast.type", ngl_ls$name, ignore.case = TRUE)[1]
  message(sprintf("  NGL base scene: BANC layer '%s' [%d], FAFB layers: st=%s, conn=%s, nb=%s",
                  if (!is.na(banc_layer_idx)) ngl_ls$name[banc_layer_idx] else "NOT FOUND",
                  banc_layer_idx,
                  if (!is.na(fafb_st_idx)) as.character(fafb_st_idx) else "NA",
                  if (!is.na(fafb_conn_idx)) as.character(fafb_conn_idx) else "NA",
                  if (!is.na(fafb_nb_idx)) as.character(fafb_nb_idx) else "NA"))
}

make_ngl_link <- function(root_id, old_fafb_id = NA, new_fafb_id = NA) {
  if (is.null(ngl_base) || is.na(banc_layer_idx)) return(NA_character_)
  tryCatch({
    sc <- ngl_base
    # BANC neuron (blue)
    sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(root_id)
    sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
    # Old FAFB match from SeaTable (red)
    if (!is.na(fafb_st_idx) && !is.na(old_fafb_id) && old_fafb_id != "") {
      sc[["layers"]][[fafb_st_idx]][["segments"]] <- as.character(old_fafb_id)
      sc[["layers"]][[fafb_st_idx]][["hiddenSegments"]] <- NULL
    }
    # New connectivity-matched FAFB neuron (green)
    if (!is.na(fafb_conn_idx) && !is.na(new_fafb_id) && new_fafb_id != "") {
      sc[["layers"]][[fafb_conn_idx]][["segments"]] <- as.character(new_fafb_id)
      sc[["layers"]][[fafb_conn_idx]][["hiddenSegments"]] <- NULL
    }
    as.character(sc)
  }, error = function(e) NA_character_)
}

for (side in sides_to_process) {
  message(sprintf("\n=== %s side discrepancies ===", side))

  align_stage <- if (nzchar(suffix_raw)) paste0("align-", suffix_raw) else "align"
  results_file <- alignment_path(align_stage, query = "banc", target = target_name,
                                  region = region_name, side = side,
                                  vq = banc.version, vt = nblast_version,
                                  ext = "csv", dir = data_dir)
  meta_file <- alignment_path("prep-banc-meta", query = "banc", target = target_name,
                               region = region_name, side = side,
                               vq = banc.version, syn = syn_source,
                               ext = "csv", dir = data_dir)
  fafb_file <- alignment_path("prep-target-meta", query = "banc", target = target_name,
                               region = region_name, side = side,
                               vq = banc.version, vt = nblast_version,
                               ext = "csv", dir = data_dir)

  if (!all(file.exists(c(results_file, meta_file, fafb_file)))) {
    message("  Skipping: missing files")
    next
  }

  results <- readr::read_csv(results_file,
    col_types = readr::cols(root_888 = "c", best_target_match = "c"), show_col_types = FALSE)
  banc_meta <- readr::read_csv(meta_file,
    col_types = readr::cols(root_id = "c", fafb_match = "c"), show_col_types = FALSE)
  fafb_meta <- readr::read_csv(fafb_file,
    col_types = readr::cols(target_id = "c"), show_col_types = FALSE)

  # Target-dataset type -> super_class majority vote
  target_type_sc <- fafb_meta %>%
    dplyr::filter(!is.na(target_cell_type), !is.na(target_super_class)) %>%
    dplyr::group_by(target_cell_type) %>%
    dplyr::summarise(new_super_class = names(sort(table(target_super_class),
                                                  decreasing = TRUE))[1],
                     .groups = "drop")

  # Join results to meta
  disc <- results %>%
    dplyr::filter(assigned_cell_type != "") %>%
    dplyr::left_join(
      banc_meta %>%
        dplyr::select(root_id, super_class, cell_type, fafb_cell_type, fafb_match) %>%
        dplyr::distinct(root_id, .keep_all = TRUE),
      by = c("root_888" = "root_id")
    ) %>%
    dplyr::left_join(target_type_sc, by = c("assigned_cell_type" = "target_cell_type")) %>%
    # Determine existing type: cell_type or fafb_cell_type (excluding auto:)
    dplyr::mutate(
      existing_type = dplyr::case_when(
        !is.na(cell_type) & cell_type != "" & !grepl("^auto:", cell_type) ~ cell_type,
        !is.na(fafb_cell_type) & fafb_cell_type != "" & !grepl("^auto:", fafb_cell_type) ~ fafb_cell_type,
        TRUE ~ NA_character_
      ),
      has_nblast_match = !is.na(best_target_match) & best_target_match != "",
      type_disagrees = !is.na(existing_type) & existing_type != assigned_cell_type,
      sc_disagrees = !is.na(super_class) & super_class != "" &
                     !is.na(new_super_class) & super_class != new_super_class
    ) %>%
    dplyr::filter(type_disagrees | sc_disagrees)

  n_sc <- sum(disc$sc_disagrees)
  n_type <- sum(disc$type_disagrees & !disc$sc_disagrees)
  message(sprintf("  Found %d discrepancies (%d super_class + %d cell_type only)",
                  nrow(disc), n_sc, n_type))

  # Generate NGL links (with FAFB neurons in separate layers)
  if (!is.null(ngl_base) && nrow(disc) > 0) {
    message("  Generating neuroglancer links...")
    disc$ngl_link <- mapply(make_ngl_link,
      root_id = disc$root_888,
      old_fafb_id = disc$fafb_match,
      new_fafb_id = disc$best_target_match,
      USE.NAMES = FALSE)
    n_links <- sum(!is.na(disc$ngl_link))
    message(sprintf("  Generated %d/%d NGL links", n_links, nrow(disc)))
  } else {
    disc$ngl_link <- NA_character_
  }

  out <- disc %>%
    dplyr::select(
      root_888,
      old_fafb_match = fafb_match,
      new_fafb_match = best_target_match,
      existing_type,
      new_assigned_cell_type = assigned_cell_type,
      old_super_class = super_class,
      new_super_class,
      sc_disagrees,
      type_disagrees,
      confidence,
      has_nblast_match,
      ngl_link
    ) %>%
    # Super_class discrepancies first, then by confidence
    dplyr::arrange(dplyr::desc(sc_disagrees), dplyr::desc(confidence))

  disc_extra <- if (nzchar(suffix_raw)) suffix_raw else NULL
  out_file <- alignment_path("validate-discrepancies", query = "banc", target = target_name,
                              region = region_name, side = side,
                              vq = banc.version, vt = nblast_version,
                              extra = disc_extra,
                              ext = "csv", dir = data_dir)
  readr::write_csv(out, out_file)
  message(sprintf("  Saved: %s", out_file))

  # Summary
  transitions <- out %>%
    dplyr::count(old_super_class, new_super_class, sort = TRUE)
  message("  Top transitions:")
  for (i in seq_len(min(5, nrow(transitions)))) {
    r <- transitions[i, ]
    message(sprintf("    %s -> %s: %d", r$old_super_class, r$new_super_class, r$n))
  }
}

})
