#' banc-al-type-changes — Detect BANC↔FAFB type mismatches for ALLN/ALPN neurons.
#'
#' Reads Codex BANC↔FAFB connectivity alignment scores for the antennal-lobe
#' local + projection neuron cutout, joins to SeaTable to compare against
#' the curated `fafb_cell_type`, emits review CSV + per-row neuroglancer
#' links, then optionally generates PNG review images.
#'
#' @section Reads:
#'   - `data/codex/alln_alpn/banc_fafb_min_syn_5_alln_alpn_..._alignment_scores.csv.gz`
#'   - SeaTable `banc_meta`
#'   - `franken_meta()` (FAFB cell_type lookup)
#'
#' @section Writes:
#'   - `data/codex/al_type_changes.csv` — mismatches with NGL URLs
#'   - `inst/images/al_review/*.png` — per-mismatch comparison images

###########################################################
### Compare FAFB connectivity matches vs seatable types
### for ALLN/ALPN (antennal lobe) neurons
###
### Reads codex connectivity alignment scores (BANC↔FAFB),
### joins to seatable to get fafb_cell_type, compares
### against the FAFB cell_type of the matched neuron.
### Outputs a CSV of mismatches with neuroglancer links,
### then generates PNG review images for each mismatch.
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: AL type change analysis (BANC↔FAFB) ###")

###########################
### Read data           ###
###########################

# Read connectivity alignment scores (no header)
# Columns: banc_id, fafb_id, score, n_synapses
conn <- readr::read_csv(
  "data/codex/alln_alpn/banc_fafb_min_syn_5_alln_alpn_norm_lr_gdt_node_alignment_scores.csv.gz",
  col_names = c("banc_id", "fafb_id", "score", "n_synapses"),
  col_types = readr::cols(banc_id = "c", fafb_id = "c", score = "d", n_synapses = "d"),
  show_col_types = FALSE
)

# Remove dummy/test rows (IDs < 15 digits are synthetic)
conn <- conn %>% dplyr::filter(nchar(banc_id) > 15, nchar(fafb_id) > 15)
message(sprintf("  Loaded %d connectivity matches", nrow(conn)))

# Keep best match per BANC neuron (highest score)
conn <- conn %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::distinct(banc_id, .keep_all = TRUE)
message(sprintf("  %d unique BANC neurons after dedup", nrow(conn)))

# Query BANC seatable
bc <- banctable_query(
  "SELECT _id, root_id, root_626, supervoxel_id, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, hemibrain_cell_type, fafb_match, malecns_cell_type, malecns_match, side, status FROM banc_meta"
) %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::mutate(
    root_626 = as.character(root_626),
    root_id = as.character(root_id),
    supervoxel_id = as.character(supervoxel_id),
    fafb_match = as.character(fafb_match)
  ) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)
message(sprintf("  SeaTable: %d neurons with root_626", nrow(bc)))

# Read FAFB metadata via franken_meta (available on O2 + locally)
message("  Loading FAFB metadata from franken_meta()...")
fm <- franken_meta()
fafb_meta <- fm %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::select(fafb_root_id = fafb_id, fafb_cell_type_meta = cell_type,
                fafb_side = side, fafb_super_class = super_class) %>%
  dplyr::mutate(fafb_root_id = as.character(fafb_root_id))
message(sprintf("  FAFB meta: %d neurons", nrow(fafb_meta)))

###########################
### Join and compare    ###
###########################

# Join connectivity matches to BANC seatable (banc_id = root_626)
df <- conn %>%
  dplyr::inner_join(bc, by = c("banc_id" = "root_626"))
message(sprintf("  %d matches joined to seatable", nrow(df)))

# Join FAFB metadata to get the connectivity-matched FAFB neuron's cell_type
df <- df %>%
  dplyr::left_join(
    fafb_meta %>% dplyr::distinct(fafb_root_id, .keep_all = TRUE),
    by = c("fafb_id" = "fafb_root_id")
  )

# The "connectivity type" is the cell_type of the matched FAFB neuron
df <- df %>%
  dplyr::rename(connectivity_type = fafb_cell_type_meta)

# Determine the "seatable type" — fafb_cell_type from BANC seatable
df <- df %>%
  dplyr::mutate(seatable_type = fafb_cell_type)

# Find mismatches: connectivity_type vs seatable fafb_cell_type
# A mismatch is where both exist and differ, OR connectivity type exists but seatable doesn't
df <- df %>%
  dplyr::mutate(
    has_conn_type = !is.na(connectivity_type) & connectivity_type != "",
    has_st_type = !is.na(seatable_type) & seatable_type != "",
    is_mismatch = dplyr::case_when(
      has_conn_type & has_st_type ~ connectivity_type != seatable_type,
      has_conn_type & !has_st_type ~ TRUE,   # new type from connectivity
      !has_conn_type & has_st_type ~ FALSE,   # connectivity match has no type — not useful
      TRUE ~ FALSE
    )
  )

mismatches <- df %>% dplyr::filter(is_mismatch)
message(sprintf("  Mismatches: %d / %d (%.1f%%)",
                nrow(mismatches), nrow(df), 100 * nrow(mismatches) / nrow(df)))

if (nrow(mismatches) == 0) {
  message("  No mismatches found. Nothing to write.")
  return(invisible())
}

###############################################
### FAFB neuron lookup for neuroglancer     ###
###############################################

# For visualising the seatable type (fafb_cell_type) in neuroglancer:
#   1. If fafb_match exists, use that specific FAFB neuron
#   2. Else, find FAFB neurons with same side + cell_type
#   3. If fafb_cell_type is blank, fall back to hemibrain_cell_type
#      and find any FAFB neuron with that type + same side → mark hemibrain_transfer

# Build FAFB lookup table: cell_type + side → root_ids
fafb_by_type_side <- fafb_meta %>%
  dplyr::filter(!is.na(fafb_cell_type_meta), fafb_cell_type_meta != "") %>%
  dplyr::group_by(fafb_cell_type_meta, fafb_side) %>%
  dplyr::summarise(fafb_ids = list(fafb_root_id), .groups = "drop")

get_fafb_ids_for_seatable_type <- function(row_fafb_match, row_seatable_type,
                                            row_hb_type, row_side) {
  hemibrain_transfer <- FALSE
  type_to_use <- row_seatable_type

  # If fafb_match exists, use it directly
  if (!is.na(row_fafb_match) && row_fafb_match != "") {
    return(list(ids = row_fafb_match, hemibrain_transfer = FALSE))
  }

  # If seatable fafb_cell_type is blank, try hemibrain_cell_type
  if (is.na(type_to_use) || type_to_use == "") {
    if (!is.na(row_hb_type) && row_hb_type != "") {
      type_to_use <- row_hb_type
      hemibrain_transfer <- TRUE
    } else {
      return(list(ids = character(0), hemibrain_transfer = FALSE))
    }
  }

  # Look up FAFB neurons by type + side
  side_to_use <- if (!is.na(row_side) && row_side != "") row_side else NA
  matches <- fafb_by_type_side %>%
    dplyr::filter(fafb_cell_type_meta == type_to_use)
  if (!is.na(side_to_use)) {
    side_matches <- matches %>% dplyr::filter(fafb_side == side_to_use)
    if (nrow(side_matches) > 0) matches <- side_matches
  }
  ids <- unlist(matches$fafb_ids)
  if (is.null(ids)) ids <- character(0)
  list(ids = ids, hemibrain_transfer = hemibrain_transfer)
}

###########################################
### Neuroglancer URLs                   ###
###########################################

message("  Building neuroglancer links...")

# Decode the base FAFB neuroglancer scene
ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4657695008227328"
ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                    return = "text", cache = TRUE)
ngl_base <- fafbseg::ngl_decode_scene(
  fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

# Find layer indices
ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))
message(sprintf("  Scene layers: %s", paste(ngl_ls$name, collapse = ", ")))

# BANC segmentation layer
banc_layer_idx <- match("v626 neurons", ngl_ls$name)
if (is.na(banc_layer_idx)) banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)
if (is.na(banc_layer_idx)) banc_layer_idx <- grep("banc|segmentation", ngl_ls$name, ignore.case = TRUE)[1]

# FAFB layers: look for fafb/flywire layers
fafb_layer_idxs <- grep("fafb|flywire", ngl_ls$name, ignore.case = TRUE)
# connectivity match layer (the specific matched neuron)
fafb_conn_idx <- if (length(fafb_layer_idxs) >= 1) fafb_layer_idxs[1] else NA_integer_
# seatable type layer (fafb_cell_type neurons)
fafb_st_idx <- if (length(fafb_layer_idxs) >= 2) fafb_layer_idxs[2] else fafb_conn_idx

message(sprintf("  BANC layer: '%s' [%d]",
                if (!is.na(banc_layer_idx)) ngl_ls$name[banc_layer_idx] else "NOT FOUND",
                banc_layer_idx))
message(sprintf("  FAFB layers: conn='%s' [%d], seatable='%s' [%d]",
                if (!is.na(fafb_conn_idx)) ngl_ls$name[fafb_conn_idx] else "NOT FOUND", fafb_conn_idx,
                if (!is.na(fafb_st_idx)) ngl_ls$name[fafb_st_idx] else "NOT FOUND", fafb_st_idx))

first_error <- NULL
ngl_urls <- character(nrow(mismatches))
hb_transfers <- logical(nrow(mismatches))

for (i in seq_len(nrow(mismatches))) {
  row <- mismatches[i, ]

  # Look up FAFB IDs for the seatable type
  st_result <- get_fafb_ids_for_seatable_type(
    row$fafb_match, row$seatable_type, row$hemibrain_cell_type, row$side
  )
  hb_transfers[i] <- st_result$hemibrain_transfer

  tryCatch({
    sc <- ngl_base

    # Set BANC neuron (use root_id for modern segmentation, fall back to banc_id=root_626)
    banc_rid <- if (!is.na(row$root_id) && row$root_id != "") row$root_id else row$banc_id
    if (!is.na(banc_layer_idx)) {
      sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
      sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
    }

    # Set connectivity-matched FAFB neuron
    if (!is.na(fafb_conn_idx)) {
      sc[["layers"]][[fafb_conn_idx]][["segments"]] <- as.character(row$fafb_id)
      sc[["layers"]][[fafb_conn_idx]][["hiddenSegments"]] <- NULL
    }

    # Set seatable-type FAFB neurons
    if (!is.na(fafb_st_idx) && fafb_st_idx != fafb_conn_idx && length(st_result$ids) > 0) {
      sc[["layers"]][[fafb_st_idx]][["segments"]] <- as.character(st_result$ids)
      sc[["layers"]][[fafb_st_idx]][["hiddenSegments"]] <- NULL
    }

    ngl_urls[i] <- as.character(sc)
  }, error = function(e) {
    if (is.null(first_error)) first_error <<- conditionMessage(e)
    ngl_urls[i] <<- NA_character_
  })
}
if (!is.null(first_error)) message(sprintf("  First NGL error: %s", first_error))
message(sprintf("  Generated %d/%d neuroglancer links",
                sum(!is.na(ngl_urls)), nrow(mismatches)))

mismatches$hemibrain_transfer <- hb_transfers
mismatches$neuroglancer_url <- ngl_urls

###########################
### Save CSV            ###
###########################

out <- mismatches %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::transmute(
    root_626 = banc_id,
    root_id,
    fafb_id,
    score,
    n_synapses,
    super_class,
    cell_class,
    cell_sub_class,
    side,
    seatable_type,
    connectivity_type,
    hemibrain_cell_type,
    hemibrain_transfer,
    fafb_match,
    neuroglancer_url
  )

out_file <- "data/codex/alln_alpn/banc_fafb_alln_alpn_type_mismatches.csv"
readr::write_csv(out, out_file)
message(sprintf("  Saved %d mismatches to %s", nrow(out), out_file))

# Summary
message("\n  === Summary ===")
message(sprintf("  Total connectivity matches: %d", nrow(df)))
message(sprintf("  Mismatches: %d", nrow(out)))
message(sprintf("  With seatable type: %d", sum(out$seatable_type != "" & !is.na(out$seatable_type))))
message(sprintf("  Using hemibrain_cell_type fallback: %d", sum(out$hemibrain_transfer)))
message(sprintf("  New types (no seatable type): %d",
                sum(is.na(out$seatable_type) | out$seatable_type == "")))

###########################################
### Review images                       ###
###########################################
message("\n=== Generating review images ===")

# Helper: download a BANC mesh and simplify
get_banc_mesh <- function(id, percent = 0.1) {
  message("    Downloading BANC mesh: ", id)
  m <- tryCatch(bancr::banc_read_neuron_meshes(id), error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  Rvcg::vcgQEdecim(m[[1]], percent = percent)
}

# Helper: download a FAFB mesh (via FlyWire) and simplify
# FAFB meshes are already in BANC space when downloaded via pre-saved OBJ files.
# If OBJ not available, download from FlyWire cloudvolume.
get_fafb_mesh <- function(id, percent = 0.1) {
  message("    Downloading FAFB mesh: ", id)
  # Try pre-saved OBJ files first (already in BANC nm space)
  if (exists("banc.nblast.fafb.obj.save.path")) {
    obj_file <- file.path(banc.nblast.fafb.obj.save.path, paste0(id, ".obj"))
    if (file.exists(obj_file)) {
      m <- tryCatch(
        nat::as.neuronlist(readobj::read.obj(obj_file, convert.rgl = TRUE)),
        error = function(e) NULL
      )
      if (!is.null(m) && length(m) > 0) {
        return(Rvcg::vcgQEdecim(m[[1]], percent = percent))
      }
    }
    # Try versioned path
    if (exists("banc.nblast.version")) {
      obj_file_v <- file.path(banc.nblast.fafb.obj.save.path, banc.nblast.version, paste0(id, ".obj"))
      if (file.exists(obj_file_v)) {
        m <- tryCatch(
          nat::as.neuronlist(readobj::read.obj(obj_file_v, convert.rgl = TRUE)),
          error = function(e) NULL
        )
        if (!is.null(m) && length(m) > 0) {
          return(Rvcg::vcgQEdecim(m[[1]], percent = percent))
        }
      }
    }
  }
  # Fall back to cloudvolume download (FlyWire space → need to transform to BANC)
  m <- tryCatch({
    fafbseg::with_segmentation("flywire31",
      fafbseg::read_cloudvolume_meshes(id))
  }, error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  mesh <- m[[1]]
  # Transform FlyWire nm → FAFB14 → JRC2018F → BANC nm
  mesh <- tryCatch({
    mesh_jrc <- nat.templatebrains::xform_brain(mesh / 1e3,
      sample = "FAFB14", reference = "JRC2018F")
    bancr::banc_to_JRC2018F(mesh_jrc, region = "brain",
      method = "tpsreg", banc.units = "nm", inverse = TRUE)
  }, error = function(e) {
    warning("  Transform failed for FAFB ", id, ": ", e$message)
    NULL
  })
  if (is.null(mesh)) return(NULL)
  Rvcg::vcgQEdecim(mesh, percent = percent)
}

# Merge multiple mesh3d objects into one
merge_meshes <- function(mesh_list) {
  mesh_list <- mesh_list[!vapply(mesh_list, is.null, logical(1))]
  if (length(mesh_list) == 0) return(NULL)
  if (length(mesh_list) == 1) return(mesh_list[[1]])
  result <- mesh_list[[1]]
  for (j in seq_along(mesh_list)[-1]) {
    m <- mesh_list[[j]]
    offset <- ncol(result$vb)
    result$vb <- cbind(result$vb, m$vb)
    result$it <- cbind(result$it, m$it + offset)
  }
  result
}

# Output directory
on_o2 <- dir.exists("/n/data1/hms/neurobio/wilson")
if (on_o2) {
  dir.images <- file.path("inst", "images", "al_type_review")
} else {
  dir.images <- file.path(path.expand("~/Downloads"), "al_type_review")
}
dir.create(dir.images, recursive = TRUE, showWarnings = FALSE)
message("  Output directory: ", dir.images)

# Simplify neuropils once (brain only for AL neurons)
banc_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_neuropil.surf), percent = 0.05)
banc_brain_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_brain_neuropil.surf), percent = 0.05)
banc_vnc_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_vnc_neuropil.surf), percent = 0.05)

# Track completed images
completed <- list.files(dir.images, pattern = "\\.png$")
completed_ids <- gsub("root_id_([0-9]+)_.*", "\\1", completed)

n_done <- 0L
n_skip <- 0L
n_fail <- 0L

for (i in seq_len(nrow(out))) {
  row <- out[i, ]
  banc_id <- row$root_626
  current_root_id <- row$root_id
  seatable_type <- row$seatable_type
  connectivity_type <- row$connectivity_type

  # Skip if already done
  if (current_root_id %in% completed_ids || banc_id %in% completed_ids) {
    n_skip <- n_skip + 1L
    next
  }

  message(sprintf("  [%d/%d] banc_id=%s seatable=%s connectivity=%s",
                  i, nrow(out), banc_id,
                  ifelse(is.na(seatable_type), "NA", seatable_type),
                  ifelse(is.na(connectivity_type), "NA", connectivity_type)))

  tryCatch({

    side <- row$side
    if (is.na(side) || side %in% c("", "NA")) side <- "unknown"
    svid <- bc$supervoxel_id[match(banc_id, bc$root_626)]

    # Download BANC neuron mesh (neuron1 = blue)
    banc_mesh <- get_banc_mesh(current_root_id)
    if (is.null(banc_mesh)) {
      warning("  No mesh for BANC neuron: ", current_root_id)
      n_fail <- n_fail + 1L
      next
    }

    # Download seatable-type FAFB neuron mesh (neuron2 = red = old/existing type)
    st_result <- get_fafb_ids_for_seatable_type(
      row$fafb_match, seatable_type, row$hemibrain_cell_type, side
    )
    if (length(st_result$ids) > 0) {
      st_mesh_id <- st_result$ids[1]
      st_mesh <- get_fafb_mesh(st_mesh_id)
    } else {
      st_mesh <- NULL
    }

    # Download connectivity-matched FAFB neuron mesh (neuron3 = green = new match)
    conn_mesh <- get_fafb_mesh(row$fafb_id)

    if (is.null(conn_mesh) && is.null(st_mesh)) {
      warning("  No FAFB meshes available for ", banc_id)
      n_fail <- n_fail + 1L
      next
    }

    # Info labels — show the three current cell type assignments
    banc_cell_type <- bc$cell_type[match(banc_id, bc$root_626)]
    if (is.na(banc_cell_type) || banc_cell_type == "") banc_cell_type <- "untyped"

    neuron1.info <- sprintf("BANC root_id:%s\ncell_type:%s\nside:%s",
                            current_root_id, banc_cell_type, side)
    st_label <- if (st_result$hemibrain_transfer) {
      sprintf("hemibrain_cell_type:%s", row$hemibrain_cell_type)
    } else {
      sprintf("fafb_cell_type:%s", ifelse(is.na(seatable_type), "untyped", seatable_type))
    }
    neuron2.info <- sprintf("SeaTable type (old)\n%s", st_label)
    neuron3.info <- sprintf("Connectivity match (new)\nfafb_id:%s\ncell_type:%s",
                            row$fafb_id,
                            ifelse(is.na(connectivity_type), "untyped", connectivity_type))

    # Construct filename
    st_clean <- gsub("[^A-Za-z0-9._-]", "", ifelse(is.na(seatable_type), "NA", seatable_type))
    conn_clean <- gsub("[^A-Za-z0-9._-]", "", ifelse(is.na(connectivity_type), "NA", connectivity_type))
    filename <- file.path(dir.images,
      sprintf("root_id_%s_supervoxel_id_%s_seatable_%s_connectivity_%s.png",
              current_root_id, nullToNA(svid), st_clean, conn_clean))

    # Plot — brain region (blue=BANC, red=seatable old, green=connectivity new)
    banc_neuron_comparison_plot(neuron1 = banc_mesh,
                                neuron2 = st_mesh,
                                neuron3 = conn_mesh,
                                neuron1.info = neuron1.info,
                                neuron2.info = neuron2.info,
                                neuron3.info = neuron3.info,
                                filename = filename,
                                banc_neuropil = banc_neuropil,
                                banc_brain_neuropil = banc_brain_neuropil,
                                banc_vnc_neuropil = banc_vnc_neuropil,
                                region = "brain")
    n_done <- n_done + 1L
    message(sprintf("    -> saved (%d done so far)", n_done))

  }, error = function(e) {
    warning(sprintf("  Failed for banc_id %s: %s", banc_id, e$message))
    n_fail <<- n_fail + 1L
  })

  # Free mesh memory
  rm(list = intersect(c("banc_mesh", "conn_mesh", "st_mesh"), ls()))
  gc()
}

message(sprintf("### Done: %d images created, %d skipped, %d failed ###", n_done, n_skip, n_fail))

###############################################################
### SeaTable update: reviewed AL type changes                ###
### accept_new == "T": update fafb_match, fafb_cell_type,    ###
### and cell_type from the connectivity-matched FAFB neuron. ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################

reviewed_csv_path <- "data/codex/alln_alpn/banc_fafb_alln_alpn_type_mismatches_reviewed.csv"
if (file.exists(reviewed_csv_path)) {
  message("=== SeaTable update: reviewed AL type changes (accept_new == T) ===")

  # Read reviewed CSV and filter to accepted
  accepted <- readr::read_csv(reviewed_csv_path,
                               col_types = readr::cols(root_626 = "c", root_id = "c",
                                                        fafb_id = "c", .default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(accept_new = trimws(accept_new)) %>%
    dplyr::filter(accept_new == "T")

  message(sprintf("  %d accepted type changes", nrow(accepted)))

  # Join to seatable to get _id and current values
  accepted <- accepted %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, `_id`,
                                           current_cell_type = cell_type,
                                           current_fafb_cell_type = fafb_cell_type,
                                           current_fafb_match = fafb_match,
                                           current_malecns_cell_type = malecns_cell_type),
                     by = "root_626") %>%
    dplyr::filter(!is.na(`_id`), `_id` != "")

  # Look up cell_type of the connectivity-matched FAFB neuron from franken_meta
  accepted <- accepted %>%
    dplyr::left_join(
      fafb_meta %>% dplyr::distinct(fafb_root_id, .keep_all = TRUE) %>%
        dplyr::select(fafb_root_id, new_fafb_cell_type = fafb_cell_type_meta),
      by = c("fafb_id" = "fafb_root_id")
    )

  # Resolve the effective new fafb_cell_type for each row
  accepted <- accepted %>%
    dplyr::mutate(
      effective_fafb_cell_type = dplyr::if_else(!is.na(new_fafb_cell_type),
                                                 new_fafb_cell_type, connectivity_type)
    )

  # --- malecns_cell_type from maleCNS SeaTable (fafb_cell_type → cell_type mapping) ---
  mcns_type_map <- tryCatch({
    banctable_query("SELECT malecns_09_id, cell_type, fafb_cell_type FROM malecns",
                    base = "cns_meta") %>%
      dplyr::filter(!is.na(fafb_cell_type), fafb_cell_type != "",
                    !is.na(cell_type), cell_type != "") %>%
      dplyr::distinct(fafb_cell_type, .keep_all = TRUE) %>%
      dplyr::select(fafb_cell_type, malecns_ct = cell_type)
  }, error = function(e) {
    message("  Could not query maleCNS SeaTable: ", e$message)
    NULL
  })

  if (!is.null(mcns_type_map) && nrow(mcns_type_map) > 0) {
    accepted <- accepted %>%
      dplyr::left_join(mcns_type_map,
                       by = c("effective_fafb_cell_type" = "fafb_cell_type"))
    n_malecns <- sum(!is.na(accepted$malecns_ct))
    message(sprintf("  malecns_cell_type resolved for %d/%d neurons via maleCNS SeaTable",
                    n_malecns, nrow(accepted)))
  } else {
    accepted$malecns_ct <- NA_character_
    message("  maleCNS type mapping unavailable; malecns_cell_type will not be updated")
  }

  # Build push data frame
  push_df <- accepted %>%
    dplyr::transmute(
      `_id`,
      fafb_match = fafb_id,
      fafb_cell_type = effective_fafb_cell_type,
      # Update cell_type if current matches old seatable_type, or was blank
      cell_type = dplyr::case_when(
        !is.na(seatable_type) & !is.na(current_cell_type) &
          current_cell_type == seatable_type ~ effective_fafb_cell_type,
        is.na(current_cell_type) | current_cell_type == "" ~ effective_fafb_cell_type,
        TRUE ~ current_cell_type
      ),
      # Only update malecns_cell_type if currently blank/NA or matches old fafb_cell_type mapping
      malecns_cell_type = dplyr::case_when(
        !is.na(malecns_ct) &
          (is.na(current_malecns_cell_type) | current_malecns_cell_type == "") ~ malecns_ct,
        !is.na(malecns_ct) ~ malecns_ct,
        TRUE ~ current_malecns_cell_type
      )
    ) %>%
    as.data.frame()

  # Replace NAs with empty string for SeaTable text columns
  push_df$fafb_match[is.na(push_df$fafb_match)] <- ""
  push_df$fafb_cell_type[is.na(push_df$fafb_cell_type)] <- ""
  push_df$cell_type[is.na(push_df$cell_type)] <- ""
  push_df$malecns_cell_type[is.na(push_df$malecns_cell_type)] <- ""

  n_ct_update <- sum(push_df$cell_type != "" &
    push_df$cell_type != accepted$current_cell_type, na.rm = TRUE)
  n_malecns_update <- sum(push_df$malecns_cell_type != "" &
    (is.na(accepted$current_malecns_cell_type) |
       accepted$current_malecns_cell_type == "" |
       push_df$malecns_cell_type != accepted$current_malecns_cell_type), na.rm = TRUE)

  message(sprintf("  Pushing update for %d neurons:", nrow(push_df)))
  message(sprintf("    fafb_match: %d non-empty", sum(push_df$fafb_match != "")))
  message(sprintf("    fafb_cell_type: %d updated", sum(push_df$fafb_cell_type != "")))
  message(sprintf("    cell_type: %d updated (where old == seatable_type or was blank)",
                  n_ct_update))
  message(sprintf("    malecns_cell_type: %d updated", n_malecns_update))

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_df,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  SeaTable update complete")

  rm(accepted, push_df); gc()
} else {
  message("  Skipping SeaTable update: reviewed CSV not found")
}

###############################################################
### SeaTable update: non-mismatch fafb_match                ###
### For connectivity matches that are NOT mismatches (not    ###
### in the reviewed CSV), set fafb_match to the best         ###
### connectivity match FAFB ID.                              ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################

message("=== SeaTable update: non-mismatch fafb_match ===")

# Non-mismatches: rows in df where connectivity type agrees (or is absent)
non_mismatches <- df %>%
  dplyr::filter(!is_mismatch) %>%
  dplyr::filter(!is.na(`_id`), `_id` != "")

push_nonmismatch <- non_mismatches %>%
  dplyr::transmute(
    `_id`,
    fafb_match = as.character(fafb_id)
  ) %>%
  dplyr::filter(!is.na(fafb_match), fafb_match != "") %>%
  as.data.frame()

message(sprintf("  %d non-mismatch neurons: setting fafb_match", nrow(push_nonmismatch)))

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_nonmismatch,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  SeaTable update complete")

rm(non_mismatches, push_nonmismatch); gc()

###############################################################
### SeaTable update: malecns_match from alignment scores     ###
### Uses banc_mcns_min_syn_5_alln_alpn alignment scores.     ###
### Only sets malecns_match if the matched maleCNS neuron's  ###
### cell_type agrees with the current malecns_cell_type      ###
### (i.e. the match does not contradict the type assignment). ###
### COMMENTED OUT — uncomment to run the seatable update.     ###
###############################################################

message("=== SeaTable update: malecns_match from BANC-maleCNS alignment ===")

mcns_align_file <- "data/codex/alln_alpn/banc_mcns_min_syn_5_alln_alpn_norm_lr_gdt_node_alignment_scores.csv.gz"
if (!file.exists(mcns_align_file)) {
  message("  Skipping malecns_match update: alignment CSV not found")
} else {

  mcns_align <- readr::read_csv(mcns_align_file,
                                 col_names = c("banc_id", "mcns_id", "score", "n_synapses"),
                                 col_types = readr::cols(banc_id = "c", mcns_id = "c",
                                                          score = "d", n_synapses = "d"),
                                 show_col_types = FALSE)
  # Remove dummy/test rows
  mcns_align <- mcns_align %>%
    dplyr::filter(nchar(banc_id) > 15, nchar(mcns_id) > 15)
  message(sprintf("  Loaded %d BANC-maleCNS alignment matches", nrow(mcns_align)))

  # Keep best match per BANC neuron
  mcns_align <- mcns_align %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::distinct(banc_id, .keep_all = TRUE)
  message(sprintf("  %d unique BANC neurons after dedup", nrow(mcns_align)))

  # Look up cell_type for each maleCNS neuron
  mcns_types <- tryCatch({
    banctable_query("SELECT malecns_09_id, cell_type FROM malecns",
                    base = "cns_meta") %>%
      dplyr::filter(!is.na(malecns_09_id), !is.na(cell_type), cell_type != "") %>%
      dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
      dplyr::distinct(malecns_09_id, .keep_all = TRUE)
  }, error = function(e) {
    message("  Could not query maleCNS SeaTable: ", e$message)
    NULL
  })

  if (!is.null(mcns_types) && nrow(mcns_types) > 0) {
    mcns_validated <- mcns_align %>%
      # Join to BANC SeaTable to get _id, malecns_cell_type, current malecns_match
      dplyr::inner_join(bc %>% dplyr::select(root_626, `_id`, malecns_cell_type,
                                              current_malecns_match = malecns_match),
                        by = c("banc_id" = "root_626")) %>%
      dplyr::filter(!is.na(`_id`), `_id` != "",
                    !is.na(malecns_cell_type), malecns_cell_type != "") %>%
      # Join to maleCNS types to get the matched neuron's cell_type
      dplyr::left_join(mcns_types %>% dplyr::rename(match_cell_type = cell_type),
                       by = c("mcns_id" = "malecns_09_id")) %>%
      # Only keep where matched cell_type agrees with malecns_cell_type
      dplyr::filter(!is.na(match_cell_type), match_cell_type == malecns_cell_type) %>%
      dplyr::distinct(banc_id, .keep_all = TRUE)

    message(sprintf("  maleCNS alignment matches with cell_type agreement: %d/%d",
                    nrow(mcns_validated), nrow(mcns_align)))

    # Only push rows where the match has actually changed
    mcns_push <- mcns_validated %>%
      dplyr::filter(is.na(current_malecns_match) | current_malecns_match == "" |
                      current_malecns_match != mcns_id) %>%
      dplyr::transmute(`_id`, malecns_match = mcns_id) %>%
      as.data.frame()

    message(sprintf("  maleCNS matches to update (changed from current): %d", nrow(mcns_push)))

    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = mcns_push,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  malecns_match push complete")

    rm(mcns_validated, mcns_push); gc()
  } else {
    message("  maleCNS types unavailable; skipping malecns_match update")
  }

  rm(mcns_align); gc()
}

###############################################################
### BANC-maleCNS mismatch file                               ###
### Compare malecns_cell_type from BANC seatable (after the  ###
### planned updates above) against the cell_type of the best ###
### BANC↔maleCNS connectivity match. Output mismatches with  ###
### neuroglancer links for manual review.                    ###
###############################################################

message("\n=== BANC-maleCNS AL mismatch analysis ===")

mcns_align_file_mm <- "data/codex/alln_alpn/banc_mcns_min_syn_5_alln_alpn_norm_lr_gdt_node_alignment_scores.csv.gz"
if (!file.exists(mcns_align_file_mm)) {
  message("  Skipping BANC-maleCNS mismatch: alignment CSV not found")
} else {

  # Read alignment scores
  mcns_conn <- readr::read_csv(mcns_align_file_mm,
                                col_names = c("banc_id", "mcns_id", "score", "n_synapses"),
                                col_types = readr::cols(banc_id = "c", mcns_id = "c",
                                                         score = "d", n_synapses = "d"),
                                show_col_types = FALSE)
  mcns_conn <- mcns_conn %>% dplyr::filter(nchar(banc_id) > 15)
  message(sprintf("  Loaded %d BANC-maleCNS connectivity matches", nrow(mcns_conn)))

  # Keep best match per BANC neuron
  mcns_conn <- mcns_conn %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::distinct(banc_id, .keep_all = TRUE)
  message(sprintf("  %d unique BANC neurons after dedup", nrow(mcns_conn)))

  # Simulate the planned seatable state: apply accept_new == "T" updates
  # to malecns_cell_type so mismatches reflect the post-update state
  bc_planned <- bc
  if (file.exists(reviewed_csv_path) && exists("mcns_type_map") &&
      !is.null(mcns_type_map) && nrow(mcns_type_map) > 0) {
    reviewed <- readr::read_csv(reviewed_csv_path,
                                 col_types = readr::cols(root_626 = "c", fafb_id = "c",
                                                          .default = "c"),
                                 show_col_types = FALSE) %>%
      dplyr::mutate(accept_new = trimws(accept_new)) %>%
      dplyr::filter(accept_new == "T") %>%
      dplyr::select(root_626, fafb_id, connectivity_type)

    # Resolve effective fafb_cell_type for accepted rows
    reviewed <- reviewed %>%
      dplyr::left_join(
        fafb_meta %>% dplyr::distinct(fafb_root_id, .keep_all = TRUE) %>%
          dplyr::select(fafb_root_id, new_fafb_ct = fafb_cell_type_meta),
        by = c("fafb_id" = "fafb_root_id")
      ) %>%
      dplyr::mutate(eff_fafb_ct = dplyr::if_else(!is.na(new_fafb_ct),
                                                   new_fafb_ct, connectivity_type)) %>%
      dplyr::left_join(mcns_type_map, by = c("eff_fafb_ct" = "fafb_cell_type")) %>%
      dplyr::filter(!is.na(malecns_ct)) %>%
      dplyr::select(root_626, planned_malecns_ct = malecns_ct)

    # Overlay planned updates onto bc
    bc_planned <- bc_planned %>%
      dplyr::left_join(reviewed, by = "root_626") %>%
      dplyr::mutate(malecns_cell_type = dplyr::if_else(
        !is.na(planned_malecns_ct), planned_malecns_ct, malecns_cell_type
      )) %>%
      dplyr::select(-planned_malecns_ct)

    message(sprintf("  Applied %d planned malecns_cell_type updates for mismatch detection",
                    nrow(reviewed)))
  }

  # Join connectivity matches to BANC seatable (planned state)
  mcns_df <- mcns_conn %>%
    dplyr::inner_join(bc_planned %>% dplyr::select(root_626, `_id`, root_id, supervoxel_id,
                                                    super_class, cell_class, cell_sub_class,
                                                    cell_type, malecns_cell_type, malecns_match,
                                                    side, status),
                      by = c("banc_id" = "root_626"))
  message(sprintf("  %d matches joined to seatable", nrow(mcns_df)))

  # Look up cell_type of the matched maleCNS neuron
  mcns_types_mm <- tryCatch({
    banctable_query("SELECT malecns_09_id, cell_type FROM malecns",
                    base = "cns_meta") %>%
      dplyr::filter(!is.na(malecns_09_id), !is.na(cell_type), cell_type != "") %>%
      dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
      dplyr::distinct(malecns_09_id, .keep_all = TRUE)
  }, error = function(e) {
    message("  Could not query maleCNS SeaTable: ", e$message)
    NULL
  })

  if (is.null(mcns_types_mm) || nrow(mcns_types_mm) == 0) {
    message("  maleCNS types unavailable; skipping BANC-maleCNS mismatch file")
  } else {
    mcns_df <- mcns_df %>%
      dplyr::left_join(mcns_types_mm %>% dplyr::rename(mcns_connectivity_type = cell_type),
                       by = c("mcns_id" = "malecns_09_id"))

    # Find mismatches: connectivity type vs seatable malecns_cell_type
    mcns_df <- mcns_df %>%
      dplyr::mutate(
        has_conn_type = !is.na(mcns_connectivity_type) & mcns_connectivity_type != "",
        has_st_type = !is.na(malecns_cell_type) & malecns_cell_type != "",
        is_mismatch_mcns = dplyr::case_when(
          has_conn_type & has_st_type ~ mcns_connectivity_type != malecns_cell_type,
          has_conn_type & !has_st_type ~ TRUE,
          !has_conn_type & has_st_type ~ FALSE,
          TRUE ~ FALSE
        )
      )

    mcns_mismatches <- mcns_df %>% dplyr::filter(is_mismatch_mcns)
    message(sprintf("  BANC-maleCNS mismatches: %d / %d (%.1f%%)",
                    nrow(mcns_mismatches), nrow(mcns_df),
                    100 * nrow(mcns_mismatches) / max(nrow(mcns_df), 1)))

    if (nrow(mcns_mismatches) > 0) {

      ###############################################
      ### maleCNS neuroglancer URLs               ###
      ###############################################

      message("  Building maleCNS neuroglancer links...")

      mcns_ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5930652703653888"
      mcns_ngl_url2 <- sub("#!middleauth+", "?", mcns_ngl_url, fixed = TRUE)
      mcns_ngl_parts <- unlist(strsplit(mcns_ngl_url2, "?", fixed = TRUE))
      mcns_ngl_json <- fafbseg::flywire_fetch(mcns_ngl_parts[2],
                                                token = bancr:::banc_token(),
                                                return = "text", cache = TRUE)
      mcns_ngl_base <- fafbseg::ngl_decode_scene(
        fafbseg::ngl_encode_url(mcns_ngl_json, baseurl = mcns_ngl_parts[1]))

      mcns_ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(mcns_ngl_base))
      message(sprintf("  Scene layers: %s", paste(mcns_ngl_ls$name, collapse = ", ")))

      # BANC segmentation layer
      mcns_banc_layer_idx <- match("v626 neurons", mcns_ngl_ls$name)
      if (is.na(mcns_banc_layer_idx))
        mcns_banc_layer_idx <- match("segmentation proofreading", mcns_ngl_ls$name)
      if (is.na(mcns_banc_layer_idx))
        mcns_banc_layer_idx <- grep("banc|segmentation", mcns_ngl_ls$name, ignore.case = TRUE)[1]

      # maleCNS layers
      mcns_layer_idxs <- grep("malecns|male.?cns", mcns_ngl_ls$name, ignore.case = TRUE)
      mcns_conn_layer_idx <- if (length(mcns_layer_idxs) >= 1) mcns_layer_idxs[1] else NA_integer_
      mcns_st_layer_idx <- if (length(mcns_layer_idxs) >= 2) mcns_layer_idxs[2] else mcns_conn_layer_idx

      message(sprintf("  BANC layer: '%s' [%d]",
                      if (!is.na(mcns_banc_layer_idx)) mcns_ngl_ls$name[mcns_banc_layer_idx] else "NOT FOUND",
                      mcns_banc_layer_idx))
      message(sprintf("  maleCNS layers: conn='%s' [%d], seatable='%s' [%d]",
                      if (!is.na(mcns_conn_layer_idx)) mcns_ngl_ls$name[mcns_conn_layer_idx] else "NOT FOUND",
                      mcns_conn_layer_idx,
                      if (!is.na(mcns_st_layer_idx)) mcns_ngl_ls$name[mcns_st_layer_idx] else "NOT FOUND",
                      mcns_st_layer_idx))

      # Build a maleCNS lookup: cell_type + side → root_ids (for seatable type visualisation)
      mcns_by_type_side <- mcns_types_mm %>%
        dplyr::left_join(
          tryCatch({
            banctable_query("SELECT malecns_09_id, side FROM malecns",
                            base = "cns_meta") %>%
              dplyr::mutate(malecns_09_id = as.character(malecns_09_id))
          }, error = function(e) {
            data.frame(malecns_09_id = character(0), side = character(0),
                       stringsAsFactors = FALSE)
          }),
          by = "malecns_09_id"
        ) %>%
        dplyr::filter(!is.na(cell_type), cell_type != "") %>%
        dplyr::group_by(cell_type, side) %>%
        dplyr::summarise(mcns_ids = list(malecns_09_id), .groups = "drop")

      get_mcns_ids_for_type <- function(row_malecns_match, row_malecns_cell_type, row_side) {
        # If malecns_match exists, use it directly
        if (!is.na(row_malecns_match) && row_malecns_match != "") {
          return(row_malecns_match)
        }
        if (is.na(row_malecns_cell_type) || row_malecns_cell_type == "") {
          return(character(0))
        }
        side_to_use <- if (!is.na(row_side) && row_side != "") row_side else NA
        matches <- mcns_by_type_side %>%
          dplyr::filter(cell_type == row_malecns_cell_type)
        if (!is.na(side_to_use)) {
          side_matches <- matches %>% dplyr::filter(side == side_to_use)
          if (nrow(side_matches) > 0) matches <- side_matches
        }
        ids <- unlist(matches$mcns_ids)
        if (is.null(ids)) character(0) else ids
      }

      mcns_first_error <- NULL
      mcns_ngl_urls <- character(nrow(mcns_mismatches))

      for (i in seq_len(nrow(mcns_mismatches))) {
        row <- mcns_mismatches[i, ]

        st_ids <- get_mcns_ids_for_type(row$malecns_match, row$malecns_cell_type, row$side)

        tryCatch({
          sc <- mcns_ngl_base

          # Set BANC neuron
          banc_rid <- if (!is.na(row$root_id) && row$root_id != "") row$root_id else row$banc_id
          if (!is.na(mcns_banc_layer_idx)) {
            sc[["layers"]][[mcns_banc_layer_idx]][["segments"]] <- as.character(banc_rid)
            sc[["layers"]][[mcns_banc_layer_idx]][["hiddenSegments"]] <- NULL
          }

          # Set connectivity-matched maleCNS neuron
          if (!is.na(mcns_conn_layer_idx)) {
            sc[["layers"]][[mcns_conn_layer_idx]][["segments"]] <- as.character(row$mcns_id)
            sc[["layers"]][[mcns_conn_layer_idx]][["hiddenSegments"]] <- NULL
          }

          # Set seatable-type maleCNS neurons
          if (!is.na(mcns_st_layer_idx) && mcns_st_layer_idx != mcns_conn_layer_idx &&
              length(st_ids) > 0) {
            sc[["layers"]][[mcns_st_layer_idx]][["segments"]] <- as.character(st_ids)
            sc[["layers"]][[mcns_st_layer_idx]][["hiddenSegments"]] <- NULL
          }

          mcns_ngl_urls[i] <- as.character(sc)
        }, error = function(e) {
          if (is.null(mcns_first_error)) mcns_first_error <<- conditionMessage(e)
          mcns_ngl_urls[i] <<- NA_character_
        })
      }
      if (!is.null(mcns_first_error))
        message(sprintf("  First NGL error: %s", mcns_first_error))
      message(sprintf("  Generated %d/%d neuroglancer links",
                      sum(!is.na(mcns_ngl_urls)), nrow(mcns_mismatches)))

      mcns_mismatches$neuroglancer_url <- mcns_ngl_urls

      # Save CSV
      mcns_out <- mcns_mismatches %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::transmute(
          root_626 = banc_id,
          root_id,
          mcns_id,
          score,
          n_synapses,
          super_class,
          cell_class,
          cell_sub_class,
          side,
          malecns_cell_type,
          mcns_connectivity_type,
          malecns_match,
          neuroglancer_url
        )

      mcns_out_file <- "data/codex/alln_alpn/banc_mcns_alln_alpn_type_mismatches.csv"
      readr::write_csv(mcns_out, mcns_out_file)
      message(sprintf("  Saved %d maleCNS mismatches to %s", nrow(mcns_out), mcns_out_file))

      # Summary
      message("\n  === BANC-maleCNS mismatch summary ===")
      message(sprintf("  Total connectivity matches: %d", nrow(mcns_df)))
      message(sprintf("  Mismatches: %d", nrow(mcns_out)))
      message(sprintf("  With seatable malecns_cell_type: %d",
                      sum(!is.na(mcns_out$malecns_cell_type) & mcns_out$malecns_cell_type != "")))
      message(sprintf("  New types (no malecns_cell_type): %d",
                      sum(is.na(mcns_out$malecns_cell_type) | mcns_out$malecns_cell_type == "")))

      rm(mcns_mismatches, mcns_out); gc()
    } else {
      message("  No BANC-maleCNS mismatches found. Nothing to write.")
    }

    rm(mcns_df, mcns_types_mm); gc()
  }

  rm(mcns_conn); gc()
}

return(invisible())
})
