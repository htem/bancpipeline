#' banc-al-review-images — Generate comparison PNGs for AL connectivity mismatches.
#'
#' For each row in `al_connectivity_mismatches.csv`, render a three-mesh
#' comparison PNG: BANC neuron (blue), curated SeaTable-type FAFB neurons
#' (red), connectivity-matched FAFB neuron (green). Uses `fafb_match`
#' where available; falls back to all FAFB neurons of `fafb_cell_type` on
#' the same side.
#'
#' @section Reads:
#'   - `data/codex/al_connectivity_mismatches.csv`
#'   - BANC meshes via `bancr::banc_read_neuron_meshes()`
#'   - FAFB meshes (pre-saved OBJ where available; `fafbseg` fallback)
#'
#' @section Writes:
#'   - `inst/images/al_review/*.png` — per-mismatch comparison images

###########################################################
### AL connectivity match review images
###
### For each entry in al_connectivity_mismatches.csv,
### create a comparison PNG showing:
###   - BANC neuron (blue)
###   - SeaTable-type FAFB neuron(s) (red)
###   - Connectivity-matched FAFB neuron (green)
###
### Uses fafb_match where available; falls back to
### all FAFB neurons of fafb_cell_type + same side.
###########################################################
source("banc/banc-startup.R")

local({

options(warn = 1)
message("### banc: AL connectivity match review images ###")

###########################
### Helper functions    ###
###########################

# Download a BANC mesh and simplify
get_banc_mesh <- function(id, percent = 0.1) {
  message("    Downloading BANC mesh: ", id)
  m <- tryCatch(bancr::banc_read_neuron_meshes(id), error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  Rvcg::vcgQEdecim(m[[1]], percent = percent)
}

# Download a FAFB mesh and simplify
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
  # Fall back to cloudvolume download (FlyWire space → transform to BANC)
  m <- tryCatch({
    fafbseg::with_segmentation("flywire31",
      fafbseg::read_cloudvolume_meshes(id))
  }, error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  mesh <- m[[1]]
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

###########################
### Read data           ###
###########################

mismatch_csv <- "data/codex/al_connectivity_mismatches.csv"
if (!file.exists(mismatch_csv)) {
  message("  Mismatch CSV not found: ", mismatch_csv)
  message("  Run banc-al-connectivity-matches.R first.")
  return(invisible())
}

out <- readr::read_csv(mismatch_csv,
                        col_types = readr::cols(root_626 = "c", root_id = "c",
                                                fafb_id = "c", fafb_match = "c",
                                                .default = "c"),
                        show_col_types = FALSE)
message(sprintf("  Loaded %d mismatches from %s", nrow(out), mismatch_csv))

if (nrow(out) == 0) {
  message("  Nothing to do.")
  return(invisible())
}

# Load FAFB metadata for type-based lookups when fafb_match is absent
message("  Loading FAFB metadata from franken_meta()...")
fm <- franken_meta()
fafb_meta <- fm %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::select(fafb_root_id = fafb_id, fafb_cell_type_meta = cell_type,
                fafb_side = side) %>%
  dplyr::mutate(fafb_root_id = as.character(fafb_root_id))

fafb_by_type_side <- fafb_meta %>%
  dplyr::filter(!is.na(fafb_cell_type_meta), fafb_cell_type_meta != "") %>%
  dplyr::group_by(fafb_cell_type_meta, fafb_side) %>%
  dplyr::summarise(fafb_ids = list(fafb_root_id), .groups = "drop")

get_fafb_ids_for_seatable_type <- function(row_fafb_match, row_seatable_type,
                                            row_hb_type, row_side) {
  hemibrain_transfer <- FALSE
  type_to_use <- row_seatable_type

  if (!is.na(row_fafb_match) && row_fafb_match != "") {
    return(list(ids = row_fafb_match, hemibrain_transfer = FALSE))
  }

  if (is.na(type_to_use) || type_to_use == "") {
    if (!is.na(row_hb_type) && row_hb_type != "") {
      type_to_use <- row_hb_type
      hemibrain_transfer <- TRUE
    } else {
      return(list(ids = character(0), hemibrain_transfer = FALSE))
    }
  }

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

###########################
### Setup output        ###
###########################

dir.images <- file.path(path.expand("~/Desktop"), "matching", "al_review")
dir.create(dir.images, recursive = TRUE, showWarnings = FALSE)
message("  Output directory: ", dir.images)

# Simplify neuropils once (brain only for AL neurons)
banc_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_neuropil.surf), percent = 0.05)
banc_brain_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_brain_neuropil.surf), percent = 0.05)
banc_vnc_neuropil <- Rvcg::vcgQEdecim(rgl::as.mesh3d(banc_vnc_neuropil.surf), percent = 0.05)

# Track completed images
completed <- list.files(dir.images, pattern = "\\.png$")
completed_ids <- gsub("root_id_([0-9]+)_.*", "\\1", completed)

###########################
### Generate images     ###
###########################

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
      st_meshes <- lapply(st_result$ids, function(id) tryCatch(get_fafb_mesh(id), error = function(e) NULL))
      st_mesh <- merge_meshes(st_meshes)
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

    # Info labels
    neuron1.info <- sprintf("BANC root_id:%s\ncell_type:%s\nside:%s",
                            current_root_id,
                            ifelse(is.na(seatable_type), "untyped", seatable_type),
                            side)
    st_label <- if (isTRUE(as.logical(row$hemibrain_transfer))) {
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
      sprintf("root_id_%s_seatable_%s_connectivity_%s.png",
              current_root_id, st_clean, conn_clean))

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
  rm(list = intersect(c("banc_mesh", "conn_mesh", "st_mesh", "st_meshes"), ls()))
  gc()
}

message(sprintf("### Done: %d images created, %d skipped, %d failed ###", n_done, n_skip, n_fail))

})
