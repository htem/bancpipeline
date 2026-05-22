#' banc-vnc-type-review-images — Generate review PNGs for VNC type changes.
#'
#' For each row in `vnc_type_changes_reveiwed.csv` with `accept_new == NA`,
#' render a three-mesh comparison PNG: BANC neuron (blue), best MANC NBLAST
#' match for the old SeaTable type (red), best MANC NBLAST match for the
#' new effective type (green). MANC meshes transformed
#' MANC → JRCVNC2018F → BANC.
#'
#' @section Reads:
#'   - `data/codex/vnc_type_changes_reveiwed.csv`
#'   - `banc_manc_v1.2.1_nblast.feather` (GCS-cached)
#'   - MANC meshes via `malevnc::read_manc_meshes()`
#'   - BANC meshes via `bancr::banc_read_neuron_meshes()`
#'
#' @section Writes:
#'   - `inst/images/vnc_review/*.png` — per-row comparison images

###########################################################
### VNC type change review images
###
### For each entry in vnc_type_changes_reveiwed.csv with
### NA accept_new, create a comparison PNG showing:
###   - BANC neuron (blue)
###   - Best MANC NBLAST match for seatable_type (red)
###   - Best MANC NBLAST match for effective_type (green)
###
### Uses NBLAST feather (from GCS) to pick single best match
### per type. Falls back to franken_meta if no NBLAST hit.
### If a type is NA, no neurons are shown for that type.
###
### Meshes are downloaded fresh and freed after each PNG.
###   MANC: read_manc_meshes -> nm->um -> MANC->JRCVNC2018F -> JRCVNC2018F->BANC (tpsreg)
###   BANC: banc_read_neuron_meshes (already in BANC nm)
###########################################################
source("banc/banc-startup.R")

local({

options(warn = 1)  # print warnings immediately for debugging
message("### banc: VNC type review images ###")

# Ensure MANC dataset is configured for read_manc_meshes
tryCatch(malevnc:::choose_malevnc_dataset("MANC"), error = function(e) NULL)

# Ensure Saalfeld lab registrations available for xform_brain
tryCatch(nat.jrcbrains::register_saalfeldlab_registrations(), error = function(e) NULL)

###########################
### Helper functions   ###
###########################

# Download a BANC mesh and simplify (no local caching)
get_banc_mesh <- function(id, percent = 0.1) {
  message("    Downloading BANC mesh: ", id)
  m <- tryCatch(bancr::banc_read_neuron_meshes(id), error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  Rvcg::vcgQEdecim(m[[1]], percent = percent)
}

# Download a MANC mesh and transform to BANC space (no local caching)
get_manc_mesh <- function(id, percent = 0.1) {
  message("    Downloading + transforming MANC mesh: ", id)
  manc.mesh <- tryCatch(read_manc_meshes(id), error = function(e) NULL)
  if (is.null(manc.mesh) || length(manc.mesh) == 0) return(NULL)

  manc.mesh.jrcvnc2018f <- tryCatch(
    nat.templatebrains::xform_brain(manc.mesh / 1e3, reference = "JRCVNC2018F", sample = "MANC"),
    error = function(e) NULL
  )
  if (is.null(manc.mesh.jrcvnc2018f)) return(NULL)

  mesh3d.banc <- tryCatch(
    bancr::banc_to_JRC2018F(manc.mesh.jrcvnc2018f, region = "vnc",
                             method = "tpsreg", banc.units = "nm", inverse = TRUE),
    error = function(e) NULL
  )
  if (is.null(mesh3d.banc)) return(NULL)

  m <- if (is.list(mesh3d.banc) && !inherits(mesh3d.banc, "mesh3d")) mesh3d.banc[[1]] else mesh3d.banc
  Rvcg::vcgQEdecim(m, percent = percent)
}

# Pick best MANC ID for a BANC neuron + cell_type from NBLAST data
# Tries both banc_id (root_626) and current_root_id
# Uses which() for NA-safe subsetting (== on NAs gives NA, not FALSE)
best_nblast_match <- function(banc_id, current_root_id, cell_type, nb_data) {
  if (is.null(nb_data) || is.na(cell_type) || cell_type == "") return(NULL)
  for (rid in unique(c(banc_id, current_root_id))) {
    idx <- which(nb_data$pt_root_id == rid & nb_data$match_cell_type == cell_type)
    if (length(idx) > 0) {
      hits <- nb_data[idx, , drop = FALSE]
      best <- which.max(hits$score)
      if (length(best) > 0) {
        id <- hits$match_id[best]
        if (!is.na(id) && id != "") return(id)
      }
    }
  }
  NULL
}

# Fallback: pick first same-side MANC ID from franken_meta
# Uses which() for NA-safe subsetting
fallback_manc_id <- function(cell_type, side, fm) {
  if (is.null(fm) || is.na(cell_type) || cell_type == "") return(NULL)
  idx <- which(fm$cell_type == cell_type)
  if (length(idx) == 0) return(NULL)
  fm_sub <- fm[idx, , drop = FALSE]
  if (side != "unknown") {
    side_idx <- which(fm_sub$side == side)
    if (length(side_idx) > 0) fm_sub <- fm_sub[side_idx, , drop = FALSE]
  }
  id <- fm_sub$manc_id[1]
  if (is.na(id) || id == "") NULL else id
}

# Type-wide fallback: get ALL same-side MANC IDs for a cell type
# from NBLAST data and/or franken_meta (when per-neuron NBLAST lookup fails)
# Uses which() for NA-safe subsetting and filters out NA IDs
all_type_manc_ids <- function(cell_type, side, nb_data, fm) {
  if (is.na(cell_type) || cell_type == "") return(character(0))
  ids <- character(0)
  if (!is.null(nb_data)) {
    idx <- which(nb_data$match_cell_type == cell_type)
    ids <- unique(nb_data$match_id[idx])
  }
  if (!is.null(fm)) {
    idx <- which(fm$cell_type == cell_type)
    ids <- unique(c(ids, fm$manc_id[idx]))
  }
  ids <- ids[!is.na(ids) & ids != ""]
  if (length(ids) == 0) return(character(0))
  # Filter to same side using franken_meta side column
  if (side != "unknown" && !is.null(fm)) {
    side_idx <- which(fm$manc_id %in% ids & fm$side == side)
    if (length(side_idx) > 0) return(unique(fm$manc_id[side_idx]))
  }
  ids
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

# Type changes needing review (accept_new is NA)
reviewed <- readr::read_csv("data/codex/vnc_type_changes_reveiwed.csv",
                             col_types = readr::cols(banc_id = "c", root_id = "c",
                                                      .default = "c"),
                             show_col_types = FALSE) %>%
  dplyr::filter(is.na(accept_new)) %>%
  dplyr::select(banc_id, root_id, seatable_type, effective_type)
message(sprintf("  %d entries with NA accept_new to process", nrow(reviewed)))

if (nrow(reviewed) == 0) {
  message("  Nothing to do.")
  return(invisible())
}

# GCS cache directory (shared across downloads)
nblast_cache_dir <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache_dir, showWarnings = FALSE, recursive = TRUE)

# Get banc.meta for side and supervoxel_id lookup
# Try SeaTable first, fall back to GCS feather
banc.meta <- tryCatch({
  message("  Querying SeaTable for banc.meta...")
  banctable_query()
}, error = function(e) {
  warning("  banctable_query failed: ", e$message, " — falling back to GCS")
  gcs_meta_path <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_746_meta.feather"
  meta_local <- file.path(nblast_cache_dir, basename(gcs_meta_path))
  if (!file.exists(meta_local)) {
    message("  Downloading banc meta from GCS...")
    system2("gsutil", c("cp", gcs_meta_path, meta_local), stdout = TRUE, stderr = TRUE)
  }
  bm <- arrow::read_feather(meta_local)
  # Normalize column names: banc_746_id -> root_id
  if ("banc_746_id" %in% names(bm) && !"root_id" %in% names(bm)) {
    bm$root_id <- as.character(bm$banc_746_id)
  }
  if (!"root_626" %in% names(bm)) bm$root_626 <- NA_character_
  message(sprintf("  GCS banc meta loaded: %d rows", nrow(bm)))
  bm
})

# Load NBLAST scores from GCS (primary source for best-match selection)
message("  Loading BANC-MANC NBLAST scores from GCS...")
gcs_nblast_path <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_manc_v1.2.1_nblast.feather"
nblast_local <- file.path(nblast_cache_dir, basename(gcs_nblast_path))

if (!file.exists(nblast_local)) {
  message("  Downloading NBLAST feather from GCS...")
  system2("gsutil", c("cp", gcs_nblast_path, nblast_local),
          stdout = TRUE, stderr = TRUE)
}

nb_data <- tryCatch({
  nb <- arrow::read_feather(nblast_local)
  nb$pt_root_id <- as.character(nb$pt_root_id)
  nb$match_id <- as.character(nb$match_id)
  nb$match_cell_type <- as.character(nb$match_cell_type)
  message(sprintf("  NBLAST data loaded: %d rows", nrow(nb)))
  nb
}, error = function(e) {
  warning("  Failed to load NBLAST data: ", e$message)
  NULL
})

# Load franken_meta as optional fallback (for entries without NBLAST hits)
fm <- tryCatch({
  f <- franken_meta()
  f <- f[!is.na(f$manc_id) & f$manc_id != "", ]
  f$manc_id <- as.character(f$manc_id)
  message(sprintf("  franken_meta loaded: %d rows", nrow(f)))
  f
}, error = function(e) {
  warning("  franken_meta unavailable: ", e$message, " — NBLAST-only mode")
  NULL
})

if (is.null(nb_data) && is.null(fm)) {
  message("  No NBLAST data and no franken_meta available. Cannot proceed.")
  return(invisible())
}

# Output directory: use Downloads when not on O2
on_o2 <- dir.exists("/n/data1/hms/neurobio/wilson")
if (on_o2) {
  dir.images <- file.path("inst", "images", "vnc_type_review")
} else {
  dir.images <- file.path(path.expand("~/Downloads"), "cell_type_review")
}
dir.create(dir.images, recursive = TRUE, showWarnings = FALSE)
message("  Output directory: ", dir.images)

# Simplify neuropils once
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

for (i in seq_len(nrow(reviewed))) {
  row <- reviewed[i, ]
  banc_id <- row$banc_id
  seatable_type <- row$seatable_type
  effective_type <- row$effective_type

  # Skip if both types are NA — nothing to compare
  if ((is.na(seatable_type) || seatable_type == "") &&
      (is.na(effective_type) || effective_type == "")) {
    n_skip <- n_skip + 1L
    next
  }

  # Skip if already done
  if (banc_id %in% completed_ids) {
    n_skip <- n_skip + 1L
    next
  }

  message(sprintf("  [%d/%d] banc_id=%s seatable=%s effective=%s",
                  i, nrow(reviewed), banc_id,
                  ifelse(is.na(seatable_type), "NA", seatable_type),
                  ifelse(is.na(effective_type), "NA", effective_type)))

  tryCatch({

    # Look up side and supervoxel_id from banc.meta
    meta_row <- banc.meta[match(banc_id, banc.meta$root_626), ]
    if (is.na(meta_row$root_id[1])) {
      meta_row <- banc.meta[match(banc_id, banc.meta$root_id), ]
    }
    if (is.na(meta_row$root_id[1])) {
      warning("  banc_id not found in banc.meta: ", banc_id)
      n_fail <- n_fail + 1L
      next
    }
    current_root_id <- meta_row$root_id[1]
    svid <- meta_row$supervoxel_id[1]
    side <- meta_row$side[1]
    if (is.na(side) || side %in% c("", "NA")) side <- "unknown"

    # Pick single best MANC neuron for seatable_type (skip if NA)
    st_chosen <- NULL
    st_all_ids <- character(0)
    if (!is.na(seatable_type) && seatable_type != "") {
      st_chosen <- best_nblast_match(banc_id, current_root_id, seatable_type, nb_data)
      if (is.null(st_chosen)) st_chosen <- fallback_manc_id(seatable_type, side, fm)
      if (is.null(st_chosen)) {
        # Type-wide fallback: show ALL same-side neurons of this type
        st_all_ids <- all_type_manc_ids(seatable_type, side, nb_data, fm)
        message(sprintf("    No per-neuron NBLAST for seatable_type '%s'; showing %d type-mates",
                        seatable_type, length(st_all_ids)))
      }
    }

    # Pick single best MANC neuron for effective_type (skip if NA)
    eff_chosen <- NULL
    eff_all_ids <- character(0)
    if (!is.na(effective_type) && effective_type != "") {
      eff_chosen <- best_nblast_match(banc_id, current_root_id, effective_type, nb_data)
      if (is.null(eff_chosen)) eff_chosen <- fallback_manc_id(effective_type, side, fm)
      if (is.null(eff_chosen)) {
        eff_all_ids <- all_type_manc_ids(effective_type, side, nb_data, fm)
        message(sprintf("    No per-neuron NBLAST for effective_type '%s'; showing %d type-mates",
                        effective_type, length(eff_all_ids)))
      }
    }

    # Download BANC neuron mesh
    banc_mesh <- get_banc_mesh(current_root_id)
    if (is.null(banc_mesh)) {
      warning("  No mesh for BANC neuron: ", current_root_id)
      n_fail <- n_fail + 1L
      next
    }

    # Download MANC mesh(es) for seatable_type
    if (!is.null(st_chosen)) {
      st_mesh <- get_manc_mesh(st_chosen)
    } else if (length(st_all_ids) > 0) {
      st_meshes <- lapply(st_all_ids, function(id) tryCatch(get_manc_mesh(id), error = function(e) NULL))
      st_mesh <- merge_meshes(st_meshes)
    } else {
      st_mesh <- NULL
    }

    # Download MANC mesh(es) for effective_type
    if (!is.null(eff_chosen)) {
      eff_mesh <- get_manc_mesh(eff_chosen)
    } else if (length(eff_all_ids) > 0) {
      eff_meshes <- lapply(eff_all_ids, function(id) tryCatch(get_manc_mesh(id), error = function(e) NULL))
      eff_mesh <- merge_meshes(eff_meshes)
    } else {
      eff_mesh <- NULL
    }

    if (is.null(st_mesh) && is.null(eff_mesh)) {
      warning("  No MANC meshes for either type: ", seatable_type, " / ", effective_type)
      n_fail <- n_fail + 1L
      next
    }

    # Info labels
    st_label <- ifelse(is.na(seatable_type), "NA", seatable_type)
    eff_label <- ifelse(is.na(effective_type), "NA", effective_type)

    neuron1.info <- sprintf("root_id:%s\nside:%s", current_root_id, side)
    if (!is.null(st_chosen)) {
      neuron2.info <- sprintf("seatable_type:%s\nMANC:%s", st_label, st_chosen)
    } else if (length(st_all_ids) > 0) {
      neuron2.info <- sprintf("seatable_type:%s\nMANC: %d neurons (whole type, %s)",
                               st_label, length(st_all_ids), side)
    } else {
      neuron2.info <- ""
    }
    if (!is.null(eff_chosen)) {
      neuron3.info <- sprintf("effective_type:%s\nMANC:%s", eff_label, eff_chosen)
    } else if (length(eff_all_ids) > 0) {
      neuron3.info <- sprintf("effective_type:%s\nMANC: %d neurons (whole type, %s)",
                               eff_label, length(eff_all_ids), side)
    } else {
      neuron3.info <- ""
    }

    # Construct filename
    st_clean <- gsub("[^A-Za-z0-9._-]", "", st_label)
    eff_clean <- gsub("[^A-Za-z0-9._-]", "", eff_label)
    filename <- file.path(dir.images,
      sprintf("root_id_%s_supervoxel_id_%s_seatable_%s_effective_%s.png",
              current_root_id, nullToNA(svid), st_clean, eff_clean))

    # Plot
    banc_neuron_comparison_plot(neuron1 = banc_mesh,
                                neuron2 = st_mesh,
                                neuron3 = eff_mesh,
                                neuron1.info = neuron1.info,
                                neuron2.info = neuron2.info,
                                neuron3.info = neuron3.info,
                                filename = filename,
                                banc_neuropil = banc_neuropil,
                                banc_brain_neuropil = banc_brain_neuropil,
                                banc_vnc_neuropil = banc_vnc_neuropil,
                                region = "vnc")
    n_done <- n_done + 1L
    message(sprintf("    -> saved (%d done so far)", n_done))

  }, error = function(e) {
    warning(sprintf("  Failed for banc_id %s: %s", banc_id, e$message))
    n_fail <<- n_fail + 1L
  })

  # Free mesh memory after each iteration
  rm(list = intersect(c("banc_mesh", "st_mesh", "eff_mesh", "st_meshes", "eff_meshes"), ls()))
  gc()
}

message(sprintf("### Done: %d created, %d skipped, %d failed ###", n_done, n_skip, n_fail))

})
