#' banc-dimorphic-density — Density plots of sexually dimorphic presynapses (combined VNC + brain).
#'
#' VNC: JRCVNC2018U space for BANC, maleCNS, MANC. Brain: JRC2018U for
#' BANC, FAFB, maleCNS. Reads the `sexually_dimorphic` column from each
#' dataset (set by `banc-vnc-sexually-dimorphic.R`) and produces faceted
#' 3×5 (VNC) / 3×6 (brain) density plots using `stat_density_2d_filled`
#' with the magma colour scale.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - `franken_meta()` (MANC + FAFB)
#'   - maleCNS SeaTable
#'   - per-dataset synapse parquets (transformed into JRCVNC2018U / JRC2018U)
#'
#' @section Writes:
#'   - VNC: `inst/images/dimorphic/vnc_{dorsal,side}_density.pdf`
#'   - Brain: `inst/images/dimorphic/brain_{dorsal,front}_density.pdf`

###########################################################
### Density plots of sexually dimorphic presynapses
### Combined VNC + Brain script
###
### VNC: JRCVNC2018U space for BANC, maleCNS, MANC
### Brain: JRC2018U space for BANC, FAFB, maleCNS
###
### Dimorphism: Reads sexually_dimorphic column
### directly from SeaTable (set by banc-vnc-sexually-dimorphic.R).
###   BANC: banc_meta.sexually_dimorphic
###   MANC: franken_meta.sexually_dimorphic (dataset=MANC)
###   FAFB: franken_meta.sexually_dimorphic (dataset=FAFB)
###   maleCNS: malecns.sexually_dimorphic (base=cns_meta)
###
### VNC output: 3x5 faceted density plots
###   Rows = dataset (BANC, maleCNS, MANC)
###   Columns = super_class (ascending, descending, sensory, ventral_nerve_cord_intrinsic, effector)
###   Two views: dorsal ("vnc") and lateral ("vnc_side")
###
### Brain output: 3x6 faceted density plots
###   Rows = dataset (BANC, FAFB, maleCNS)
###   Columns = super_class (ascending, central_brain_intrinsic, descending, effector, sensory, visual)
###   Two views: dorsal ("brain") and anterior ("brain_front")
###
### Template: JRCVNC2018U / JRC2018U mesh in light grey
### Density: stat_density_2d_filled with magma colour scale
###
### Coordinate transforms:
###   BANC nm -> banc_to_JRC2018F(vnc/brain) -> JRCVNC2018F/JRC2018F -> JRCVNC2018U/JRC2018U
###   MANC um -> MANC -> JRCVNC2018F -> JRCVNC2018U (on O2: via BANC nm)
###   FAFB nm -> FAFB14 -> JRC2018F -> JRC2018U
###   maleCNS 8nm voxels -> nm -> JRCFIB2022M -> BANC nm -> JRCVNC2018U/JRC2018U
###
### Usage: Rscript banc/utilities/banc-dimorphic-density.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: sexually dimorphic synapse density plots (VNC + brain) ###")
t_start <- Sys.time()

library(ggplot2)
library(dplyr)
library(nat.flybrains)
library(nat.templatebrains)
library(nat.ggplot)

# Register JRC template brain registrations (needed for xform_brain)
nat.jrcbrains::register_saalfeldlab_registrations()

# Detect environment: O2 cluster has this path
on_o2 <- dir.exists("/n/data1/hms/neurobio/wilson")
banc.connectome.save.path <- '/n/data1/hms/neurobio/wilson/connectomes/'
  
# GCS bucket for synapse data (requires gsutil authentication).
# Files are cached locally on first download; subsequent runs reuse cache.
gcs_bucket <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
local_cache <- file.path(tempdir(), "dimorphic_density_cache")
dir.create(local_cache, showWarnings = FALSE, recursive = TRUE)

###########################################################
### Shared helper functions
###########################################################

# Helper: download from GCS if not already cached, return local path
gcs_cache <- function(gcs_path, cache_dir = local_cache) {
  local_file <- file.path(cache_dir, basename(gcs_path))
  if (!file.exists(local_file)) {
    message(sprintf("  Caching %s ...", basename(gcs_path)))
    result <- system2("gsutil", c("cp", gcs_path, local_file),
                      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
      unlink(local_file)  # remove partial download
      stop(sprintf("gsutil cp failed for %s:\n%s", basename(gcs_path),
                   paste(result, collapse = "\n")))
    }
    if (!file.exists(local_file)) {
      stop(sprintf("gsutil cp did not create %s", local_file))
    }
  } else {
    message(sprintf("  Using cached %s", basename(local_file)))
  }
  local_file
}

# Robust wrapper for arrow::open_dataset that handles HPC filesystems
# (Lustre, GPFS) where arrow's URI scheme detection fails with
# "Unrecognized filesystem type in URI: file:///_".
# Falls back to reading the full file as an Arrow Table.
open_parquet <- function(path) {
  tryCatch(
    arrow::open_dataset(path, format = "parquet"),
    error = function(e) {
      if (grepl("filesystem", e$message, ignore.case = TRUE)) {
        message("  open_dataset failed on filesystem; using read_parquet fallback")
        arrow::read_parquet(path, as_data_frame = FALSE)
      } else {
        stop(e)
      }
    }
  )
}

# Project 3D XYZ to 2D using a rotation matrix.
# Returns data.frame with X, Y columns (Z is depth, used for rendering order).
project_points <- function(xyz, rotation_matrix) {
  projected <- t(rotation_matrix[1:3, 1:3] %*% t(as.matrix(xyz)))
  data.frame(X = projected[, 1], Y = projected[, 2], Z = projected[, 3])
}

# Fix super_class based on cell_type naming conventions:
#   AN* = ascending, DN* = descending, SN* = sensory, SA[letter] = sensory ascending
# This prevents ascending/descending/sensory neurons from leaking into
# ventral_nerve_cord_intrinsic.
fix_super_class <- function(df, ct_col = "cell_type") {
  if (!ct_col %in% colnames(df)) return(df)
  ct <- df[[ct_col]]
  fixed <- dplyr::case_when(
    grepl("^AN", ct) ~ "ascending",
    grepl("^DN", ct) ~ "descending",
    grepl("^SN", ct) ~ "sensory",
    grepl("^SA[A-Za-z]", ct) ~ "ascending",
    TRUE ~ df$super_class
  )
  n_changed <- sum(fixed != df$super_class, na.rm = TRUE)
  if (n_changed > 0) {
    changes <- df[which(fixed != df$super_class), c(ct_col, "super_class")]
    changes$new_super_class <- fixed[which(fixed != df$super_class)]
    message(sprintf("  Fixed %d super_class assignments by cell_type prefix:", n_changed))
    for (i in seq_len(min(nrow(changes), 20))) {
      message(sprintf("    %s: %s -> %s", changes[[ct_col]][i],
                       changes$super_class[i], changes$new_super_class[i]))
    }
    if (nrow(changes) > 20) message(sprintf("    ... and %d more", nrow(changes) - 20))
  }
  df$super_class <- fixed
  df
}

# Remap super_class for VNC faceting: anything containing "sensory" -> "sensory",
# motor/visceral_circulatory -> "effector", then filter to super_classes.
remap_super_class <- function(sc) {
  super_classes <- c("ascending", "descending", "sensory",
                     "ventral_nerve_cord_intrinsic", "effector")
  dplyr::case_when(
    grepl("sensory", sc, ignore.case = TRUE) ~ "sensory",
    grepl("motor|visceral_circulatory", sc, ignore.case = TRUE) ~ "effector",
    sc %in% super_classes ~ sc,
    is.na(sc) | sc == "" ~ "ventral_nerve_cord_intrinsic",
    TRUE ~ sc
  )
}

# Remap super_class for brain faceting:
#   sensory (any super_class containing "sensory")
#   visual (any super_class containing "optic" or "visual")
#   central_brain_intrinsic (any super_class containing "central")
#   ascending, descending (exact)
#   effector (motor or visceral_circulatory)
remap_brain_super_class <- function(sc) {
  dplyr::case_when(
    grepl("sensory", sc, ignore.case = TRUE) ~ "sensory",
    grepl("optic|visual", sc, ignore.case = TRUE) ~ "visual",
    grepl("central", sc, ignore.case = TRUE) ~ "central_brain_intrinsic",
    grepl("motor|visceral_circulatory", sc, ignore.case = TRUE) ~ "effector",
    sc == "ascending" ~ "ascending",
    sc == "descending" ~ "descending",
    TRUE ~ NA_character_
  )
}

# Add dimorphism group: "dimorphic" or "sex-specific"
add_dimorphism_group <- function(df) {
  df$dimorphism_group <- dplyr::if_else(
    df$dimorphism == "dimorphic", "dimorphic", "sex-specific")
  df
}

###########################################################
### Rotation matrices (all 4 views)
###########################################################

# VNC dorsal view: 270-degree clockwise rotation so VNC long axis is vertical.
# Maps: display_X = -old_Y, display_Y = old_X, depth = old_Z
jrc_vnc_mat <- matrix(c(
  0,-1, 0, 0,
  1, 0, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
), 4, 4, byrow = TRUE)

# VNC lateral view: Z->display plane, 270-degree clockwise rotation
# Maps: display_X = -old_Y, display_Y = old_Z, depth = old_X
jrc_vnc_side_mat <- matrix(c(
  0,-1, 0, 0,
  0, 0, 1, 0,
  1, 0, 0, 0,
  0, 0, 0, 1
), 4, 4, byrow = TRUE)

# Brain dorsal view (top-down): X = left-right, Y = anterior-posterior
jrc_brain_dor_mat <- matrix(c(
  1, 0, 0, 0,
  0,-1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
), 4, 4, byrow = TRUE)

# Brain anterior view (front-on): X = left-right, Y = dorsal-ventral (dorsal up)
jrc_brain_ant_mat <- matrix(c(
  1, 0, 0, 0,
  0, 0,-1, 0,
  0, 1, 0, 0,
  0, 0, 0, 1
), 4, 4, byrow = TRUE)

###########################################################
### Shared density plot function
###########################################################

# Core density plot function.
# For each facet (dataset x super_class), renders:
#   1. Template mesh as light grey polygon background
#   2. Filled 2D density contours (magma colour scale)
#   3. Template outline as thin border
# normalize_by: NULL (ndensity per panel), or "dataset" (subsample to smallest dataset)
make_density_plot <- function(df, rotation_matrix, template_2d, normalize_by = NULL, neuron_counts = NULL) {

  # Replicate the template outline for every facet panel
  facets <- expand.grid(
    super_class = levels(df$super_class),
    dataset = levels(df$dataset),
    stringsAsFactors = FALSE
  )
  template_all <- do.call(rbind, lapply(seq_len(nrow(facets)), function(i) {
    d <- template_2d
    d$super_class <- factor(facets$super_class[i], levels = levels(df$super_class))
    d$dataset <- factor(facets$dataset[i], levels = levels(df$dataset))
    d$group <- d$group + (i - 1) * max(template_2d$group, na.rm = TRUE)
    d
  }))

  # Project synapse coordinates to 2D
  syn_proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  syn_2d <- data.frame(px = syn_proj$X, py = syn_proj$Y,
                        super_class = df$super_class,
                        dataset = df$dataset)

  if (!is.null(normalize_by) && normalize_by == "dataset") {
    # Subsample each dataset to the size of the smallest, so density
    # magnitudes are directly comparable across datasets.
    set.seed(42)
    ds_counts <- table(syn_2d$dataset[, drop = TRUE])
    min_n <- min(ds_counts)
    syn_2d <- do.call(rbind, lapply(names(ds_counts), function(ds) {
      sub <- syn_2d[syn_2d$dataset == ds, ]
      if (nrow(sub) > min_n) sub[sample(nrow(sub), min_n), ] else sub
    }))
  }

  p <- ggplot() +
    geom_polygon(data = template_all,
                 aes(x = X, y = Y, group = group),
                 fill = "grey80", colour = NA, alpha = 0.1) +
    stat_density_2d_filled(data = syn_2d,
                            aes(x = px, y = py),
                            h = c(5, 5), bins = 20, n = 200,
                            contour_var = "ndensity",
                            breaks = seq(0.001, 1, length.out = 20),
                            alpha = 0.8) +
    scale_fill_viridis_d(option = "magma", name = "Normalized\ndensity",
                          guide = guide_legend(override.aes = list(alpha = 1)))

  # Add neuron count annotations below each panel's template mesh
  if (!is.null(neuron_counts)) {
    x_mid <- mean(range(template_2d$X, na.rm = TRUE))
    y_bottom <- min(template_2d$Y, na.rm = TRUE)
    count_df <- neuron_counts %>%
      dplyr::filter(super_class %in% levels(df$super_class),
                    dataset %in% levels(df$dataset)) %>%
      dplyr::mutate(
        super_class = factor(super_class, levels = levels(df$super_class)),
        dataset = factor(dataset, levels = levels(df$dataset)),
        label = paste0("n=", n_neurons),
        x = x_mid,
        y = y_bottom
      )
    p <- p + geom_text(data = count_df,
                       aes(x = x, y = y, label = label),
                       hjust = 0.5, vjust = 1.5,
                       size = 3.5, colour = "grey30",
                       inherit.aes = FALSE)
  }

  p + facet_grid(dataset ~ super_class) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.spacing = unit(0.5, "cm"),
          legend.position = "right")
}

###########################################################
### Density subtraction plot helper
###########################################################

# Compute BANC - other dataset density subtraction on a 2D grid.
# Each dataset's density is normalized by its total synapse count before subtraction.
# Returns a ggplot with diverging red (BANC higher) / blue (other higher) scale.
#
# df: data.frame with X, Y, Z, super_class, dataset columns (all dimorphism merged)
# banc_name: name of BANC dataset level
# other_name: name of comparison dataset level
# rotation_matrix, template_2d: as for make_density_plot
# grid_n: number of bins per axis for the 2D histogram
# zlim: symmetric colour limit; if NULL, computed from data
# Returns: list(plot = ggplot, zlim = numeric) so caller can unify scale
make_subtraction_plot <- function(df, banc_name, other_name,
                                  rotation_matrix, template_2d,
                                  super_classes, grid_n = 100, zlim = NULL) {

  # Project to 2D
  syn_proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  df$px <- syn_proj$X
  df$py <- syn_proj$Y

  # Shared grid extent
  x_range <- range(template_2d$X, na.rm = TRUE)
  y_range <- range(template_2d$Y, na.rm = TRUE)
  x_breaks <- seq(x_range[1], x_range[2], length.out = grid_n + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = grid_n + 1)
  x_mids <- (x_breaks[-1] + x_breaks[-(grid_n + 1)]) / 2
  y_mids <- (y_breaks[-1] + y_breaks[-(grid_n + 1)]) / 2

  # Build per-dataset, per-super_class normalized density grids
  results <- list()
  for (sc in super_classes) {
    for (ds in c(banc_name, other_name)) {
      sub <- df[df$dataset == ds & df$super_class == sc, ]
      n_total <- sum(df$dataset == ds)  # total synapses in this dataset (all super_classes)
      if (n_total == 0 || nrow(sub) == 0) {
        grid_vals <- matrix(0, nrow = grid_n, ncol = grid_n)
      } else {
        h <- hist2d_fast(sub$px, sub$py, x_breaks, y_breaks)
        grid_vals <- h / n_total  # normalize by dataset total
      }
      results[[paste(sc, ds, sep = "|")]] <- grid_vals
    }
  }

  # Compute subtraction and assemble into long-form data.frame
  diff_df <- do.call(rbind, lapply(super_classes, function(sc) {
    banc_grid <- results[[paste(sc, banc_name, sep = "|")]]
    other_grid <- results[[paste(sc, other_name, sep = "|")]]
    diff_grid <- banc_grid - other_grid
    expand.grid(x = x_mids, y = y_mids) %>%
      dplyr::mutate(diff = as.vector(diff_grid),
                    super_class = sc)
  }))
  diff_df$super_class <- factor(diff_df$super_class, levels = super_classes)

  # Determine symmetric colour limits
  if (is.null(zlim)) {
    zlim <- max(abs(diff_df$diff), na.rm = TRUE)
  }

  # Template outline for each facet
  facets <- data.frame(super_class = super_classes, stringsAsFactors = FALSE)
  template_all <- do.call(rbind, lapply(seq_len(nrow(facets)), function(i) {
    d <- template_2d
    d$super_class <- factor(facets$super_class[i], levels = super_classes)
    d$group <- d$group + (i - 1) * max(template_2d$group, na.rm = TRUE)
    d
  }))

  # Drop near-zero bins for cleaner plot
  diff_df_plot <- diff_df[abs(diff_df$diff) > zlim * 0.01, ]

  p <- ggplot() +
    geom_polygon(data = template_all,
                 aes(x = X, y = Y, group = group),
                 fill = "grey90", colour = "grey70", linewidth = 0.3, alpha = 0.3) +
    geom_tile(data = diff_df_plot,
              aes(x = x, y = y, fill = diff), alpha = 0.85) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-zlim, zlim),
                         name = paste0(banc_name, " \u2212 ", other_name)) +
    facet_wrap(~ super_class, nrow = 1) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.spacing = unit(0.5, "cm"),
          legend.position = "right")

  list(plot = p, zlim = zlim)
}

# Fast 2D histogram: bin px, py into grid defined by x_breaks, y_breaks
# Returns matrix of counts (nrow = length(x_breaks)-1, ncol = length(y_breaks)-1)
hist2d_fast <- function(px, py, x_breaks, y_breaks) {
  xbin <- findInterval(px, x_breaks, rightmost.closed = TRUE)
  ybin <- findInterval(py, y_breaks, rightmost.closed = TRUE)
  nx <- length(x_breaks) - 1
  ny <- length(y_breaks) - 1
  # Drop out-of-range
  valid <- xbin >= 1 & xbin <= nx & ybin >= 1 & ybin <= ny
  tab <- table(factor(xbin[valid], levels = seq_len(nx)),
               factor(ybin[valid], levels = seq_len(ny)))
  as.matrix(tab)
}

###########################################################
### Load dimorphism from SeaTable (shared across VNC + brain)
###########################################################

message("\n=== Loading dimorphism from SeaTable ===")

# 1. BANC (banc_meta) — sexually_dimorphic column
gcs_banc_fallback <- FALSE
bc_base <- tryCatch({
  bq <- banctable_query_cached(paste0(
    "SELECT root_id, root_626, super_class, cell_type, sexually_dimorphic ",
    "FROM banc_meta")) %>%
    dplyr::mutate(root_626 = as.character(root_626)) %>%
    dplyr::filter(sexually_dimorphic %in% c("dimorphic", "female-specific")) %>%
    dplyr::rename(dimorphism = sexually_dimorphic) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  if (nrow(bq) == 0) stop("No dimorphic BANC neurons found")
  bq <- fix_super_class(bq, "cell_type")
  bq$cell_type <- NULL
  bq
}, error = function(e) {
  warning("SeaTable unavailable: ", e$message, " -- using GCS banc_746_meta")
  gcs_banc_fallback <<- TRUE
  gcs_meta_file <- gcs_cache(paste0(gcs_bucket, "/banc/banc_746_meta.feather"))
  bm <- arrow::read_feather(gcs_meta_file)
  bm$root_id <- as.character(bm$banc_746_id)
  bm$root_626 <- as.character(bm$banc_746_id)
  if (!"sexually_dimorphic" %in% names(bm)) bm$sexually_dimorphic <- NA_character_
  bm %>%
    dplyr::filter(sexually_dimorphic %in% c("dimorphic", "female-specific")) %>%
    dplyr::rename(dimorphism = sexually_dimorphic) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
})
message(sprintf("  BANC: %d dimorphic neurons", nrow(bc_base)))

# 2. MANC (franken_meta, dataset=MANC) — sexually_dimorphic column
fm <- tryCatch({
  res <- banctable_query_cached(paste0(
    "SELECT manc_id, cell_type, super_class, cell_class, cell_sub_class, ",
    "sexually_dimorphic, dataset FROM franken_meta"),
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  if (!"manc_id" %in% colnames(res)) {
    id_col <- intersect(c("manc_121_id", "bodyid"), colnames(res))[1]
    if (!is.na(id_col)) res$manc_id <- as.character(res[[id_col]])
  }
  if (!"dataset" %in% colnames(res)) res$dataset <- "MANC"
  if (!"sexually_dimorphic" %in% colnames(res)) res$sexually_dimorphic <- NA_character_
  res %>%
    dplyr::filter(grepl("MANC", dataset, ignore.case = TRUE), !is.na(manc_id)) %>%
    dplyr::mutate(manc_id = as.character(manc_id)) %>%
    dplyr::distinct(manc_id, .keep_all = TRUE)
}, error = function(e) {
  warning("Could not query franken_meta for MANC: ", e$message)
  NULL
})
if (!is.null(fm)) message(sprintf("  MANC: %d neurons (%d dimorphic)",
  nrow(fm), sum(fm$sexually_dimorphic %in% c("dimorphic", "male-specific"), na.rm = TRUE)))

# 3. FAFB (franken_meta, dataset=FAFB/FlyWire) — sexually_dimorphic column
fafb.meta <- tryCatch({
  res <- banctable_query_cached(paste0(
    "SELECT fafb_id, cell_type, super_class, cell_class, cell_sub_class, ",
    "sexually_dimorphic, dataset FROM franken_meta"),
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  res %>%
    dplyr::filter(grepl("FAFB|FlyWire", dataset, ignore.case = TRUE),
                  !is.na(fafb_id)) %>%
    dplyr::mutate(fafb_id = as.character(fafb_id)) %>%
    dplyr::distinct(fafb_id, .keep_all = TRUE)
}, error = function(e) {
  warning("Could not query franken_meta for FAFB: ", e$message)
  NULL
})
if (!is.null(fafb.meta)) message(sprintf("  FAFB: %d neurons (%d dimorphic)",
  nrow(fafb.meta), sum(fafb.meta$sexually_dimorphic %in% c("dimorphic", "female-specific"), na.rm = TRUE)))

# 4. maleCNS (malecns table, cns_meta base) — sexually_dimorphic column
malecns.meta <- tryCatch({
  res <- banctable_query_cached(paste0(
    "SELECT malecns_09_id, cell_type, manc_cell_type, super_class, ",
    "cell_class, cell_sub_class, sexually_dimorphic FROM malecns"),
    base = "cns_meta")
  if (is.null(res) || nrow(res) == 0) stop("No data returned")
  res %>%
    dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
    dplyr::rename(dimorphism = sexually_dimorphic)
}, error = function(e) {
  warning("Could not query malecns: ", e$message)
  # Fallback to GCS feather
  tryCatch({
    if (on_o2) {
      arrow::read_feather(file.path('/n/data1/hms/neurobio/wilson/malecns',
                                     "malecns_09_meta.feather"))
    } else {
      malecns_meta_file <- gcs_cache(paste0(gcs_bucket, "/malecns/malecns_09_meta.feather"))
      arrow::read_feather(malecns_meta_file)
    }
  }, error = function(e2) {
    warning("Could not load maleCNS meta: ", e2$message)
    NULL
  })
})
if (!is.null(malecns.meta)) message(sprintf("  maleCNS: %d neurons (%d dimorphic)",
  nrow(malecns.meta), sum(malecns.meta$sexually_dimorphic %in% c("dimorphic", "male-specific"), na.rm = TRUE)))

###########################################################
###########################################################
###                                                     ###
###                  VNC SECTION                         ###
###                                                     ###
###########################################################
###########################################################

message("\n\n========================================")
message("=== VNC DENSITY PLOTS ==================")
message("========================================")

###########################
### VNC transform funcs ###
###########################

# Transform BANC coordinates (nm) to JRCVNC2018U (um).
# Two-step: BANC -> JRC2018F (VNC region) -> JRCVNC2018U.
banc_to_jrcvnc2018u <- function(xyz) {
  xyz_jrc <- bancr::banc_to_JRC2018F(xyz, region = "vnc",
                                       method = "tpsreg",
                                       banc.units = "nm",
                                       inverse = FALSE)
  nat.templatebrains::xform_brain(xyz_jrc,
                                   sample = "JRCVNC2018F",
                                   reference = "JRCVNC2018U")
}

# Transform MANC coordinates (um) to JRCVNC2018U (um).
# Two-step: MANC -> JRCVNC2018F -> JRCVNC2018U.
manc_to_jrcvnc2018u <- function(xyz) {
  xyz_jrc <- nat.templatebrains::xform_brain(xyz,
                                              sample = "MANC",
                                              reference = "JRCVNC2018F")
  nat.templatebrains::xform_brain(xyz_jrc,
                                   sample = "JRCVNC2018F",
                                   reference = "JRCVNC2018U")
}

# Shared helper: transform maleCNS nm coords to BANC nm via Python navis.
# Processes in chunks to avoid memory issues with large point sets.
malecns_to_banc_nm <- function(xyz_nm, chunk_size = 500000L) {
  library(reticulate)
  navis <- import("navis", convert = FALSE)
  flybrains <- import("flybrains")
  np <- import("numpy", convert = FALSE)

  n_pts <- nrow(xyz_nm)
  n_chunks <- ceiling(n_pts / chunk_size)
  xyz_banc <- matrix(NA_real_, nrow = n_pts, ncol = 3)

  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1) * chunk_size + 1):min(i * chunk_size, n_pts)
    message(sprintf("  maleCNS navis chunk %d/%d (%d points)", i, n_chunks, length(idx)))
    coords_py <- np$array(r_to_py(xyz_nm[idx, , drop = FALSE]))
    coords_banc_py <- navis$xform_brain(coords_py,
                                         source = "JRCFIB2022M",
                                         target = "BANC")
    xyz_banc[idx, ] <- py_to_r(coords_banc_py)
  }

  xyz_banc
}

# Transform maleCNS coordinates (nm) to JRCVNC2018U via BANC.
malecns_to_jrcvnc2018u <- function(xyz_nm, chunk_size = 500000L) {
  xyz_banc <- malecns_to_banc_nm(xyz_nm, chunk_size = chunk_size)
  banc_to_jrcvnc2018u(xyz_banc)
}

###########################
### VNC template surface ##
###########################

jrc_surf <- nat.flybrains::JRCVNC2018U.surf

###########################
### VNC settings        ###
###########################

vnc_plot.dir <- file.path("images", "dimorphic_density")
dir.create(vnc_plot.dir, recursive = TRUE, showWarnings = FALSE)

super_classes <- c("ascending", "descending", "sensory", "ventral_nerve_cord_intrinsic", "effector")

###########################
### BANC (VNC)          ###
###########################

message("\n=== BANC (VNC) ===")

# Apply VNC-specific remap to shared bc_base
bc_vnc <- bc_base %>%
  dplyr::mutate(super_class = remap_super_class(super_class)) %>%
  dplyr::filter(super_class %in% super_classes) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

message(sprintf("  %d dimorphic BANC neurons for VNC", nrow(bc_vnc)))

# Track neuron counts per super_class and dimorphism for facet labels
vnc_banc_neuron_counts <- bc_vnc %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "BANC")

# Read BANC synapse data with lazy loading
banc_vnc_ids <- na.omit(unique(as.character(bc_vnc$root_626)))
if (on_o2) {
  banc_syn_file <- Sys.glob(file.path(banc.connectome.save.path,
                                       "banc/banc_*_synapses.parquet"))[1]
  banc_vnc.pre <- open_parquet(banc_syn_file) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::filter(pre %in% banc_vnc_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_vnc %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("pre" = "root_626"))
} else {
  banc_syn_file <- gcs_cache(paste0(gcs_bucket, "/banc/banc_746_synapses.parquet"))
  banc_vnc.pre <- open_parquet(banc_syn_file) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::filter(pre %in% banc_vnc_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_vnc %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("pre" = "root_626"))
}
gc()
message(sprintf("  %d presynapses from dimorphic neurons", nrow(banc_vnc.pre)))

# Transform presynapse XYZ from BANC nm to JRCVNC2018U um
xyz_banc <- nat::xyzmatrix(banc_vnc.pre)
xyz_u <- banc_to_jrcvnc2018u(xyz_banc)

# Keep only points that fall inside the JRCVNC2018U template
inside <- nat::pointsinside(xyz_u, jrc_surf)
inside[is.na(inside)] <- FALSE

banc_vnc.density <- data.frame(X = xyz_u[inside, 1],
                                Y = xyz_u[inside, 2],
                                Z = xyz_u[inside, 3],
                                super_class = banc_vnc.pre$super_class[inside],
                                dimorphism = banc_vnc.pre$dimorphism[inside],
                                dataset = "BANC")
message(sprintf("  %d presynapses within JRCVNC2018U", sum(inside)))
rm(banc_vnc.pre, xyz_banc, xyz_u, inside); gc()

# BANC VNC postsynapses
message("  Loading BANC postsynapses...")
if (on_o2) {
  banc_vnc.post <- open_parquet(banc_syn_file) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::filter(post %in% banc_vnc_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_vnc %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("post" = "root_626"))
} else {
  banc_vnc.post <- open_parquet(banc_syn_file) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::filter(post %in% banc_vnc_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_vnc %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("post" = "root_626"))
}
gc()
message(sprintf("  %d postsynapses to dimorphic neurons", nrow(banc_vnc.post)))

xyz_banc_post <- nat::xyzmatrix(banc_vnc.post)
xyz_u_post <- banc_to_jrcvnc2018u(xyz_banc_post)
inside_post <- nat::pointsinside(xyz_u_post, jrc_surf)
inside_post[is.na(inside_post)] <- FALSE

banc_vnc.density.post <- data.frame(X = xyz_u_post[inside_post, 1],
                                     Y = xyz_u_post[inside_post, 2],
                                     Z = xyz_u_post[inside_post, 3],
                                     super_class = banc_vnc.post$super_class[inside_post],
                                     dimorphism = banc_vnc.post$dimorphism[inside_post],
                                     dataset = "BANC")
message(sprintf("  %d postsynapses within JRCVNC2018U", sum(inside_post)))
rm(banc_vnc.post, xyz_banc_post, xyz_u_post, inside_post); gc()

###########################
### MANC               ###
###########################

message("\n=== MANC ===")

if (is.null(fm)) {
  message("  Skipping MANC (metadata unavailable)")
  manc.dim <- data.frame(manc_id = character(0), super_class = character(0),
                          dimorphism = character(0))
} else {
  manc.dim <- fm %>%
    dplyr::filter(sexually_dimorphic %in% c("dimorphic", "male-specific")) %>%
    dplyr::rename(dimorphism = sexually_dimorphic)
  manc.dim <- fix_super_class(manc.dim, "cell_type")
  manc.dim <- manc.dim %>%
    dplyr::mutate(super_class = remap_super_class(super_class)) %>%
    dplyr::filter(super_class %in% super_classes) %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::select(manc_id, super_class, dimorphism)
}
message(sprintf("  %d dimorphic MANC neurons", nrow(manc.dim)))

vnc_manc_neuron_counts <- manc.dim %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "MANC")

# Skip synapse download/transform if no dimorphic neurons found
if (nrow(manc.dim) == 0) {
  message("  No dimorphic MANC neurons in target super_classes; skipping")
  manc_vnc.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                                  super_class = character(0), dimorphism = character(0),
                                  dataset = character(0))
  manc_vnc.density.post <- manc_vnc.density
} else {

# Read MANC synapse data with lazy loading.
# On O2: pre-transformed to BANC nm space. Off O2: native MANC um from GCS.
if (on_o2) {
  manc_syn_file <- file.path(banc.connectome.save.path, "manc/manc_121_synapses.feather")
  manc.pre <- arrow::read_feather(manc_syn_file) %>%
    dplyr::rename(pre_id = pre, post_id = post) %>%
    dplyr::filter(prepost == 0, pre_id %in% manc.dim$manc_id) %>%
    dplyr::collect() %>%
    dplyr::inner_join(manc.dim, by = c("pre_id" = "manc_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(manc.pre)))
  # O2 MANC synapses are already in BANC nm space
  xyz_u_manc <- banc_to_jrcvnc2018u(nat::xyzmatrix(manc.pre))
} else {
  manc_syn_file <- gcs_cache(paste0(gcs_bucket, "/manc/manc_121_synapses.parquet"))
  manc_ids_str <- as.character(manc.dim$manc_id)

  manc.pre <- open_parquet(manc_syn_file) %>%
    dplyr::filter(prepost == 0L, pre %in% manc_ids_str) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(manc.dim %>% dplyr::mutate(manc_id = as.character(manc_id)),
                      by = c("pre" = "manc_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(manc.pre)))
  # GCS MANC synapses are in native MANC um space
  xyz_u_manc <- manc_to_jrcvnc2018u(nat::xyzmatrix(manc.pre))
}

# Filter to points inside JRCVNC2018U template
inside_manc <- nat::pointsinside(xyz_u_manc, jrc_surf)
inside_manc[is.na(inside_manc)] <- FALSE

manc_vnc.density <- data.frame(X = xyz_u_manc[inside_manc, 1],
                                Y = xyz_u_manc[inside_manc, 2],
                                Z = xyz_u_manc[inside_manc, 3],
                                super_class = manc.pre$super_class[inside_manc],
                                dimorphism = manc.pre$dimorphism[inside_manc],
                                dataset = "MANC")
message(sprintf("  %d presynapses within JRCVNC2018U", sum(inside_manc)))
rm(manc.pre, xyz_u_manc, inside_manc); gc()

# MANC postsynapses: input synapses to dimorphic neurons
message("  Loading MANC postsynapses...")
if (on_o2) {
  manc.post <- arrow::read_feather(manc_syn_file) %>%
    dplyr::rename(pre_id = pre, post_id = post) %>%
    dplyr::filter(prepost == 1, post_id %in% manc.dim$manc_id) %>%
    dplyr::collect() %>%
    dplyr::inner_join(manc.dim, by = c("post_id" = "manc_id"))
  gc()
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(manc.post)))
  xyz_u_manc_post <- banc_to_jrcvnc2018u(nat::xyzmatrix(manc.post))
} else {
  manc.post <- open_parquet(manc_syn_file) %>%
    dplyr::filter(prepost == 1L, post %in% manc_ids_str) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(manc.dim %>% dplyr::mutate(manc_id = as.character(manc_id)),
                      by = c("post" = "manc_id"))
  gc()
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(manc.post)))
  xyz_u_manc_post <- manc_to_jrcvnc2018u(nat::xyzmatrix(manc.post))
}

inside_manc_post <- nat::pointsinside(xyz_u_manc_post, jrc_surf)
inside_manc_post[is.na(inside_manc_post)] <- FALSE

manc_vnc.density.post <- data.frame(X = xyz_u_manc_post[inside_manc_post, 1],
                                     Y = xyz_u_manc_post[inside_manc_post, 2],
                                     Z = xyz_u_manc_post[inside_manc_post, 3],
                                     super_class = manc.post$super_class[inside_manc_post],
                                     dimorphism = manc.post$dimorphism[inside_manc_post],
                                     dataset = "MANC")
message(sprintf("  %d postsynapses within JRCVNC2018U", sum(inside_manc_post)))
rm(manc.post, xyz_u_manc_post, inside_manc_post); gc()

} # end else (manc.dim > 0)

###########################
### maleCNS (VNC)       ###
###########################

message("\n=== maleCNS (VNC) ===")

# Use dimorphism column directly from SeaTable
if (is.null(malecns.meta)) {
  malecns_vnc.dim <- data.frame(malecns_09_id = character(0), super_class = character(0),
                                 dimorphism = character(0))
} else {
  malecns_vnc.dim <- malecns.meta %>%
    dplyr::filter(dimorphism %in% c("dimorphic", "male-specific"))
  mcns_ct_col <- if ("manc_cell_type" %in% colnames(malecns_vnc.dim)) "manc_cell_type" else "cell_type"
  malecns_vnc.dim <- fix_super_class(malecns_vnc.dim, mcns_ct_col)
  malecns_vnc.dim <- malecns_vnc.dim %>%
    dplyr::mutate(super_class = remap_super_class(super_class)) %>%
    dplyr::filter(super_class %in% super_classes) %>%
    dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
    dplyr::select(malecns_09_id, super_class, dimorphism)
}
message(sprintf("  %d dimorphic maleCNS neurons", nrow(malecns_vnc.dim)))

vnc_malecns_neuron_counts <- malecns_vnc.dim %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "maleCNS")

# Read maleCNS synapse data with lazy loading.
# maleCNS coordinates are in 8nm voxel space (JRCFIB2022M).
# Transform chain: 8nm voxels -> nm -> JRCFIB2022M -> BANC nm -> JRCVNC2018U.
malecns_vnc_ids_str <- as.character(malecns_vnc.dim$malecns_09_id)
if (on_o2) {
  malecns.path <- '/n/data1/hms/neurobio/wilson/connectomes/malecns'
  malecns_syn_file <- file.path(malecns.path, "malecns_09_synapses.parquet")
  malecns_vnc.pre <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(pre %in% malecns_vnc_ids_str) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_vnc.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("pre" = "malecns_09_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(malecns_vnc.pre)))
} else {
  malecns_syn_file <- gcs_cache(paste0(gcs_bucket, "/malecns/malecns_09_synapses.parquet"))
  malecns_vnc.pre <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(pre %in% malecns_vnc_ids_str) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_vnc.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("pre" = "malecns_09_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(malecns_vnc.pre)))
}

# Convert from 8nm voxel coordinates to nm, then transform to JRCVNC2018U
xyz_nm <- nat::xyzmatrix(malecns_vnc.pre) * 8
xyz_u_malecns <- malecns_to_jrcvnc2018u(xyz_nm)

# Filter to points inside JRCVNC2018U template
inside_malecns <- nat::pointsinside(xyz_u_malecns, jrc_surf)
inside_malecns[is.na(inside_malecns)] <- FALSE

malecns_vnc.density <- data.frame(X = xyz_u_malecns[inside_malecns, 1],
                                   Y = xyz_u_malecns[inside_malecns, 2],
                                   Z = xyz_u_malecns[inside_malecns, 3],
                                   super_class = malecns_vnc.pre$super_class[inside_malecns],
                                   dimorphism = malecns_vnc.pre$dimorphism[inside_malecns],
                                   dataset = "maleCNS")
message(sprintf("  %d presynapses within JRCVNC2018U", sum(inside_malecns)))
rm(malecns_vnc.pre, xyz_nm, xyz_u_malecns, inside_malecns); gc()

# maleCNS VNC postsynapses
message("  Loading maleCNS postsynapses...")
if (on_o2) {
  malecns_vnc.post <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(post %in% malecns_vnc_ids_str) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_vnc.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("post" = "malecns_09_id"))
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(malecns_vnc.post)))
} else {
  malecns_vnc.post <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(post %in% malecns_vnc_ids_str) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_vnc.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("post" = "malecns_09_id"))
  gc()
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(malecns_vnc.post)))
}

xyz_nm_post <- nat::xyzmatrix(malecns_vnc.post) * 8
xyz_u_malecns_post <- malecns_to_jrcvnc2018u(xyz_nm_post)

inside_malecns_post <- nat::pointsinside(xyz_u_malecns_post, jrc_surf)
inside_malecns_post[is.na(inside_malecns_post)] <- FALSE

malecns_vnc.density.post <- data.frame(X = xyz_u_malecns_post[inside_malecns_post, 1],
                                        Y = xyz_u_malecns_post[inside_malecns_post, 2],
                                        Z = xyz_u_malecns_post[inside_malecns_post, 3],
                                        super_class = malecns_vnc.post$super_class[inside_malecns_post],
                                        dimorphism = malecns_vnc.post$dimorphism[inside_malecns_post],
                                        dataset = "maleCNS")
message(sprintf("  %d postsynapses within JRCVNC2018U", sum(inside_malecns_post)))
rm(malecns_vnc.post, xyz_nm_post, xyz_u_malecns_post, inside_malecns_post); gc()

###########################
### VNC combine + save  ###
###########################

message("\n=== Combining VNC datasets ===")

vnc_neuron_counts <- dplyr::bind_rows(vnc_banc_neuron_counts, vnc_manc_neuron_counts, vnc_malecns_neuron_counts)
message("\n  VNC neuron counts per dataset x super_class x dimorphism_group:")
for (dg in c("dimorphic", "sex-specific")) {
  message(sprintf("  --- %s ---", dg))
  for (sc in super_classes) {
    counts_str <- paste(sapply(c("BANC", "maleCNS", "MANC"), function(ds) {
      n <- vnc_neuron_counts$n_neurons[vnc_neuron_counts$dataset == ds &
                                         vnc_neuron_counts$super_class == sc &
                                         vnc_neuron_counts$dimorphism_group == dg]
      sprintf("%s=%d", ds, if (length(n)) n else 0L)
    }), collapse = ", ")
    message(sprintf("    %s: %s", sc, counts_str))
  }
}

# Presynapses
vnc_all.density <- dplyr::bind_rows(banc_vnc.density, manc_vnc.density, malecns_vnc.density) %>%
  add_dimorphism_group()
vnc_all.density$dataset <- factor(vnc_all.density$dataset, levels = c("BANC", "maleCNS", "MANC"))
vnc_all.density$super_class <- factor(vnc_all.density$super_class, levels = super_classes)

vnc_feather_path <- file.path(vnc_plot.dir, "dimorphic_presynapses_jrcvnc2018u.feather")
arrow::write_feather(vnc_all.density, vnc_feather_path)
message(sprintf("  %d total presynapse points across all datasets", nrow(vnc_all.density)))
message(sprintf("  Saved: %s", vnc_feather_path))

# Postsynapses
vnc_all.density.post <- dplyr::bind_rows(banc_vnc.density.post, manc_vnc.density.post, malecns_vnc.density.post) %>%
  add_dimorphism_group()
vnc_all.density.post$dataset <- factor(vnc_all.density.post$dataset, levels = c("BANC", "maleCNS", "MANC"))
vnc_all.density.post$super_class <- factor(vnc_all.density.post$super_class, levels = super_classes)

vnc_feather_path_post <- file.path(vnc_plot.dir, "dimorphic_postsynapses_jrcvnc2018u.feather")
arrow::write_feather(vnc_all.density.post, vnc_feather_path_post)
message(sprintf("  %d total postsynapse points across all datasets", nrow(vnc_all.density.post)))
message(sprintf("  Saved: %s", vnc_feather_path_post))

# Log dimorphism group counts
for (dg in c("dimorphic", "sex-specific")) {
  n_pre <- sum(vnc_all.density$dimorphism_group == dg)
  n_post <- sum(vnc_all.density.post$dimorphism_group == dg)
  message(sprintf("  %s: %d presynapses, %d postsynapses", dg, n_pre, n_post))
  for (ds in c("BANC", "maleCNS", "MANC")) {
    n <- sum(vnc_all.density$dimorphism_group == dg & vnc_all.density$dataset == ds)
    message(sprintf("    %s: %d presynapses", ds, n))
  }
}

###########################
### VNC Plotting        ###
###########################

message("\n=== VNC Plotting ===")

template_2d_vnc <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_mat)
template_2d_side <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_side_mat)

for (dg in c("dimorphic", "sex-specific")) {
  dg_label <- gsub("-", "_", dg)  # for filenames

  for (syn_type in c("pre", "post")) {
    df_all <- if (syn_type == "pre") vnc_all.density else vnc_all.density.post
    df_sub <- df_all %>% dplyr::filter(dimorphism_group == dg)
    # Re-level factors to drop unused levels
    df_sub$dataset <- factor(df_sub$dataset, levels = levels(df_all$dataset))
    df_sub$super_class <- factor(df_sub$super_class, levels = super_classes)

    syn_label <- if (syn_type == "pre") "presynapses" else "postsynapses"

    if (nrow(df_sub) == 0) {
      message(sprintf("  Skipping %s %s: no data", dg, syn_label))
      next
    }

    # Neuron counts for this dimorphism group
    counts_sub <- vnc_neuron_counts %>% dplyr::filter(dimorphism_group == dg)

    for (view in c("vnc", "vnc_side")) {
      rot_mat <- if (view == "vnc") jrc_vnc_mat else jrc_vnc_side_mat
      tmpl_2d <- if (view == "vnc") template_2d_vnc else template_2d_side
      view_label <- if (view == "vnc") "dorsal" else "lateral"

      # Unnormalized plot
      fname <- sprintf("%s_%s_%s.pdf", dg_label, syn_label, view)
      message(sprintf("  Plotting %s %s %s...", dg, syn_label, view_label))
      g <- make_density_plot(df_sub, rot_mat, tmpl_2d, neuron_counts = counts_sub)
      ggsave(file.path(vnc_plot.dir, fname), g, width = 16, height = 12,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname))

      # Normalized by dataset
      fname_norm <- sprintf("%s_%s_%s_normalized.pdf", dg_label, syn_label, view)
      message(sprintf("  Plotting %s %s %s (normalized)...", dg, syn_label, view_label))
      g_norm <- make_density_plot(df_sub, rot_mat, tmpl_2d, normalize_by = "dataset",
                                  neuron_counts = counts_sub)
      ggsave(file.path(vnc_plot.dir, fname_norm), g_norm, width = 16, height = 12,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname_norm))
    }
  }
}

###########################
### VNC Subtraction     ###
###########################

message("\n=== VNC Density Subtraction Plots ===")
message("  (all dimorphism labels merged; densities normalized by dataset total)")

vnc_comparisons <- c("maleCNS", "MANC")

for (syn_type in c("pre", "post")) {
  df_all <- if (syn_type == "pre") vnc_all.density else vnc_all.density.post
  syn_label <- if (syn_type == "pre") "presynapses" else "postsynapses"

  for (view in c("vnc", "vnc_side")) {
    rot_mat <- if (view == "vnc") jrc_vnc_mat else jrc_vnc_side_mat
    tmpl_2d <- if (view == "vnc") template_2d_vnc else template_2d_side
    view_label <- if (view == "vnc") "dorsal" else "lateral"

    # First pass: compute zlim across all comparisons for a shared colour scale
    sub_results <- list()
    global_zlim <- 0
    for (other in vnc_comparisons) {
      df_pair <- df_all %>% dplyr::filter(dataset %in% c("BANC", other))
      if (nrow(df_pair[df_pair$dataset == "BANC", ]) == 0 ||
          nrow(df_pair[df_pair$dataset == other, ]) == 0) next
      res <- make_subtraction_plot(df_pair, "BANC", other, rot_mat, tmpl_2d,
                                    super_classes, grid_n = 100, zlim = NULL)
      sub_results[[other]] <- res
      global_zlim <- max(global_zlim, res$zlim)
    }

    # Second pass: replot with shared zlim (both directions)
    for (other in names(sub_results)) {
      df_pair <- df_all %>% dplyr::filter(dataset %in% c("BANC", other))

      # BANC minus other
      res <- make_subtraction_plot(df_pair, "BANC", other, rot_mat, tmpl_2d,
                                    super_classes, grid_n = 100, zlim = global_zlim)
      fname <- sprintf("subtraction_BANC_minus_%s_%s_%s.pdf",
                        gsub(" ", "", other), syn_label, view)
      message(sprintf("  Plotting subtraction BANC - %s %s %s...", other, syn_label, view_label))
      ggsave(file.path(vnc_plot.dir, fname), res$plot, width = 16, height = 5,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname))

      # Other minus BANC
      res_rev <- make_subtraction_plot(df_pair, other, "BANC", rot_mat, tmpl_2d,
                                        super_classes, grid_n = 100, zlim = global_zlim)
      fname_rev <- sprintf("subtraction_%s_minus_BANC_%s_%s.pdf",
                            gsub(" ", "", other), syn_label, view)
      message(sprintf("  Plotting subtraction %s - BANC %s %s...", other, syn_label, view_label))
      ggsave(file.path(vnc_plot.dir, fname_rev), res_rev$plot, width = 16, height = 5,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname_rev))
    }
  }
}

# Clean up VNC-specific objects
rm(banc_vnc.density, banc_vnc.density.post,
   manc_vnc.density, manc_vnc.density.post,
   malecns_vnc.density, malecns_vnc.density.post,
   vnc_all.density, vnc_all.density.post,
   bc_vnc, manc.dim); gc()

message(sprintf("\n### VNC density plots complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

###########################################################
###########################################################
###                                                     ###
###                  BRAIN SECTION                       ###
###                                                     ###
###########################################################
###########################################################

message("\n\n========================================")
message("=== BRAIN DENSITY PLOTS ================")
message("========================================")

t_brain_start <- Sys.time()

###########################
### Brain transform funcs #
###########################

# BANC brain nm -> JRC2018U um
banc_to_jrc2018u <- function(xyz) {
  xyz_jrc <- bancr::banc_to_JRC2018F(xyz, region = "brain",
                                       method = "tpsreg",
                                       banc.units = "nm",
                                       inverse = FALSE)
  nat.templatebrains::xform_brain(xyz_jrc,
                                   sample = "JRC2018F",
                                   reference = "JRC2018U")
}

# FAFB14 nm -> JRC2018U um
fafb_to_jrc2018u <- function(xyz) {
  xyz_jrc <- nat.templatebrains::xform_brain(xyz,
                                              sample = "FAFB14",
                                              reference = "JRC2018F")
  nat.templatebrains::xform_brain(xyz_jrc,
                                   sample = "JRC2018F",
                                   reference = "JRC2018U")
}

# Transform maleCNS coordinates (nm) to JRC2018U via BANC brain.
malecns_to_jrc2018u <- function(xyz_nm, chunk_size = 500000L) {
  xyz_banc <- malecns_to_banc_nm(xyz_nm, chunk_size = chunk_size)
  banc_to_jrc2018u(xyz_banc)
}

###########################
### Brain template      ###
###########################

jrc_brain_surf <- tryCatch({
  nat.flybrains::JRC2018U.surf
}, error = function(e) {
  tryCatch({
    get("JFRC2018U.surf", envir = asNamespace("nat.flybrains"))
  }, error = function(e2) {
    message("  Warning: JRC2018U.surf not found, trying JRC2018F.surf")
    tryCatch(nat.flybrains::JRC2018F.surf, error = function(e3) NULL)
  })
})
if (is.null(jrc_brain_surf)) {
  stop("Could not load brain template surface. Check nat.flybrains/nat.jrcbrains installation.")
}
message("  Brain template surface loaded")

###########################
### Brain settings      ###
###########################

brain_plot.dir <- file.path("images", "dimorphic_density_brain")
dir.create(brain_plot.dir, recursive = TRUE, showWarnings = FALSE)

brain_super_classes <- c("sensory", "visual", "central_brain_intrinsic",
                         "ascending", "descending", "effector")

###########################
### BANC (brain)        ###
###########################

message("\n=== BANC (brain) ===")

# Apply brain-specific remap to shared bc_base
bc_brain <- bc_base %>%
  dplyr::mutate(super_class = remap_brain_super_class(super_class)) %>%
  dplyr::filter(!is.na(super_class)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("  %d dimorphic BANC neurons with brain super_classes", nrow(bc_brain)))

brain_banc_neuron_counts <- bc_brain %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "BANC")

# Read BANC synapse data -- presynapses
message("  Loading BANC presynapses...")
banc_brain_ids <- na.omit(unique(as.character(bc_brain$root_626)))
if (on_o2) {
  banc_syn_file <- Sys.glob(file.path(banc.connectome.save.path,
                                       "banc/banc_*_synapses.parquet"))[1]
  banc_brain.pre <- open_parquet(banc_syn_file) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::filter(pre %in% banc_brain_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_brain %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("pre" = "root_626"))
} else {
  banc_syn_file <- gcs_cache(paste0(gcs_bucket, "/banc/banc_746_synapses.parquet"))
  banc_brain.pre <- open_parquet(banc_syn_file) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::filter(pre %in% banc_brain_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_brain %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("pre" = "root_626"))
}
gc()
message(sprintf("  %d presynapses from dimorphic neurons", nrow(banc_brain.pre)))

# Transform to JRC2018U brain space
xyz_banc <- nat::xyzmatrix(banc_brain.pre)
xyz_u <- banc_to_jrc2018u(xyz_banc)

# Keep only points inside brain template
inside <- nat::pointsinside(xyz_u, jrc_brain_surf)
inside[is.na(inside)] <- FALSE

banc_brain.density <- data.frame(X = xyz_u[inside, 1], Y = xyz_u[inside, 2], Z = xyz_u[inside, 3],
                                  super_class = banc_brain.pre$super_class[inside],
                                  dimorphism = banc_brain.pre$dimorphism[inside],
                                  dataset = "BANC")
message(sprintf("  %d presynapses within JRC2018U brain", sum(inside)))
rm(banc_brain.pre, xyz_banc, xyz_u, inside); gc()

# BANC brain postsynapses
message("  Loading BANC postsynapses...")
if (on_o2) {
  banc_brain.post <- open_parquet(banc_syn_file) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::filter(post %in% banc_brain_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_brain %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("post" = "root_626"))
} else {
  banc_brain.post <- open_parquet(banc_syn_file) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::filter(post %in% banc_brain_ids) %>%
    dplyr::collect() %>%
    dplyr::rename(X = x, Y = y, Z = z) %>%
    dplyr::inner_join(bc_brain %>% dplyr::select(root_626, super_class, dimorphism) %>%
                        dplyr::mutate(root_626 = as.character(root_626)),
                      by = c("post" = "root_626"))
}
gc()
message(sprintf("  %d postsynapses to dimorphic neurons", nrow(banc_brain.post)))

xyz_banc_post <- nat::xyzmatrix(banc_brain.post)
xyz_u_post <- banc_to_jrc2018u(xyz_banc_post)
inside_post <- nat::pointsinside(xyz_u_post, jrc_brain_surf)
inside_post[is.na(inside_post)] <- FALSE

banc_brain.density.post <- data.frame(X = xyz_u_post[inside_post, 1], Y = xyz_u_post[inside_post, 2],
                                       Z = xyz_u_post[inside_post, 3],
                                       super_class = banc_brain.post$super_class[inside_post],
                                       dimorphism = banc_brain.post$dimorphism[inside_post],
                                       dataset = "BANC")
message(sprintf("  %d postsynapses within JRC2018U brain", sum(inside_post)))
rm(banc_brain.post, xyz_banc_post, xyz_u_post, inside_post); gc()

###########################
### FAFB (brain)        ###
###########################

message("\n=== FAFB (brain) ===")

if (is.null(fafb.meta)) {
  message("  Skipping FAFB (metadata unavailable)")
  fafb.dim <- data.frame(fafb_id = character(0), super_class = character(0),
                           dimorphism = character(0))
  fafb_neuron_counts <- data.frame(super_class = character(0), dimorphism_group = character(0),
                                    n_neurons = integer(0), dataset = character(0))
  fafb.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                               super_class = character(0), dimorphism = character(0),
                               dataset = character(0))
  fafb.density.post <- fafb.density
} else {
  fafb.dim <- fafb.meta %>%
    dplyr::filter(sexually_dimorphic %in% c("dimorphic", "female-specific")) %>%
    dplyr::rename(dimorphism = sexually_dimorphic)
  fafb.dim <- fix_super_class(fafb.dim, "cell_type")
  fafb.dim <- fafb.dim %>%
    dplyr::mutate(super_class = remap_brain_super_class(super_class)) %>%
    dplyr::filter(!is.na(super_class)) %>%
    dplyr::distinct(fafb_id, .keep_all = TRUE) %>%
    dplyr::select(fafb_id, super_class, dimorphism)
}
message(sprintf("  %d dimorphic FAFB neurons with brain super_classes", nrow(fafb.dim)))

fafb_neuron_counts <- fafb.dim %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "FAFB")

if (nrow(fafb.dim) == 0) {
  message("  No dimorphic FAFB neurons in brain super_classes; skipping")
  fafb.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                               super_class = character(0), dimorphism = character(0),
                               dataset = character(0))
  fafb.density.post <- fafb.density
} else {

# Load FAFB synapse data
fafb_syn_available <- FALSE
if (on_o2) {
  fafb_syn_file <- file.path(banc.connectome.save.path,
                              "fafb/fafb_783_synapses.parquet")
  if (file.exists(fafb_syn_file)) fafb_syn_available <- TRUE
}
if (!fafb_syn_available) {
  # Try GCS
  fafb_syn_file <- tryCatch({
    gcs_cache(paste0(gcs_bucket, "/fafb/fafb_783_synapses.parquet"))
  }, error = function(e) NULL)
  if (!is.null(fafb_syn_file) && file.exists(fafb_syn_file)) {
    fafb_syn_available <- TRUE
  }
}

if (!fafb_syn_available) {
  message("  FAFB synapse data not available; skipping synapse density")
  fafb.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                               super_class = character(0), dimorphism = character(0),
                               dataset = character(0))
  fafb.density.post <- fafb.density
} else {
  message(sprintf("  Loading FAFB synapse data from %s", basename(fafb_syn_file)))
  fafb_ids_str <- as.character(fafb.dim$fafb_id)

  # FAFB presynapses
  fafb.pre <- open_parquet(fafb_syn_file) %>%
    dplyr::select(dplyr::any_of(c("pre", "pre_pt_root_id", "x", "y", "z",
                                   "pre_x", "pre_y", "pre_z"))) %>%
    dplyr::collect()

  # Detect column naming (pre or pre_pt_root_id)
  pre_col <- intersect(c("pre", "pre_pt_root_id"), colnames(fafb.pre))[1]
  xyz_cols <- if (all(c("x", "y", "z") %in% colnames(fafb.pre))) c("x", "y", "z") else
              if (all(c("pre_x", "pre_y", "pre_z") %in% colnames(fafb.pre))) c("pre_x", "pre_y", "pre_z") else NULL

  if (is.na(pre_col) || is.null(xyz_cols)) {
    message("  FAFB synapse column format not recognized; skipping")
    fafb.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                                 super_class = character(0), dimorphism = character(0),
                                 dataset = character(0))
    fafb.density.post <- fafb.density
  } else {
    fafb.pre <- fafb.pre %>%
      dplyr::filter(as.character(.data[[pre_col]]) %in% fafb_ids_str) %>%
      dplyr::rename(X = !!xyz_cols[1], Y = !!xyz_cols[2], Z = !!xyz_cols[3]) %>%
      dplyr::inner_join(fafb.dim, by = stats::setNames("fafb_id", pre_col))
    gc()
    message(sprintf("  %d presynapses from dimorphic FAFB neurons", nrow(fafb.pre)))

    # Transform FAFB14 nm -> JRC2018U
    xyz_fafb <- nat::xyzmatrix(fafb.pre)
    xyz_u_fafb <- fafb_to_jrc2018u(xyz_fafb)

    inside_fafb <- nat::pointsinside(xyz_u_fafb, jrc_brain_surf)
    inside_fafb[is.na(inside_fafb)] <- FALSE

    fafb.density <- data.frame(X = xyz_u_fafb[inside_fafb, 1], Y = xyz_u_fafb[inside_fafb, 2],
                                Z = xyz_u_fafb[inside_fafb, 3],
                                super_class = fafb.pre$super_class[inside_fafb],
                                dimorphism = fafb.pre$dimorphism[inside_fafb],
                                dataset = "FAFB")
    message(sprintf("  %d presynapses within JRC2018U brain", sum(inside_fafb)))
    rm(fafb.pre, xyz_fafb, xyz_u_fafb, inside_fafb); gc()

    # FAFB postsynapses
    message("  Loading FAFB postsynapses...")
    post_col <- intersect(c("post", "post_pt_root_id"), colnames(
      open_parquet(fafb_syn_file) %>%
        head(1) %>% dplyr::collect()
    ))[1]
    post_xyz_cols <- if (all(c("x", "y", "z") %in% colnames(
      open_parquet(fafb_syn_file) %>%
        head(1) %>% dplyr::collect()
    ))) c("x", "y", "z") else c("post_x", "post_y", "post_z")

    fafb.post <- open_parquet(fafb_syn_file) %>%
      dplyr::collect() %>%
      dplyr::filter(as.character(.data[[post_col]]) %in% fafb_ids_str) %>%
      dplyr::rename(X = !!post_xyz_cols[1], Y = !!post_xyz_cols[2], Z = !!post_xyz_cols[3]) %>%
      dplyr::inner_join(fafb.dim, by = stats::setNames("fafb_id", post_col))
    gc()
    message(sprintf("  %d postsynapses to dimorphic FAFB neurons", nrow(fafb.post)))

    xyz_fafb_post <- nat::xyzmatrix(fafb.post)
    xyz_u_fafb_post <- fafb_to_jrc2018u(xyz_fafb_post)
    inside_fafb_post <- nat::pointsinside(xyz_u_fafb_post, jrc_brain_surf)
    inside_fafb_post[is.na(inside_fafb_post)] <- FALSE

    fafb.density.post <- data.frame(X = xyz_u_fafb_post[inside_fafb_post, 1],
                                     Y = xyz_u_fafb_post[inside_fafb_post, 2],
                                     Z = xyz_u_fafb_post[inside_fafb_post, 3],
                                     super_class = fafb.post$super_class[inside_fafb_post],
                                     dimorphism = fafb.post$dimorphism[inside_fafb_post],
                                     dataset = "FAFB")
    message(sprintf("  %d postsynapses within JRC2018U brain", sum(inside_fafb_post)))
    rm(fafb.post, xyz_fafb_post, xyz_u_fafb_post, inside_fafb_post); gc()
  }
}
} # end fafb.dim > 0

###########################
### maleCNS (brain)     ###
###########################

message("\n=== maleCNS (brain) ===")

if (is.null(malecns.meta)) {
  message("  Skipping maleCNS (metadata unavailable)")
  malecns_brain.dim <- data.frame(malecns_09_id = character(0), super_class = character(0),
                                   dimorphism = character(0))
  malecns_brain.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                                       super_class = character(0), dimorphism = character(0),
                                       dataset = character(0))
  malecns_brain.density.post <- malecns_brain.density
  brain_malecns_neuron_counts <- data.frame(super_class = character(0), dimorphism_group = character(0),
                                             n_neurons = integer(0), dataset = character(0))
} else {

malecns_brain.dim <- malecns.meta %>%
  dplyr::filter(dimorphism %in% c("dimorphic", "male-specific"))
mcns_brain_ct_col <- if ("manc_cell_type" %in% colnames(malecns_brain.dim)) "manc_cell_type" else "cell_type"
malecns_brain.dim <- fix_super_class(malecns_brain.dim, mcns_brain_ct_col)
malecns_brain.dim <- malecns_brain.dim %>%
  dplyr::mutate(super_class = remap_brain_super_class(super_class)) %>%
  dplyr::filter(!is.na(super_class)) %>%
  dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
  dplyr::select(malecns_09_id, super_class, dimorphism)
message(sprintf("  %d dimorphic maleCNS neurons with brain super_classes", nrow(malecns_brain.dim)))

brain_malecns_neuron_counts <- malecns_brain.dim %>%
  dplyr::mutate(dimorphism_group = dplyr::if_else(
    dimorphism == "dimorphic", "dimorphic", "sex-specific")) %>%
  dplyr::count(super_class, dimorphism_group, name = "n_neurons") %>%
  dplyr::mutate(dataset = "maleCNS")

if (nrow(malecns_brain.dim) == 0) {
  message("  No dimorphic maleCNS neurons; skipping")
  malecns_brain.density <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0),
                                       super_class = character(0), dimorphism = character(0),
                                       dataset = character(0))
  malecns_brain.density.post <- malecns_brain.density
} else {

# maleCNS synapse data
if (on_o2) {
  malecns.path <- '/n/data1/hms/neurobio/wilson/connectomes/malecns'
  malecns_syn_file <- file.path(malecns.path, "malecns_09_synapses.parquet")
  malecns_brain_ids_str <- as.character(malecns_brain.dim$malecns_09_id)
  malecns_brain.pre <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(pre %in% malecns_brain_ids_str) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_brain.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("pre" = "malecns_09_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(malecns_brain.pre)))
} else {
  malecns_syn_file <- gcs_cache(paste0(gcs_bucket, "/malecns/malecns_09_synapses.parquet"))
  malecns_brain_ids_str <- as.character(malecns_brain.dim$malecns_09_id)
  malecns_brain.pre <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(pre %in% malecns_brain_ids_str) %>%
    dplyr::select(pre, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_brain.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("pre" = "malecns_09_id"))
  gc()
  message(sprintf("  %d presynapses from dimorphic neurons", nrow(malecns_brain.pre)))
}

# maleCNS: 8nm voxels -> nm -> JRC2018U brain
xyz_nm <- nat::xyzmatrix(malecns_brain.pre) * 8
xyz_u_malecns <- malecns_to_jrc2018u(xyz_nm)

inside_malecns <- nat::pointsinside(xyz_u_malecns, jrc_brain_surf)
inside_malecns[is.na(inside_malecns)] <- FALSE

malecns_brain.density <- data.frame(X = xyz_u_malecns[inside_malecns, 1],
                                     Y = xyz_u_malecns[inside_malecns, 2],
                                     Z = xyz_u_malecns[inside_malecns, 3],
                                     super_class = malecns_brain.pre$super_class[inside_malecns],
                                     dimorphism = malecns_brain.pre$dimorphism[inside_malecns],
                                     dataset = "maleCNS")
message(sprintf("  %d presynapses within JRC2018U brain", sum(inside_malecns)))
rm(malecns_brain.pre, xyz_nm, xyz_u_malecns, inside_malecns); gc()

# maleCNS brain postsynapses
message("  Loading maleCNS postsynapses...")
if (on_o2) {
  malecns_brain.post <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(post %in% malecns_brain_ids_str) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_brain.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("post" = "malecns_09_id"))
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(malecns_brain.post)))
} else {
  malecns_brain.post <- open_parquet(malecns_syn_file) %>%
    dplyr::filter(post %in% malecns_brain_ids_str) %>%
    dplyr::select(post, x, y, z) %>%
    dplyr::collect() %>%
    dplyr::inner_join(malecns_brain.dim %>% dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
                      by = c("post" = "malecns_09_id"))
  gc()
  message(sprintf("  %d postsynapses to dimorphic neurons", nrow(malecns_brain.post)))
}

xyz_nm_post <- nat::xyzmatrix(malecns_brain.post) * 8
xyz_u_malecns_post <- malecns_to_jrc2018u(xyz_nm_post)

inside_malecns_post <- nat::pointsinside(xyz_u_malecns_post, jrc_brain_surf)
inside_malecns_post[is.na(inside_malecns_post)] <- FALSE

malecns_brain.density.post <- data.frame(X = xyz_u_malecns_post[inside_malecns_post, 1],
                                          Y = xyz_u_malecns_post[inside_malecns_post, 2],
                                          Z = xyz_u_malecns_post[inside_malecns_post, 3],
                                          super_class = malecns_brain.post$super_class[inside_malecns_post],
                                          dimorphism = malecns_brain.post$dimorphism[inside_malecns_post],
                                          dataset = "maleCNS")
message(sprintf("  %d postsynapses within JRC2018U brain", sum(inside_malecns_post)))
rm(malecns_brain.post, xyz_nm_post, xyz_u_malecns_post, inside_malecns_post); gc()

} # end malecns_brain.dim > 0
} # end malecns.meta not null

###########################
### Brain combine + save ##
###########################

message("\n=== Combining brain datasets ===")

brain_datasets <- c("BANC", "FAFB", "maleCNS")

brain_neuron_counts <- dplyr::bind_rows(brain_banc_neuron_counts, fafb_neuron_counts, brain_malecns_neuron_counts)
message("\n  Brain neuron counts per dataset x super_class x dimorphism_group:")
for (dg in c("dimorphic", "sex-specific")) {
  message(sprintf("  --- %s ---", dg))
  for (sc in brain_super_classes) {
    counts_str <- paste(sapply(brain_datasets, function(ds) {
      n <- brain_neuron_counts$n_neurons[brain_neuron_counts$dataset == ds &
                                           brain_neuron_counts$super_class == sc &
                                           brain_neuron_counts$dimorphism_group == dg]
      sprintf("%s=%d", ds, if (length(n)) n else 0L)
    }), collapse = ", ")
    message(sprintf("    %s: %s", sc, counts_str))
  }
}

# Presynapses
brain_all.density <- dplyr::bind_rows(banc_brain.density, fafb.density, malecns_brain.density) %>%
  add_dimorphism_group()
brain_all.density$dataset <- factor(brain_all.density$dataset, levels = brain_datasets)
brain_all.density$super_class <- factor(brain_all.density$super_class, levels = brain_super_classes)

brain_feather_path <- file.path(brain_plot.dir, "dimorphic_presynapses_jrc2018u_brain.feather")
arrow::write_feather(brain_all.density, brain_feather_path)
message(sprintf("  %d total presynapse points across all datasets", nrow(brain_all.density)))

# Postsynapses
brain_all.density.post <- dplyr::bind_rows(banc_brain.density.post, fafb.density.post, malecns_brain.density.post) %>%
  add_dimorphism_group()
brain_all.density.post$dataset <- factor(brain_all.density.post$dataset, levels = brain_datasets)
brain_all.density.post$super_class <- factor(brain_all.density.post$super_class, levels = brain_super_classes)

brain_feather_path_post <- file.path(brain_plot.dir, "dimorphic_postsynapses_jrc2018u_brain.feather")
arrow::write_feather(brain_all.density.post, brain_feather_path_post)
message(sprintf("  %d total postsynapse points across all datasets", nrow(brain_all.density.post)))

###########################
### Brain Plotting      ###
###########################

message("\n=== Brain Plotting ===")

template_2d_dor <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_brain_surf), rotation_matrix = jrc_brain_dor_mat)
template_2d_ant <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_brain_surf), rotation_matrix = jrc_brain_ant_mat)

# Generate plots
for (dg in c("dimorphic", "sex-specific")) {
  dg_label <- gsub("-", "_", dg)

  for (syn_type in c("pre", "post")) {
    df_all <- if (syn_type == "pre") brain_all.density else brain_all.density.post
    df_sub <- df_all %>% dplyr::filter(dimorphism_group == dg)
    df_sub$dataset <- factor(df_sub$dataset, levels = levels(df_all$dataset))
    df_sub$super_class <- factor(df_sub$super_class, levels = brain_super_classes)

    syn_label <- if (syn_type == "pre") "presynapses" else "postsynapses"

    if (nrow(df_sub) == 0) {
      message(sprintf("  Skipping %s %s: no data", dg, syn_label))
      next
    }

    counts_sub <- brain_neuron_counts %>% dplyr::filter(dimorphism_group == dg)

    for (view in c("brain", "brain_front")) {
      rot_mat <- if (view == "brain") jrc_brain_dor_mat else jrc_brain_ant_mat
      tmpl_2d <- if (view == "brain") template_2d_dor else template_2d_ant
      view_label <- if (view == "brain") "dorsal" else "anterior"

      fname <- sprintf("brain_%s_%s_%s.pdf", dg_label, syn_label, view)
      message(sprintf("  Plotting %s %s %s...", dg, syn_label, view_label))
      g <- make_density_plot(df_sub, rot_mat, tmpl_2d, neuron_counts = counts_sub)
      ggsave(file.path(brain_plot.dir, fname), g, width = 19, height = 12,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname))

      fname_norm <- sprintf("brain_%s_%s_%s_normalized.pdf", dg_label, syn_label, view)
      message(sprintf("  Plotting %s %s %s (normalized)...", dg, syn_label, view_label))
      g_norm <- make_density_plot(df_sub, rot_mat, tmpl_2d, normalize_by = "dataset",
                                  neuron_counts = counts_sub)
      ggsave(file.path(brain_plot.dir, fname_norm), g_norm, width = 19, height = 12,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname_norm))
    }
  }
}

###########################
### Brain Subtraction   ###
###########################

message("\n=== Brain Density Subtraction Plots ===")

brain_comparisons <- c("FAFB", "maleCNS")

for (syn_type in c("pre", "post")) {
  df_all <- if (syn_type == "pre") brain_all.density else brain_all.density.post
  syn_label <- if (syn_type == "pre") "presynapses" else "postsynapses"

  for (view in c("brain", "brain_front")) {
    rot_mat <- if (view == "brain") jrc_brain_dor_mat else jrc_brain_ant_mat
    tmpl_2d <- if (view == "brain") template_2d_dor else template_2d_ant
    view_label <- if (view == "brain") "dorsal" else "anterior"

    # First pass: compute shared zlim
    sub_results <- list()
    global_zlim <- 0
    for (other in brain_comparisons) {
      df_pair <- df_all %>% dplyr::filter(dataset %in% c("BANC", other))
      if (nrow(df_pair[df_pair$dataset == "BANC", ]) == 0 ||
          nrow(df_pair[df_pair$dataset == other, ]) == 0) next
      res <- make_subtraction_plot(df_pair, "BANC", other, rot_mat, tmpl_2d,
                                    brain_super_classes, grid_n = 100, zlim = NULL)
      sub_results[[other]] <- res
      global_zlim <- max(global_zlim, res$zlim)
    }

    # Second pass: replot with shared zlim (both directions)
    for (other in names(sub_results)) {
      df_pair <- df_all %>% dplyr::filter(dataset %in% c("BANC", other))

      # BANC minus other
      res <- make_subtraction_plot(df_pair, "BANC", other, rot_mat, tmpl_2d,
                                    brain_super_classes, grid_n = 100, zlim = global_zlim)
      fname <- sprintf("brain_subtraction_BANC_minus_%s_%s_%s.pdf",
                        gsub(" ", "", other), syn_label, view)
      message(sprintf("  Plotting subtraction BANC - %s %s %s...", other, syn_label, view_label))
      ggsave(file.path(brain_plot.dir, fname), res$plot, width = 19, height = 5,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname))

      # Other minus BANC
      res_rev <- make_subtraction_plot(df_pair, other, "BANC", rot_mat, tmpl_2d,
                                        brain_super_classes, grid_n = 100, zlim = global_zlim)
      fname_rev <- sprintf("brain_subtraction_%s_minus_BANC_%s_%s.pdf",
                            gsub(" ", "", other), syn_label, view)
      message(sprintf("  Plotting subtraction %s - BANC %s %s...", other, syn_label, view_label))
      ggsave(file.path(brain_plot.dir, fname_rev), res_rev$plot, width = 19, height = 5,
             dpi = 300, bg = "white")
      message(sprintf("  Saved: %s", fname_rev))
    }
  }
}

message(sprintf("\n### banc: brain density plots complete [%s] ###",
                format(round(difftime(Sys.time(), t_brain_start, units = "mins"), 1))))

message(sprintf("\n### banc: all density plots complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})

