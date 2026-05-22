#' banc-visualise-clusters — Render per-cluster synapse-density neuroanatomy PNGs.
#'
#' For each spectral cluster produces full-CNS dorsal, brain frontal, and
#' VNC ventral synapse-density views via `nat.ggplot` + `bancr` neuropil
#' surfaces. Also emits a super_class composition stacked bar PDF and a
#' top-3 neuropil bar chart per cluster. Saved to both bancpipeline and
#' BANC-project image trees.
#'
#' @section Reads:
#'   - `data/cns_network/spectral_clustering_..._{v2|v3}.csv`
#'   - Synapse data: `<banc.save.path>/v<ver>/synapses_v2_human_readable.csv`
#'     (or earlier `synapses_250226_human_readable.csv` fallback)
#'
#' @section Writes:
#'   - `inst/images/cns_network/{cluster}_{main,brain,vnc}_neuroanatomy.png`
#'   - `inst/images/cns_network/{cluster}_super_class.pdf`
#'   - mirrored to `BANC-project/images/cns_network/`
#'
#' @section CLI:
#'   --source {v2,v3}              default v3
#'   --min-connection-strength N   default 1
#'   --cluster-count K             default 13
#'
#' @section Paper:
#'   Methods §"Spectral clustering".

###########################################################
### Visualise Spectral Clusters
###
### Produces synapse density plots for each spectral cluster
### with brain/VNC split views using nat.ggplot.
###
### For each cluster saves:
###   {cluster}_neuropil_bar.png       (top-3 neuropil bar chart)
###   {cluster}_main_neuroanatomy.png  (full CNS dorsal)
###   {cluster}_brain_neuroanatomy.png (brain frontal)
###   {cluster}_vnc_neuroanatomy.png   (VNC ventral)
###
### Input:
###   - data/cns_network/spectral_clustering_..._banc_version_{version}_..._{v2|v3}.csv
###   - Synapse data (v{banc.version} or earlier from banc.save.path)
###
### Output (saved to BOTH bancpipeline and BANC-project):
###   - {inst,/n/.../BANC-project}/images/cns_network/{cluster}_*.png
###   - {inst,/n/.../BANC-project}/images/cns_network/{cluster}_super_class.pdf
###
### Usage: Rscript banc/clustering/banc-visualise-clusters.R [--source v2|v3] \
###          [--min-connection-strength N] [--cluster-count K]
###########################################################
source("banc/banc-startup.R")
library(nat.ggplot)

local({
  message("### banc: Visualise spectral clusters ###")
  t_start <- Sys.time()

  # CLI: --source v2|v3 (default v3); dataset_version follows banc.version
  .cli_args <- commandArgs(trailingOnly = TRUE)
  i <- which(.cli_args == "--source")
  syn_source <- if (length(i) == 1 && length(.cli_args) >= i + 1) tolower(.cli_args[i + 1]) else "v3"
  if (!syn_source %in% c("v2", "v3"))
    stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", syn_source))

  .parse_int_flag <- function(flag, default) {
    j <- which(.cli_args == flag)
    if (length(j) == 1 && length(.cli_args) >= j + 1) return(as.integer(.cli_args[j + 1]))
    default
  }
  min_connection_strength <- .parse_int_flag("--min-connection-strength", 1L)
  cluster_count <- .parse_int_flag("--cluster-count", 14L)

  # ---------------------------------------------------------------
  # 1. Configuration
  # ---------------------------------------------------------------
  # Plots save to BOTH bancpipeline (working copy alongside the script) and
  # BANC-project (canonical publication location). Mirrors the spectral CSV
  # behaviour in banc-spectral-clustering.R.
  output_dirs <- c(
    "inst/images/cns_network",
    "/n/data1/hms/neurobio/wilson/banc/BANC-project/images/cns_network"
  )
  for (od in output_dirs) {
    if (!dir.exists(od)) dir.create(od, recursive = TRUE)
  }

  .save_all <- function(plot, basename, ...) {
    for (od in output_dirs) {
      ggplot2::ggsave(plot = plot, filename = file.path(od, basename), ...)
    }
  }

  dataset_version <- as.integer(banc.version)
  # New: spectral CSV filename now carries _v2 / _v3 suffix.
  clustering_file <- file.path(
    "data/cns_network",
    sprintf("spectral_clustering_min_connection_strength_%d_banc_version_%d_cluster_count_%d_cluster_seed_10_embedding_seed_3_%s.csv",
            min_connection_strength, dataset_version, cluster_count, syn_source)
  )

  if (!file.exists(clustering_file)) {
    stop("Clustering CSV not found: ", clustering_file,
         "\nRun banc-spectral-clustering.py or .R first.")
  }

  # ---------------------------------------------------------------
  # 2. Load clustering results
  # ---------------------------------------------------------------
  message("Loading clustering results...")
  clustering <- readr::read_csv(clustering_file,
                                col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(spectral_cluster = as.integer(spectral_cluster))

  cluster_ids <- unique(clustering$root_id)
  message(sprintf("  %d neurons in %d clusters",
                  length(cluster_ids),
                  length(unique(clustering$spectral_cluster))))

  # ---------------------------------------------------------------
  # 3. Load synapse data
  # ---------------------------------------------------------------
  message("Loading synapse data...")

  # Try current-version synapses first, fall back to older revisions
  syn_path <- NULL
  syn_candidates <- c(
    file.path(banc.save.path, paste0("v", banc.version),
              "synapses_v2_human_readable.csv"),
    file.path(banc.save.path, "v850", "synapses_v2_human_readable.csv"),
    file.path(banc.save.path, "synapses_250226_human_readable.csv"),
    bancsynapses
  )
  for (sp in syn_candidates) {
    if (file.exists(sp)) { syn_path <- sp; break }
  }
  if (is.null(syn_path)) {
    stop("Could not find synapse data. Checked: ",
         paste(syn_candidates, collapse = ", "))
  }
  message(sprintf("  Using: %s", syn_path))

  column_names <- c("id", "pre_x", "pre_y", "pre_z",
                     "post_x", "post_y", "post_z",
                     "ctr_x", "ctr_y", "ctr_z", "size",
                     "pre_supervoxel_id", "pre_root_id",
                     "post_supervoxel_id", "post_root_id")
  desired_columns <- c("id", "size", "pre_root_id", "post_root_id",
                        "ctr_x", "ctr_y", "ctr_z")
  col_types <- readr::cols(
    id = readr::col_character(),
    size = readr::col_double(),
    pre_root_id = readr::col_character(),
    post_root_id = readr::col_character(),
    ctr_x = readr::col_double(),
    ctr_y = readr::col_double(),
    ctr_z = readr::col_double(),
    .default = readr::col_double()
  )

  banc.syns <- vroom::vroom(
    syn_path,
    col_names = column_names,
    col_select = dplyr::all_of(desired_columns),
    col_types = col_types,
    skip = 1
  ) %>%
    dplyr::rename(X = ctr_x, Y = ctr_y, Z = ctr_z) %>%
    dplyr::filter(
      pre_root_id %in% cluster_ids | post_root_id %in% cluster_ids,
      pre_root_id != post_root_id,
      size > 5
    )
  message(sprintf("  %d synapses loaded (filtered to cluster neurons, size > 5)",
                  nrow(banc.syns)))

  # ---------------------------------------------------------------
  # 4. Template plots
  # ---------------------------------------------------------------
  message("Setting up plot templates...")

  g.anat <- ggplot2::ggplot() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::guides(fill = "none", color = "none") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0, size = 8,
                                          face = "bold", colour = "black"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.spacing = ggplot2::unit(0, "cm"),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "")

  # Neuropil surface base layers
  g.base.main <- g.anat +
    geom_neuron(x = banc_neuropil.surf,
                cols = c("grey60", "grey30"),
                rotation_matrix = bancr:::banc_rotation_matrices[["main"]],
                alpha = 0.1)
  g.base.brain <- g.anat +
    geom_neuron(x = banc_brain_neuropil.surf,
                cols = c("grey60", "grey30"),
                rotation_matrix = bancr:::banc_rotation_matrices[["front"]],
                alpha = 0.1)
  g.base.vnc <- g.anat +
    geom_neuron(x = banc_vnc_neuropil.surf,
                cols = c("grey60", "grey30"),
                rotation_matrix = bancr:::banc_rotation_matrices[["vnc"]],
                alpha = 0.1)

  # ---------------------------------------------------------------
  # 5. Per-cluster synapse density plots
  # ---------------------------------------------------------------
  clusters <- sort(unique(clustering$spectral_cluster))
  message(sprintf("Generating plots for %d clusters...", length(clusters)))

  for (cl in clusters) {
    try({
      message(sprintf("  Cluster %d...", cl))

      # Get root_ids for this cluster
      cl_ids <- clustering$root_id[clustering$spectral_cluster == cl]

      # Filter synapses to this cluster's neurons
      cl_syns <- banc.syns %>%
        dplyr::filter(
          pre_root_id %in% cl_ids | post_root_id %in% cl_ids
        )

      # Sample up to 1e5 synapses
      if (nrow(cl_syns) > 1e5) {
        cl_syns <- cl_syns %>% dplyr::slice_sample(n = 1e5)
      }

      if (nrow(cl_syns) < 100) {
        message(sprintf("    Skipping cluster %d: only %d synapses", cl, nrow(cl_syns)))
        next
      }

      # Get synapse positions and filter to inside neuropil
      pts <- nat::xyzmatrix(cl_syns)
      inside <- nat::pointsinside(pts, banc_neuropil.surf)
      pts <- pts[inside, ]

      if (nrow(pts) < 100) {
        message(sprintf("    Skipping cluster %d: too few points inside neuropil", cl))
        next
      }

      # --- Neuropil distribution bar chart ---
      message(sprintf("    Classifying %d synapses into neuropils...", nrow(pts)))
      chunk_size <- 10000
      neurons.np <- as.data.frame(pts) %>%
        dplyr::mutate(id = dplyr::row_number(),
                      neuropil = NA, region = NA, side = NA,
                      chunk = ceiling(dplyr::row_number() / chunk_size)) %>%
        dplyr::group_by(chunk) %>%
        dplyr::group_split()
      neurons.np <- purrr::map(neurons.np, function(chunk) {
        data <- tryCatch(pointsinside_banc(chunk), error = function(e) chunk)
        data <- tryCatch(pointsnearby_banc(data), error = function(e) data)
        data
      })
      neurons.np.df <- dplyr::bind_rows(neurons.np) %>%
        dplyr::mutate(
          neuropil = gsub("ITO_optic_|ITO_midbrain_|COURT_vnc_|_L|_R|_right|_left|MANC_.*|\\,.*", "", neuropil),
          neuropil = gsub("MB_.*", "MB", neuropil)
        ) %>%
        dplyr::arrange(region, neuropil)

      top_neuropils <- neurons.np.df %>%
        dplyr::count(neuropil, sort = TRUE) %>%
        dplyr::slice_head(n = 3) %>%
        dplyr::pull(neuropil)

      data_plot <- neurons.np.df %>%
        dplyr::mutate(
          neuropil_plot = ifelse(neuropil %in% top_neuropils, neuropil, "other"),
          region_plot   = ifelse(neuropil %in% top_neuropils, region, "other")
        ) %>%
        dplyr::count(neuropil_plot, region_plot) %>%
        dplyr::mutate(percent = n / sum(n) * 100)

      g.np <- ggplot2::ggplot(data_plot,
                              ggplot2::aes(x = reorder(neuropil_plot, n),
                                           y = percent, fill = region_plot)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = "", y = "", fill = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
      if (exists("paper.cols")) {
        g.np <- g.np + ggplot2::scale_fill_manual(values = paper.cols)
      }

      .save_all(g.np, sprintf("%d_neuropil_bar.png", cl),
                width = 2, height = 2, dpi = 300)

      # --- Main (full CNS) view ---
      x_main <- as.data.frame(
        t(bancr:::banc_rotation_matrices[["main"]][, 1:3] %*% t(as.matrix(pts)))
      )
      x_main <- x_main[, -4, drop = FALSE]
      colnames(x_main) <- c("X", "Y", "Z")

      g.main <- g.base.main +
        ggplot2::stat_density_2d(
          data = x_main,
          ggplot2::aes(x = X, y = Y, fill = ggplot2::after_stat(level)),
          n = 100, geom = "polygon", alpha = 0.5
        ) +
        ggplot2::scale_fill_viridis_c(option = "C") +
        ggplot2::theme_void() +
        ggplot2::guides(fill = "none", color = "none")

      .save_all(g.main, sprintf("%d_main_neuroanatomy.png", cl),
                width = 10, height = 10, dpi = 300)

      # --- Brain view ---
      pts_brain <- banc_decapitate(pts, invert = TRUE, OmitFailures = TRUE)
      if (isTRUE(NROW(pts_brain) > 50)) {
        x_brain <- as.data.frame(
          t(bancr:::banc_rotation_matrices[["front"]][, 1:3] %*% t(as.matrix(pts_brain)))
        )
        x_brain <- x_brain[, -4, drop = FALSE]
        colnames(x_brain) <- c("X", "Y", "Z")

        g.brain <- g.base.brain +
          ggplot2::stat_density_2d(
            data = x_brain,
            ggplot2::aes(x = X, y = Y, fill = ggplot2::after_stat(level)),
            n = 100, geom = "polygon", alpha = 0.5
          ) +
          ggplot2::scale_fill_viridis_c(option = "C") +
          ggplot2::theme_void() +
          ggplot2::guides(fill = "none", color = "none")

        .save_all(g.brain, sprintf("%d_brain_neuroanatomy.png", cl),
                  width = 10, height = 10, dpi = 300)
      }

      # --- VNC view ---
      pts_vnc <- banc_decapitate(pts, invert = FALSE, OmitFailures = TRUE)
      if (isTRUE(NROW(pts_vnc) > 50)) {
        x_vnc <- as.data.frame(
          t(bancr:::banc_rotation_matrices[["vnc"]][, 1:3] %*% t(as.matrix(pts_vnc)))
        )
        x_vnc <- x_vnc[, -4, drop = FALSE]
        colnames(x_vnc) <- c("X", "Y", "Z")

        g.vnc <- g.base.vnc +
          ggplot2::stat_density_2d(
            data = x_vnc,
            ggplot2::aes(x = X, y = Y, fill = ggplot2::after_stat(level)),
            n = 100, geom = "polygon", alpha = 0.5
          ) +
          ggplot2::scale_fill_viridis_c(option = "C") +
          ggplot2::theme_void() +
          ggplot2::guides(fill = "none", color = "none")

        .save_all(g.vnc, sprintf("%d_vnc_neuroanatomy.png", cl),
                  width = 10, height = 10, dpi = 300)
      }

      message(sprintf("    Cluster %d done", cl))
    })
  }

  # ---------------------------------------------------------------
  # 6. Per-cluster super class composition (one PDF each)
  # ---------------------------------------------------------------
  message("Generating per-cluster super_class composition plots...")

  # Load meta for super_class
  data_candidates <- c(
    file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
              paste0("banc_", dataset_version)),
    file.path("lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data",
              paste0("banc_", dataset_version))
  )
  data_dir <- NULL
  for (d in data_candidates) {
    if (dir.exists(d)) { data_dir <- d; break }
  }
  meta <- arrow::read_feather(
    file.path(data_dir, paste0("banc_", dataset_version, "_meta.feather"))
  )
  id_col <- paste0("root_", dataset_version)
  meta$root_id <- meta[[id_col]]

  # Join super_class to clustering results
  cl_sc <- clustering %>%
    dplyr::left_join(
      meta %>% dplyr::select(root_id, super_class) %>%
        dplyr::distinct(root_id, .keep_all = TRUE),
      by = "root_id"
    ) %>%
    dplyr::mutate(
      super_class = ifelse(is.na(super_class) | super_class == "",
                           "unknown", super_class)
    )

  for (cl in clusters) {
    try({
      cl_data <- cl_sc %>% dplyr::filter(spectral_cluster == cl)
      sc_counts <- cl_data %>% dplyr::count(super_class, sort = TRUE)
      n_total <- nrow(cl_data)

      # Build cluster title
      if ("cns_network" %in% names(clustering)) {
        cl_label <- clustering %>%
          dplyr::filter(spectral_cluster == cl) %>%
          dplyr::pull(cns_network) %>%
          unique() %>%
          paste(collapse = "/")
        cl_title <- sprintf("%s (cluster %d)", cl_label, cl)
      } else {
        cl_title <- sprintf("Cluster %d", cl)
      }

      g.sc <- ggplot2::ggplot(sc_counts,
                               ggplot2::aes(x = cl_title, y = n, fill = super_class)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_text(
          ggplot2::aes(x = cl_title, y = n_total),
          label = n_total, inherit.aes = FALSE,
          vjust = -0.3, size = 4
        ) +
        ggplot2::labs(x = "", y = "Number of neurons", fill = "Super class") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(size = 10, face = "bold"),
          legend.position = "right"
        )

      .save_all(g.sc, sprintf("%d_super_class.pdf", cl),
                width = 5, height = 5)
    })
  }
  message(sprintf("  Saved %d per-cluster super_class PDFs", length(clusters)))

  message(sprintf("### banc: Cluster visualisation complete [%s] ###",
                  format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
})
