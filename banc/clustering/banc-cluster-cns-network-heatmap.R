#' banc-cluster-cns-network-heatmap — Cross-tabulate new spectral clusters vs old cns_network.
#'
#' Cross-tabulates each NEW `spectral_cluster` from the latest spectral
#' clustering CSV against the OLD `cns_network` labels snapshotted from
#' SeaTable. Emits two heatmaps per source (counts + row-normalised) and
#' saves to both bancpipeline and BANC-project.
#'
#' @section Reads:
#'   - `data/cns_network/spectral_clustering_..._{v2|v3}.csv`
#'   - `data/cns_network/seatable_cns_network_snapshot.csv`
#'
#' @section Writes:
#'   - `inst/images/cns_network/cluster_vs_old_{counts,frac}_<src>.pdf` (+ BANC-project mirror)
#'
#' @section CLI:
#'   --source {v2,v3}              default v3
#'   --min-connection-strength N   default 1
#'   --cluster-count K             default 13
#'
#' @section Paper:
#'   Methods §"Spectral clustering" (qualitative diagnostic).

###########################################################
### Heatmap: new spectral clusters vs. old SeaTable cns_network
###
### Cross-tabulates each NEW spectral_cluster from the latest
### spectral clustering CSV against the OLD cns_network labels
### snapshotted from SeaTable (saved by
### banc-spectral-clustering.R as
### data/cns_network/seatable_cns_network_snapshot.csv).
###
### Two heatmaps are produced per source (v2|v3):
###   - counts: raw neuron counts per (cluster, old label)
###   - row-normalised: fraction of each cluster going to each label
###
### Saved to BOTH bancpipeline (working copy) and BANC-project.
###
### Usage:
###   Rscript banc/clustering/banc-cluster-cns-network-heatmap.R \
###     [--source v2|v3] [--min-connection-strength N] [--cluster-count K]
###########################################################
source("banc/banc-startup.R")

local({
  message("### banc: cluster x cns_network heatmap ###")

  .cli_args <- commandArgs(trailingOnly = TRUE)
  .parse_str_flag <- function(flag, default) {
    j <- which(.cli_args == flag)
    if (length(j) == 1 && length(.cli_args) >= j + 1) return(tolower(.cli_args[j + 1]))
    default
  }
  .parse_int_flag <- function(flag, default) {
    j <- which(.cli_args == flag)
    if (length(j) == 1 && length(.cli_args) >= j + 1) return(as.integer(.cli_args[j + 1]))
    default
  }

  syn_source <- .parse_str_flag("--source", "v3")
  if (!syn_source %in% c("v2", "v3"))
    stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", syn_source))
  min_connection_strength <- .parse_int_flag("--min-connection-strength", 1L)
  cluster_count <- .parse_int_flag("--cluster-count", 14L)
  dataset_version <- as.integer(banc.version)

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

  clustering_file <- file.path(
    "data/cns_network",
    sprintf("spectral_clustering_min_connection_strength_%d_banc_version_%d_cluster_count_%d_cluster_seed_10_embedding_seed_3_%s.csv",
            min_connection_strength, dataset_version, cluster_count, syn_source)
  )
  snapshot_file <- "data/cns_network/seatable_cns_network_snapshot.csv"
  for (f in c(clustering_file, snapshot_file)) {
    if (!file.exists(f)) stop("Required input not found: ", f)
  }

  clustering <- readr::read_csv(clustering_file,
                                col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(spectral_cluster = as.integer(spectral_cluster)) %>%
    dplyr::select(root_id, spectral_cluster)

  # The snapshot uses root_{version} as its ID column, mirror back to root_id
  id_col <- paste0("root_", dataset_version)
  snapshot <- readr::read_csv(snapshot_file,
                              col_types = readr::cols(.default = "c"))
  if (!id_col %in% names(snapshot) && "root_id" %in% names(snapshot)) {
    snapshot <- snapshot %>% dplyr::rename(!!id_col := root_id)
  }
  snapshot <- snapshot %>%
    dplyr::rename(root_id = !!id_col) %>%
    dplyr::filter(!is.na(cns_network), cns_network != "")

  joined <- clustering %>%
    dplyr::inner_join(snapshot %>% dplyr::select(root_id, old_cns_network = cns_network),
                      by = "root_id")
  message(sprintf("  Joined: %d neurons (clustering %d, snapshot %d)",
                  nrow(joined), nrow(clustering), nrow(snapshot)))

  counts <- joined %>%
    dplyr::count(spectral_cluster, old_cns_network, name = "n") %>%
    dplyr::group_by(spectral_cluster) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()

  # Order old labels by total neurons assigned (largest first)
  old_levels <- counts %>%
    dplyr::group_by(old_cns_network) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(old_cns_network)

  counts <- counts %>%
    dplyr::mutate(
      old_cns_network = factor(old_cns_network, levels = old_levels),
      spectral_cluster = factor(spectral_cluster,
                                levels = sort(unique(spectral_cluster)))
    )

  base_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  g_counts <- ggplot2::ggplot(counts,
                              ggplot2::aes(x = old_cns_network,
                                           y = spectral_cluster, fill = n)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::geom_text(ggplot2::aes(label = n), size = 2.5) +
    ggplot2::scale_fill_viridis_c(option = "C", trans = "log10",
                                  na.value = "white") +
    ggplot2::labs(
      title = sprintf("New spectral_cluster x old cns_network (counts) | source=%s",
                      syn_source),
      x = "old cns_network (SeaTable snapshot)",
      y = "new spectral_cluster",
      fill = "neurons"
    ) +
    base_theme

  g_frac <- ggplot2::ggplot(counts,
                            ggplot2::aes(x = old_cns_network,
                                         y = spectral_cluster, fill = frac)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", frac)), size = 2.5) +
    ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1),
                                  na.value = "white") +
    ggplot2::labs(
      title = sprintf("New spectral_cluster x old cns_network (row-normalised) | source=%s",
                      syn_source),
      x = "old cns_network (SeaTable snapshot)",
      y = "new spectral_cluster",
      fill = "fraction of cluster"
    ) +
    base_theme

  n_cols <- length(old_levels)
  n_rows <- length(unique(counts$spectral_cluster))
  w <- max(6, 0.45 * n_cols + 2)
  h <- max(4, 0.45 * n_rows + 2)

  base <- sprintf(
    "spectral_cluster_vs_cns_network_heatmap_min%d_banc%d_k%d_%s",
    min_connection_strength, dataset_version, cluster_count, syn_source)

  .save_all(g_counts, paste0(base, "_counts.pdf"), width = w, height = h)
  .save_all(g_counts, paste0(base, "_counts.png"), width = w, height = h, dpi = 200)
  .save_all(g_frac,   paste0(base, "_frac.pdf"),   width = w, height = h)
  .save_all(g_frac,   paste0(base, "_frac.png"),   width = w, height = h, dpi = 200)

  message("### heatmap done ###")
})
