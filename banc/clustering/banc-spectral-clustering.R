#' banc-spectral-clustering — R port of the BANC CNS-network spectral clustering.
#'
#' Column-normalised symmetrised adjacency → normalised Laplacian → KMeans
#' on bottom-k eigenvectors → UMAP. Cluster labels assigned by majority vote
#' against existing SeaTable `cns_network` annotations.
#'
#' @section Reads:
#'   - `banc_<ver>_meta.feather`, `banc_<ver>_edgelist_simple_<src>.feather`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `data/cns_network/spectral_clustering_..._<src>.csv`
#'
#' @section CLI:
#'   --source {v2,v3} --min-connection-strength N --cluster-count K --cluster-seed N --embedding-seed N
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   BANC-project/R/figures/panels_cns_networks.R, panels_cns_network_diagram.R,
#'   panels_cns_network_analyses.R (read `.banc_spectral_csv` for Fig. 6);
#'   R/annotations/banc-spectral-cluster-update.R (pushes cluster labels to SeaTable).
#'
#' @section Schema:
#'   banc_888_cns_network_spectral_clustering_v2.md.
#'
#' @section Paper:
#'   Methods §"Spectral clustering".

###########################################################
### Spectral Clustering of BANC CNS Network
###
### R port of banc-spectral-clustering.py
###
### Performs spectral clustering on the neuron-level
### connectivity graph for central brain, VNC, neck
### connective, and visual neurons.
###
### Algorithm:
###   1. Load BANC meta + simple edgelist from feather (source=v2 or v3)
###   2. Quality filter (exclude glia, trachea, unproofread)
###   3. Assign clustering_set (central brain, VNC, neck, visual)
###   4. Filter edgelist, prune to strongly connected component
###   5. Build column-normalized, symmetrized adjacency matrix
###   6. Spectral clustering (normalized Laplacian, KMeans)
###   7. UMAP of spectral embedding
###   8. Assign cns_network labels from seatable
###   9. Save CSV (output suffix carries the source: _v2 or _v3)
###
### Parameters:
###   min_connection_strength = 1 (paper Methods §"Spectral clustering"; override via --min-connection-strength)
###   cluster_count = 13 (paper Methods §"Spectral clustering"; override via --cluster-count)
###   cluster_seed = 10
###   embedding_seed = 3
###   UMAP: n_neighbors=100, metric=cosine, min_dist=0
###
### Input:
###   - banc_{version}_meta.feather
###   - banc_{version}_edgelist_simple_{v2|v3}.feather
###
### Output:
###   - data/cns_network/spectral_clustering_..._banc_version_{version}_..._{v2|v3}.csv
###   - data/cns_network/seatable_cns_network_snapshot.csv
###
### Usage: Rscript banc/clustering/banc-spectral-clustering.R [--source v2|v3] [--min-connection-strength N]
###########################################################
source("banc/banc-startup.R")

local({
  # Source selection
  .parse_source <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    i <- which(args == "--source")
    if (length(i) == 1 && length(args) >= i + 1) return(tolower(args[i + 1]))
    env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
    if (!is.na(env) && nzchar(env)) return(tolower(env))
    if (exists("banc.synapse.source.default")) return(tolower(banc.synapse.source.default))
    "v3"
  }
  syn_source <- .parse_source()
  if (!syn_source %in% c("v2", "v3"))
    stop(sprintf("--source must be 'v2' or 'v3' (got '%s')", syn_source))

  .parse_int_flag <- function(flag, default) {
    args <- commandArgs(trailingOnly = TRUE)
    i <- which(args == flag)
    if (length(i) == 1 && length(args) >= i + 1) return(as.integer(args[i + 1]))
    default
  }

  message(sprintf("### banc: Spectral clustering v%s | source=%s ###",
                  banc.version, syn_source))
  t_start <- Sys.time()

  # ---------------------------------------------------------------
  # Parameters
  # ---------------------------------------------------------------
  min_connection_strength <- .parse_int_flag("--min-connection-strength", 1L)
  cluster_count <- .parse_int_flag("--cluster-count", 13L)
  cluster_seed <- 10
  embedding_seed <- 3
  umap_n_neighbors <- 100
  dataset_version <- as.integer(banc.version)
  id_col <- paste0("root_", dataset_version)
  # Save to bancpipeline AND BANC-project. The latter is the canonical
  # publication location; the former keeps the working copy alongside the
  # script for fast local iteration.
  output_dirs <- c(
    "data/cns_network",
    "/n/data1/hms/neurobio/wilson/banc/BANC-project/data/cns_network"
  )

  # ---------------------------------------------------------------
  # 1. Load data
  # ---------------------------------------------------------------
  message("Loading data...")

  # Auto-detect data directory
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
  if (is.null(data_dir)) stop("Could not find BANC ", dataset_version, " data directory")

  meta <- arrow::read_feather(
    file.path(data_dir, paste0("banc_", dataset_version, "_meta.feather"))
  )
  edgelist <- arrow::read_feather(
    file.path(data_dir, paste0("banc_", dataset_version,
                                "_edgelist_simple_", syn_source, ".feather"))
  )

  # Use the current version's root_{version} column as canonical ID, renamed to root_id
  meta$root_id <- meta[[id_col]]
  meta <- meta[!is.na(meta$root_id) & meta$root_id != "", ]
  message(sprintf("  Meta: %d neurons, Edgelist: %d edges",
                  nrow(meta), nrow(edgelist)))

  # ---------------------------------------------------------------
  # 2. Quality filter
  # ---------------------------------------------------------------
  message("Applying quality filter...")
  meta <- meta %>%
    dplyr::filter(
      !grepl("glia|trachea|not_a_neuron|merge|debris",
             super_class, ignore.case = TRUE),
      !grepl("GLIA|TRACHEA|NOT_A_NEURON|DEBRIS|MERGE|DELETE", status),
      proofread == "TRUE" | roughly_proofread == "TRUE"
    )
  message(sprintf("  After quality filter: %d neurons", nrow(meta)))

  # ---------------------------------------------------------------
  # 3. Assign clustering_set
  # ---------------------------------------------------------------
  # Identify neck_connective neurons via super_class (ascending/descending) —
  # backwards-compatible with the now-deprecated region=="neck_connective" filter.
  clustering_regions <- c("central_brain", "ventral_nerve_cord")
  excluded_sc <- c("sensory", "motor", "efferent", "afferent", "visceral")

  meta <- meta %>%
    dplyr::mutate(
      clustering_set = dplyr::case_when(
        super_class %in% c("visual_centrifugal", "visual_projection") ~ "visual",
        grepl("ascending|descending", super_class) ~ "neck_connective",
        region %in% clustering_regions &
          !grepl(paste(excluded_sc, collapse = "|"), super_class) ~ region,
        TRUE ~ NA_character_
      )
    )

  clustering_values <- c("central_brain", "neck_connective",
                         "ventral_nerve_cord", "visual")
  cluster_ids <- meta %>%
    dplyr::filter(clustering_set %in% clustering_values) %>%
    dplyr::pull(root_id) %>%
    unique()
  message(sprintf("  Neurons in clustering set: %d", length(cluster_ids)))

  # ---------------------------------------------------------------
  # 4. Filter edgelist + prune
  # ---------------------------------------------------------------
  message("Filtering edgelist...")
  el <- edgelist %>%
    dplyr::filter(
      count >= min_connection_strength,
      pre != post,
      pre %in% cluster_ids,
      post %in% cluster_ids
    )
  message(sprintf("  Edges after filter: %d", nrow(el)))

  # Iterative pruning to strongly connected component
  message("Pruning to strongly connected component...")
  prev_n <- -1
  while (TRUE) {
    both <- intersect(unique(el$pre), unique(el$post))
    if (length(both) == prev_n) break
    prev_n <- length(both)
    el <- el %>% dplyr::filter(pre %in% both, post %in% both)
    message(sprintf("  Pruning: %d neurons", prev_n))
  }

  neuron_ids <- sort(both)
  n <- length(neuron_ids)
  message(sprintf("  Final: %d neurons, %d edges", n, nrow(el)))

  # ---------------------------------------------------------------
  # 5. Build sparse adjacency matrix
  # ---------------------------------------------------------------
  message("Building adjacency matrix...")
  id_map <- setNames(seq_along(neuron_ids), neuron_ids)

  # Aggregate weights per pre-post pair
  el_agg <- el %>%
    dplyr::group_by(pre, post) %>%
    dplyr::summarise(weight = sum(count), .groups = "drop")

  i_idx <- id_map[el_agg$pre]
  j_idx <- id_map[el_agg$post]

  adj <- Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = el_agg$weight,
    dims = c(n, n),
    dimnames = list(neuron_ids, neuron_ids)
  )

  # Column-normalize
  col_sums <- Matrix::colSums(adj)
  col_sums[col_sums == 0] <- 1
  adj_norm <- adj %*% Matrix::Diagonal(n, 1 / col_sums)

  # Symmetrize
  adj_sym <- 0.5 * (adj_norm + Matrix::t(adj_norm))
  message("  Matrix built and symmetrized")

  # ---------------------------------------------------------------
  # 6. Spectral clustering
  # ---------------------------------------------------------------
  message(sprintf("Computing spectral clustering (k=%d)...", cluster_count))

  # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
  deg <- Matrix::rowSums(adj_sym)
  deg[deg == 0] <- 1
  deg_inv_sqrt <- 1 / sqrt(deg)
  D_inv_sqrt <- Matrix::Diagonal(n, deg_inv_sqrt)
  lap <- Matrix::Diagonal(n) - D_inv_sqrt %*% adj_sym %*% D_inv_sqrt

  # Bottom-k eigenvectors
  eig <- RSpectra::eigs_sym(lap, k = cluster_count, which = "SM")
  eigvec <- eig$vectors

  # Row-normalize
  row_norms <- sqrt(rowSums(eigvec^2))
  row_norms[row_norms == 0] <- 1
  embedding <- eigvec / row_norms

  # KMeans
  set.seed(cluster_seed)
  km <- kmeans(embedding, centers = cluster_count, nstart = 25, iter.max = 100)
  labels <- km$cluster  # 1-indexed

  for (cl in sort(unique(labels))) {
    message(sprintf("  Cluster %d: %d neurons", cl, sum(labels == cl)))
  }

  # ---------------------------------------------------------------
  # 7. UMAP embedding
  # ---------------------------------------------------------------
  message("Computing UMAP embedding...")
  umap_result <- uwot::umap(
    embedding,
    n_neighbors = umap_n_neighbors,
    metric = "cosine",
    min_dist = 0,
    n_components = 2,
    seed = embedding_seed,
    n_threads = 1
  )

  # ---------------------------------------------------------------
  # 8. Build result + assign cns_network labels
  # ---------------------------------------------------------------
  result <- data.frame(
    root_id = neuron_ids,
    spectral_cluster = labels,
    umap_x = umap_result[, 1],
    umap_y = umap_result[, 2],
    stringsAsFactors = FALSE
  )

  # Join supervoxel_id and position from meta
  meta_cols <- meta %>%
    dplyr::select(root_id, supervoxel_id, position) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  result <- result %>%
    dplyr::left_join(meta_cols, by = "root_id")

  # Get cns_network from seatable and save snapshot for Python
  message("Querying seatable for cns_network labels...")
  st <- tryCatch({
    banctable_query(paste0("SELECT ", id_col, ", cns_network FROM banc"))
  }, error = function(e) {
    message("  Seatable query failed, using cns_network from meta feather")
    meta %>% dplyr::select(root_id, cns_network)
  })

  # Normalise column names
  if (id_col %in% names(st)) {
    names(st)[names(st) == id_col] <- "root_id"
  }
  st <- st %>%
    dplyr::filter(!is.na(cns_network), cns_network != "") %>%
    dplyr::distinct(root_id, .keep_all = TRUE)

  # Save seatable snapshot to bancpipeline (working copy used by downstream
  # scripts and the heatmap). The snapshot is the OLD cns_network values used
  # for the cluster->label majority vote in step 9, so it pairs with the
  # spectral CSV emitted below.
  primary_dir <- output_dirs[1]
  if (!dir.exists(primary_dir)) dir.create(primary_dir, recursive = TRUE)
  snapshot_path <- file.path(primary_dir, "seatable_cns_network_snapshot.csv")
  readr::write_csv(
    st %>% dplyr::rename(!!id_col := root_id),
    snapshot_path
  )
  message(sprintf("  Saved seatable snapshot: %s (%d entries)", snapshot_path, nrow(st)))

  # Assign cns_network per cluster by plain per-cluster majority vote on the
  # OLD SeaTable cns_network labels (neurons missing a label are discounted).
  # Same OLD label may be assigned to multiple clusters when the network
  # genuinely splits a region (e.g., leg VNC into left/right). Clusters with
  # zero overlap with any OLD label keep their numeric cluster string.
  result <- result %>%
    dplyr::left_join(
      st %>% dplyr::select(root_id, existing_cns = cns_network),
      by = "root_id"
    )
  result$cns_network <- as.character(result$spectral_cluster)

  cluster_ids <- sort(unique(result$spectral_cluster))
  old_labels <- sort(unique(result$existing_cns[
    !is.na(result$existing_cns) & result$existing_cns != ""
  ]))

  if (length(old_labels) == 0L) {
    message("  No old cns_network labels available — using cluster numbers")
  } else {
    for (cl in cluster_ids) {
      mask <- result$spectral_cluster == cl
      ex <- result$existing_cns[mask]
      ex <- ex[!is.na(ex) & ex != ""]
      if (length(ex) == 0L) {
        message(sprintf("  Cluster %d -> '%s' (no labelled neurons — using cluster number)",
                        cl, as.character(cl)))
        next
      }
      tab <- sort(table(ex), decreasing = TRUE)
      lbl <- names(tab)[1]
      result$cns_network[mask] <- lbl
      message(sprintf("  Cluster %d -> '%s' (%d/%d neurons, %.1f%% of labelled)",
                      cl, lbl, tab[1], sum(mask),
                      100 * tab[1] / length(ex)))
    }
  }

  result <- result %>% dplyr::select(-existing_cns)

  # ---------------------------------------------------------------
  # 9. Save output
  # ---------------------------------------------------------------
  output_basename <- sprintf(
    "spectral_clustering_min_connection_strength_%d_banc_version_%d_cluster_count_%d_cluster_seed_%d_embedding_seed_%d_%s.csv",
    min_connection_strength, dataset_version, cluster_count,
    cluster_seed, embedding_seed, syn_source)
  # Dated copy preserves historical runs alongside the canonical (overwritten)
  # filename. Date is the calendar date the script started (UTC-naïve, fine for
  # day-granularity ordering).
  dated_basename <- sub("\\.csv$",
                        sprintf("_%s.csv", format(Sys.Date(), "%Y-%m-%d")),
                        output_basename)

  output_df <- result %>%
    dplyr::select(root_id, supervoxel_id, position,
                  spectral_cluster, umap_x, umap_y, cns_network)

  for (od in output_dirs) {
    if (!dir.exists(od)) dir.create(od, recursive = TRUE)
    for (.bn in c(output_basename, dated_basename)) {
      out_path <- file.path(od, .bn)
      readr::write_csv(output_df, out_path)
      message(sprintf("\nSaved: %s (%d neurons, %d clusters)",
                      out_path, nrow(output_df),
                      length(unique(output_df$spectral_cluster))))
    }
  }

  message(sprintf("### banc: Spectral clustering complete [%s] ###",
                  format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
})
