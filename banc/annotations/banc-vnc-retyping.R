#' banc-vnc-retyping — Cluster + retype female-specific efferent VNC neurons.
#'
#' Clusters ALL female-specific efferent VNC neurons (typed + untyped)
#' using NBLAST + upstream connectivity. Already-typed neurons act as
#' anchors that define cluster boundaries; clusters with a dominant
#' existing type inherit that name; clusters without get new `MNadfXX`
#' names. Only untyped neurons are assigned.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - BANC mesh / skeleton stores; CAVE upstream partners
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `cell_type` (untyped neurons only)
#'   - diagnostic plots under `inst/images/vnc_retyping/`

###########################################################
### VNC retyping: female-specific efferent neurons
###
### Clusters ALL female-specific efferent VNC neurons
### (typed + untyped) using NBLAST + upstream connectivity.
### Already-typed neurons act as reference anchors that
### define cluster boundaries. Clusters with a dominant
### existing type inherit that name; clusters without get
### new MNadfXX names. Only untyped neurons are assigned.
###
###   1. Download neuron meshes, convert to dotprops
###   2. Run all-by-all NBLAST (symmetric)
###   3. Fetch upstream connectivity, build type-level matrix
###   4. Compute cosine similarity on connectivity profiles
###   5. Combine NBLAST + connectivity into joint distance
###   6. Hierarchical clustering with typed anchors
###   7. Assign types: existing type or MNadfXX
###
### Output: cluster assignments + diagnostic plots
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: VNC retyping — female-specific efferent neurons ###")
t_start <- Sys.time()

###########################
### 1. Query SeaTable   ###
###########################

bc <- banctable_query("SELECT _id, root_id, root_626, supervoxel_id, super_class, flow, cell_type, side, region, sexually_dimorphic, status FROM banc_meta") %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::mutate(root_626 = as.character(root_626),
                root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE) %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON|DEBRIS", status, ignore.case = TRUE) | is.na(status))

# ALL female-specific efferent VNC neurons (typed + untyped)
all_fem_eff <- bc %>%
  dplyr::filter(sexually_dimorphic == "female-specific",
                flow == "efferent",
                grepl("ventral_nerve_cord", region)) %>%
  dplyr::mutate(is_typed = !is.na(cell_type) & cell_type != "")

n_typed <- sum(all_fem_eff$is_typed)
n_untyped <- sum(!all_fem_eff$is_typed)
message(sprintf("  %d female-specific efferent VNC neurons (%d typed, %d untyped)",
                nrow(all_fem_eff), n_typed, n_untyped))

if (n_untyped == 0) {
  message("  All neurons already typed. Nothing to do.")
  return(invisible())
}

if (nrow(all_fem_eff) < 3) {
  message("  Too few neurons for clustering. Stopping.")
  return(invisible())
}

all_ids <- all_fem_eff$root_id

###########################
### 2. Neuron meshes    ###
###########################

message("  Downloading neuron meshes...")

dps <- tryCatch({
  meshes <- banc_read_neuron_meshes(all_ids, OmitFailures = TRUE)
  nat::dotprops(meshes / 1e3, k = 20, resample = FALSE, OmitFailures = TRUE)
}, error = function(e) {
  warning("  Mesh download failed: ", e$message)
  NULL
})

if (is.null(dps) || length(dps) == 0) {
  message("  No dotprops available. Stopping.")
  return(invisible())
}

message(sprintf("  %d dotprops available for NBLAST", length(dps)))

# Filter to only those with dotprops
valid_ids <- names(dps)
all_fem_eff <- all_fem_eff %>% dplyr::filter(root_id %in% valid_ids)
n_typed <- sum(all_fem_eff$is_typed)
n_untyped <- sum(!all_fem_eff$is_typed)

message(sprintf("  After mesh filter: %d typed, %d untyped", n_typed, n_untyped))

if (n_untyped == 0) {
  message("  All neurons with meshes are already typed. Stopping.")
  return(invisible())
}

# Mirror left-side neurons to right so NBLAST is side-invariant
left_ids <- all_fem_eff$root_id[which(all_fem_eff$side == "left")]
left_ids <- intersect(left_ids, names(dps))
if (length(left_ids) > 0) {
  message(sprintf("  Mirroring %d left-side neurons to right...", length(left_ids)))
  dps[left_ids] <- nat::nlapply(dps[left_ids], function(x) {
    banc_mirror(x, banc.units = "um")
  }, OmitFailures = TRUE)
}

###########################
### 3. All-by-all NBLAST ###
###########################

message("  Running all-by-all NBLAST...")

nb_raw <- nat.nblast::nblast(query = dps, target = dps,
                              UseAlpha = TRUE,
                              normalised = TRUE,
                              smat = nat.nblast::smat_alpha.fcwb)

# Symmetric average: (M + M') / 2
nb_sym <- (nb_raw + t(nb_raw)) / 2

# Convert to distance
nblast_dist <- as.dist(1 - nb_sym)

message(sprintf("  NBLAST done: %dx%d matrix", nrow(nb_sym), ncol(nb_sym)))

###########################
### 4. Upstream connectivity ###
###########################

message("  Fetching upstream connectivity...")

valid_ids <- all_fem_eff$root_id
upstream_list <- list()
for (i in seq_along(valid_ids)) {
  rid <- valid_ids[i]
  if (i %% 10 == 0 || i == 1) message(sprintf("    [%d/%d] %s", i, length(valid_ids), rid))

  ups <- tryCatch(
    banc_partner_summary(rid, partners = "inputs", threshold = 1),
    error = function(e) {
      warning(sprintf("    Failed upstream for %s: %s", rid, e$message))
      NULL
    }
  )

  if (!is.null(ups) && nrow(ups) > 0) {
    ups$target_id <- rid
    upstream_list[[rid]] <- ups
  }
}

upstream_df <- dplyr::bind_rows(upstream_list)
message(sprintf("  %d upstream connections for %d neurons",
                nrow(upstream_df), length(unique(upstream_df$target_id))))

# Map upstream partner IDs to cell types from seatable
partner_ids <- unique(as.character(upstream_df$partner))
partner_types <- bc %>%
  dplyr::filter(root_id %in% partner_ids) %>%
  dplyr::select(root_id, partner_type = cell_type) %>%
  dplyr::mutate(partner_type = dplyr::if_else(
    is.na(partner_type) | partner_type == "", "untyped", partner_type
  ))

upstream_df <- upstream_df %>%
  dplyr::mutate(partner = as.character(partner)) %>%
  dplyr::left_join(partner_types, by = c("partner" = "root_id")) %>%
  dplyr::mutate(partner_type = dplyr::if_else(
    is.na(partner_type), "untyped", partner_type
  ))

# Build neuron x upstream_type matrix (synapse-weighted)
conn_matrix <- upstream_df %>%
  dplyr::group_by(target_id, partner_type) %>%
  dplyr::summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = partner_type, values_from = weight,
                     values_fill = 0) %>%
  tibble::column_to_rownames("target_id") %>%
  as.matrix()

# Ensure all neurons are in the matrix
missing <- setdiff(valid_ids, rownames(conn_matrix))
if (length(missing) > 0) {
  empty_mat <- matrix(0, nrow = length(missing), ncol = ncol(conn_matrix),
                      dimnames = list(missing, colnames(conn_matrix)))
  conn_matrix <- rbind(conn_matrix, empty_mat)
}
conn_matrix <- conn_matrix[valid_ids, , drop = FALSE]

message(sprintf("  Connectivity matrix: %d neurons x %d upstream types",
                nrow(conn_matrix), ncol(conn_matrix)))

###########################
### 5. Cosine similarity ###
###########################

cosine_sim <- function(mat) {
  norms <- sqrt(rowSums(mat^2))
  norms[norms == 0] <- 1
  mat_norm <- mat / norms
  sim <- tcrossprod(mat_norm)
  diag(sim) <- 1
  sim
}

conn_sim <- cosine_sim(conn_matrix)
conn_dist <- as.dist(1 - conn_sim)

message("  Cosine similarity computed")

###########################
### 6. Combined clustering ###
###   with typed anchors  ###
###########################

message("  Clustering...")

combined_dist <- (nblast_dist + conn_dist) / 2

hc <- hclust(combined_dist, method = "ward.D2")

# --- Determine optimal k using typed neurons as guides ---
# Strategy: sweep k values and pick the one that best preserves
# typed neurons of the same cell_type within the same cluster,
# while also having good silhouette scores.
#
# "Anchor purity" = for each existing type with >=2 neurons in the
# clustering, what fraction are in the same cluster?

typed_df <- all_fem_eff %>%
  dplyr::filter(is_typed) %>%
  dplyr::select(root_id, cell_type)

# Only consider types with >= 2 representatives (can't measure cohesion for singletons)
type_counts <- table(typed_df$cell_type)
multi_types <- names(type_counts[type_counts >= 2])

max_k <- min(floor(nrow(all_fem_eff) / 2), 50)
if (max_k < 2) max_k <- 2

k_metrics <- data.frame(k = 2:max_k, silhouette = NA_real_, anchor_purity = NA_real_)

for (ki in seq_len(nrow(k_metrics))) {
  k <- k_metrics$k[ki]
  cl <- cutree(hc, k = k)

  # Silhouette
  if (length(unique(cl)) >= 2) {
    sil <- cluster::silhouette(cl, combined_dist)
    k_metrics$silhouette[ki] <- mean(sil[, 3])
  }

  # Anchor purity: for each multi-type, what fraction of its neurons
  # are in the most common cluster?
  if (length(multi_types) > 0) {
    purities <- sapply(multi_types, function(tp) {
      ids <- typed_df$root_id[typed_df$cell_type == tp]
      ids <- intersect(ids, names(cl))
      if (length(ids) < 2) return(1)
      cl_ids <- cl[ids]
      max(table(cl_ids)) / length(cl_ids)
    })
    k_metrics$anchor_purity[ki] <- mean(purities)
  } else {
    k_metrics$anchor_purity[ki] <- 1
  }
}

# Combined score: balance silhouette with anchor purity
# Anchor purity should be high (don't split known types across clusters)
k_metrics$combined <- k_metrics$silhouette * k_metrics$anchor_purity

best_k <- k_metrics$k[which.max(k_metrics$combined)]
best_row <- k_metrics[k_metrics$k == best_k, ]
message(sprintf("  Best k: %d (silhouette=%.3f, anchor_purity=%.3f, combined=%.3f)",
                best_k, best_row$silhouette, best_row$anchor_purity, best_row$combined))

# Cut tree at best k
clusters <- cutree(hc, k = best_k)

###########################
### 7. Assign type names ###
###########################

# For each cluster, check if it has a dominant existing cell_type
cluster_types <- data.frame(
  root_id = names(clusters),
  cluster = unname(clusters),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(all_fem_eff %>% dplyr::select(root_id, cell_type, is_typed),
                   by = "root_id")

# Per-cluster type assignment
cluster_labels <- cluster_types %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_typed = sum(is_typed),
    n_untyped = sum(!is_typed),
    # Dominant type: most common non-empty cell_type in the cluster
    dominant_type = {
      typed_in_cluster <- cell_type[is_typed]
      if (length(typed_in_cluster) == 0) NA_character_
      else names(sort(table(typed_in_cluster), decreasing = TRUE))[1]
    },
    dominant_count = {
      typed_in_cluster <- cell_type[is_typed]
      if (length(typed_in_cluster) == 0) 0L
      else as.integer(max(table(typed_in_cluster)))
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # A cluster inherits the existing type if typed neurons form a majority
    # or if it's the only type in the cluster
    use_existing = n_typed > 0 & (dominant_count >= n_typed * 0.5)
  )

# Assign MNadfXX numbers only to clusters without a dominant existing type
novel_clusters <- cluster_labels %>%
  dplyr::filter(!use_existing) %>%
  dplyr::arrange(cluster)

pad_width <- max(2, nchar(as.character(nrow(novel_clusters))))
novel_clusters$assigned_type <- sprintf("MNadf%0*d", pad_width, seq_len(nrow(novel_clusters)))

# Merge: existing-type clusters get the dominant type, novel clusters get MNadfXX
cluster_labels <- cluster_labels %>%
  dplyr::left_join(novel_clusters %>% dplyr::select(cluster, assigned_type),
                   by = "cluster") %>%
  dplyr::mutate(
    final_type = dplyr::if_else(use_existing, dominant_type, assigned_type)
  )

message(sprintf("  %d clusters with existing types, %d novel (MNadf) types",
                sum(cluster_labels$use_existing), sum(!cluster_labels$use_existing)))

# Build result: only untyped neurons get assignments
result <- cluster_types %>%
  dplyr::filter(!is_typed) %>%
  dplyr::left_join(cluster_labels %>% dplyr::select(cluster, final_type),
                   by = "cluster") %>%
  dplyr::left_join(all_fem_eff %>% dplyr::select(root_id, `_id`, root_626,
                                                   supervoxel_id, side, super_class),
                   by = "root_id") %>%
  dplyr::select(`_id`, root_id, root_626, supervoxel_id, side, super_class,
                cluster, new_cell_type = final_type)

# Summary
message("  Cluster assignments:")
for (i in seq_len(nrow(cluster_labels))) {
  cl <- cluster_labels[i, ]
  message(sprintf("    Cluster %d -> %s: %d typed + %d untyped = %d total%s",
                  cl$cluster, cl$final_type, cl$n_typed, cl$n_untyped, cl$n_total,
                  if (cl$use_existing) sprintf(" (anchor: %s)", cl$dominant_type) else " (novel)"))
}

###########################
### 8. Diagnostic plots ###
###########################

# Prepare label info for dendrogram
label_df <- cluster_types %>%
  dplyr::left_join(cluster_labels %>% dplyr::select(cluster, final_type, use_existing),
                   by = "cluster")

# Plot 1: Dendrogram with typed neurons highlighted
png("inst/images/vnc_retyping_dendrogram.png", width = 1400, height = 700, res = 150)
dend <- as.dendrogram(hc)
dend_order <- order.dendrogram(dend)
ordered_ids <- labels(dend)

# Color by cluster
n_clusters <- best_k
cluster_colors <- rainbow(n_clusters, s = 0.7, v = 0.9)
leaf_clusters <- clusters[ordered_ids]
leaf_colors <- cluster_colors[leaf_clusters]

# Bold/shape for typed vs untyped
leaf_typed <- all_fem_eff$is_typed[match(ordered_ids, all_fem_eff$root_id)]
leaf_colors[leaf_typed] <- adjustcolor(leaf_colors[leaf_typed], red.f = 0.7,
                                        green.f = 0.7, blue.f = 0.7)

dendextend::labels_colors(dend) <- leaf_colors
plot(dend,
     main = sprintf("Female-specific efferent VNC (n=%d: %d typed + %d untyped, k=%d)",
                     nrow(all_fem_eff), n_typed, n_untyped, best_k),
     ylab = "Combined distance (NBLAST + connectivity)",
     leaflab = "none")
rect.hclust(hc, k = best_k, border = cluster_colors)

# Add legend for cluster -> type mapping
type_labels <- cluster_labels$final_type
legend("topright", legend = sprintf("C%d: %s (n=%d)",
                                     cluster_labels$cluster,
                                     type_labels,
                                     cluster_labels$n_total),
       fill = cluster_colors[cluster_labels$cluster],
       cex = 0.5, ncol = 2, bg = "white")
dev.off()
message("  Saved: vnc_retyping_dendrogram.png")

# Plot 2: k optimization (silhouette + anchor purity)
png("inst/images/vnc_retyping_k_optimization.png", width = 900, height = 500, res = 150)
par(mar = c(5, 4, 4, 4) + 0.1)
plot(k_metrics$k, k_metrics$silhouette, type = "b", pch = 19, col = "steelblue",
     xlab = "Number of clusters (k)", ylab = "Silhouette width",
     main = "Cluster optimization: silhouette + anchor purity",
     ylim = c(0, max(k_metrics$silhouette, na.rm = TRUE) * 1.1))
par(new = TRUE)
plot(k_metrics$k, k_metrics$anchor_purity, type = "b", pch = 17, col = "darkgreen",
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1.1))
axis(4, col = "darkgreen", col.axis = "darkgreen")
mtext("Anchor purity", side = 4, line = 3, col = "darkgreen")
abline(v = best_k, col = "red", lty = 2)
legend("bottomleft", legend = c("Silhouette", "Anchor purity", sprintf("Best k=%d", best_k)),
       col = c("steelblue", "darkgreen", "red"),
       pch = c(19, 17, NA), lty = c(1, 1, 2), cex = 0.8)
dev.off()
message("  Saved: vnc_retyping_k_optimization.png")

# Plot 3: NBLAST similarity heatmap
cluster_order <- order(clusters)
nb_ordered <- nb_sym[cluster_order, cluster_order]
typed_ordered <- all_fem_eff$is_typed[match(rownames(nb_ordered), all_fem_eff$root_id)]

png("inst/images/vnc_retyping_nblast_heatmap.png", width = 800, height = 800, res = 150)
image(nb_ordered, col = hcl.colors(50, "RdYlBu"),
      main = "NBLAST similarity (ordered by cluster)",
      xaxt = "n", yaxt = "n")
cumsum_clusters <- cumsum(table(clusters[cluster_order]))
boundaries <- cumsum_clusters / nrow(nb_ordered)
abline(v = boundaries, col = "black", lwd = 0.5)
abline(h = boundaries, col = "black", lwd = 0.5)
# Mark typed neurons on axes
typed_pos <- which(typed_ordered) / length(typed_ordered)
points(typed_pos, rep(-0.02, length(typed_pos)), pch = "|", col = "black", cex = 0.5, xpd = TRUE)
points(rep(-0.02, length(typed_pos)), typed_pos, pch = "-", col = "black", cex = 0.5, xpd = TRUE)
dev.off()
message("  Saved: vnc_retyping_nblast_heatmap.png")

# Plot 4: Connectivity cosine similarity heatmap
conn_ordered <- conn_sim[cluster_order, cluster_order]

png("inst/images/vnc_retyping_connectivity_heatmap.png", width = 800, height = 800, res = 150)
image(conn_ordered, col = hcl.colors(50, "RdYlBu"),
      main = "Upstream connectivity cosine similarity (ordered by cluster)",
      xaxt = "n", yaxt = "n")
abline(v = boundaries, col = "black", lwd = 0.5)
abline(h = boundaries, col = "black", lwd = 0.5)
points(typed_pos, rep(-0.02, length(typed_pos)), pch = "|", col = "black", cex = 0.5, xpd = TRUE)
points(rep(-0.02, length(typed_pos)), typed_pos, pch = "-", col = "black", cex = 0.5, xpd = TRUE)
dev.off()
message("  Saved: vnc_retyping_connectivity_heatmap.png")

###########################
### 9. Save results     ###
###########################

csv_path <- "data/codex/vnc_retyping_female_specific_efferent.csv"
readr::write_csv(result, csv_path)
message(sprintf("  Saved: %s (%d untyped neurons, %d types assigned)",
                csv_path, nrow(result), length(unique(result$new_cell_type))))

# Full cluster info including typed neurons (for reference)
full_cluster_csv <- cluster_types %>%
  dplyr::left_join(cluster_labels %>% dplyr::select(cluster, final_type, use_existing),
                   by = "cluster") %>%
  dplyr::left_join(all_fem_eff %>% dplyr::select(root_id, `_id`, root_626,
                                                   supervoxel_id, side, super_class),
                   by = "root_id") %>%
  dplyr::select(`_id`, root_id, root_626, side, super_class,
                existing_cell_type = cell_type, is_typed,
                cluster, assigned_type = final_type)

readr::write_csv(full_cluster_csv,
                 "data/codex/vnc_retyping_full_clustering.csv")
message("  Saved: vnc_retyping_full_clustering.csv (all neurons incl. typed anchors)")

# Save R objects for further analysis
rds_path <- "data/codex/vnc_retyping_objects.rds"
saveRDS(list(
  all_fem_eff = all_fem_eff,
  nblast_sim = nb_sym,
  conn_sim = conn_sim,
  conn_matrix = conn_matrix,
  combined_dist = combined_dist,
  hclust = hc,
  clusters = clusters,
  cluster_labels = cluster_labels,
  k_metrics = k_metrics,
  best_k = best_k
), rds_path)
message(sprintf("  Saved: %s", rds_path))

# SeaTable update (commented out — review first)
push_types <- result %>%
  dplyr::select(`_id`, cell_type = new_cell_type) %>%
  as.data.frame()

message(sprintf("\n  Ready to push %d cell_type assignments to SeaTable.", nrow(push_types)))
message("  Review the plots and CSV first. To push:")
message("  banctable_update_rows(base='banc_meta', table='banc_meta', df=push_types, append_allowed=FALSE, chunksize=200)")

message(sprintf("### Done [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
