# `banc/clustering/` — CNS-network spectral clustering (k = 13)

Connectivity-based spectral clustering of the BANC CNS network. Builds a column-normalised symmetrised adjacency matrix from the edgelist, computes the bottom-k eigenvectors of its normalised Laplacian, KMeans on the spectral embedding, and assigns human-readable labels (`CNS_01`..`CNS_13`) by majority vote against existing SeaTable annotations. Output is `banc_888_cns_network_spectral_clustering_v2.csv` in the public release. Paper Methods §"Spectral clustering".

Sits after `banc/metrics/banc-calculate-connectivity.R` (edgelist input) and before `banc/share/banc-data.R` (which packages the CSV into the release).

## Scripts

| File | Purpose |
|---|---|
| `banc-spectral-clustering.R` | R implementation: spectral clustering + UMAP + majority-vote labels against SeaTable `cns_network`. |
| `banc-spectral-clustering.py` | Python port for SLURM / off-O2 runs. Accepts `--seatable-csv` to consume the R-generated SeaTable snapshot. |
| `banc-spectral-clustering-run.R` | Dispatcher: invoke the Python script from R with `--data-dir` auto-detection. |
| `banc-update-seatable-clusters.R` | Push final `cns_cluster` + `cns_network` columns to SeaTable, joined on `root_888`. |
| `banc-visualise-clusters.R` | Per-cluster synapse-density PNGs (CNS dorsal, brain frontal, VNC ventral) + super_class composition stacked-bar PDF. |
| `banc-cluster-cns-network-heatmap.R` | Cross-tabulation heatmaps comparing new spectral clusters vs old SeaTable `cns_network` labels. |
| `install-deps.sh` | One-off: pip-install Python dependencies for the spectral clustering script. |

## Executable for users?

Partial. The Python script accepts `--source gcs` and runs anywhere; the R variant + SeaTable push need lab credentials.

See [`banc/README.md`](../README.md), top-level [`README.md`](../../README.md), and [`docs/data-products.md`](../../docs/data-products.md).
