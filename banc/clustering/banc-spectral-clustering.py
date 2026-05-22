#!/usr/bin/env python3
"""
Spectral Clustering of BANC CNS Network

Port of python/spectral_clustering/banc_spectral_clustering_final.ipynb
for use on O2 (SLURM) and locally.

Algorithm:
  1. Load BANC meta + simple edgelist from feather (source=v2 or v3)
  2. Quality filter (exclude glia, trachea, unproofread, etc.)
  3. Assign clustering_set (central brain, VNC, neck, visual)
  4. Filter edgelist, prune to strongly connected component
  5. Build column-normalized, symmetrized adjacency matrix
  6. Spectral clustering (normalized Laplacian, bottom-k eigenvectors, KMeans)
  7. UMAP of spectral embedding
  8. Assign cns_network labels (from local meta or GCS via --gcs)
  9. Save CSV (output suffix carries the source: _v2 or _v3)

Parameters (paper Methods §"Spectral clustering"):
  min_connection_strength = 1
  cluster_count = 13
  cluster_seed = 10
  embedding_seed = 3
  UMAP: n_neighbors=100, metric=cosine, min_dist=0, n_components=2

Usage:
  python banc-spectral-clustering.py --data-dir path/to/banc_888/ --source v3
  python banc-spectral-clustering.py  # auto-detects O2 vs local, source=v3
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize
import umap
import plotly.express as px


# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Spectral clustering of BANC CNS connectivity network"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Directory containing banc_{version}_meta.feather and "
             "banc_{version}_edgelist_simple_{v2|v3}.feather. "
             "Auto-detects O2 vs local if not specified."
    )
    parser.add_argument(
        "--source", type=str, default=None, choices=[None, "v2", "v3"],
        help="Synapse source for the edgelist suffix. "
             "Defaults to $BANC_SYN_SOURCE or 'v3'."
    )
    parser.add_argument(
        "--output-dir", type=str,
        default=os.path.join(
            os.environ.get(
                "BANC_PROJECT_DATA_DIR",
                "/n/data1/hms/neurobio/wilson/banc/BANC-project/data"),
            "cns_network"),
        help=("Output directory for clustering CSV. Default: "
              "$BANC_PROJECT_DATA_DIR/cns_network, falling back to "
              "/n/data1/.../BANC-project/data/cns_network. Set "
              "BANC_PROJECT_DATA_DIR to redirect to a different repo checkout "
              "(e.g. on another machine).")
    )
    parser.add_argument(
        "--gcs", action="store_true", default=False,
        help="Fetch cns_network labels from GCS meta feather "
             "(gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_<ver>/) "
             "via gsutil. If not set, uses cns_network from the local meta."
    )
    parser.add_argument(
        "--min-connection-strength", type=int, default=1,
        help="Minimum synapse count per connection (default: 1, paper Methods)"
    )
    parser.add_argument(
        "--cluster-count", type=int, default=13,
        help="Number of spectral clusters (default: 13, paper Methods)"
    )
    parser.add_argument(
        "--cluster-seed", type=int, default=10,
        help="Random seed for KMeans (default: 10)"
    )
    parser.add_argument(
        "--embedding-seed", type=int, default=3,
        help="Random seed for UMAP (default: 3)"
    )
    parser.add_argument(
        "--banc-version", type=int, default=888,
        help="BANC dataset version (default: 888)"
    )
    args = parser.parse_args()
    if args.source is None:
        args.source = os.environ.get("BANC_SYN_SOURCE", "v3").lower()
    if args.source not in ("v2", "v3"):
        parser.error(f"--source must be 'v2' or 'v3' (got {args.source!r})")
    return args


def detect_data_dir(version):
    """Auto-detect data directory: O2 cluster vs local."""
    candidates = [
        # O2 HPC
        f"/n/data1/hms/neurobio/wilson/connectomes/banc/banc_{version}",
        # Local (relative to repo root)
        f"lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_{version}",
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    print(f"ERROR: Could not find data directory. Tried: {candidates}",
          file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------
# Spectral clustering
# ---------------------------------------------------------------

def spectral_clustering_from_adj(adj_sym, num_clusters, random_state=None):
    """
    Spectral clustering on a symmetric adjacency matrix.

    1. Compute normalized Laplacian
    2. Find bottom-k eigenvectors (smallest eigenvalues)
    3. Row-normalize the eigenvector embedding
    4. KMeans clustering

    Returns (labels, embedding) where labels are 1-indexed.
    """
    # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
    lap = laplacian(adj_sym, normed=True)

    # Bottom-k eigenvectors
    eigval, eigvec = eigsh(lap, k=num_clusters, which="SM")

    # Row-normalize
    embedding = normalize(eigvec)

    # KMeans
    km = KMeans(n_clusters=num_clusters, random_state=random_state, n_init=10)
    labels = km.fit_predict(embedding) + 1  # 1-indexed

    return labels, embedding


def perform_umap(embedding, n_neighbors=100, random_state=3):
    """UMAP projection of spectral embedding to 2D."""
    model = umap.UMAP(
        n_neighbors=n_neighbors,
        metric="cosine",
        min_dist=0,
        n_components=2,
        random_state=random_state,
        n_jobs=1,
    )
    return model.fit_transform(embedding)


# ---------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------

def main():
    args = parse_args()
    t_start = time.time()

    # Resolve data directory
    data_dir = args.data_dir or detect_data_dir(args.banc_version)
    print(f"Data directory: {data_dir}")

    # --- 1. Load data ---
    print("Loading meta and edgelist...")
    meta_path = os.path.join(data_dir, f"banc_{args.banc_version}_meta.feather")
    el_path = os.path.join(
        data_dir,
        f"banc_{args.banc_version}_edgelist_simple_{args.source}.feather",
    )
    meta = pd.read_feather(meta_path)
    edgelist = pd.read_feather(el_path)
    print(f"  Meta: {meta.shape[0]} neurons, Edgelist: {edgelist.shape[0]} edges")

    # Use root_{banc.version} as the canonical ID column
    id_col = f"root_{args.banc_version}"
    meta = meta.rename(columns={id_col: "root_id"})
    meta = meta[meta["root_id"].notna() & (meta["root_id"] != "")].copy()

    # --- 2. Quality filter ---
    print("Applying quality filter...")
    quality_mask = (
        ~meta["super_class"].str.contains(
            "glia|trachea|not_a_neuron|merge|debris",
            case=False, na=False
        )
        & ~meta["status"].str.contains(
            "GLIA|TRACHEA|NOT_A_NEURON|DEBRIS|MERGE|DELETE",
            na=False
        )
        & (
            (meta["proofread"].astype(str) == "TRUE")
            | (meta["roughly_proofread"].astype(str) == "TRUE")
        )
    )
    meta = meta[quality_mask].copy()
    print(f"  After quality filter: {meta.shape[0]} neurons")

    # --- 3. Assign clustering_set ---
    # Identify neck_connective neurons via super_class (ascending/descending) —
    # backwards-compatible with the now-deprecated region == "neck_connective" filter.
    clustering_regions = ["central_brain", "ventral_nerve_cord"]
    excluded_sc = ["sensory", "motor", "efferent", "afferent", "visceral"]
    excluded_pattern = "|".join(excluded_sc)

    meta["clustering_set"] = None

    # Visual projection / centrifugal
    visual_mask = meta["super_class"].isin(["visual_centrifugal", "visual_projection"])
    meta.loc[visual_mask, "clustering_set"] = "visual"

    # Ascending/descending neurons → neck_connective bucket
    neck_mask = (
        meta["super_class"].str.contains("ascending|descending", na=False)
        & meta["clustering_set"].isna()
    )
    meta.loc[neck_mask, "clustering_set"] = "neck_connective"

    # Region-based, excluding sensory/motor/etc.
    for region in clustering_regions:
        region_mask = (
            (meta["region"] == region)
            & ~meta["super_class"].str.contains(excluded_pattern, na=False)
            & meta["clustering_set"].isna()
        )
        meta.loc[region_mask, "clustering_set"] = region

    clustering_values = ["central_brain", "neck_connective",
                         "ventral_nerve_cord", "visual"]
    cluster_ids = set(
        meta.loc[meta["clustering_set"].isin(clustering_values), "root_id"].values
    )
    print(f"  Neurons in clustering set: {len(cluster_ids)}")

    # --- 4. Filter edgelist + prune ---
    print("Filtering edgelist...")
    el = edgelist[
        (edgelist["count"] >= args.min_connection_strength)
        & (edgelist["pre"] != edgelist["post"])
        & (edgelist["pre"].isin(cluster_ids))
        & (edgelist["post"].isin(cluster_ids))
    ].copy()
    print(f"  Edges after filter: {el.shape[0]}")

    # Iterative pruning to strongly connected component
    print("Pruning to strongly connected component...")
    prev_n = -1
    while True:
        both = set(el["pre"].values) & set(el["post"].values)
        if len(both) == prev_n:
            break
        prev_n = len(both)
        el = el[el["pre"].isin(both) & el["post"].isin(both)]
        print(f"  Pruning: {prev_n} neurons")

    neuron_ids = sorted(both)
    n = len(neuron_ids)
    print(f"  Final: {n} neurons, {el.shape[0]} edges")

    # --- 5. Build sparse adjacency matrix ---
    print("Building adjacency matrix...")
    id_map = {nid: i for i, nid in enumerate(neuron_ids)}

    # Aggregate weights per pre-post pair
    el_agg = el.groupby(["pre", "post"], as_index=False)["count"].sum()
    el_agg.rename(columns={"count": "weight"}, inplace=True)

    row_idx = el_agg["pre"].map(id_map).values
    col_idx = el_agg["post"].map(id_map).values
    weights = el_agg["weight"].values.astype(np.float32)

    adj = coo_matrix((weights, (row_idx, col_idx)), shape=(n, n))
    adj_csc = csc_matrix(adj)

    # Column-normalize
    col_sums = np.array(adj_csc.sum(axis=0)).flatten()
    col_sums[col_sums == 0] = 1
    adj_norm = adj_csc.multiply(1.0 / col_sums)

    # Symmetrize: 0.5 * (A_norm + A_norm^T)
    adj_sym = 0.5 * (adj_norm + adj_norm.T)
    adj_sym = adj_sym.astype(np.float32)
    print("  Matrix built and symmetrized")

    # --- 6. Spectral clustering ---
    print(f"Computing spectral clustering (k={args.cluster_count})...")
    labels, embedding = spectral_clustering_from_adj(
        adj_sym, args.cluster_count, random_state=args.cluster_seed
    )
    for cl in sorted(set(labels)):
        print(f"  Cluster {cl}: {(labels == cl).sum()} neurons")

    # --- 7. UMAP ---
    print("Computing UMAP embedding...")
    umap_coords = perform_umap(
        embedding,
        n_neighbors=100,
        random_state=args.embedding_seed,
    )

    # --- 8. Build result dataframe ---
    result = pd.DataFrame({
        "root_id": neuron_ids,
        "spectral_cluster": labels,
        "umap_x": umap_coords[:, 0],
        "umap_y": umap_coords[:, 1],
    })

    # Join supervoxel_id and position from meta
    meta_cols = meta[["root_id", "supervoxel_id", "position"]].drop_duplicates(
        subset="root_id"
    )
    result = result.merge(meta_cols, on="root_id", how="left")

    # --- 9. Assign cns_network labels ---
    print("Assigning cns_network labels...")
    if args.gcs:
        # Fetch latest meta for cns_network labels.
        # On O2 the versioned save.path already has it; otherwise download
        # from GCS via gsutil.
        o2_meta = os.path.join(
            f"/n/data1/hms/neurobio/wilson/connectomes/banc/banc_{args.banc_version}",
            f"banc_{args.banc_version}_meta.feather",
        )
        if os.path.isfile(o2_meta):
            print(f"  Reading labels from O2: {o2_meta}")
            gcs_meta = pd.read_feather(o2_meta)
        else:
            gcs_path = (
                f"gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/"
                f"banc_{args.banc_version}/banc_{args.banc_version}_meta.feather"
            )
            local_tmp = os.path.join(
                tempfile.gettempdir(),
                f"banc_{args.banc_version}_meta_gcs.feather",
            )
            print(f"  Downloading {gcs_path} ...")
            subprocess.run(["gsutil", "cp", gcs_path, local_tmp], check=True)
            gcs_meta = pd.read_feather(local_tmp)
        id_col = f"root_{args.banc_version}"
        if id_col in gcs_meta.columns:
            gcs_meta = gcs_meta.rename(columns={id_col: "root_id"})
        st_sub = gcs_meta[["root_id", "cns_network"]].dropna(subset=["cns_network"])
        st_sub = st_sub[st_sub["cns_network"] != ""]
        st_map = dict(zip(st_sub["root_id"], st_sub["cns_network"]))
        print(f"  GCS labels: {len(st_map)} neurons with cns_network")
    else:
        # Use cns_network from the locally loaded meta feather
        st_sub = meta[["root_id", "cns_network"]].dropna(subset=["cns_network"])
        st_sub = st_sub[st_sub["cns_network"] != ""]
        st_map = dict(zip(st_sub["root_id"], st_sub["cns_network"]))

    result["existing_cns_network"] = result["root_id"].map(st_map)

    # Assign cns_network per cluster by Hungarian (bijective) matching.
    # Pair each new spectral_cluster with at most one OLD cns_network label so
    # total neuron overlap on the assignment is maximised. Independent per-
    # cluster majority votes let two new clusters both claim e.g. "right
    # olfactory" while no cluster ends up labelled "left olfactory" —
    # Hungarian forbids that. Clusters left unmatched (all labels claimed
    # elsewhere, or zero overlap with every label) keep their numeric string.
    result["cns_network"] = result["spectral_cluster"].astype(str)
    cluster_ids = sorted(result["spectral_cluster"].unique())
    valid_existing = result["existing_cns_network"].dropna()
    valid_existing = valid_existing[valid_existing != ""]
    old_labels = sorted(valid_existing.unique())

    if len(old_labels) == 0:
        print("  No old cns_network labels available — using cluster numbers")
    else:
        K = len(cluster_ids)
        L = len(old_labels)
        label_to_col = {lbl: j for j, lbl in enumerate(old_labels)}
        # K x L overlap matrix (rows = new clusters, cols = old labels)
        M = np.zeros((K, L), dtype=np.int64)
        for i, cl in enumerate(cluster_ids):
            mask = result["spectral_cluster"] == cl
            ex = result.loc[mask, "existing_cns_network"].dropna()
            ex = ex[ex != ""]
            if len(ex) > 0:
                vc = ex.value_counts()
                for lbl, n in vc.items():
                    if lbl in label_to_col:
                        M[i, label_to_col[lbl]] = n
        # Pad to square so any unmatchable row can take a zero-cost padding col
        N = max(K, L)
        M_sq = np.zeros((N, N), dtype=np.int64)
        M_sq[:K, :L] = M
        row_ind, col_ind = linear_sum_assignment(M_sq, maximize=True)
        # row_ind is 0..N-1 in order; col_ind[i] is the column assigned to row i
        for i, cl in enumerate(cluster_ids):
            j = col_ind[i]
            mask = result["spectral_cluster"] == cl
            if j < L and M[i, j] > 0:
                lbl = old_labels[j]
                result.loc[mask, "cns_network"] = lbl
                print(f"  Cluster {cl} -> '{lbl}' "
                      f"({M[i, j]}/{mask.sum()} neurons)")
            else:
                print(f"  Cluster {cl} -> '{cl}' "
                      f"(no available label — using cluster number)")

    result.drop(columns=["existing_cns_network"], inplace=True)

    # --- 10. Save output ---
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(
        args.output_dir,
        f"spectral_clustering_min_connection_strength_"
        f"{args.min_connection_strength}_banc_version_{args.banc_version}"
        f"_cluster_count_{args.cluster_count}_cluster_seed_"
        f"{args.cluster_seed}_embedding_seed_{args.embedding_seed}"
        f"_{args.source}.csv"
    )

    # Final column order
    result = result[["root_id", "supervoxel_id", "position",
                     "spectral_cluster", "umap_x", "umap_y", "cns_network"]]
    result.to_csv(output_file, index=False)

    elapsed = time.time() - t_start
    print(f"\nSaved: {output_file}")
    print(f"  {result.shape[0]} neurons, "
          f"{result['spectral_cluster'].nunique()} clusters")
    print(f"  Total time: {elapsed / 60:.1f} minutes")

    # --- 11. Interactive UMAP plot ---
    print("Generating interactive UMAP plot...")
    plot_df = result.copy()
    plot_df["spectral_cluster"] = plot_df["spectral_cluster"].astype(str)

    fig = px.scatter(
        plot_df,
        x="umap_x", y="umap_y",
        color="cns_network",
        hover_data=["root_id", "spectral_cluster", "cns_network"],
        title=f"BANC v{args.banc_version} Spectral Clustering UMAP "
              f"(k={args.cluster_count})",
        labels={"umap_x": "UMAP 1", "umap_y": "UMAP 2",
                "cns_network": "CNS Network"},
        width=1000, height=800,
    )
    fig.update_traces(marker=dict(size=2, opacity=0.6))
    fig.update_layout(legend=dict(itemsizing="constant"))

    html_file = os.path.join(
        args.output_dir,
        f"spectral_clustering_umap_banc_version_{args.banc_version}_{args.source}.html"
    )
    fig.write_html(html_file)
    print(f"  UMAP plot saved: {html_file}")


if __name__ == "__main__":
    main()
