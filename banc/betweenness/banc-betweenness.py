#!/usr/bin/env python3
"""
Source-Target Betweenness Centrality for BANC

Based on: python/betweeness/Results_afferent_to_efferent.ipynb
          (Zuoyu Zhang & Tatsuo Okubo, 2025/10/31)

Computes two betweenness centrality measures on the BANC connectivity
graph, filtered to neurons only (using meta feather):

  1. Afferent-to-efferent: sources = sensory neurons,
     targets = motor + visceral_circulatory neurons
  2. All-to-all: standard betweenness (no source/target restriction)

Both computations are long-running (hours). Use --mode to run one at a
time, or omit to run both sequentially.

Usage:
  python banc-betweenness.py                          # both modes
  python banc-betweenness.py --mode afferent_efferent # just aff->eff
  python banc-betweenness.py --mode all_to_all        # just all-to-all
  python banc-betweenness.py --data-dir path/to/banc_888/ --source v3
"""

import argparse
import os
import sys
import time

import igraph as ig
import numpy as np
import pandas as pd


# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------

AFFERENT_SUPER_CLASSES = ["sensory", "sensory_ascending", "sensory_descending"]
EFFERENT_SUPER_CLASSES = [
    "motor",
    "visceral_circulatory",
    "ascending_visceral_circulatory",
]

# Super classes to exclude from the neuron graph (non-neurons)
NON_NEURON_PATTERNS = [
    "glia", "trachea", "not_a_neuron", "merge", "orpha", "tadpole",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Betweenness centrality on BANC connectivity graph"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Directory containing banc_{version}_meta.feather and "
             "banc_{version}_edgelist_simple_{v2|v3}.feather. "
             "Auto-detects O2 vs local if not specified."
    )
    parser.add_argument(
        "--output-dir", type=str,
        default=os.path.join(
            os.environ.get(
                "BANC_PROJECT_DATA_DIR",
                "/n/data1/hms/neurobio/wilson/banc/BANC-project/data"),
            "betweeness"),
        help=("Output directory for betweenness CSVs. Default: "
              "$BANC_PROJECT_DATA_DIR/betweeness, falling back to "
              "/n/data1/.../BANC-project/data/betweeness. Set "
              "BANC_PROJECT_DATA_DIR to redirect to a different repo checkout.")
    )
    parser.add_argument(
        "--banc-version", type=int, default=888,
        help="BANC dataset version (default: 888)"
    )
    parser.add_argument(
        "--source", type=str, default=None, choices=[None, "v2", "v3"],
        help="Synapse source for the edgelist suffix. "
             "Defaults to $BANC_SYN_SOURCE or 'v3'."
    )
    parser.add_argument(
        "--mode", type=str, default="both",
        choices=["afferent_efferent", "all_to_all", "both"],
        help="Which betweenness to compute (default: both)"
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
        f"/n/data1/hms/neurobio/wilson/connectomes/banc/banc_{version}",
        f"lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_{version}",
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    print(f"ERROR: Could not find data directory. Tried: {candidates}",
          file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    args = parse_args()
    t_start = time.time()

    data_dir = args.data_dir or detect_data_dir(args.banc_version)
    print(f"Data directory: {data_dir}")

    id_col = f"root_{args.banc_version}"

    # --- Load data ---
    print("Loading meta and edgelist...")
    meta = pd.read_feather(
        os.path.join(data_dir, f"banc_{args.banc_version}_meta.feather")
    )
    edgelist = pd.read_feather(
        os.path.join(data_dir,
                     f"banc_{args.banc_version}_edgelist_simple_{args.source}.feather")
    )
    print(f"  Meta: {meta.shape[0]} rows, Edgelist: {edgelist.shape[0]} edges")

    # --- Identify neurons ---
    # Keep only proofread or roughly_proofread neurons, exclude non-neuron super_classes.
    proofread_mask = (
        (meta["proofread"].astype(str) == "TRUE")
        | (meta["roughly_proofread"].astype(str) == "TRUE")
    )
    pattern = "|".join(NON_NEURON_PATTERNS)
    superclass_mask = ~meta["super_class"].str.contains(pattern, case=False, na=True)
    neuron_mask = proofread_mask & superclass_mask
    meta_neurons = meta[neuron_mask].copy()
    meta_neurons["root_id"] = meta_neurons[id_col].astype(str)
    neuron_ids = set(meta_neurons["root_id"].dropna().unique())
    print(f"  Neurons in meta: {len(neuron_ids)}")

    # --- Filter edgelist to neuron-only graph ---
    edgelist["pre"] = edgelist["pre"].astype(str)
    edgelist["post"] = edgelist["post"].astype(str)
    el = edgelist[
        edgelist["pre"].isin(neuron_ids) & edgelist["post"].isin(neuron_ids)
    ].copy()
    print(f"  Edges after neuron filter: {el.shape[0]}")

    # --- Build node mapping ---
    unique_nodes = sorted(
        set(el["pre"].unique()) | set(el["post"].unique())
    )
    node_to_index = {nid: i for i, nid in enumerate(unique_nodes)}
    n = len(unique_nodes)
    print(f"  Unique nodes in graph: {n}")

    # --- Build igraph directed graph ---
    print("Building igraph graph...")
    edges = list(zip(
        el["pre"].map(node_to_index).values,
        el["post"].map(node_to_index).values,
    ))
    g = ig.Graph(n, edges=edges, directed=True)
    print(f"  Nodes: {g.vcount()}, Edges: {g.ecount()}")

    # --- Build root_id -> meta lookup for output ---
    meta_lookup = (
        meta_neurons[["root_id", "super_class", "cell_type"]]
        .drop_duplicates(subset="root_id")
        .set_index("root_id")
    )

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Afferent-to-efferent betweenness ---
    if args.mode in ("afferent_efferent", "both"):
        print("\n--- Afferent-to-efferent betweenness ---")

        afferent_roots = set(
            meta_neurons.loc[
                meta_neurons["super_class"].isin(AFFERENT_SUPER_CLASSES), "root_id"
            ]
        )
        efferent_roots = set(
            meta_neurons.loc[
                meta_neurons["super_class"].isin(EFFERENT_SUPER_CLASSES), "root_id"
            ]
        )

        afferent_idx = [node_to_index[r] for r in afferent_roots if r in node_to_index]
        efferent_idx = [node_to_index[r] for r in efferent_roots if r in node_to_index]

        print(f"  Afferents in graph: {len(afferent_idx)}")
        print(f"  Efferents in graph: {len(efferent_idx)}")

        t0 = time.time()
        print("  Computing betweenness (this may take hours)...")
        bc = g.betweenness(
            vertices=None,
            directed=True,
            weights=None,
            sources=afferent_idx,
            targets=efferent_idx,
        )
        print(f"  Done in {(time.time() - t0) / 60:.1f} minutes")

        result = pd.DataFrame({
            "vertex_id": range(n),
            id_col: unique_nodes,
            "betweenness": bc,
        })
        result = result.join(meta_lookup, on=id_col, how="left")
        result = result[["vertex_id", id_col, "super_class", "cell_type", "betweenness"]]

        outfile = os.path.join(
            args.output_dir,
            f"betweenness_afferent_to_efferent_banc_{args.banc_version}_{args.source}.csv",
        )
        result.to_csv(outfile, index=False)
        print(f"  Saved: {outfile} ({result.shape[0]} rows)")
        dated = outfile.replace(".csv", f"_{time.strftime('%Y-%m-%d')}.csv")
        result.to_csv(dated, index=False)
        print(f"  Saved (dated): {dated}")

    # --- All-to-all betweenness ---
    if args.mode in ("all_to_all", "both"):
        print("\n--- All-to-all betweenness ---")

        t0 = time.time()
        print("  Computing betweenness (this may take hours)...")
        bc = g.betweenness(
            vertices=None,
            directed=True,
            weights=None,
        )
        print(f"  Done in {(time.time() - t0) / 60:.1f} minutes")

        result = pd.DataFrame({
            "vertex_id": range(n),
            id_col: unique_nodes,
            "betweenness": bc,
        })
        result = result.join(meta_lookup, on=id_col, how="left")
        result = result[["vertex_id", id_col, "super_class", "cell_type", "betweenness"]]

        outfile = os.path.join(
            args.output_dir,
            f"betweenness_all_to_all_banc_{args.banc_version}_{args.source}.csv",
        )
        result.to_csv(outfile, index=False)
        print(f"  Saved: {outfile} ({result.shape[0]} rows)")
        dated = outfile.replace(".csv", f"_{time.strftime('%Y-%m-%d')}.csv")
        result.to_csv(dated, index=False)
        print(f"  Saved (dated): {dated}")

    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed / 60:.1f} minutes")


if __name__ == "__main__":
    main()
