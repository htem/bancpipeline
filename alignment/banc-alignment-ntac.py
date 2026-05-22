"""
Optic Lobe Cell Type Assignment: NTAC Comparison

Runs NTAC (Neuronal Type Assignment from Connectivity) on the BANC optic
lobe graph using the same seed/holdout split as the iterative alignment,
for fair comparison.

NTAC operates within a single connectome graph — it does NOT use FAFB data.
It propagates type labels through the BANC graph based on connectivity patterns.

Reference: Schwartzman et al. (2025) Nature Communications, doi:10.1038/s41467-025-68044-1

Usage:
  python alignment/banc-alignment-ntac.py --side right
  python alignment/banc-alignment-ntac.py --side both --no-holdout --output-suffix full
"""

import argparse
import functools
import os
import sys
import time
import numpy as np
import pandas as pd
import pyarrow.feather as pf
import scipy.sparse as sp
from ntac import Ntac

# Shared stratified-holdout splitter (single source of truth)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alignment_splits import build_target_type_vocab, stratified_holdout_split

print = functools.partial(print, flush=True)


# Module-level: which `root_<NNN>` column to read from the seeds CSV.
# Auto-detected by detect_root_col() in main() based on the seeds CSV header
# (e.g. "root_850" if data was prepped at v850, ROOT_COL if v888).
ROOT_COL = "root_id"


def detect_root_col(seeds_path):
    """Read the seeds CSV header and return the first column matching
    `root_<digits>`. Falls back to "root_id". Lets the script work
    against either v850 or v888 prep without re-prep."""
    cols = list(pd.read_csv(seeds_path, nrows=0).columns)
    for c in cols:
        if c.startswith("root_") and c[5:].isdigit():
            return c
    return "root_id"


def main():
    parser = argparse.ArgumentParser(description="NTAC comparison for optic lobe")
    parser.add_argument("--side", default="right", choices=["right", "left", "both"])
    parser.add_argument("--data-dir", default="data/optic_lobe")
    parser.add_argument("--file-prefix", default="optic")
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--output-suffix", default=None)
    parser.add_argument("--holdout-fraction", type=float, default=None,
                        help="Fraction of typed neurons to hold out (random split). "
                             "Mutually exclusive with --stratified-holdout.")
    parser.add_argument("--stratified-holdout", type=float, default=None,
                        help="Stratified holdout (same fraction per cell_type, "
                             "identical split to banc-alignment-run.py at the "
                             "same fraction). Mutually exclusive with --holdout-fraction.")
    parser.add_argument("--no-holdout", action="store_true",
                        help="Use ALL typed neurons as seeds, no holdout. "
                             "For production runs where we want predictions for untyped neurons.")
    args = parser.parse_args()

    side = args.side
    data_dir = args.data_dir
    fp = args.file_prefix

    if sum(x is not None and x is not False
           for x in (args.holdout_fraction, args.stratified_holdout,
                     args.no_holdout if args.no_holdout else None)) > 1:
        parser.error("Pass at most one of --holdout-fraction / --stratified-holdout / --no-holdout")

    # Load data (same as alignment pipeline)
    print(f"=== Loading data for {side} side (prefix={fp}) ===")
    seeds_path = f"{data_dir}/banc_{fp}_{side}_seeds.csv"
    global ROOT_COL
    ROOT_COL = detect_root_col(seeds_path)
    print(f"  Auto-detected seeds root column: {ROOT_COL}")

    banc_meta = pd.read_csv(f"{data_dir}/banc_{fp}_{side}_meta.csv",
                            dtype={"root_id": str}, low_memory=False)
    target_meta = pd.read_csv(f"{data_dir}/fafb_{fp}_{side}_meta.csv",
                            dtype={"target_id": str}, low_memory=False)
    banc_el = pf.read_feather(f"{data_dir}/banc_{fp}_{side}_edgelist.feather")
    seeds = pd.read_csv(seeds_path, dtype={ROOT_COL: str})

    target_type_vocab = build_target_type_vocab(target_meta)

    # Filter to proofread
    proofread = banc_meta["proofread"].fillna("").astype(str).str.upper() == "TRUE"
    roughly = banc_meta["roughly_proofread"].fillna("").astype(str).str.upper() == "TRUE"
    banc_meta = banc_meta[proofread | roughly].reset_index(drop=True)

    pool_ids = banc_meta["root_id"].astype(str).unique()
    n_pool = len(pool_ids)
    id_to_idx = {rid: i for i, rid in enumerate(pool_ids)}
    print(f"  Pool: {n_pool} neurons")

    # Build adjacency (pool neurons only, symmetric)
    banc_el = banc_el.copy()
    banc_el["pre"] = banc_el["pre"].astype(str)
    banc_el["post"] = banc_el["post"].astype(str)
    pool_set = set(pool_ids)
    banc_el = banc_el[banc_el["pre"].isin(pool_set) & banc_el["post"].isin(pool_set)]

    rows = banc_el["pre"].map(id_to_idx).dropna().astype(int).values
    cols = banc_el["post"].map(id_to_idx).dropna().astype(int).values
    vals_count = banc_el["count"].values[:len(rows)].astype(np.float32)
    adj = sp.csr_array((vals_count, (rows, cols)), shape=(n_pool, n_pool))
    adj = adj + adj.T
    print(f"  Adjacency: {adj.nnz} non-zeros")

    # Labels array and holdout bookkeeping
    labels = np.full(n_pool, "?", dtype=object)
    # Track which pool rows are in holdout — written to the output CSV so the
    # eval script does not have to replicate the split logic.
    is_holdout_arr = np.zeros(n_pool, dtype=bool)

    if args.stratified_holdout is not None:
        seed_list, holdout_list = stratified_holdout_split(
            seeds, id_to_idx, target_type_vocab,
            args.stratified_holdout, seed=42)
        for idx, ct, rid in seed_list:
            labels[idx] = ct
        holdout_indices = np.array([e[0] for e in holdout_list], dtype=np.int64)
        holdout_true_types = np.array([e[1] for e in holdout_list], dtype=object)
        is_holdout_arr[holdout_indices] = True
        n_seeded = (labels != "?").sum()
        print(f"  Stratified holdout: {n_seeded} seeds "
              f"({100*(1-args.stratified_holdout):.0f}% of typed), "
              f"{len(holdout_indices)} holdout "
              f"({100*args.stratified_holdout:.0f}%), seed=42")
    elif args.no_holdout:
        # All typed neurons become seeds
        for _, row in seeds.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in id_to_idx and pd.notna(ct) and ct in set(target_type_vocab):
                labels[id_to_idx[rid]] = ct
        holdout_indices = np.array([], dtype=np.int64)
        holdout_true_types = np.array([], dtype=object)
        n_seeded = (labels != "?").sum()
        print(f"  No holdout: {n_seeded} seeds "
              f"(all typed neurons in FAFB vocab)")
    else:
        # Random split (legacy behaviour — kept for backwards compat)
        if args.holdout_fraction is None:
            args.holdout_fraction = 0.2
        all_typed = []
        for _, row in seeds.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in id_to_idx and pd.notna(ct) and ct in set(target_type_vocab):
                all_typed.append((id_to_idx[rid], ct))
        typed_indices = np.array([t[0] for t in all_typed])
        typed_types = np.array([t[1] for t in all_typed])
        n_typed = len(all_typed)
        np.random.seed(42)
        perm = np.random.permutation(n_typed)
        n_holdout_split = int(args.holdout_fraction * n_typed)
        holdout_idx = perm[:n_holdout_split]
        seed_idx = perm[n_holdout_split:]
        for i in seed_idx:
            labels[typed_indices[i]] = typed_types[i]
        holdout_indices = typed_indices[holdout_idx]
        holdout_true_types = typed_types[holdout_idx]
        is_holdout_arr[holdout_indices] = True
        n_seeded = (labels != "?").sum()
        print(f"  Random holdout: {n_seeded} seeds "
              f"({100*(1-args.holdout_fraction):.0f}% of {n_typed} typed), "
              f"{len(holdout_indices)} holdout "
              f"({100*args.holdout_fraction:.0f}%), seed=42")

    # Run NTAC
    print(f"\n=== Running NTAC (seeded, {args.max_iter} iterations) ===")
    t0 = time.time()
    model = Ntac(adj, labels=labels, lr=0.3, topk=1, verbose=False)

    best_acc = 0
    best_iter = 0
    best_partition = None

    for i in range(args.max_iter):
        model.step()
        partition = model.get_partition()

        # Count assigned
        n_assigned = sum(1 for p in partition if pd.notna(p) and str(p) != "?")
        elapsed = time.time() - t0

        if len(holdout_indices) > 0:
            # Evaluate holdout
            pred = partition[holdout_indices]
            pred = np.array([str(p) if pd.notna(p) else "?" for p in pred])
            valid = pred != "?"
            correct = (pred[valid] == holdout_true_types[valid]).sum() if valid.any() else 0
            n_valid = valid.sum()
            acc = 100 * correct / n_valid if n_valid > 0 else 0

            mi1_mask = holdout_true_types == "Mi1"
            mi1_pred = pred[mi1_mask]
            mi1_acc = 100 * (mi1_pred == "Mi1").sum() / mi1_mask.sum() if mi1_mask.sum() > 0 else 0

            print(f"  Iter {i+1:3d}: holdout={acc:.1f}%, Mi1={mi1_acc:.1f}%, "
                  f"assigned={n_assigned}/{n_pool}, {elapsed:.1f}s")

            if acc > best_acc:
                best_acc = acc
                best_iter = i + 1
                best_partition = partition.copy()
        else:
            # No holdout — just track progress
            print(f"  Iter {i+1:3d}: assigned={n_assigned}/{n_pool}, {elapsed:.1f}s")
            best_partition = partition.copy()
            best_iter = i + 1

    if len(holdout_indices) > 0:
        print(f"\n  Best holdout: {best_acc:.1f}% at iteration {best_iter}")
    print(f"  Total time: {time.time() - t0:.1f}s")

    # Save results (is_holdout column lets the eval script restrict holdout
    # accuracy to the same neurons the paired align run used)
    results = pd.DataFrame({
        ROOT_COL: pool_ids,
        "ntac_cell_type": best_partition,
        "is_holdout": is_holdout_arr,
    })
    suffix = f"_{args.output_suffix}" if args.output_suffix else ""
    out_file = f"{data_dir}/banc_{fp}_{side}_ntac{suffix}.csv"
    results.to_csv(out_file, index=False)
    print(f"  Saved: {out_file}")


if __name__ == "__main__":
    main()
