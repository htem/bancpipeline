"""
Optic Lobe Alignment: Super_class Discrepancy Detection

Identifies BANC neurons where the connectivity-assigned cell type implies
a different super_class than what is in SeaTable. Excludes cases where
the original super_class is NA or blank.

Usage:
  python alignment/banc-alignment-discrepancies.py --side right
  python alignment/banc-alignment-discrepancies.py --side left
  python alignment/banc-alignment-discrepancies.py  # both sides
"""

import argparse
import functools
import pandas as pd
import numpy as np

print = functools.partial(print, flush=True)


def build_discrepancies(side, data_dir="data/optic_lobe"):
    """Build discrepancies CSV for one side."""
    print(f"\n=== Building discrepancies for {side} side ===")

    results = pd.read_csv(f"{data_dir}/banc_optic_{side}_alignment_results.csv",
                          dtype={"root_888": str, "best_target_match": str})
    banc_meta = pd.read_csv(f"{data_dir}/banc_optic_{side}_meta.csv",
                            dtype={"root_id": str}, low_memory=False)
    target_meta = pd.read_csv(f"{data_dir}/fafb_optic_{side}_meta.csv",
                            dtype={"target_id": str})

    # FAFB type -> super_class (majority vote)
    target_type_sc = target_meta.dropna(subset=["target_cell_type", "target_super_class"])
    target_type_sc = target_type_sc.groupby("target_cell_type")["target_super_class"].agg(
        lambda x: x.mode().iloc[0] if len(x) > 0 else None).to_dict()

    # Join results to BANC meta
    # Get ngl_link where available; use empty string for missing
    ngl_col = "ngl_link" if "ngl_link" in banc_meta.columns else None
    cols_to_get = ["root_id", "super_class", "fafb_cell_type", "fafb_match"]
    if ngl_col:
        cols_to_get.append(ngl_col)
    banc_cols = banc_meta[cols_to_get].copy()
    if not ngl_col:
        banc_cols["ngl_link"] = ""
    banc_cols["fafb_match"] = banc_cols["fafb_match"].astype(str).replace("nan", "")
    banc_cols = banc_cols.drop_duplicates("root_id")

    merged = results.merge(banc_cols, left_on="root_888", right_on="root_id", how="left")

    # Map assigned type to FAFB super_class
    merged["new_super_class"] = merged["assigned_cell_type"].map(target_type_sc)

    # Filter to typed neurons with non-blank old super_class
    typed = merged[
        (merged["assigned_cell_type"] != "") &
        (merged["super_class"].notna()) &
        (merged["super_class"] != "")
    ].copy()

    # Find discrepancies: old super_class != new super_class
    disc = typed[
        (typed["new_super_class"].notna()) &
        (typed["super_class"] != typed["new_super_class"])
    ].copy()

    disc["has_nblast_match"] = disc["best_target_match"].fillna("").str.len() > 0

    out = disc[[
        "root_888",
        "fafb_match",       # old FAFB match from SeaTable
        "best_target_match",  # new FAFB match from connectivity
        "fafb_cell_type",   # old FAFB cell type from SeaTable
        "assigned_cell_type",  # new connectivity-assigned type
        "super_class",      # old super_class
        "new_super_class",  # new super_class
        "confidence",
        "has_nblast_match",
    ]].rename(columns={
        "fafb_match": "old_fafb_match",
        "best_target_match": "new_fafb_match",
        "fafb_cell_type": "old_fafb_cell_type",
        "assigned_cell_type": "new_assigned_cell_type",
        "super_class": "old_super_class",
    }).sort_values("confidence", ascending=False)

    # Add ngl_link from meta (where available)
    out = out.merge(
        banc_cols[["root_id", "ngl_link"]].rename(columns={"root_id": "root_888"}),
        on="root_888", how="left", suffixes=("", "_meta")
    )
    out["ngl_link"] = out["ngl_link"].fillna("")

    out_file = f"{data_dir}/banc_optic_{side}_discrepancies.csv"
    out.to_csv(out_file, index=False)

    # Summary
    n_typed = len(typed)
    n_disc = len(out)
    n_no_match = (~out["has_nblast_match"]).sum()
    print(f"  Typed with super_class: {n_typed}")
    print(f"  Super_class discrepancies: {n_disc} ({100*n_disc/max(n_typed,1):.1f}%)")
    print(f"  Without FAFB match: {n_no_match}")
    print(f"  Mean confidence (discrepancies): {out['confidence'].mean():.4f}")
    print(f"  Saved: {out_file}")

    # Top discrepancy transitions
    if n_disc > 0:
        transitions = out.groupby(["old_super_class", "new_super_class"]).size().reset_index(name="n")
        transitions = transitions.sort_values("n", ascending=False)
        print("  Top super_class transitions:")
        for _, r in transitions.head(10).iterrows():
            print(f"    {r['old_super_class']} -> {r['new_super_class']}: {r['n']}")


def main():
    parser = argparse.ArgumentParser(description="Detect super_class discrepancies")
    parser.add_argument("--side", default=None, choices=["right", "left"],
                        help="Side to process (default: both)")
    parser.add_argument("--data-dir", default="data/optic_lobe")
    args = parser.parse_args()

    sides = [args.side] if args.side else ["right", "left"]
    for side in sides:
        try:
            build_discrepancies(side, args.data_dir)
        except FileNotFoundError as e:
            print(f"  Skipping {side}: {e}")


if __name__ == "__main__":
    main()
