"""stratified_holdout_split invariants.

Holdout logic is shared between align.py and NTAC — if the fraction-per-type
contract drifts, anchor/holdout sets diverge silently and holdout accuracy
stops being comparable across methods.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _make_seeds(types, root_col="root_850"):
    """Build a seeds_df with a root_<NNN> column and cell_type."""
    return pd.DataFrame({
        root_col:    [f"r{i}" for i in range(len(types))],
        "cell_type": types,
    })


def test_stratified_split_fraction_per_type(splits):
    """Holdout count per type = max(1, round(frac * n))."""
    types = ["A"] * 10 + ["B"] * 4    # 10 A's + 4 B's
    seeds = _make_seeds(types)
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    vocab = ["A", "B"]

    _, holdout = splits.stratified_holdout_split(
        seeds, pool_id_to_idx, vocab, fraction=0.5, seed=42
    )
    counts = {t: 0 for t in vocab}
    for _, ct, _ in holdout:
        counts[ct] += 1
    assert counts["A"] == 5       # round(0.5 * 10) = 5
    assert counts["B"] == 2       # round(0.5 * 4) = 2


def test_stratified_split_min_one_per_type(splits):
    """Tiny fractions still put >=1 neuron per type in holdout."""
    seeds = _make_seeds(["A"] * 20 + ["B"] * 2)
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    _, holdout = splits.stratified_holdout_split(
        seeds, pool_id_to_idx, ["A", "B"], fraction=0.01, seed=42
    )
    counts = {t: 0 for t in ["A", "B"]}
    for _, ct, _ in holdout:
        counts[ct] += 1
    assert counts["A"] >= 1
    assert counts["B"] >= 1


def test_stratified_split_disjoint_sets(splits):
    """Seed ∩ holdout = ∅ (by pool_idx)."""
    seeds = _make_seeds(["A"] * 10 + ["B"] * 10 + ["C"] * 5)
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    seed_entries, holdout_entries = splits.stratified_holdout_split(
        seeds, pool_id_to_idx, ["A", "B", "C"], fraction=0.3, seed=42
    )
    seed_idxs = {e[0] for e in seed_entries}
    holdout_idxs = {e[0] for e in holdout_entries}
    assert not (seed_idxs & holdout_idxs)


def test_stratified_split_deterministic(splits):
    """Same seed -> identical split."""
    seeds = _make_seeds(["A"] * 20 + ["B"] * 20)
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    a = splits.stratified_holdout_split(seeds, pool_id_to_idx, ["A", "B"],
                                         fraction=0.4, seed=42)
    b = splits.stratified_holdout_split(seeds, pool_id_to_idx, ["A", "B"],
                                         fraction=0.4, seed=42)
    assert a == b


def test_stratified_split_respects_vocab(splits):
    """Neurons whose cell_type is outside fafb_type_vocab are silently
    dropped — matches align.py's type_to_idx filter."""
    seeds = _make_seeds(["A", "A", "OutOfVocab", "B"])
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    seed_entries, holdout_entries = splits.stratified_holdout_split(
        seeds, pool_id_to_idx, ["A", "B"], fraction=0.5, seed=42
    )
    all_ct = {e[1] for e in seed_entries} | {e[1] for e in holdout_entries}
    assert all_ct <= {"A", "B"}


def test_stratified_split_fraction_one_all_holdout(splits):
    """fraction=1.0 puts every typed neuron in holdout, none as seed."""
    seeds = _make_seeds(["A"] * 5 + ["B"] * 5)
    pool_id_to_idx = {rid: i for i, rid in enumerate(seeds["root_850"])}
    seed_entries, holdout_entries = splits.stratified_holdout_split(
        seeds, pool_id_to_idx, ["A", "B"], fraction=1.0, seed=42
    )
    assert seed_entries == []
    assert len(holdout_entries) == 10
