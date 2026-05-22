"""Profile-builder correctness tests.

Adjacency and type matrices are the input to centroid + profile computation.
An axis swap or a stale index map would silently poison every downstream
score. Tests use toy graphs with hand-computed expectations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_build_adjacency_pool_first_then_boundary(align):
    """Pool neurons must occupy indices 0..n_pool-1 in input order; boundary
    neurons land at n_pool..n_all-1 in sorted order."""
    el_df = pd.DataFrame({
        "pre":  ["A", "X", "B", "Y"],
        "post": ["B", "A", "Y", "X"],
        "norm": [0.1, 0.2, 0.3, 0.4],
    })
    adj, id_to_idx, n_pool = align.build_adjacency(el_df, ["A", "B"])
    assert n_pool == 2
    assert id_to_idx["A"] == 0
    assert id_to_idx["B"] == 1
    # Boundary neurons sorted
    assert id_to_idx["X"] == 2
    assert id_to_idx["Y"] == 3
    assert adj.shape == (4, 4)


def test_build_adjacency_weight_placement(align):
    """Edge (pre, post, norm) must land at adj[pre_idx, post_idx]."""
    el_df = pd.DataFrame({
        "pre":  ["A", "B"],
        "post": ["B", "A"],
        "norm": [0.7, 0.3],
    })
    adj, id_to_idx, _ = align.build_adjacency(el_df, ["A", "B"])
    a, b = id_to_idx["A"], id_to_idx["B"]
    dense = adj.toarray()
    assert dense[a, b] == pytest.approx(0.7)
    assert dense[b, a] == pytest.approx(0.3)
    # Self-edges unset
    assert dense[a, a] == 0
    assert dense[b, b] == 0


def test_build_adjacency_deduplicates_pool(align):
    """Duplicates in neuron_ids collapse — function only emits each unique id
    once in the pool."""
    el_df = pd.DataFrame({
        "pre": ["A"], "post": ["A"], "norm": [1.0],
    })
    _, id_to_idx, n_pool = align.build_adjacency(el_df, ["A", "A"])
    assert n_pool == 1
    assert id_to_idx["A"] == 0


def test_build_type_matrix_onehot(align):
    """type_matrix[i, t] = 1 iff neuron i has cell_type mapped to type index t."""
    meta = pd.DataFrame({
        "root_id":   ["A", "B", "C"],
        "cell_type": ["T5a", "Mi1", ""],    # empty string -> unassigned
    })
    id_to_idx = {"A": 0, "B": 1, "C": 2}
    type_to_idx = {"T5a": 0, "Mi1": 1}
    m = align.build_type_matrix(meta, "root_id", id_to_idx, type_to_idx, n_all=3)
    dense = m.toarray()
    assert dense[0, 0] == 1.0 and dense[0, 1] == 0.0     # A -> T5a
    assert dense[1, 1] == 1.0 and dense[1, 0] == 0.0     # B -> Mi1
    assert dense[2].sum() == 0                           # C -> unassigned


def test_build_type_matrix_respects_n_all(align):
    """Neurons with id_to_idx index >= n_all are dropped (boundary overflow)."""
    meta = pd.DataFrame({"root_id": ["A", "B"], "cell_type": ["T5a", "Mi1"]})
    id_to_idx = {"A": 0, "B": 5}  # B beyond n_all
    type_to_idx = {"T5a": 0, "Mi1": 1}
    m = align.build_type_matrix(meta, "root_id", id_to_idx, type_to_idx, n_all=3)
    assert m.shape == (3, 2)
    assert m.toarray()[0, 0] == 1.0
    # B never lands
    assert m.toarray()[:, 1].sum() == 0
