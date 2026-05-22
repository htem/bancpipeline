"""Metric correctness tests for banc-alignment-run.py.

Hand-rolled fixtures with pencil-and-paper expected values. These are the
scoring primitives the optimisation loop rests on — silent regressions
here translate directly to silently-wrong type assignments in production.
"""
from __future__ import annotations

import numpy as np
import pytest


# ---- cosine_batch -----------------------------------------------------------

def test_cosine_batch_identity_rows(align):
    """Identical profile and centroid vectors score 1.0."""
    profiles = np.array([[1.0, 2.0, 3.0],
                         [0.5, 0.5, 0.0]], dtype=np.float32)
    centroids = profiles.copy()
    scores = align.cosine_batch(profiles, centroids)
    np.testing.assert_allclose(np.diag(scores), 1.0, atol=1e-6)


def test_cosine_batch_orthogonal_zero(align):
    """Disjoint-support vectors score 0."""
    profiles = np.array([[1.0, 0.0]], dtype=np.float32)
    centroids = np.array([[0.0, 1.0]], dtype=np.float32)
    assert align.cosine_batch(profiles, centroids)[0, 0] == pytest.approx(0.0)


def test_cosine_batch_handcalc(align):
    """One-off cosine matches pencil-and-paper."""
    p = np.array([[3.0, 4.0]], dtype=np.float32)       # ||p|| = 5
    c = np.array([[4.0, 3.0]], dtype=np.float32)       # ||c|| = 5
    # dot = 12 + 12 = 24  -> cos = 24/25
    assert align.cosine_batch(p, c)[0, 0] == pytest.approx(24 / 25, rel=1e-6)


def test_cosine_batch_zero_vector_safe(align):
    """Zero row must not produce NaN (norms_p[norms_p==0] = 1 path)."""
    profiles = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
    centroids = np.array([[1.0, 1.0, 1.0]], dtype=np.float32)
    score = align.cosine_batch(profiles, centroids)
    assert np.isfinite(score).all()
    assert score[0, 0] == pytest.approx(0.0)


# ---- weighted_jaccard_batch -------------------------------------------------

def test_weighted_jaccard_identity(align):
    """Identical profiles and centroids -> 1.0."""
    profiles = np.array([[2.0, 3.0, 0.0],
                         [1.0, 0.0, 4.0]], dtype=np.float32)
    centroids = profiles.copy()
    scores = align.weighted_jaccard_batch(profiles, centroids)
    np.testing.assert_allclose(np.diag(scores), 1.0, atol=1e-6)


def test_weighted_jaccard_handcalc(align):
    """sum(min) / sum(max) on a 3-d fixture."""
    p = np.array([[1.0, 3.0, 0.0]], dtype=np.float32)
    c = np.array([[2.0, 3.0, 1.0]], dtype=np.float32)
    # min = (1, 3, 0) -> 4;  max = (2, 3, 1) -> 6;  jac = 4/6
    assert align.weighted_jaccard_batch(p, c)[0, 0] == pytest.approx(4 / 6, rel=1e-6)


def test_weighted_jaccard_disjoint_supports(align):
    """Disjoint supports -> 0."""
    p = np.array([[1.0, 0.0]], dtype=np.float32)
    c = np.array([[0.0, 1.0]], dtype=np.float32)
    # min sum = 0; max sum = 2; -> 0/2 = 0
    assert align.weighted_jaccard_batch(p, c)[0, 0] == pytest.approx(0.0)


def test_weighted_jaccard_both_zero_safe(align):
    """If a (neuron, type) pair has total-zero support, maxs==0 path must not
    produce NaN — the implementation replaces maxs==0 with 1 to emit 0."""
    p = np.array([[0.0, 0.0]], dtype=np.float32)
    c = np.array([[0.0, 0.0]], dtype=np.float32)
    assert np.isfinite(align.weighted_jaccard_batch(p, c)).all()


def test_weighted_jaccard_chunking_consistent(align):
    """Result must be identical whether inputs fit in one chunk or cross the
    hardcoded chunk_size=200 boundary. Guards against index-off-by-one bugs
    in chunked loop."""
    rng = np.random.RandomState(0)
    profiles = rng.rand(250, 5).astype(np.float32)
    centroids = rng.rand(10, 5).astype(np.float32)
    full = align.weighted_jaccard_batch(profiles, centroids)
    # recompute in two halves and stitch
    top = align.weighted_jaccard_batch(profiles[:100], centroids)
    bot = align.weighted_jaccard_batch(profiles[100:], centroids)
    np.testing.assert_allclose(full[:100], top, atol=1e-6)
    np.testing.assert_allclose(full[100:], bot, atol=1e-6)
