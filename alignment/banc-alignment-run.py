"""
Optic Lobe Cross-Dataset Matching: Iterative Connectivity Alignment

Assigns FAFB cell_type labels to BANC optic lobe neurons by iteratively
comparing type-level connectivity profiles between the two datasets.

Algorithm:
  1. FAFB types are ground truth (fixed centroids)
  2. BANC neurons initialized from SeaTable seeds + soft NBLAST probabilities
  3. Each iteration: build type-connectivity profiles, score against FAFB
     centroids via cosine similarity, blend with NBLAST regularizer,
     apply soft reassignment with temperature annealing
  4. Early stopping: keep best-holdout-iteration assignments

Features:
  - Cosine similarity scoring
  - Soft assignments with temperature annealing (reduces oscillation)
  - Type capacity constraints: soft Sinkhorn-like scaling prevents
    over-assignment to rare types (FAFB count = max capacity)
  - NBLAST hard constraint: types with no NBLAST match > threshold excluded
  - NBLAST regularizer: adaptive alpha blends connectivity + NBLAST scores
  - Per-iteration metrics: holdout (overall, per-super_class, Mi1),
    NT consistency, type count correlation (BANC vs FAFB)

Input (from banc-alignment-prep.R --region optic-lobe):
  data/optic_lobe/banc_optic_{side}_meta.csv
  data/optic_lobe/fafb_optic_{side}_meta.csv
  data/optic_lobe/banc_optic_{side}_edgelist.feather
  data/optic_lobe/fafb_optic_{side}_edgelist.feather
  data/optic_lobe/banc_fafb_optic_{side}_nblast.csv
  data/optic_lobe/banc_optic_{side}_seeds.csv

Output:
  data/optic_lobe/banc_optic_{side}_alignment_results.csv
"""

import argparse
import functools
import os
import time
import numpy as np
import pandas as pd
import pyarrow.feather as pf
import scipy.sparse as sp
from collections import Counter
from scipy.stats import spearmanr

# Import alignment_path() from the local mirror in alignment/. The script is
# typically invoked from the repo root (cwd), so add this file's dir to sys.path.
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alignment_paths import alignment_path  # noqa: E402

# Force unbuffered output so progress is visible in real time
print = functools.partial(print, flush=True)


# Maps the legacy --file-prefix flag to the v2 grammar's `region` tag.
_PREFIX_TO_REGION = {"optic": "optic-lobe", "brain": "whole-brain"}


def _region_for(file_prefix):
    """Return v2 region tag for a legacy file_prefix. Unknown prefixes are
    used verbatim (caller must ensure they're a valid region tag)."""
    return _PREFIX_TO_REGION.get(file_prefix, file_prefix)


# Module-level: which `root_<NNN>` column to read from the seeds CSV.
# Auto-detected by detect_root_col() at startup based on the prep's seeds CSV
# header (e.g. "root_850" if data was prepped at v850, ROOT_COL if v888).
# Defaults to "root_id" until detection runs.
ROOT_COL = "root_id"


def detect_root_col(seeds_path):
    """Read the seeds CSV header and return the first column matching
    `root_<digits>`. Falls back to "root_id". Lets the script work
    against either v850 or v888 prep without re-prep, by aligning the
    column name with whatever the prep actually wrote."""
    cols = list(pd.read_csv(seeds_path, nrows=0).columns)
    for c in cols:
        if c.startswith("root_") and c[5:].isdigit():
            return c
    return "root_id"


def load_data(side, data_dir="data/optic_lobe", file_prefix="optic",
              target_name="fafb", banc_version="888", nblast_version="783",
              syn_source="v3"):
    """Load all prepared data for the given side."""
    print(f"=== Loading data for {side} side (prefix={file_prefix}) ===")

    region = _region_for(file_prefix)

    seeds_path = alignment_path("prep-seeds", query="banc", target=target_name,
                                region=region, side=side,
                                vq=banc_version, vt=nblast_version,
                                ext="csv", dir=data_dir)
    global ROOT_COL
    ROOT_COL = detect_root_col(seeds_path)
    print(f"  Auto-detected seeds root column: {ROOT_COL}")

    banc_meta = pd.read_csv(
        alignment_path("prep-banc-meta", query="banc", target=target_name,
                       region=region, side=side,
                       vq=banc_version, syn=syn_source,
                       ext="csv", dir=data_dir),
        dtype={"root_id": str, "supervoxel_id": str})
    target_meta = pd.read_csv(
        alignment_path("prep-target-meta", query="banc", target=target_name,
                       region=region, side=side,
                       vq=banc_version, vt=nblast_version,
                       ext="csv", dir=data_dir),
        dtype={"target_id": str})
    banc_el = pf.read_feather(
        alignment_path("prep-banc-edges", query="banc", target=target_name,
                       region=region, side=side,
                       vq=banc_version, syn=syn_source,
                       ext="feather", dir=data_dir))
    target_el = pf.read_feather(
        alignment_path("prep-target-edges", query="banc", target=target_name,
                       region=region, side=side,
                       vq=banc_version, vt=nblast_version,
                       ext="feather", dir=data_dir))
    nblast = pd.read_csv(
        alignment_path("prep-nblast", query="banc", target=target_name,
                       region=region, side=side,
                       vq=banc_version, vt=nblast_version,
                       ext="csv", dir=data_dir),
        dtype={"banc_id": str, "target_id": str})
    seeds = pd.read_csv(seeds_path, dtype={ROOT_COL: str})

    # Filter BANC to proofread / roughly_proofread only
    proofread = banc_meta["proofread"].fillna("").astype(str).str.upper() == "TRUE"
    roughly = banc_meta["roughly_proofread"].fillna("").astype(str).str.upper() == "TRUE"
    banc_pre = len(banc_meta)
    banc_meta = banc_meta[proofread | roughly].reset_index(drop=True)
    print(f"  BANC: {len(banc_meta)} proofread neurons (of {banc_pre}), {len(banc_el)} edges")
    print(f"  target: {len(target_meta)} neurons, {len(target_el)} edges")
    print(f"  NBLAST: {len(nblast)} scores for {nblast['banc_id'].nunique()} BANC neurons")
    print(f"  Seeds: {(seeds['tier']==1).sum()} T1, {(seeds['tier']==2).sum()} T2, "
          f"{(seeds['tier']==3).sum()} T3")

    return banc_meta, target_meta, banc_el, target_el, nblast, seeds


def build_type_index(target_meta):
    """Build type vocabulary from target-dataset typed neurons."""
    typed = target_meta.dropna(subset=["target_cell_type"])
    typed = typed[typed["target_cell_type"] != ""]
    types = sorted(typed["target_cell_type"].unique())
    type_to_idx = {t: i for i, t in enumerate(types)}
    print(f"  Type vocabulary: {len(types)} FAFB types")
    return types, type_to_idx


def build_adjacency(el_df, neuron_ids):
    """Build sparse adjacency matrix from edgelist.

    Uses normalized weights (norm column). Neurons not in neuron_ids
    are included as extra rows/cols (boundary neurons) to preserve
    cross-region connectivity signal, but only the first len(neuron_ids)
    rows correspond to the optic pool.

    Returns:
        adj: CSR matrix (n_all x n_all)
        id_to_idx: dict mapping neuron ID -> matrix index
        n_pool: number of optic pool neurons (first n_pool indices)
    """
    el_df = el_df.copy()
    el_df["pre"] = el_df["pre"].astype(str)
    el_df["post"] = el_df["post"].astype(str)

    # Pool neurons get indices 0..n_pool-1 (deduplicate, preserving order)
    seen = set()
    unique_ids = []
    for nid in neuron_ids:
        if nid not in seen:
            seen.add(nid)
            unique_ids.append(nid)
    neuron_ids = unique_ids
    pool_set = set(neuron_ids)
    id_to_idx = {nid: i for i, nid in enumerate(neuron_ids)}
    n_pool = len(neuron_ids)

    # Boundary neurons (in edgelist but not in pool) get indices n_pool..
    boundary_ids = sorted(
        set(el_df["pre"]).union(set(el_df["post"])) - pool_set
    )
    for nid in boundary_ids:
        id_to_idx[nid] = len(id_to_idx)
    n_all = len(id_to_idx)

    # Build sparse matrix
    rows = el_df["pre"].map(id_to_idx)
    cols = el_df["post"].map(id_to_idx)
    valid = rows.notna() & cols.notna()
    rows = rows[valid].astype(int).values
    cols = cols[valid].astype(int).values
    vals = el_df.loc[valid, "norm"].values.astype(np.float32)

    adj = sp.csr_matrix((vals, (rows, cols)), shape=(n_all, n_all))
    print(f"  Adjacency: {n_all} neurons ({n_pool} pool + {n_all - n_pool} boundary), "
          f"{adj.nnz} non-zeros")
    return adj, id_to_idx, n_pool


def build_type_matrix(meta_df, id_col, id_to_idx, type_to_idx, n_all):
    """Build sparse neuron-to-type assignment matrix (n_all x n_types).
    Used for FAFB centroid computation (hard assignments, fixed).
    """
    n_types = len(type_to_idx)
    fids = meta_df[id_col].astype(str).values
    cts = meta_df["target_cell_type"].values

    row_indices = np.array([id_to_idx.get(str(fid), -1) for fid in fids], dtype=np.int64)
    col_indices = np.array([type_to_idx.get(ct, -1) if pd.notna(ct) and ct != "" else -1
                            for ct in cts], dtype=np.int64)

    valid = (row_indices >= 0) & (col_indices >= 0) & (row_indices < n_all)
    rows = row_indices[valid]
    cols = col_indices[valid]

    return sp.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(n_all, n_types))


def build_target_centroids(target_adj, target_id_to_idx, target_meta, types, type_to_idx, n_target_pool):
    """Compute FAFB type centroids: average type-connectivity profile per type.
    Fully vectorized via sparse matrix operations.

    Returns:
        input_centroids: (n_types x n_types) dense array
        output_centroids: (n_types x n_types) dense array
    """
    n_types = len(types)
    n_all = target_adj.shape[0]

    target_type_mat = build_type_matrix(target_meta, "target_id", target_id_to_idx, type_to_idx, n_all)

    input_profiles = target_adj.T.tocsr() @ target_type_mat   # (n_all x n_types)
    output_profiles = target_adj @ target_type_mat             # (n_all x n_types)

    input_centroids = np.asarray((target_type_mat.T @ input_profiles).todense()).astype(np.float32)
    output_centroids = np.asarray((target_type_mat.T @ output_profiles).todense()).astype(np.float32)

    type_counts = np.asarray(target_type_mat.sum(axis=0)).flatten().astype(np.float32)
    type_counts[type_counts == 0] = 1
    input_centroids /= type_counts[:, None]
    output_centroids /= type_counts[:, None]

    # 2-hop centroids: neighbors' type profiles
    input_2hop = target_adj.T.tocsr() @ input_profiles   # (n_all x n_types)
    output_2hop = target_adj @ output_profiles            # (n_all x n_types)
    input_2hop_centroids = np.asarray((target_type_mat.T @ input_2hop).todense()).astype(np.float32)
    output_2hop_centroids = np.asarray((target_type_mat.T @ output_2hop).todense()).astype(np.float32)
    input_2hop_centroids /= type_counts[:, None]
    output_2hop_centroids /= type_counts[:, None]

    print(f"  FAFB centroids: {n_types} types, "
          f"1-hop density: {(input_centroids > 0).mean():.3f}, "
          f"2-hop density: {(input_2hop_centroids > 0).mean():.3f}")

    return input_centroids, output_centroids, input_2hop_centroids, output_2hop_centroids


def cosine_batch(profiles, centroids):
    """Compute cosine similarity between each row of profiles and centroids.

    Args:
        profiles: (n_neurons x d) dense array (non-negative)
        centroids: (n_types x d) dense array (non-negative)

    Returns:
        scores: (n_neurons x n_types) dense array, values in [0, 1]
    """
    norms_p = np.linalg.norm(profiles, axis=1, keepdims=True)
    norms_c = np.linalg.norm(centroids, axis=1, keepdims=True)
    norms_p[norms_p == 0] = 1  # zero vectors stay zero after normalization
    norms_c[norms_c == 0] = 1
    return (profiles / norms_p) @ (centroids / norms_c).T


def weighted_jaccard_batch(profiles, centroids):
    """Compute weighted Jaccard similarity between each row of profiles and centroids.

    weighted_jaccard(A, B) = sum(min(a_i, b_i)) / sum(max(a_i, b_i))

    Args:
        profiles: (n_neurons x d) dense array
        centroids: (n_types x d) dense array

    Returns:
        scores: (n_neurons x n_types) dense array
    """
    n_neurons = profiles.shape[0]
    n_types = centroids.shape[0]
    scores = np.zeros((n_neurons, n_types), dtype=np.float32)

    chunk_size = 200
    for start in range(0, n_neurons, chunk_size):
        end = min(start + chunk_size, n_neurons)
        chunk = profiles[start:end]
        c3 = chunk[:, None, :]
        t3 = centroids[None, :, :]
        mins = np.minimum(c3, t3).sum(axis=2)
        maxs = np.maximum(c3, t3).sum(axis=2)
        maxs[maxs == 0] = 1
        scores[start:end] = mins / maxs

    return scores


def build_nblast_type_scores(nblast, target_meta, banc_id_to_idx, type_to_idx,
                             n_banc_pool, threshold=0.1):
    """Precompute per-neuron per-type max NBLAST scores and allowed type mask.

    For each BANC neuron, finds the maximum NBLAST score to any FAFB neuron
    of each type. For neurons WITH any NBLAST data, only types with at least
    one match scoring >= threshold are allowed; all other types (including
    those with zero NBLAST hits) are blocked. Neurons with NO NBLAST data
    at all are unconstrained (all types allowed).

    Returns:
        nblast_type_scores: (n_pool x n_types) float32
        nblast_allowed: (n_pool x n_types) bool
        has_nblast: (n_pool,) bool
    """
    n_types = len(type_to_idx)
    nblast_type_scores = np.zeros((n_banc_pool, n_types), dtype=np.float32)
    nblast_type_covered = np.zeros((n_banc_pool, n_types), dtype=bool)
    has_nblast = np.zeros(n_banc_pool, dtype=bool)

    # FAFB ID -> type index lookup
    target_type_lookup = {}
    for fid, ct in zip(target_meta["target_id"].astype(str), target_meta["target_cell_type"]):
        if pd.notna(ct) and ct in type_to_idx:
            target_type_lookup[fid] = type_to_idx[ct]

    # Map NBLAST rows to matrix indices (vectorized via pandas)
    banc_idx = nblast["banc_id"].astype(str).map(banc_id_to_idx)
    type_idx = nblast["target_id"].astype(str).map(target_type_lookup)
    nb_scores = nblast["nblast_score"].values.astype(np.float32)

    # Mark has_nblast for all neurons with any NBLAST entry
    valid_banc = banc_idx.notna()
    all_banc_nb = banc_idx[valid_banc].astype(int).values
    all_banc_nb = all_banc_nb[all_banc_nb < n_banc_pool]
    has_nblast[np.unique(all_banc_nb)] = True

    # Filter to entries where both BANC and FAFB type are resolved
    valid = valid_banc & type_idx.notna()
    b_idx = banc_idx[valid].astype(int).values
    t_idx = type_idx[valid].astype(int).values
    scores = nb_scores[valid.values]

    # Further filter to pool neurons
    pool_mask = b_idx < n_banc_pool
    b_idx = b_idx[pool_mask]
    t_idx = t_idx[pool_mask]
    scores = scores[pool_mask]

    # Max aggregation and coverage
    np.maximum.at(nblast_type_scores, (b_idx, t_idx), scores)
    nblast_type_covered[b_idx, t_idx] = True

    # Allowed type mask:
    # - Neurons WITH NBLAST data: only types with score >= threshold are allowed.
    #   Types with no NBLAST hits at all are blocked (same as below-threshold).
    # - Neurons WITHOUT any NBLAST data: all types allowed (unconstrained).
    nblast_allowed = nblast_type_scores >= threshold         # strict: need positive evidence
    nblast_allowed[~has_nblast] = True                       # unconstrained if no NBLAST at all

    n_with = has_nblast.sum()
    avg_allowed = nblast_allowed[has_nblast].mean() * n_types if n_with > 0 else n_types
    print(f"  NBLAST type scores: {n_with}/{n_banc_pool} neurons with data")
    print(f"  Avg allowed types per NBLAST neuron: {avg_allowed:.0f}/{n_types}")

    return nblast_type_scores, nblast_allowed, has_nblast


def load_forbidden_matches(path, banc_id_to_idx, type_to_idx, n_banc_pool):
    """Load reviewed (root_id, cell_type) pairs that are forbidden assignments.

    The CSV must have two columns: ``root_id`` (BANC neuron, matching the
    meta CSV's root_id used as align.py's pool key) and ``cell_type`` (FAFB
    cell type to veto). Out-of-pool ids and unknown types are dropped with
    a count log.

    Returns:
        forbidden_neuron_arr (int32[N]): pool indices of forbidden pairs
        forbidden_type_arr   (int32[N]): type indices of forbidden pairs
        forbidden_pairs_set  (set[(int,int)]): for fast post-hoc verification
        forbidden_per_neuron (dict[int, set[int]]): for greedy-loop checks
    """
    print(f"\n=== Loading forbidden matches: {path} ===")
    df = pd.read_csv(path, dtype={"root_id": str})
    if "root_id" not in df.columns or "cell_type" not in df.columns:
        raise ValueError(
            f"Forbidden matches CSV must have columns 'root_id' and 'cell_type'; "
            f"got {df.columns.tolist()}")
    n_in = len(df)
    df = df.dropna(subset=["root_id", "cell_type"])
    df = df.drop_duplicates(subset=["root_id", "cell_type"])

    neuron_idx, type_idx = [], []
    n_oop = 0
    n_unk_type = 0
    for rid, ct in zip(df["root_id"].astype(str), df["cell_type"].astype(str)):
        if rid not in banc_id_to_idx:
            n_oop += 1
            continue
        ni = banc_id_to_idx[rid]
        if ni >= n_banc_pool:
            n_oop += 1
            continue
        if ct not in type_to_idx:
            n_unk_type += 1
            continue
        neuron_idx.append(ni)
        type_idx.append(type_to_idx[ct])

    forbidden_neuron_arr = np.array(neuron_idx, dtype=np.int32)
    forbidden_type_arr = np.array(type_idx, dtype=np.int32)
    forbidden_pairs_set = set(zip(neuron_idx, type_idx))
    forbidden_per_neuron = {}
    for ni, ti in zip(neuron_idx, type_idx):
        forbidden_per_neuron.setdefault(int(ni), set()).add(int(ti))

    n_unique_neurons = len(forbidden_per_neuron)
    print(f"  Read {n_in} rows -> kept {len(forbidden_neuron_arr)} pairs "
          f"({n_unique_neurons} unique neurons)")
    if n_oop or n_unk_type:
        print(f"  Skipped: {n_oop} out-of-pool ids, {n_unk_type} unknown types")
    return (forbidden_neuron_arr, forbidden_type_arr,
            forbidden_pairs_set, forbidden_per_neuron)


def _build_results_df(hard_assignments, scores, banc_pool_ids, types,
                      n_banc_pool, seeds, anchor_mask, nblast, target_meta,
                      verbose=True):
    # Per-neuron assigned-type score.
    assigned_scores = np.zeros(n_banc_pool, dtype=np.float32)
    for i in range(n_banc_pool):
        if hard_assignments[i] >= 0:
            assigned_scores[i] = scores[i, hard_assignments[i]]
    # Runner-up: best score EXCLUDING the assigned type.
    scores_copy = scores.copy()
    for i in range(n_banc_pool):
        if hard_assignments[i] >= 0:
            scores_copy[i, hard_assignments[i]] = -1e9
    runner_up_scores = scores_copy.max(axis=1)
    runner_up_idx = scores_copy.argmax(axis=1)

    seeds_tier = seeds.set_index(ROOT_COL)["tier"].to_dict()

    rows = []
    for i in range(n_banc_pool):
        rid = banc_pool_ids[i]
        tier_val = int(seeds_tier.get(rid, 3))
        if hard_assignments[i] >= 0:
            assigned_type = types[hard_assignments[i]]
            ru_type = types[runner_up_idx[i]] if runner_up_scores[i] > -1e8 else ""
            ru_score = float(runner_up_scores[i]) if runner_up_scores[i] > -1e8 else 0.0
        else:
            assigned_type = ""
            ru_type = ""
            ru_score = 0.0
        rows.append({
            ROOT_COL: rid,
            "assigned_cell_type": assigned_type,
            "best_score": float(assigned_scores[i]),
            "runner_up_type": ru_type,
            "runner_up_score": ru_score,
            "confidence": float(assigned_scores[i]) - ru_score,
            "tier": tier_val,
            "is_anchor": bool(anchor_mask[i]),
        })
    df = pd.DataFrame(rows)

    if len(nblast) > 0:
        if verbose:
            print("  Finding best FAFB matches...")
        target_type_lookup = target_meta.set_index("target_id")["target_cell_type"].to_dict()
        nblast_wt = nblast.copy()
        nblast_wt["target_type"] = nblast_wt["target_id"].map(target_type_lookup)
        nblast_wt = nblast_wt.dropna(subset=["target_type"])
        idx_max = nblast_wt.groupby(["banc_id", "target_type"])["nblast_score"].idxmax()
        best_matches = nblast_wt.loc[idx_max, ["banc_id", "target_type", "target_id"]].copy()
        best_matches_lookup = best_matches.set_index(["banc_id", "target_type"])["target_id"].to_dict()
        df["best_target_match"] = [
            best_matches_lookup.get((rid, atype), "")
            for rid, atype in zip(df[ROOT_COL], df["assigned_cell_type"])
        ]
    else:
        df["best_target_match"] = ""
    return df


def run_alignment(banc_meta, target_meta, banc_el, target_el, nblast, seeds,
                  data_dir="data/optic_lobe",
                  manual_labels=None,
                  random_holdout=None,
                  stratified_holdout=None,
                  max_iter=30, convergence_pct=1.0,
                  tau_start=2.0, tau_end=0.5,
                  alpha_start=0.3, alpha_end=0.8,
                  nblast_threshold=0.25,
                  capacity_scale=1.5,
                  metric="cosine",
                  prior_weight=0.0,
                  hop2_weight=0.5,
                  ensemble_top_k=30,
                  ensemble_blend=0.5,
                  nt_weight=0.5,
                  ind_weight=0.0,
                  ind_metric="auto",
                  chunked=True,
                  conf_gate=True,
                  conf_gate_threshold=0.1,
                  bilateral=False,
                  forbidden_matches=None,
                  apply_soma_rule=True,
                  checkpoint_path=None,
                  checkpoint_every=0,
                  momentum=0.0,
                  skip_stage1=False,
                  # Path-grammar context (v2). Passed through so capacity_file
                  # resolution uses the same alignment_path() params the prep
                  # script wrote with.
                  region="optic-lobe",
                  side="right",
                  target_name="fafb",
                  banc_version="888",
                  nblast_version="783",
                  syn_source="v3"):
    """Run iterative connectivity alignment with soft assignments.

    Args:
        manual_labels: Use known cell_type labels as near-fixed anchors.
            None/False = don't use (holdout stays hidden).
            True = use ALL holdout labels as near-fixed anchors.
            List of type names = only use those specific types (e.g. ["Mi1"]).
            Manual-label neurons count toward type capacity.
        capacity_scale: Allow each type up to capacity_scale * target_count neurons.
        metric: Similarity metric — "cosine", "jaccard", or "ensemble".
        prior_weight: Weight for type frequency log-prior.
        hop2_weight: Weight for 2-hop connectivity profiles (0=disable, 1=equal to 1-hop).
        ensemble_top_k: Number of cosine candidates for Jaccard re-ranking in ensemble mode.
        ensemble_blend: Weight of cosine vs Jaccard in ensemble (0=pure Jaccard, 1=pure cosine).
        ind_weight: Weight for individual FAFB profile scoring (0=disabled).
            Blends centroid and best individual scores per type:
            (1-w)*centroid + w*best_individual. Values 0.1-0.3 recommended.
        ind_metric: Metric for individual profile scoring.
            "auto"=match main metric, "cosine", "jaccard", "ensemble".
        chunked: Compute profiles and scores in neuron-chunks to reduce peak memory.
            Essential for whole-brain runs. No accuracy difference.
        conf_gate: Only include confident neurons in type profiles. Neurons with
            max soft_prob < conf_gate_threshold contribute reduced weight, preventing
            noisy early assignments from polluting neighbors' profiles.
        conf_gate_threshold: Min soft_prob confidence to get full weight (default 0.1).
        bilateral: Enforce bilateral consistency — penalize left-right type count
            imbalance. Requires 'side' column in banc_meta. For whole-brain bilateral runs.
        forbidden_matches: Optional path to a CSV with columns
            (root_id, cell_type) listing reviewed false-positive pairs.
            When provided, those (neuron, type) assignments are vetoed via
            score masking. Anchored neurons win over forbidden pairs (the
            forbidden entry is dropped with a warning if it conflicts).
        apply_soma_rule: If True (default), enforce a soma-presence rule:
            BANC neurons that have a nucleus (nucleus_id present) cannot be
            assigned to FAFB types whose super_class contains "sensory" or
            to L1-L5 / R7 / R8 (these neurons' somata sit outside the CNS
            in their source ganglion). Anchors are exempt — if a neuron is
            already anchored to a soma-forbidden type, the rule is skipped
            for that neuron. Requires `nucleus_id` column in banc_meta CSV.
    """

    print("\n=== Building data structures ===")

    # Type vocabulary from FAFB
    types, type_to_idx = build_type_index(target_meta)
    n_types = len(types)

    # FAFB adjacency and centroids (fixed)
    # Deduplicate pool IDs (bilateral mode can have duplicates)
    target_pool_ids = target_meta["target_id"].astype(str).unique()
    target_adj, target_id_to_idx, n_target_pool = build_adjacency(target_el, target_pool_ids)
    # hop2_weight passed as parameter (default 0.5)

    print("  Computing FAFB type centroids (1-hop + 2-hop)...")
    input_1h, output_1h, input_2h, output_2h = build_target_centroids(
        target_adj, target_id_to_idx, target_meta, types, type_to_idx, n_target_pool)
    target_centroids = np.hstack([input_1h, output_1h,
                                hop2_weight * input_2h, hop2_weight * output_2h])

    # Individual FAFB neuron profiles materialised LATER, only for the FAFB
    # neurons that actually appear in ind_target_arr (the per-BANC top-k NBLAST
    # matches). At WB scale (n_target_pool ~156k, n_types ~8800) materialising
    # the full (n_target_pool, 4*n_types) dense array is ~22 GB and OOMs; the
    # NBLAST-reachable subset is typically 3-5x smaller.
    target_ind_profiles = None
    target_ind_profiles_normed = None
    target_full_to_local = None
    n_target_all = target_adj.shape[0]
    # Always need target_type_mat for centroid computation
    target_type_mat = build_type_matrix(target_meta, "target_id", target_id_to_idx, type_to_idx, n_target_all)
    if ind_weight > 0:
        print("  Will compute individual FAFB neuron profiles after ind index built")
    else:
        print("  Skipping individual FAFB profiles (ind_weight=0)")

    # FAFB type counts for type-count correlation metric and capacity constraints
    target_type_counts = np.zeros(n_types, dtype=int)
    for ct in target_meta["target_cell_type"]:
        if pd.notna(ct) and ct in type_to_idx:
            target_type_counts[type_to_idx[ct]] += 1

    # Log-prior from type frequency (zero-mean so common types get bonus)
    log_prior = np.log(target_type_counts.astype(np.float32) + 1)
    log_prior -= log_prior.mean()
    if prior_weight > 0:
        print(f"  Type frequency prior: weight={prior_weight}, "
              f"range [{log_prior.min():.2f}, {log_prior.max():.2f}]")

    # Type capacity: per-type from FAFB L-R variability (if available), else fixed multiplier
    capacity_file = alignment_path("prep-capacity", query="banc", target=target_name,
                                    region=region, side=side,
                                    vq=banc_version, vt=nblast_version,
                                    ext="csv", dir=data_dir)
    type_capacity_left = None
    type_capacity_right = None
    if os.path.exists(capacity_file):
        cap_df = pd.read_csv(capacity_file)
        cap_lookup = dict(zip(cap_df["cell_type"], cap_df["capacity"]))
        type_capacity = np.array([cap_lookup.get(t, max(target_type_counts[i] * capacity_scale, 1))
                                  for i, t in enumerate(types)], dtype=np.float32)
        print(f"  Type capacity: per-type from FAFB L-R variability "
              f"(range {type_capacity.min():.0f}-{type_capacity.max():.0f}, "
              f"mean ratio {type_capacity.mean() / max(target_type_counts.mean(), 1):.2f}x)")
        # Per-side hard caps for bilateral mode (cap_left/cap_right columns,
        # added by prep.R; absent in older capacity files).
        if "cap_left" in cap_df.columns and "cap_right" in cap_df.columns:
            cap_l_lookup = dict(zip(cap_df["cell_type"], cap_df["cap_left"]))
            cap_r_lookup = dict(zip(cap_df["cell_type"], cap_df["cap_right"]))
            type_capacity_left = np.array(
                [cap_l_lookup.get(t, type_capacity[i]) for i, t in enumerate(types)],
                dtype=np.float32)
            type_capacity_right = np.array(
                [cap_r_lookup.get(t, type_capacity[i]) for i, t in enumerate(types)],
                dtype=np.float32)
            print(f"  Per-side caps loaded: "
                  f"left range {type_capacity_left.min():.0f}-{type_capacity_left.max():.0f}, "
                  f"right range {type_capacity_right.min():.0f}-{type_capacity_right.max():.0f}")
    else:
        type_capacity = (target_type_counts * capacity_scale).astype(np.float32)
        type_capacity[type_capacity < 1] = 1
        print(f"  Type capacity: {capacity_scale}x FAFB counts "
              f"(range {type_capacity.min():.0f}-{type_capacity.max():.0f})")

    # BANC adjacency
    banc_pool_ids = banc_meta["root_id"].astype(str).unique()
    banc_adj, banc_id_to_idx, n_banc_pool = build_adjacency(banc_el, banc_pool_ids)
    n_banc_all = banc_adj.shape[0]
    banc_adj_t = banc_adj.T.tocsr()

    # NBLAST type scores and constraints
    print("  Computing NBLAST type scores...")
    nblast_type_scores, nblast_allowed, has_nblast = build_nblast_type_scores(
        nblast, target_meta, banc_id_to_idx, type_to_idx, n_banc_pool, nblast_threshold)

    # Forbidden (root_id, cell_type) pairs from manual review (optional).
    # Loaded here because we need banc_id_to_idx and type_to_idx; anchor
    # conflicts are resolved later, after anchor_mask is built.
    forbidden_neuron_arr = np.zeros(0, dtype=np.int32)
    forbidden_type_arr = np.zeros(0, dtype=np.int32)
    forbidden_pairs_set = set()
    forbidden_per_neuron = {}
    if forbidden_matches is not None:
        (forbidden_neuron_arr, forbidden_type_arr,
         forbidden_pairs_set, forbidden_per_neuron) = load_forbidden_matches(
            forbidden_matches, banc_id_to_idx, type_to_idx, n_banc_pool)

    # ---------------------------------------------------------------
    # Initialize soft assignment probabilities (n_pool x n_types)
    # ---------------------------------------------------------------
    print("\n=== Initializing BANC assignments (soft) ===")
    soft_probs = np.zeros((n_banc_pool, n_types), dtype=np.float32)
    anchor_mask = np.zeros(n_banc_pool, dtype=bool)

    # Holdout modes: random, stratified, or normal (tier-based)
    random_holdout_indices = None
    random_holdout_true_types = None
    holdout_fraction = random_holdout or stratified_holdout
    if stratified_holdout is not None:
        from alignment_splits import stratified_holdout_split
        seed_list, holdout_list = stratified_holdout_split(
            seeds, banc_id_to_idx, type_to_idx.keys(),
            stratified_holdout, seed=42)
        n_t1 = 0
        for idx, ct, rid in seed_list:
            if idx < n_banc_pool:
                soft_probs[idx] = 0
                soft_probs[idx, type_to_idx[ct]] = 1.0
                anchor_mask[idx] = True
                n_t1 += 1
        random_holdout_indices = np.array([x[0] for x in holdout_list], dtype=np.int32)
        random_holdout_true_types = np.array([x[1] for x in holdout_list])
        # Re-group for the "no seeds" diagnostic
        typed_by_type = {}
        for _, row in seeds.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in banc_id_to_idx and pd.notna(ct) and ct in type_to_idx:
                typed_by_type.setdefault(ct, []).append(rid)
        n_types_with_no_seed = sum(
            1 for ct, neurons in typed_by_type.items()
            if len(neurons) <= max(1, round(stratified_holdout * len(neurons)))
        )
        print(f"  Stratified holdout: {n_t1} seeds, {len(holdout_list)} holdout "
              f"({len(typed_by_type)} types, {n_types_with_no_seed} with no seeds), seed=42")
    elif random_holdout is not None:
        # Random: uniform permutation (unchanged)
        all_typed = []
        for _, row in seeds.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in banc_id_to_idx and pd.notna(ct) and ct in type_to_idx:
                all_typed.append((banc_id_to_idx[rid], ct, rid))
        rng = np.random.RandomState(42)
        n_typed = len(all_typed)
        perm = rng.permutation(n_typed)
        n_holdout = int(random_holdout * n_typed)
        holdout_perm = perm[:n_holdout]
        seed_perm = perm[n_holdout:]
        n_t1 = 0
        for i in seed_perm:
            idx, ct, rid = all_typed[i]
            if idx < n_banc_pool:
                soft_probs[idx] = 0
                soft_probs[idx, type_to_idx[ct]] = 1.0
                anchor_mask[idx] = True
                n_t1 += 1
        random_holdout_indices = np.array([all_typed[i][0] for i in holdout_perm], dtype=np.int32)
        random_holdout_true_types = np.array([all_typed[i][1] for i in holdout_perm])
        print(f"  Random holdout: {n_t1} seeds ({100*(1-random_holdout):.0f}%), "
              f"{len(holdout_perm)} holdout ({100*random_holdout:.0f}%), seed=42")
    else:
        # Normal mode: Tier 1 anchors (non-holdout)
        tier1 = seeds[(seeds["tier"] == 1) & (seeds["is_holdout"] == False)]
        n_t1 = 0
        for _, row in tier1.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in banc_id_to_idx and ct in type_to_idx:
                idx = banc_id_to_idx[rid]
                if idx < n_banc_pool:
                    soft_probs[idx] = 0
                    soft_probs[idx, type_to_idx[ct]] = 1.0
                    anchor_mask[idx] = True
                    n_t1 += 1

    # Manual labels: promote holdout neurons to near-fixed anchors
    manual_label_mask = np.zeros(n_banc_pool, dtype=bool)
    n_manual = 0
    if manual_labels:
        holdout_seeds = seeds[seeds["is_holdout"] == True]
        for _, row in holdout_seeds.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if not (rid in banc_id_to_idx and pd.notna(ct) and ct in type_to_idx):
                continue
            # If manual_labels is a list, only include specified types
            if isinstance(manual_labels, list) and ct not in manual_labels:
                continue
            idx = banc_id_to_idx[rid]
            if idx < n_banc_pool:
                soft_probs[idx] = 0
                soft_probs[idx, type_to_idx[ct]] = 1.0
                anchor_mask[idx] = True
                manual_label_mask[idx] = True
                n_manual += 1
        print(f"  Manual labels: {n_manual} neurons promoted to anchors "
              f"({'all types' if manual_labels is True else ', '.join(manual_labels)})")

    # Precompute anchor indices and their type indices for fast restore
    anchor_indices = np.where(anchor_mask)[0]
    anchor_type_indices = soft_probs[anchor_mask].argmax(axis=1)

    # Resolve forbidden vs anchor conflicts: anchor wins. A forbidden pair
    # whose (neuron, type) matches an anchor's actual type is contradictory
    # by construction — drop the forbidden entry and warn so it can be
    # cleaned up upstream in SeaTable.
    if len(forbidden_neuron_arr) > 0:
        anchor_type_for_idx = {int(ai): int(ti)
                               for ai, ti in zip(anchor_indices, anchor_type_indices)}
        keep = np.ones(len(forbidden_neuron_arr), dtype=bool)
        for k in range(len(forbidden_neuron_arr)):
            ni_i = int(forbidden_neuron_arr[k])
            ti_i = int(forbidden_type_arr[k])
            if anchor_type_for_idx.get(ni_i, -1) == ti_i:
                keep[k] = False
                forbidden_pairs_set.discard((ni_i, ti_i))
                if ni_i in forbidden_per_neuron:
                    forbidden_per_neuron[ni_i].discard(ti_i)
                    if not forbidden_per_neuron[ni_i]:
                        del forbidden_per_neuron[ni_i]
        n_dropped = int((~keep).sum())
        if n_dropped > 0:
            forbidden_neuron_arr = forbidden_neuron_arr[keep]
            forbidden_type_arr = forbidden_type_arr[keep]
            print(f"  Forbidden conflicts: {n_dropped} pairs matched anchor/manual GT "
                  f"(anchor wins, dropped)")
        print(f"  Forbidden mask active: {len(forbidden_neuron_arr)} pairs over "
              f"{len(forbidden_per_neuron)} neurons")

    # ---------------------------------------------------------------
    # Soma-presence rule: certain FAFB types must only be assigned to
    # BANC neurons that LACK a nucleus_id (their somata sit outside the
    # CNS — photoreceptors in the retina, lamina monopolar cells in the
    # lamina, sensory afferents in peripheral ganglia). Anchors exempt.
    # The mask is rectangular (rows=neurons-with-nucleus, cols=no-soma
    # types) so we keep it as 1-D row/col index arrays and apply via
    # np.ix_ at the same hooks as the forbidden mask.
    # ---------------------------------------------------------------
    NO_SOMA_EXPLICIT_TYPES = ("R7", "R8", "L1", "L2", "L3", "L4", "L5")
    soma_rule_row_idx = np.zeros(0, dtype=np.int32)
    soma_rule_col_idx = np.zeros(0, dtype=np.int32)
    soma_rule_col_set = set()
    if apply_soma_rule:
        if "nucleus_id" not in banc_meta.columns:
            print(f"\n  Soma rule: SKIPPED (banc_meta has no 'nucleus_id' column)")
        else:
            print(f"\n=== Soma-presence rule ===")
            # has_nucleus aligned to banc_pool_ids order
            nucleus_lookup = dict(zip(
                banc_meta["root_id"].astype(str),
                banc_meta["nucleus_id"].fillna("").astype(str)))
            has_nucleus = np.zeros(n_banc_pool, dtype=bool)
            for i in range(n_banc_pool):
                nv = nucleus_lookup.get(banc_pool_ids[i], "").strip()
                # Treat NA, empty, "0", "nan" all as "no nucleus"
                has_nucleus[i] = nv not in ("", "0", "nan", "NA", "<NA>", "None")
            n_with_nuc = int(has_nucleus.sum())
            n_without_nuc = n_banc_pool - n_with_nuc
            print(f"  BANC pool: {n_with_nuc} with nucleus, "
                  f"{n_without_nuc} without nucleus")

            # Column mask: types in explicit list OR FAFB super_class contains
            # "sensory". A type is flagged if ANY FAFB neuron of that type has
            # a sensory super_class (types should be coherent on super_class).
            soma_required_types = set()
            if "target_super_class" in target_meta.columns:
                target_ct = target_meta["target_cell_type"].astype(str).values
                target_sc = target_meta["target_super_class"].fillna("").astype(str) \
                    .str.lower().values
                is_sens = np.array(["sensory" in s for s in target_sc], dtype=bool)
                for ct in np.unique(target_ct[is_sens]):
                    if ct in type_to_idx:
                        soma_required_types.add(ct)
            else:
                print(f"  WARNING: target_meta has no 'target_super_class' column; "
                      f"only explicit types will be flagged")
            n_sensory_types = len(soma_required_types)
            for ct in NO_SOMA_EXPLICIT_TYPES:
                if ct in type_to_idx:
                    soma_required_types.add(ct)
            n_total_types = len(soma_required_types)
            print(f"  No-soma types: {n_total_types} "
                  f"({n_sensory_types} from sensory super_class, "
                  f"{n_total_types - n_sensory_types} added from "
                  f"L1-L5/R7/R8 not already covered)")

            soma_rule_col_idx = np.array(
                sorted(type_to_idx[ct] for ct in soma_required_types),
                dtype=np.int32)
            soma_rule_col_set = set(int(c) for c in soma_rule_col_idx)
            soma_rule_row_idx = np.where(has_nucleus & ~anchor_mask)[0] \
                .astype(np.int32)
            print(f"  Soma rule active: {len(soma_rule_row_idx)} BANC neurons "
                  f"(have nucleus, not anchor) x {len(soma_rule_col_idx)} types")

            # Anchor exceptions: anchors with nucleus already assigned to a
            # no-soma type. The rule is skipped for them (anchor wins). Surface
            # them so they can be cleaned up upstream if unintended.
            if len(soma_rule_col_idx) > 0 and len(anchor_indices) > 0:
                n_anchor_violations = 0
                examples = []
                for ai, ti in zip(anchor_indices, anchor_type_indices):
                    if has_nucleus[ai] and int(ti) in soma_rule_col_set:
                        n_anchor_violations += 1
                        if len(examples) < 5:
                            examples.append((banc_pool_ids[ai], types[int(ti)]))
                if n_anchor_violations > 0:
                    print(f"  Anchor exceptions: {n_anchor_violations} anchors "
                          f"have a nucleus AND a no-soma type (kept as-is)")
                    for rid, tn in examples:
                        print(f"    {rid} -> {tn}")

    # Tier 2: initialize with soft NBLAST probs where available
    tier2 = seeds[seeds["tier"] == 2]
    n_t2 = 0
    for _, row in tier2.iterrows():
        rid = str(row[ROOT_COL])
        if rid not in banc_id_to_idx:
            continue
        idx = banc_id_to_idx[rid]
        if idx >= n_banc_pool or anchor_mask[idx]:
            continue
        # Soft init from NBLAST type scores
        if has_nblast[idx]:
            nb_scores = nblast_type_scores[idx].copy()
            nb_scores[nb_scores < 0] = 0
            if nb_scores.sum() > 0:
                soft_probs[idx] = nb_scores / nb_scores.sum()
                n_t2 += 1
                continue
        # Fallback: hard assignment from seed
        ct = row["cell_type"]
        if pd.notna(ct) and ct in type_to_idx:
            soft_probs[idx, type_to_idx[ct]] = 1.0
            n_t2 += 1

    n_assigned = (soft_probs.max(axis=1) > 0).sum()
    print(f"  Initialized: {n_t1} anchors, {n_t2} NBLAST-soft, "
          f"{n_banc_pool - n_assigned} unassigned ({n_assigned}/{n_banc_pool} total)")

    # ---------------------------------------------------------------
    # NT consistency checker
    # ---------------------------------------------------------------
    banc_nt_map = {}
    for _, row in banc_meta.iterrows():
        rid = str(row["root_id"])
        nt = row.get("top_nt", None)
        ct = row.get("cell_type", "")
        if pd.notna(ct) and ct in ("R1-6", "R7", "R8"):
            nt = "histamine"
        if rid in banc_id_to_idx and pd.notna(nt) and str(nt).strip():
            banc_nt_map[banc_id_to_idx[rid]] = str(nt).lower().strip()

    target_type_nt = {}
    type_nt_votes = {}
    for _, row in target_meta.iterrows():
        ct = row.get("target_cell_type", None)
        nt = row.get("target_top_nt", None)
        if pd.notna(ct) and ct != "" and ct in ("R1-6", "R7", "R8"):
            nt = "histamine"
        if pd.notna(ct) and ct != "" and pd.notna(nt) and str(nt).strip():
            type_nt_votes.setdefault(ct, []).append(str(nt).lower().strip())
    for ct, votes in type_nt_votes.items():
        target_type_nt[ct] = Counter(votes).most_common(1)[0][0]

    # Soft NT compatibility scores: (n_pool x n_types)
    # Match → confidence score, mismatch → 1-confidence, no data → 1.0 (neutral)
    invalid_nt = {"", "na", "nan", "none", "unclear", "unknown"}
    banc_nt_conf = {}  # idx -> (nt_string, confidence_float)
    for _, row in banc_meta.iterrows():
        rid = str(row["root_id"])
        if rid not in banc_id_to_idx:
            continue
        idx = banc_id_to_idx[rid]
        if idx >= n_banc_pool:
            continue
        nt = row.get("top_nt", None)
        ct = row.get("cell_type", "")
        if pd.notna(ct) and ct in ("R1-6", "R7", "R8"):
            nt = "histamine"
        if pd.isna(nt) or str(nt).strip().lower() in invalid_nt:
            continue
        # Get confidence: try top_nt_conf, neurotransmitter_score, fall back to 0.7
        conf = row.get("top_nt_conf", None)
        if pd.isna(conf):
            conf = row.get("neurotransmitter_score", None)
        if pd.isna(conf):
            conf = 0.7  # default moderate confidence
        else:
            conf = float(conf)
            conf = max(0.1, min(conf, 1.0))  # clamp to [0.1, 1.0]
        banc_nt_conf[idx] = (str(nt).lower().strip(), conf)

    # Build soft NT score matrix, controlled by nt_weight:
    # nt_weight=0 → all 1.0 (NT ignored)
    # nt_weight=0.5 → moderate soft scoring
    # nt_weight=1.0 → hard constraint (match=1, mismatch→0)
    nt_compat = np.ones((n_banc_pool, n_types), dtype=np.float32)
    n_nt_scored = 0
    if nt_weight > 0:
        for idx, (nt, conf) in banc_nt_conf.items():
            for tidx in range(n_types):
                tname = types[tidx]
                if tname not in target_type_nt:
                    continue  # no FAFB consensus → neutral (1.0)
                if target_type_nt[tname] == nt:
                    nt_compat[idx, tidx] = 1.0  # match: no penalty
                else:
                    # mismatch: scale penalty by confidence and nt_weight
                    # nt_weight=0 → 1.0 (no penalty), nt_weight=1,conf=1 → 0.01
                    nt_compat[idx, tidx] = max(1.0 - conf * nt_weight, 0.01)
            n_nt_scored += 1
    print(f"  NT soft scoring: {n_nt_scored}/{n_banc_pool} neurons, nt_weight={nt_weight}")

    def eval_nt(hard_assignments):
        n_check = n_match = 0
        for idx, tidx in enumerate(hard_assignments):
            if tidx < 0 or idx not in banc_nt_map:
                continue
            if types[tidx] not in target_type_nt:
                continue
            n_check += 1
            if banc_nt_map[idx] == target_type_nt[types[tidx]]:
                n_match += 1
        return (100 * n_match / n_check, n_check) if n_check else (0.0, 0)

    # ---------------------------------------------------------------
    # Holdout validation
    # ---------------------------------------------------------------
    sc_lookup = banc_meta.set_index("root_id")["super_class"].to_dict()

    if random_holdout_indices is not None:
        # Random/stratified holdout mode: use the pre-computed split
        holdout_indices = random_holdout_indices
        holdout_true_types = np.array([type_to_idx[t] for t in random_holdout_true_types])
        holdout_super_classes = np.array([
            sc_lookup.get(banc_pool_ids[i], "unknown") for i in holdout_indices])
        holdout_super_classes = np.where(pd.isna(holdout_super_classes), "unknown", holdout_super_classes)
    else:
        # Normal mode: holdout = intrinsic neurons with known cell_type
        holdout = seeds[seeds["is_holdout"] == True].copy()
        holdout_indices = []
        holdout_true_types = []
        holdout_super_classes = []

        for _, row in holdout.iterrows():
            rid = str(row[ROOT_COL])
            ct = row["cell_type"]
            if rid in banc_id_to_idx and pd.notna(ct) and ct in type_to_idx:
                idx = banc_id_to_idx[rid]
                if idx < n_banc_pool and not manual_label_mask[idx]:
                    holdout_indices.append(idx)
                    holdout_true_types.append(type_to_idx[ct])
                    sc = sc_lookup.get(rid, "unknown")
                    holdout_super_classes.append(sc if pd.notna(sc) else "unknown")

    holdout_indices = np.array(holdout_indices, dtype=np.int32)
    holdout_true_types = np.array(holdout_true_types, dtype=np.int32)
    holdout_super_classes = np.array(holdout_super_classes)
    unique_scs = sorted(set(holdout_super_classes))
    print(f"  Holdout: {len(holdout_indices)} neurons, super_classes: {', '.join(unique_scs)}")

    # Mi1 holdout subset
    mi1_idx = type_to_idx.get("Mi1", -1)
    mi1_mask = holdout_true_types == mi1_idx if mi1_idx >= 0 else np.zeros(len(holdout_indices), dtype=bool)
    print(f"  Mi1 holdout: {mi1_mask.sum()} neurons")

    def eval_holdout(hard_assignments):
        if len(holdout_indices) == 0:
            return 0.0, 0
        pred = hard_assignments[holdout_indices]
        valid = pred >= 0
        if valid.sum() == 0:
            return 0.0, 0
        correct = (pred[valid] == holdout_true_types[valid]).sum()
        # Denominator = total holdout, not just assigned. Treats unassigned as
        # wrong so that best-iter selection isn't biased toward low-coverage
        # iterations where only high-confidence neurons are assigned.
        return 100 * correct / len(holdout_indices), int(len(holdout_indices))

    def eval_holdout_by_sc(hard_assignments):
        results = {}
        for sc in unique_scs:
            m = holdout_super_classes == sc
            if m.sum() == 0:
                continue
            pred = hard_assignments[holdout_indices[m]]
            true = holdout_true_types[m]
            valid = pred >= 0
            if valid.sum() == 0:
                results[sc] = (0.0, int(m.sum()))
            else:
                correct = (pred[valid] == true[valid]).sum()
                results[sc] = (100 * correct / m.sum(), int(m.sum()))
        return results

    def eval_mi1(hard_assignments, debug_label=None):
        if mi1_mask.sum() == 0:
            return 0.0, 0
        pred = hard_assignments[holdout_indices[mi1_mask]]
        true = holdout_true_types[mi1_mask]
        valid = pred >= 0
        if debug_label is not None:
            from collections import Counter
            unassigned = int((pred < 0).sum())
            wrong = pred[valid][pred[valid] != true[valid]]
            wrong_types = Counter(types[t] for t in wrong)
            print(f"    [DEBUG eval_mi1 {debug_label}] total Mi1 ho={len(pred)}, unassigned={unassigned}, "
                  f"correct={int((pred[valid]==true[valid]).sum())}, "
                  f"wrong={len(wrong)}, top wrong types: {wrong_types.most_common(8)}")
        if valid.sum() == 0:
            return 0.0, 0
        return 100 * (pred[valid] == true[valid]).sum() / valid.sum(), int(valid.sum())

    def eval_type_counts(hard_assignments):
        """Spearman correlation between BANC and FAFB type count distributions."""
        banc_counts = np.bincount(hard_assignments[hard_assignments >= 0], minlength=n_types)
        mask = (target_type_counts > 0) | (banc_counts > 0)
        if mask.sum() < 2:
            return 0.0
        rho, _ = spearmanr(target_type_counts[mask], banc_counts[mask])
        return rho

    # ---------------------------------------------------------------
    # Build lookups for neuron-level matching
    # ---------------------------------------------------------------
    # FAFB pool index: fafb_id -> position intarget_pool_ids
    fafb_pool_idx = {fid: i for i, fid in enumerate(target_pool_ids)}

    # FAFB neuron -> cell_type lookup
    fafb_type_for_id = {}
    for fid, ct in zip(target_meta["target_id"].astype(str), target_meta["target_cell_type"]):
        if pd.notna(ct) and ct != "":
            fafb_type_for_id[fid] = ct

    # Per-BANC neuron NBLAST matches: banc_id -> [(fafb_id, score), ...]
    print("  Building neuron-level NBLAST lookup...")
    nblast_banc_target_lookup = {}
    for bid, fid, sc in zip(nblast["banc_id"].astype(str),
                            nblast["target_id"].astype(str),
                            nblast["nblast_score"]):
        nblast_banc_target_lookup.setdefault(bid, []).append((fid, float(sc)))
    print(f"  NBLAST lookup: {len(nblast_banc_target_lookup)} BANC neurons, "
          f"{sum(len(v) for v in nblast_banc_target_lookup.values())} total pairs")

    # Pre-build vectorized index arrays for individual profile scoring
    # Keep only top-k NBLAST matches per neuron to control memory/speed
    ind_top_k_nblast = 20  # keep top-20 NBLAST matches per neuron
    n_ind_pairs = 0
    if ind_weight > 0:
        print(f"  Building individual scoring index (top-{ind_top_k_nblast} NBLAST per neuron)...")
        _nb_banc, _nb_fafb, _nb_type = [], [], []
        for bid, pairs in nblast_banc_target_lookup.items():
            if bid not in banc_id_to_idx:
                continue
            bi = banc_id_to_idx[bid]
            if bi >= n_banc_pool or anchor_mask[bi]:
                continue
            # Resolve and filter valid FAFB matches, keep top-k by score
            valid = []
            for fid, sc in pairs:
                if fid not in fafb_pool_idx:
                    continue
                tname = fafb_type_for_id.get(fid)
                if tname is None or tname not in type_to_idx:
                    continue
                valid.append((sc, fafb_pool_idx[fid], type_to_idx[tname]))
            valid.sort(key=lambda x: -x[0])
            for _, fi, ti in valid[:ind_top_k_nblast]:
                _nb_banc.append(bi)
                _nb_fafb.append(fi)
                _nb_type.append(ti)
        ind_banc_arr = np.array(_nb_banc, dtype=np.int32)
        ind_target_arr = np.array(_nb_fafb, dtype=np.int32)
        ind_type_arr = np.array(_nb_type, dtype=np.int32)
        n_ind_pairs = len(ind_banc_arr)
        n_ind_neurons = len(set(_nb_banc))
        print(f"  Individual scoring: {n_ind_pairs} pairs, {n_ind_neurons} BANC neurons "
              f"(filtered from {sum(len(v) for v in nblast_banc_target_lookup.values())} total)")
        del _nb_banc, _nb_fafb, _nb_type

        # Materialise FAFB individual profiles only for the unique FAFB indices
        # in ind_target_arr. 1-hops stay sparse (need full matrix for 2-hop), then
        # row-subsetted on densification. 2-hops are computed for the subset
        # only via row-subsetted sparse @ sparse multiply.
        if n_ind_pairs > 0:
            print("  Computing individual FAFB neuron profiles (subset)...")
            fafb_used_idx = np.unique(ind_target_arr).astype(np.int32)
            n_fafb_used = len(fafb_used_idx)
            target_adj_T_csr = target_adj.T.tocsr()
            target_ind_input = target_adj_T_csr @ target_type_mat
            target_ind_output = target_adj @ target_type_mat
            target_ind_input_2h = target_adj_T_csr[fafb_used_idx] @ target_ind_input
            target_ind_output_2h = target_adj[fafb_used_idx] @ target_ind_output
            target_ind_profiles = np.hstack([
                np.asarray(target_ind_input[fafb_used_idx].todense()),
                np.asarray(target_ind_output[fafb_used_idx].todense()),
                hop2_weight * np.asarray(target_ind_input_2h.todense()),
                hop2_weight * np.asarray(target_ind_output_2h.todense()),
            ]).astype(np.float32)
            target_ind_norms = np.linalg.norm(target_ind_profiles, axis=1, keepdims=True)
            target_ind_norms[target_ind_norms == 0] = 1
            target_ind_profiles_normed = target_ind_profiles / target_ind_norms
            print(f"  FAFB individual profiles: {target_ind_profiles.shape}, "
                  f"{(target_ind_norms.flatten() > 0).sum()}/{n_fafb_used} non-zero "
                  f"(materialised {n_fafb_used}/{n_target_pool} FAFB pool, "
                  f"{target_ind_profiles.nbytes / 1024**3:.2f} GB)")
            # Remap ind_target_arr to local indices for fancy indexing in the loop
            target_full_to_local = np.full(n_target_pool, -1, dtype=np.int32)
            target_full_to_local[fafb_used_idx] = np.arange(n_fafb_used, dtype=np.int32)
            ind_target_arr = target_full_to_local[ind_target_arr]
            assert (ind_target_arr >= 0).all(), "ind_target_arr remap produced -1"
            del target_ind_input, target_ind_output, target_ind_input_2h, target_ind_output_2h
            del target_adj_T_csr

    # ---------------------------------------------------------------
    # Initial evaluation
    # ---------------------------------------------------------------
    hard_init = np.full(n_banc_pool, -1, dtype=np.int32)
    has_init = soft_probs.max(axis=1) > 0
    hard_init[has_init] = soft_probs[has_init].argmax(axis=1).astype(np.int32)

    acc, n_eval = eval_holdout(hard_init)
    nt_acc, n_nt = eval_nt(hard_init)
    mi1_acc, n_mi1 = eval_mi1(hard_init, debug_label="hard_init")
    tc_rho = eval_type_counts(hard_init)
    print(f"  Init holdout: {acc:.1f}% ({n_eval}), Mi1: {mi1_acc:.1f}% ({n_mi1}), "
          f"NT: {nt_acc:.1f}% ({n_nt}), type_rho: {tc_rho:.3f}")

    # ---------------------------------------------------------------
    # Bilateral consistency setup
    # ---------------------------------------------------------------
    if bilateral:
        banc_side_arr = banc_meta["side"].fillna("").astype(str).str.lower().values
        # Map to pool indices (banc_meta may have been filtered/reordered)
        pool_side = np.array([banc_side_arr[i] if i < len(banc_side_arr) else ""
                              for i in range(n_banc_pool)])
        is_left = pool_side == "left"
        is_right = pool_side == "right"
        # FAFB expected L/R ratio per type from capacity file
        cap_file = alignment_path("prep-capacity", query="banc", target=target_name,
                                   region=region, side=side,
                                   vq=banc_version, vt=nblast_version,
                                   ext="csv", dir=data_dir)
        if os.path.exists(cap_file):
            cap_df = pd.read_csv(cap_file)
            target_lr = {}
            for _, r in cap_df.iterrows():
                ct = r["cell_type"]
                if ct in type_to_idx:
                    fl = r.get("fafb_left", 0) or 0
                    fr = r.get("fafb_right", 0) or 0
                    target_lr[type_to_idx[ct]] = (float(fl), float(fr))
            target_left_arr = np.array([target_lr.get(i, (1, 1))[0] for i in range(n_types)])
            fafb_right_arr = np.array([target_lr.get(i, (1, 1))[1] for i in range(n_types)])
            target_expected_left_ratio = target_left_arr / np.maximum(target_left_arr + fafb_right_arr, 1)
        else:
            target_expected_left_ratio = np.full(n_types, 0.5)
        print(f"  Bilateral consistency: {is_left.sum()} left, {is_right.sum()} right neurons")

    # Per-side hard cap mode is active iff bilateral and prep.R provided cap_left/cap_right
    use_per_side_cap = bilateral and (type_capacity_left is not None) and (type_capacity_right is not None)
    if use_per_side_cap:
        print(f"  Per-side hard caps: ENABLED (separate type_slots_left/right)")
    elif bilateral:
        print(f"  Per-side hard caps: disabled (no cap_left/cap_right in capacity file)")

    # ---------------------------------------------------------------
    # Iterative alignment
    # ---------------------------------------------------------------
    best_acc = acc
    best_composite = -np.inf  # for no-holdout: composite = type_rho + nt_acc/100
    best_tc_rho_at_best = 0.0
    best_nt_at_best = 0.0
    best_hard = hard_init.copy()
    best_scores = np.zeros((n_banc_pool, n_types), dtype=np.float32)
    best_iter = 0
    prev_hard = hard_init.copy()

    print(f"\n=== Iterating (max {max_iter}, convergence <{convergence_pct}%) ===")
    print(f"  tau: {tau_start} -> {tau_end}, alpha: {alpha_start} -> {alpha_end}, "
          f"nb_thresh: {nblast_threshold}, capacity: {capacity_scale}x")

    for iteration in range(max_iter):
        t0 = time.time()

        # Schedule: temperature and alpha
        progress = iteration / max(max_iter - 1, 1)
        tau = tau_start * (tau_end / tau_start) ** progress
        alpha = alpha_start + (alpha_end - alpha_start) * progress

        # Build type matrix from soft probs (n_all x n_types)
        type_mat = np.zeros((n_banc_all, n_types), dtype=np.float32)
        if conf_gate:
            # Confidence gating: scale contribution by max soft_prob
            # Neurons below threshold get reduced weight in type profiles
            conf = soft_probs.max(axis=1)
            gate = np.where(conf >= conf_gate_threshold, 1.0,
                            conf / max(conf_gate_threshold, 1e-6)).astype(np.float32)
            type_mat[:n_banc_pool] = soft_probs * gate[:, None]
        else:
            type_mat[:n_banc_pool] = soft_probs

        # 1-hop profiles via sparse @ dense (always computed fully)
        input_1hop = banc_adj_t @ type_mat   # (n_all x n_types)
        output_1hop = banc_adj @ type_mat     # (n_all x n_types)

        if chunked:
            # Chunked scoring: compute profiles + scores per neuron-chunk
            # Avoids materializing full (n_pool x 4*n_types) profiles array
            # 2-hop computed on-the-fly per chunk from full 1-hop intermediates
            conn_scores = np.zeros((n_banc_pool, n_types), dtype=np.float32)
            pchunk = 5000
            # Pre-normalize centroids for cosine (reused across chunks)
            cnorms = np.linalg.norm(target_centroids, axis=1, keepdims=True)
            cnorms[cnorms == 0] = 1
            centroids_normed = target_centroids / cnorms

            for pcs in range(0, n_banc_pool, pchunk):
                pce = min(pcs + pchunk, n_banc_pool)
                # Build chunk profiles
                if hop2_weight > 0:
                    chunk_in2 = banc_adj_t[pcs:pce] @ input_1hop
                    chunk_out2 = banc_adj[pcs:pce] @ output_1hop
                    cp = np.hstack([
                        input_1hop[pcs:pce], output_1hop[pcs:pce],
                        hop2_weight * chunk_in2, hop2_weight * chunk_out2
                    ]).astype(np.float32)
                else:
                    cp = np.hstack([
                        input_1hop[pcs:pce], output_1hop[pcs:pce]
                    ]).astype(np.float32)

                # Cosine scores for this chunk
                pnorms = np.linalg.norm(cp, axis=1, keepdims=True)
                pnorms[pnorms == 0] = 1
                chunk_cos = (cp / pnorms) @ centroids_normed.T

                if metric == "cosine":
                    conn_scores[pcs:pce] = chunk_cos
                elif metric == "jaccard":
                    chunk_size_j = min(200, pce - pcs)
                    for js in range(0, pce - pcs, chunk_size_j):
                        je = min(js + chunk_size_j, pce - pcs)
                        c3 = cp[js:je, None, :]
                        t3 = target_centroids[None, :, :]
                        mins = np.minimum(c3, t3).sum(axis=2)
                        maxs = np.maximum(c3, t3).sum(axis=2)
                        maxs[maxs == 0] = 1
                        conn_scores[pcs + js:pcs + je] = mins / maxs
                elif metric == "ensemble":
                    conn_scores[pcs:pce] = chunk_cos.copy()
                    top_k = ensemble_top_k
                    cos_w, jac_w = ensemble_blend, 1.0 - ensemble_blend
                    top_idx = np.argpartition(-chunk_cos, min(top_k, n_types - 1), axis=1)[:, :top_k]
                    # Sub-chunk the Jaccard re-ranking to control memory
                    # (pchunk, top_k, dim) array can be huge at high dim
                    dim = cp.shape[1]
                    jac_sub = max(1, min(pce - pcs, int(2e9 / (top_k * dim * 4))))  # ~2GB limit
                    for js in range(0, pce - pcs, jac_sub):
                        je = min(js + jac_sub, pce - pcs)
                        gathered = target_centroids[top_idx[js:je]]
                        p = cp[js:je, None, :]
                        mins = np.minimum(p, gathered).sum(axis=2)
                        maxs = np.maximum(p, gathered).sum(axis=2)
                        maxs[maxs == 0] = 1
                        jaccard_k = mins / maxs
                        for i in range(je - js):
                            for k in range(top_k):
                                tidx = top_idx[js + i, k]
                                conn_scores[pcs + js + i, tidx] = cos_w * chunk_cos[js + i, tidx] + jac_w * jaccard_k[i, k]

            # Save banc_profiles for ind_weight (only current chunk needed below,
            # but individual scoring needs the full array — recompute if needed)
            banc_profiles = None  # signal that we used chunked mode
        else:
            # Non-chunked: materialize full profiles (original path, higher memory)
            input_2hop = banc_adj_t @ input_1hop
            output_2hop = banc_adj @ output_1hop
            banc_profiles = np.hstack([
                input_1hop[:n_banc_pool], output_1hop[:n_banc_pool],
                hop2_weight * input_2hop[:n_banc_pool],
                hop2_weight * output_2hop[:n_banc_pool]
            ]).astype(np.float32)

            if metric == "cosine":
                conn_scores = cosine_batch(banc_profiles, target_centroids)
            elif metric == "jaccard":
                conn_scores = weighted_jaccard_batch(banc_profiles, target_centroids)
            elif metric == "ensemble":
                cosine_scores = cosine_batch(banc_profiles, target_centroids)
                conn_scores = cosine_scores.copy()
                top_k = ensemble_top_k
                cos_w, jac_w = ensemble_blend, 1.0 - ensemble_blend
                for cs in range(0, n_banc_pool, 2000):
                    ce = min(cs + 2000, n_banc_pool)
                    chunk_profiles = banc_profiles[cs:ce]
                    chunk_cosine = cosine_scores[cs:ce]
                    top_idx = np.argpartition(-chunk_cosine, min(top_k, n_types - 1), axis=1)[:, :top_k]
                    gathered = target_centroids[top_idx]
                    p = chunk_profiles[:, None, :]
                    mins = np.minimum(p, gathered).sum(axis=2)
                    maxs = np.maximum(p, gathered).sum(axis=2)
                    maxs[maxs == 0] = 1
                    jaccard_k = mins / maxs
                    for i in range(ce - cs):
                        for k in range(top_k):
                            tidx = top_idx[i, k]
                            conn_scores[cs + i, tidx] = cos_w * chunk_cosine[i, tidx] + jac_w * jaccard_k[i, k]

        # Individual FAFB profile scoring (vectorized, chunked)
        # For BANC neurons with NBLAST data, score against individual FAFB
        # neurons' connectivity profiles. Blends centroid and individual scores:
        # (1 - ind_weight) * centroid_score + ind_weight * best_individual_score
        if ind_weight > 0 and n_ind_pairs > 0:
            if banc_profiles is not None:
                # Non-chunked: use full profiles directly
                banc_norms = np.linalg.norm(banc_profiles, axis=1, keepdims=True)
                banc_norms[banc_norms == 0] = 1
                banc_normed = banc_profiles / banc_norms
            else:
                # Chunked mode: compute profiles only for unique neurons in ind index
                unique_banc_ind = np.unique(ind_banc_arr)
                _lp = [input_1hop[unique_banc_ind], output_1hop[unique_banc_ind]]
                if hop2_weight > 0:
                    _lp.append(hop2_weight * (banc_adj_t[unique_banc_ind] @ input_1hop))
                    _lp.append(hop2_weight * (banc_adj[unique_banc_ind] @ output_1hop))
                local_profiles = np.hstack(_lp).astype(np.float32)
                local_norms = np.linalg.norm(local_profiles, axis=1, keepdims=True)
                local_norms[local_norms == 0] = 1
                local_normed = local_profiles / local_norms
                # Remap ind_banc_arr to local indices for the scoring loop
                _b2l = np.full(n_banc_pool, -1, dtype=np.int32)
                _b2l[unique_banc_ind] = np.arange(len(unique_banc_ind), dtype=np.int32)
                _local_banc = _b2l[ind_banc_arr]
                # Build full-pool normalized array for indexed access
                banc_normed = np.zeros((n_banc_pool, local_normed.shape[1]), dtype=np.float32)
                banc_normed[unique_banc_ind] = local_normed
                banc_profiles = np.zeros_like(banc_normed)
                banc_profiles[unique_banc_ind] = local_profiles
                if iteration == 0:
                    print(f"    Ind profiles rebuilt for {len(unique_banc_ind)} neurons (chunked mode)")

            ind_type_scores = np.zeros((n_banc_pool, n_types), dtype=np.float32)
            ind_chunk = 100000
            for cs in range(0, n_ind_pairs, ind_chunk):
                ce = min(cs + ind_chunk, n_ind_pairs)
                bi = ind_banc_arr[cs:ce]
                fi = ind_target_arr[cs:ce]
                ti = ind_type_arr[cs:ce]

                # Resolve individual metric
                _im = metric if ind_metric == "auto" else ind_metric

                # Cosine: row-wise dot product (always computed — fast)
                chunk_cos = np.einsum('ij,ij->i',
                                      banc_normed[bi],
                                      target_ind_profiles_normed[fi])

                if _im in ("jaccard", "ensemble"):
                    chunk_b = banc_profiles[bi]
                    chunk_f = target_ind_profiles[fi]
                    chunk_mins = np.minimum(chunk_b, chunk_f).sum(axis=1)
                    chunk_maxs = np.maximum(chunk_b, chunk_f).sum(axis=1)
                    chunk_maxs[chunk_maxs == 0] = 1
                    chunk_jac = (chunk_mins / chunk_maxs).astype(np.float32)
                    if _im == "jaccard":
                        chunk_scores = chunk_jac
                    else:  # ensemble
                        chunk_scores = (ensemble_blend * chunk_cos +
                                        (1 - ensemble_blend) * chunk_jac)
                else:
                    chunk_scores = chunk_cos

                np.maximum.at(ind_type_scores, (bi, ti), chunk_scores)

            # Weighted blend where individual data exists
            has_ind = ind_type_scores > 0
            conn_scores[has_ind] = ((1 - ind_weight) * conn_scores[has_ind] +
                                    ind_weight * ind_type_scores[has_ind])

            if iteration == 0:
                n_with_ind = (ind_type_scores.max(axis=1) > 0).sum()
                print(f"    Individual profile scoring: {n_with_ind} neurons, "
                      f"{n_ind_pairs} pairs, weight={ind_weight}")

        # Bilateral consistency: penalize over-represented side per type
        if bilateral and iteration > 0:
            left_counts = np.bincount(prev_hard[is_left & (prev_hard >= 0)], minlength=n_types).astype(np.float32)
            right_counts = np.bincount(prev_hard[is_right & (prev_hard >= 0)], minlength=n_types).astype(np.float32)
            total_counts = left_counts + right_counts
            with np.errstate(invalid="ignore", divide="ignore"):
                current_left_ratio = np.where(total_counts > 5, left_counts / total_counts, 0.5)
            # Deviation from FAFB expected ratio
            lr_dev = current_left_ratio - target_expected_left_ratio  # positive = too many left
            # Penalty: reduce scores for over-represented side (scaled by deviation)
            bilateral_strength = 0.1 * progress  # ramp up over iterations
            left_penalty = np.clip(lr_dev * bilateral_strength, 0, None)    # positive when left over-represented
            right_penalty = np.clip(-lr_dev * bilateral_strength, 0, None)  # positive when right over-represented
            conn_scores[is_left] -= left_penalty[None, :]
            conn_scores[is_right] -= right_penalty[None, :]

        # Blend with NBLAST regularizer (only for neurons with NBLAST data)
        final_scores = np.where(
            has_nblast[:, None],
            alpha * conn_scores + (1 - alpha) * nblast_type_scores,
            conn_scores
        )

        # Type frequency prior: boost common types
        if prior_weight > 0:
            final_scores += prior_weight * log_prior[None, :]

        # NBLAST hard constraint: disallow types with low NBLAST evidence
        final_scores[~nblast_allowed] = -1e9

        # Soft NT scoring: multiply by compatibility (match=boost, mismatch=penalize, no data=neutral)
        final_scores *= nt_compat

        # Forbidden (root_id, cell_type) pairs: hard veto from review.
        # Apply after every score component so neither blending nor priors
        # can resurrect a rejected pair. Stage 1 already skips types with
        # score <= -1e8 (line below sets exactly that), and Stage 2 reads
        # from soft_probs which we scrub after softmax.
        if len(forbidden_neuron_arr) > 0:
            final_scores[forbidden_neuron_arr, forbidden_type_arr] = -1e9

        # Soma-presence rule: same hard veto, but rectangular (rows x cols).
        # Mask all (neuron-with-nucleus, no-soma-type) pairs to -1e9 so they
        # cannot win Stage 1 NBLAST matching or Stage 2 capacity allocation.
        if len(soma_rule_row_idx) > 0 and len(soma_rule_col_idx) > 0:
            final_scores[np.ix_(soma_rule_row_idx, soma_rule_col_idx)] = -1e9

        # Which neurons to update (not anchors, must have some signal)
        # Check which neurons have non-zero connectivity profiles
        if banc_profiles is not None:
            has_profile = banc_profiles.sum(axis=1) > 0
        else:
            # Chunked mode: check from 1-hop intermediates
            has_profile = ((np.asarray(input_1hop[:n_banc_pool]).sum(axis=1) +
                            np.asarray(output_1hop[:n_banc_pool]).sum(axis=1)) > 0)
            if hasattr(has_profile, 'A1'):  # sparse matrix .sum returns matrix
                has_profile = np.asarray(has_profile).flatten()
        update_mask = ~anchor_mask & (has_profile | has_nblast)

        # Softmax with temperature for updated neurons. Optional classic
        # momentum on the softmax output: new_soft = β·prev_soft + (1-β)·softmax.
        # Distinct from the capacity-feedback block below, which blends *hard*
        # one-hot assignments into soft_probs. Together they damp two different
        # noise sources: momentum smooths softmax fluctuation across iters;
        # feedback enforces capacity-aware assignments.
        if update_mask.any():
            scores_upd = final_scores[update_mask]
            scaled = scores_upd / tau
            scaled -= scaled.max(axis=1, keepdims=True)  # numerical stability
            exp_s = np.exp(scaled)
            row_sums = exp_s.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1
            new_soft = (exp_s / row_sums).astype(np.float32)
            if momentum > 0:
                soft_probs[update_mask] = momentum * soft_probs[update_mask] + (1 - momentum) * new_soft
            else:
                soft_probs[update_mask] = new_soft

        # Scrub forbidden pairs from soft_probs so the feedback blend cannot
        # reintroduce a vetoed type via the constrained-soft mixture, and so
        # Stage 2's argsort+break loop cleanly skips them.
        if len(forbidden_neuron_arr) > 0:
            soft_probs[forbidden_neuron_arr, forbidden_type_arr] = 0

        # Same scrub for the soma rule (rectangular mask).
        if len(soma_rule_row_idx) > 0 and len(soma_rule_col_idx) > 0:
            soft_probs[np.ix_(soma_rule_row_idx, soma_rule_col_idx)] = 0

        # =============================================================
        # Neuron-level greedy assignment
        # Stage 1: Neurons with NBLAST — match to individual FAFB neurons
        #   composite = type_conn_score(banc_i, type(fafb_j)) + nblast(banc_i, fafb_j)
        #   Each FAFB neuron can only be claimed once.
        # Stage 2: Neurons without NBLAST — type-level capacity assignment
        # =============================================================
        hard_assignments = np.full(n_banc_pool, -1, dtype=np.int32)
        target_claimed = np.zeros(n_target_pool, dtype=bool)  # track claimed FAFB neurons
        type_slots = type_capacity.copy()
        if use_per_side_cap:
            type_slots_left = type_capacity_left.copy()
            type_slots_right = type_capacity_right.copy()

        # Pre-fill anchors (anchors always get their type, even if it overshoots cap)
        hard_assignments[anchor_indices] = anchor_type_indices
        for ai, tidx in zip(anchor_indices, anchor_type_indices):
            type_slots[tidx] -= 1
            if use_per_side_cap:
                if is_left[ai]:
                    type_slots_left[tidx] -= 1
                elif is_right[ai]:
                    type_slots_right[tidx] -= 1
        type_slots = np.maximum(type_slots, 0)
        if use_per_side_cap:
            type_slots_left = np.maximum(type_slots_left, 0)
            type_slots_right = np.maximum(type_slots_right, 0)
        # Claim FAFB neurons for anchors (best NBLAST match of their type)
        for ai, ati in zip(anchor_indices, anchor_type_indices):
            rid = banc_pool_ids[ai]
            tname = types[ati]
            if rid in nblast_banc_target_lookup:
                candidates = nblast_banc_target_lookup[rid]
                best_fid = None
                best_sc = -1
                for fid, sc in candidates:
                    if fafb_type_for_id.get(fid) == tname and not target_claimed[fafb_pool_idx.get(fid, -1)] if fid in fafb_pool_idx else False:
                        if sc > best_sc:
                            best_sc = sc
                            best_fid = fid
                if best_fid is not None and best_fid in fafb_pool_idx:
                    target_claimed[fafb_pool_idx[best_fid]] = True

        # Stage 1: neurons WITH NBLAST data — neuron-level matching
        # Build scored pairs: (composite_score, banc_pool_idx, fafb_pool_idx, type_idx)
        # Stage 1 (BANC↔FAFB neuron-level bipartite matching) is skippable.
        # When --skip-stage1, all has_nblast neurons fall through to Stage 2's
        # per-neuron preference-ordered greedy (the v1 mechanism).
        stage1_neurons = (np.array([], dtype=np.int64) if skip_stage1
                          else np.where(has_nblast & ~anchor_mask &
                                        (soft_probs.max(axis=1) > 0))[0])
        if len(stage1_neurons) > 0 and iteration == 0:
            print(f"    Stage 1: {len(stage1_neurons)} neurons with NBLAST for neuron-level matching")
        elif skip_stage1 and iteration == 0:
            print(f"    Stage 1: SKIPPED (--skip-stage1) — all has_nblast neurons go to Stage 2")

        pairs = []  # (score, banc_idx, fafb_pool_idx, type_idx)
        for ni in stage1_neurons:
            rid = banc_pool_ids[ni]
            if rid not in nblast_banc_target_lookup:
                continue
            type_score_vec = final_scores[ni]  # (n_types,) type-level scores
            for fid, nb_score in nblast_banc_target_lookup[rid]:
                if fid not in fafb_pool_idx:
                    continue
                fi = fafb_pool_idx[fid]
                tname = fafb_type_for_id.get(fid)
                if tname is None or tname not in type_to_idx:
                    continue
                tidx = type_to_idx[tname]
                if type_score_vec[tidx] <= -1e8:  # filtered by constraint
                    continue
                composite = alpha * type_score_vec[tidx] + (1 - alpha) * nb_score
                pairs.append((composite, ni, fi, tidx))

        # Sort by composite score descending, assign greedily
        pairs.sort(key=lambda x: -x[0])
        for score_val, ni, fi, tidx in pairs:
            if hard_assignments[ni] >= 0:
                continue  # already assigned
            if target_claimed[fi]:
                continue  # FAFB neuron already taken
            if type_slots[tidx] <= 0:
                continue  # type full (pooled)
            if use_per_side_cap:
                if is_left[ni] and type_slots_left[tidx] <= 0:
                    continue  # left side full for this type
                if is_right[ni] and type_slots_right[tidx] <= 0:
                    continue  # right side full for this type
            hard_assignments[ni] = tidx
            target_claimed[fi] = True
            type_slots[tidx] -= 1
            if use_per_side_cap:
                if is_left[ni]:
                    type_slots_left[tidx] -= 1
                elif is_right[ni]:
                    type_slots_right[tidx] -= 1

        n_stage1 = (hard_assignments[stage1_neurons] >= 0).sum() if len(stage1_neurons) > 0 else 0

        # Stage 2: remaining neurons — type-level capacity assignment
        unassigned = np.where((hard_assignments < 0) & ~anchor_mask &
                              (soft_probs.max(axis=1) > 0))[0]
        if len(unassigned) > 0:
            best_scores_arr = soft_probs[unassigned].max(axis=1)
            sort_order = np.argsort(-best_scores_arr)
            for idx_pos in sort_order:
                ni = unassigned[idx_pos]
                type_order = np.argsort(-soft_probs[ni])
                for tidx in type_order:
                    if soft_probs[ni, tidx] <= 0:
                        break
                    if type_slots[tidx] <= 0:
                        continue
                    if use_per_side_cap:
                        if is_left[ni] and type_slots_left[tidx] <= 0:
                            continue
                        if is_right[ni] and type_slots_right[tidx] <= 0:
                            continue
                    hard_assignments[ni] = tidx
                    type_slots[tidx] -= 1
                    if use_per_side_cap:
                        if is_left[ni]:
                            type_slots_left[tidx] -= 1
                        elif is_right[ni]:
                            type_slots_right[tidx] -= 1
                    break

        # Feedback: blend hard constrained + original soft for next iteration
        # Higher blend = stickier assignments = less oscillation
        soft_probs_constrained = np.zeros_like(soft_probs)
        assigned_mask_c = hard_assignments >= 0
        soft_probs_constrained[assigned_mask_c,
                               hard_assignments[assigned_mask_c]] = 1.0
        # Anneal blend: start soft (0.3), end sticky (0.8)
        feedback_blend = 0.3 + 0.5 * progress
        soft_probs = np.where(
            assigned_mask_c[:, None],
            feedback_blend * soft_probs_constrained + (1 - feedback_blend) * soft_probs,
            soft_probs
        )

        # Convergence
        changed = (hard_assignments != prev_hard) & update_mask
        n_changed = int(changed.sum())
        n_reassignable = int(update_mask.sum())
        pct_changed = 100 * n_changed / max(n_reassignable, 1)
        n_assigned = int((hard_assignments >= 0).sum())
        elapsed = time.time() - t0

        # Evaluate
        acc, n_eval = eval_holdout(hard_assignments)
        nt_acc, n_nt = eval_nt(hard_assignments)
        # Only debug-print Mi1 misclassifications at iter 1 (the collapse moment)
        mi1_acc, n_mi1 = eval_mi1(hard_assignments,
                                   debug_label=f"iter{iteration+1}" if iteration < 2 else None)
        tc_rho = eval_type_counts(hard_assignments)

        print(f"  Iter {iteration+1:2d}: {n_changed:6d} changed ({pct_changed:5.1f}%), "
              f"{n_assigned}/{n_banc_pool} assigned, "
              f"holdout: {acc:.1f}%, Mi1: {mi1_acc:.1f}%, NT: {nt_acc:.1f}%, "
              f"type_rho: {tc_rho:.3f}, s1={n_stage1}, tau={tau:.3f}, a={alpha:.2f}, {elapsed:.1f}s")

        # Detailed reporting every 5 iterations
        if (iteration + 1) % 5 == 0 or iteration == 0:
            sc_results = eval_holdout_by_sc(hard_assignments)
            for sc, (sc_acc, sc_n) in sorted(sc_results.items()):
                print(f"         {sc}: {sc_acc:.1f}% ({sc_n})")
            # Top type count divergences
            banc_counts = np.bincount(hard_assignments[hard_assignments >= 0], minlength=n_types)
            diffs = banc_counts - target_type_counts
            top_div = np.argsort(-np.abs(diffs))[:10]
            print("         Top type count divergences (BANC - FAFB):")
            for tidx in top_div:
                if np.abs(diffs[tidx]) > 0:
                    print(f"           {types[tidx]:20s}: BANC={banc_counts[tidx]:5d}, "
                          f"FAFB={target_type_counts[tidx]:5d}, diff={diffs[tidx]:+5d}")

        # Track best iteration (early stopping)
        # When no holdout (acc always 0), use composite of type_rho + NT consistency
        # as criterion. Both metrics in [0,1]; equal-weight sum keeps the best
        # iteration good on global type distribution AND per-neuron NT match.
        if len(holdout_indices) == 0:
            composite = tc_rho + nt_acc / 100.0
            if composite > best_composite:
                best_composite = composite
                best_acc = tc_rho  # keep best_acc semantics: dominant criterion
                best_tc_rho_at_best = tc_rho
                best_nt_at_best = nt_acc
                best_hard = hard_assignments.copy()
                best_scores = final_scores.copy()
                best_iter = iteration + 1
        elif acc > best_acc:
            best_acc = acc
            best_hard = hard_assignments.copy()
            best_scores = final_scores.copy()
            best_iter = iteration + 1

        prev_hard = hard_assignments.copy()

        # Periodic checkpoint of best-so-far assignments. Lets a walltime kill
        # or manual scancel leave a usable result — same schema as the final
        # CSV so the push script can read it directly.
        if checkpoint_path and checkpoint_every > 0 \
                and (iteration + 1) % checkpoint_every == 0:
            ck_df = _build_results_df(
                best_hard, best_scores, banc_pool_ids, types,
                n_banc_pool, seeds, anchor_mask, nblast, target_meta,
                verbose=False)
            ck_df["checkpoint_iter"] = iteration + 1
            ck_df["best_iter"] = best_iter
            ck_df["best_acc"] = best_acc
            tmp = checkpoint_path + ".tmp"
            ck_df.to_csv(tmp, index=False)
            os.replace(tmp, checkpoint_path)
            print(f"  Checkpoint @ iter {iteration + 1}: {checkpoint_path} "
                  f"(best={best_acc:.1f}% @ iter {best_iter})")

        if pct_changed < convergence_pct and iteration >= 2:
            print(f"  Converged at iteration {iteration + 1}")
            break

    # Use best-iteration assignments
    if len(holdout_indices) == 0:
        print(f"\n  Best (no-holdout, composite=type_rho+NT/100): {best_composite:.3f} at iteration {best_iter}")
        print(f"    type_rho={best_tc_rho_at_best:.3f}, NT={best_nt_at_best:.1f}%")
    else:
        print(f"\n  Best holdout: {best_acc:.1f}% at iteration {best_iter}")
    if best_iter != iteration + 1:
        print(f"  Using best-iteration assignments (iter {best_iter}) instead of last (iter {iteration + 1})")
    hard_assignments = best_hard
    scores = best_scores

    # Post-hoc verification: the per-iteration mask should already prevent
    # forbidden assignments, but double-check the restored best_hard and
    # clear any violations rather than silently shipping them.
    if forbidden_pairs_set:
        n_violations = 0
        for i in range(n_banc_pool):
            t = int(hard_assignments[i])
            if t >= 0 and (int(i), t) in forbidden_pairs_set:
                hard_assignments[i] = -1
                n_violations += 1
        if n_violations > 0:
            print(f"  WARNING: {n_violations} forbidden assignments cleared post-hoc")
        else:
            print(f"  Forbidden mask: 0 violations in best iteration")

    # Same post-hoc check for the soma rule.
    if len(soma_rule_row_idx) > 0 and len(soma_rule_col_idx) > 0:
        n_soma_violations = 0
        for ni in soma_rule_row_idx:
            t = int(hard_assignments[int(ni)])
            if t >= 0 and t in soma_rule_col_set:
                hard_assignments[int(ni)] = -1
                n_soma_violations += 1
        if n_soma_violations > 0:
            print(f"  WARNING: {n_soma_violations} soma-rule violations cleared post-hoc")
        else:
            print(f"  Soma rule: 0 violations in best iteration")

    # ---------------------------------------------------------------
    # Build output
    # ---------------------------------------------------------------
    print("\n=== Building output ===")
    # For each neuron, get the score of its ASSIGNED type (not necessarily the argmax)
    # This allows negative confidence for anchored neurons whose GT type
    # scores worse than the connectivity-preferred type.
    results_df = _build_results_df(
        hard_assignments, scores, banc_pool_ids, types,
        n_banc_pool, seeds, anchor_mask, nblast, target_meta,
        verbose=True)

    # Summary
    n_typed = (results_df["assigned_cell_type"] != "").sum()
    n_unique = results_df.loc[results_df["assigned_cell_type"] != "", "assigned_cell_type"].nunique()
    n_with_match = (results_df["best_target_match"] != "").sum()
    print(f"\n=== Results ===")
    print(f"  Typed: {n_typed}/{n_banc_pool} ({100*n_typed/n_banc_pool:.1f}%)")
    print(f"  Unique types assigned: {n_unique}/{n_types}")
    print(f"  With FAFB match: {n_with_match}")
    typed_conf = results_df.loc[results_df["assigned_cell_type"] != "", "confidence"]
    print(f"  Mean confidence: {typed_conf.mean():.4f}")

    metrics = {
        "best_holdout": best_acc,
        "best_iter": best_iter,
        "mi1_at_best": eval_mi1(best_hard)[0],
        "nt_at_best": eval_nt(best_hard)[0],
        "type_rho_at_best": eval_type_counts(best_hard),
        "n_typed": int(n_typed),
        "n_unique_types": int(n_unique),
    }

    return results_df, metrics


def main():
    parser = argparse.ArgumentParser(
        description="Optic lobe cross-dataset connectivity alignment")
    parser.add_argument("--side", default="right", choices=["right", "left", "both"])
    parser.add_argument("--data-dir", default="data/optic_lobe")
    parser.add_argument("--file-prefix", default="optic",
                        help="File naming prefix (default: optic). Use 'brain' for whole-brain.")
    # v2 grammar context: must match the alignment_path() params the prep
    # script used so reads find the prep outputs. Defaults align with paper-run.
    parser.add_argument("--target-name", default="fafb",
                        help="Target dataset name (default: fafb).")
    parser.add_argument("--banc-version", default="888",
                        help="BANC version pin used by prep (default: 888).")
    parser.add_argument("--nblast-version", default="783",
                        help="NBLAST version pin used by prep (default: 783).")
    parser.add_argument("--syn-source", default="v3",
                        help="BANC synapse-table version used by prep (default: v3).")
    parser.add_argument("--manual-labels", default=None, nargs="?", const="all",
                        help="Use manual labels as anchors. No arg or 'all'=all types, "
                             "or comma-separated list e.g. 'Mi1,Tm1'.")
    parser.add_argument("--random-holdout", type=float, default=None,
                        help="Override holdout: randomly hold out this fraction of ALL typed neurons "
                             "(e.g. 0.2 = 20%% holdout). Seed=42 for reproducibility. "
                             "Use for fair comparison with NTAC.")
    parser.add_argument("--stratified-holdout", type=float, default=None,
                        help="Like --random-holdout but stratified by cell_type: same fraction "
                             "held out per type, floor at 1 holdout neuron per type. "
                             "Ensures all types represented in both seed and holdout sets.")
    parser.add_argument("--output-suffix", default=None,
                        help="Suffix for output filename, e.g. 'no_manual_labels' → "
                             "banc_optic_right_alignment_no_manual_labels.csv")
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--convergence", type=float, default=1.0,
                        help="Convergence threshold %% (default: 1.0)")
    parser.add_argument("--tau-start", type=float, default=2.0,
                        help="Initial softmax temperature (default: 2.0)")
    parser.add_argument("--tau-end", type=float, default=0.5,
                        help="Final softmax temperature (default: 0.5)")
    parser.add_argument("--alpha-start", type=float, default=0.3,
                        help="Initial connectivity weight in NBLAST blend (default: 0.3)")
    parser.add_argument("--alpha-end", type=float, default=0.8,
                        help="Final connectivity weight (default: 0.8)")
    parser.add_argument("--nblast-threshold", type=float, default=0.25,
                        help="Min NBLAST score to allow type (default: 0.25)")
    parser.add_argument("--metric", default="cosine", choices=["cosine", "jaccard", "ensemble"],
                        help="Similarity metric: cosine, jaccard, or ensemble (default: cosine)")
    parser.add_argument("--prior-weight", type=float, default=0.0,
                        help="Type frequency log-prior weight (default: 0, try 0.3-1.0)")
    parser.add_argument("--hop2-weight", type=float, default=0.5,
                        help="Weight for 2-hop profiles (default: 0.5, 0=disable)")
    parser.add_argument("--ensemble-top-k", type=int, default=30,
                        help="Top-K cosine candidates for Jaccard re-ranking (default: 30)")
    parser.add_argument("--ensemble-blend", type=float, default=0.5,
                        help="Cosine weight in ensemble (default: 0.5, 1=pure cosine)")
    parser.add_argument("--nt-weight", type=float, default=0.5,
                        help="NT constraint softness: 0=ignore NT, 0.5=soft, 1.0=near-hard (default: 0.5)")
    parser.add_argument("--ind-weight", type=float, default=0.0,
                        help="Weight for individual FAFB profile scoring (0=disabled, default: 0). "
                             "Blends centroid and best individual scores per type.")
    parser.add_argument("--ind-metric", default="auto",
                        choices=["auto", "cosine", "jaccard", "ensemble"],
                        help="Metric for individual profile scoring (default: auto=match main metric)")
    parser.add_argument("--capacity-scale", type=float, default=1.5,
                        help="Type capacity as multiple of FAFB count (default: 1.5)")
    parser.add_argument("--no-chunked", action="store_true",
                        help="Disable chunked profile computation (uses more memory)")
    parser.add_argument("--no-conf-gate", action="store_true",
                        help="Disable confidence gating on type profiles")
    parser.add_argument("--conf-gate-threshold", type=float, default=0.1,
                        help="Confidence gating threshold (default: 0.1)")
    parser.add_argument("--bilateral", action="store_true",
                        help="Enable bilateral L/R consistency constraint")
    parser.add_argument("--forbidden-matches", default=None,
                        help="Optional CSV with columns (root_id, cell_type) listing "
                             "reviewed false-positive pairs to veto. See "
                             "banc-alignment-false-positives.R.")
    parser.add_argument("--momentum", type=float, default=0.0,
                        help="Classic softmax momentum: new_soft = β·prev + (1-β)·softmax. "
                             "0=off (current behaviour), 0.3-0.7 typical. Distinct from the "
                             "capacity-feedback block which mixes hard one-hot back in.")
    parser.add_argument("--skip-stage1", action="store_true",
                        help="Skip the BANC↔FAFB bipartite matching stage. "
                             "All has_nblast neurons fall through to Stage 2 per-neuron greedy "
                             "(matches v1 behaviour; useful for diagnosing iteration regressions).")
    parser.add_argument("--no-soma-rule", action="store_true",
                        help="Disable the soma-presence rule (default ON). "
                             "When on, BANC neurons with a nucleus_id cannot be "
                             "assigned to FAFB types whose super_class contains "
                             "'sensory' or to L1-L5/R7/R8. Anchors exempt.")
    parser.add_argument("--log-file", default=None,
                        help="Write output to log file in addition to stdout")
    parser.add_argument("--checkpoint-every", type=int, default=0,
                        help="If >0, write best-so-far assignments every N "
                             "iterations to <out_file>.checkpoint.csv. Lets a "
                             "walltime kill leave usable output.")
    args = parser.parse_args()

    # Set up log file if requested
    if args.log_file:
        import sys
        class Tee:
            def __init__(self, *streams):
                self.streams = streams
            def write(self, data):
                for s in self.streams:
                    s.write(data)
                    s.flush()
            def flush(self):
                for s in self.streams:
                    s.flush()
        log_fh = open(args.log_file, "w")
        sys.stdout = Tee(sys.__stdout__, log_fh)

    t_start = time.time()
    banc_meta, target_meta, banc_el, target_el, nblast, seeds = load_data(
        args.side, args.data_dir, args.file_prefix,
        target_name=args.target_name,
        banc_version=args.banc_version,
        nblast_version=args.nblast_version,
        syn_source=args.syn_source)

    # Parse manual_labels arg
    ml = args.manual_labels
    if ml is None:
        manual_labels = None
    elif ml.lower() == "all":
        manual_labels = True
    else:
        manual_labels = [t.strip() for t in ml.split(",")]

    align_stage = f"align-{args.output_suffix}" if args.output_suffix else "align"
    out_file = alignment_path(align_stage, query="banc", target=args.target_name,
                              region=_region_for(args.file_prefix), side=args.side,
                              vq=args.banc_version, vt=args.nblast_version,
                              ext="csv", dir=args.data_dir)
    checkpoint_path = out_file.replace(".csv", "_checkpoint.csv") \
        if args.checkpoint_every > 0 else None

    results, metrics = run_alignment(
        banc_meta, target_meta, banc_el, target_el, nblast, seeds,
        data_dir=args.data_dir,
        manual_labels=manual_labels,
        random_holdout=args.random_holdout,
        stratified_holdout=args.stratified_holdout,
        max_iter=args.max_iter, convergence_pct=args.convergence,
        tau_start=args.tau_start, tau_end=args.tau_end,
        alpha_start=args.alpha_start, alpha_end=args.alpha_end,
        nblast_threshold=args.nblast_threshold,
        capacity_scale=args.capacity_scale,
        metric=args.metric,
        prior_weight=args.prior_weight,
        hop2_weight=args.hop2_weight,
        ensemble_top_k=args.ensemble_top_k,
        ensemble_blend=args.ensemble_blend,
        nt_weight=args.nt_weight,
        ind_weight=args.ind_weight,
        ind_metric=args.ind_metric,
        chunked=not args.no_chunked,
        conf_gate=not args.no_conf_gate,
        conf_gate_threshold=args.conf_gate_threshold,
        bilateral=args.bilateral,
        forbidden_matches=args.forbidden_matches,
        apply_soma_rule=not args.no_soma_rule,
        checkpoint_path=checkpoint_path,
        checkpoint_every=args.checkpoint_every,
        momentum=args.momentum,
        skip_stage1=args.skip_stage1,
        region=_region_for(args.file_prefix),
        side=args.side,
        target_name=args.target_name,
        banc_version=args.banc_version,
        nblast_version=args.nblast_version,
        syn_source=args.syn_source)

    results.to_csv(out_file, index=False)
    print(f"\nSaved: {out_file}")
    # Final results superseded the checkpoint; remove it to avoid confusion.
    if checkpoint_path and os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)
    print(f"Total time: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
