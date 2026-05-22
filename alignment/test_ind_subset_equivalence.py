"""Verify the deferred-subset materialisation of FAFB individual profiles
produces the same dot-product / Jaccard results as the original full-pool
materialisation, on synthetic data sized like the test path.

Runs locally (no SLURM, no real data) — just exercises the matrix algebra.
"""
import numpy as np
import scipy.sparse as sp


def full_materialisation(fafb_adj, fafb_type_mat, hop2_weight, n_fafb_pool):
    """Original code path."""
    fafb_ind_input = fafb_adj.T.tocsr() @ fafb_type_mat
    fafb_ind_output = fafb_adj @ fafb_type_mat
    fafb_ind_input_2h = fafb_adj.T.tocsr() @ fafb_ind_input
    fafb_ind_output_2h = fafb_adj @ fafb_ind_output
    profiles = np.hstack([
        np.asarray(fafb_ind_input[:n_fafb_pool].todense()),
        np.asarray(fafb_ind_output[:n_fafb_pool].todense()),
        hop2_weight * np.asarray(fafb_ind_input_2h[:n_fafb_pool].todense()),
        hop2_weight * np.asarray(fafb_ind_output_2h[:n_fafb_pool].todense()),
    ]).astype(np.float32)
    norms = np.linalg.norm(profiles, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return profiles, profiles / norms


def subset_materialisation(fafb_adj, fafb_type_mat, hop2_weight, n_fafb_pool, fafb_used_idx):
    """New code path (matches the production change)."""
    fafb_adj_T_csr = fafb_adj.T.tocsr()
    fafb_ind_input = fafb_adj_T_csr @ fafb_type_mat
    fafb_ind_output = fafb_adj @ fafb_type_mat
    fafb_ind_input_2h = fafb_adj_T_csr[fafb_used_idx] @ fafb_ind_input
    fafb_ind_output_2h = fafb_adj[fafb_used_idx] @ fafb_ind_output
    profiles = np.hstack([
        np.asarray(fafb_ind_input[fafb_used_idx].todense()),
        np.asarray(fafb_ind_output[fafb_used_idx].todense()),
        hop2_weight * np.asarray(fafb_ind_input_2h.todense()),
        hop2_weight * np.asarray(fafb_ind_output_2h.todense()),
    ]).astype(np.float32)
    norms = np.linalg.norm(profiles, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return profiles, profiles / norms


def main():
    rng = np.random.default_rng(42)
    n_fafb_all = 1500
    n_fafb_pool = 1200
    n_types = 80

    # Sparse adjacency (~0.5% density, weights 1-5)
    nnz = 30000
    rows = rng.integers(0, n_fafb_all, size=nnz)
    cols = rng.integers(0, n_fafb_all, size=nnz)
    data = rng.integers(1, 6, size=nnz).astype(np.float32)
    fafb_adj = sp.coo_matrix((data, (rows, cols)),
                             shape=(n_fafb_all, n_fafb_all)).tocsr()

    # Type matrix (each pool neuron assigned one type; non-pool neurons unassigned)
    type_assign = rng.integers(0, n_types, size=n_fafb_pool)
    rows_t = np.arange(n_fafb_pool)
    cols_t = type_assign
    data_t = np.ones(n_fafb_pool, dtype=np.float32)
    fafb_type_mat = sp.coo_matrix((data_t, (rows_t, cols_t)),
                                  shape=(n_fafb_all, n_types)).tocsr()

    hop2_weight = 1.0

    # Pretend ind_fafb_arr references ~30% of the pool
    fafb_used_idx = np.sort(
        rng.choice(n_fafb_pool, size=int(0.3 * n_fafb_pool), replace=False)
    ).astype(np.int32)
    n_fafb_used = len(fafb_used_idx)

    full_profiles, full_normed = full_materialisation(
        fafb_adj, fafb_type_mat, hop2_weight, n_fafb_pool)
    subset_profiles, subset_normed = subset_materialisation(
        fafb_adj, fafb_type_mat, hop2_weight, n_fafb_pool, fafb_used_idx)

    # The subset must equal the rows of full corresponding to fafb_used_idx
    expected = full_profiles[fafb_used_idx]
    expected_n = full_normed[fafb_used_idx]

    err_p = np.abs(subset_profiles - expected).max()
    err_n = np.abs(subset_normed - expected_n).max()
    print(f"  full shape: {full_profiles.shape} ({full_profiles.nbytes / 1024**2:.1f} MB)")
    print(f"  subset shape: {subset_profiles.shape} ({subset_profiles.nbytes / 1024**2:.1f} MB)")
    print(f"  max abs diff (raw profiles): {err_p:.2e}")
    print(f"  max abs diff (normed):       {err_n:.2e}")
    assert err_p < 1e-4, f"raw profile mismatch: {err_p}"
    assert err_n < 1e-5, f"normed profile mismatch: {err_n}"

    # Also check that simulating a chunk's pairwise score is identical
    # under the remap (full[fi] vs subset[fafb_full_to_local[fi]])
    fafb_full_to_local = np.full(n_fafb_pool, -1, dtype=np.int32)
    fafb_full_to_local[fafb_used_idx] = np.arange(n_fafb_used, dtype=np.int32)

    # Random query "BANC" profile
    banc_profile = rng.standard_normal(full_profiles.shape[1]).astype(np.float32)
    banc_normed = banc_profile / np.linalg.norm(banc_profile)

    # Fake 1000 (banc, fafb) pairs all targeting fafb_used_idx subset
    fi_full = rng.choice(fafb_used_idx, size=1000)
    fi_local = fafb_full_to_local[fi_full]

    score_full = full_normed[fi_full] @ banc_normed
    score_subset = subset_normed[fi_local] @ banc_normed
    err_score = np.abs(score_full - score_subset).max()
    print(f"  max abs diff (1000 random pair scores): {err_score:.2e}")
    assert err_score < 1e-5, f"score mismatch: {err_score}"

    print("PASS: subset materialisation matches full materialisation.")


if __name__ == "__main__":
    main()
