"""Shared stratified holdout split used by banc-alignment-run.py and
banc-alignment-ntac.py.

The algorithm must stay in sync across callers: at a given fraction, both
methods must hold out the exact same neurons so the eval script can compare
them on an identical validation set. This module is the single source of
truth — do not duplicate the algorithm elsewhere.
"""

import numpy as np
import pandas as pd


def _detect_root_col(seeds_df):
    """Find the `root_<digits>` column (e.g. root_850, root_888) in the seeds
    DataFrame. Lets this module work against prep outputs from either BANC
    version. Falls back to 'root_id'."""
    for c in seeds_df.columns:
        if c.startswith("root_") and c[5:].isdigit():
            return c
    return "root_id"


def build_target_type_vocab(target_meta):
    """Sorted list of FAFB cell types. Matches align.build_type_index."""
    typed = target_meta.dropna(subset=["cell_type"])
    typed = typed[typed["cell_type"] != ""]
    return sorted(typed["cell_type"].unique())


def stratified_holdout_split(seeds_df, pool_id_to_idx, target_type_vocab,
                             fraction, seed=42):
    """Split typed BANC neurons into (seed_entries, holdout_entries).

    Each entry is a tuple ``(pool_idx, cell_type, root_<NNN>)`` where the
    third field comes from whichever ``root_<digits>`` column is present in
    ``seeds_df`` (e.g. root_850 for v850 prep, root_888 for v888 prep).

    Parameters
    ----------
    seeds_df : pandas.DataFrame
        Must have some ``root_<digits>`` column plus ``cell_type`` (typically
        the banc_<prefix>_<side>_seeds.csv file).
    pool_id_to_idx : dict[str, int]
        Maps ``str(root_<NNN>)`` to pool index. Neurons absent from this map
        are ignored.
    target_type_vocab : iterable[str]
        Valid cell types. BANC neurons whose type is outside this set are
        ignored (matches align.py's ``type_to_idx`` filter).
    fraction : float in (0, 1]
        Fraction held out per type. At 1.0, all typed neurons go to holdout.
    seed : int
        ``np.random.RandomState`` seed (default 42, matching align.py).

    Notes
    -----
    Algorithm details that must not change without updating align.py's
    expectations:

    1. ``all_typed`` is built by iterating ``seeds_df`` in file order.
    2. Types are iterated in ``sorted()`` order; the RNG state advances
       across types.
    3. ``n_ho = max(1, round(fraction * n))`` — at least one neuron per type
       is held out unless the type has zero members.
    """
    root_col = _detect_root_col(seeds_df)
    type_set = set(target_type_vocab)
    all_typed = []
    for _, row in seeds_df.iterrows():
        rid = str(row[root_col])
        ct = row["cell_type"]
        if rid in pool_id_to_idx and pd.notna(ct) and ct in type_set:
            all_typed.append((pool_id_to_idx[rid], ct, rid))

    rng = np.random.RandomState(seed)

    typed_by_type = {}
    for entry in all_typed:
        typed_by_type.setdefault(entry[1], []).append(entry)

    seed_list, holdout_list = [], []
    for ct, neurons in sorted(typed_by_type.items()):
        n = len(neurons)
        n_ho = max(1, round(fraction * n))
        n_ho = min(n_ho, n)
        perm = rng.permutation(n)
        for i in perm[:n_ho]:
            holdout_list.append(neurons[i])
        for i in perm[n_ho:]:
            seed_list.append(neurons[i])
    return seed_list, holdout_list
