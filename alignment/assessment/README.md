# alignment/assessment

Test suite for the alignment code (Python aligner + R data-source helpers).
Kept scoped to invariants that don't depend on sweep outcomes, so these are
safe to run and extend while the whole-brain sweep is still settling.

## Running

**Python (pytest):**

```bash
cd /home/ab714/bancpipeline
source o2/o2_env.sh         # picks up module-loaded Python
pytest alignment/assessment
```

Or a single file:

```bash
pytest alignment/assessment/test_metrics.py -v
```

**R (testthat):**

```bash
cd /home/ab714/bancpipeline
Rscript -e 'testthat::test_dir("alignment/assessment")'
```

## What's covered

| File | Target | Type |
|------|--------|------|
| `test_metrics.py` | `cosine_batch`, `weighted_jaccard_batch` | numeric correctness, zero-safety, chunking consistency |
| `test_profile_builder.py` | `build_adjacency`, `build_type_matrix` | pool/boundary index layout, edge weight placement, n_all clamping |
| `test_splits.py` | `stratified_holdout_split` | per-type fractions, min-one-per-type, disjointness, determinism, vocab filter, fraction=1.0 |
| `test-data-sources.R` | `parse_alignment_data_args` | all CLI forms, env var fallback, validation errors |

## What's NOT covered yet

Deliberately deferred until sweep results tell us which paths matter most:

- `run_alignment` end-to-end smoke test (large fixture needed)
- Greedy capacity assignment cascade (priority; awaiting sweep)
- Forbidden-pair masking path
- Soma rule
- NBLAST-type-scores threshold path
- `resolve_alignment_paths` (needs filesystem mock)
- NTAC-side tests (script has no `--forbidden-matches`, no `--seed`, no
  stratified holdout yet — add after the methods are added)

See `project_alignment_test_plan.md` in memory for the full priority list.

## Fixtures

`conftest.py` loads `banc-alignment-run.py` via `importlib.util` because
the hyphen in the filename blocks plain `import`. The module is cached at
session scope, so the one-time import cost is paid once.
