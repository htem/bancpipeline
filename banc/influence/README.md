# `banc/influence/` — all-to-all influence scoring

Computes all-to-all influence scores (steady-state activation between every source / target neuron pair) via the `influencer` R package, which wraps a Python PETSc / SLEPc backend. Output is sharded partitioned parquets at `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/influence/all_to_all/`. Paper Methods §"Influence".

Sits after `banc/metrics/banc-calculate-connectivity.R` (edgelist input) and before `banc/share/banc-data.R` (which references the influence parquets in the release manifest).

## Scripts

| File | Purpose |
|---|---|
| `banc-build-influence.R` | Compute all-to-all influence scores. Resumable; supports SLURM-array sharding via `BANC_INFLUENCE_SHARD_{IDX,TOTAL}`. |
| `banc-aggregate-influence.R` | Aggregate the sharded all-to-all parquet into sensory-subclass→all and all→effector-subclass parquets; pushes both to GCS. |
| `banc-sync-influence.R` | Validate completeness + spot-check + rsync influence shards to GCS. Run after all build-shard jobs finish. |
| `HOW_TO_RUN.txt` | Operator notes — SLURM-array submission recipe + chunk math. |

## Executable for users?

No — needs HMS-O2 (PETSc / SLEPc build, hundreds of GB scratch, SLURM array submission rights).

See [`banc/README.md`](../README.md), [`o2/production/o2_banc_influence.sh`](../../o2/production/), top-level [`README.md`](../../README.md), and [`docs/data-products.md`](../../docs/data-products.md).
