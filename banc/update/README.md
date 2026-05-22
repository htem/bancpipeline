# `banc/update/` — SeaTable read/write

The integration boundary with the SeaTable curation database. Refreshes BANC root_ids from CAVE, snapshots SeaTable to disk, and pushes computed columns (cell_type cascade, status flags, side assignments, NBLAST + PNG matches, NT predictions, per-neuron metrics, DCV density) back to `banc_meta`. The top-level driver `banc-update-seatable.R` chains the others.

Sits at the head of the pipeline — `banc-ids.R` defines which neurons exist for the current run; the rest push computed columns back to SeaTable after the per-metric / per-NBLAST scripts finish.

## Scripts

| File | Purpose |
|---|---|
| `banc-ids.R` | Discover BANC neurons + rebuild the canonical `banc_ids.csv`. Joins SeaTable to CAVE proofreading-notes / status / nuclei. |
| `banc-ids888.R` | Standalone root_888 backfill for SeaTable banc_meta (independent of `banc-ids.R`). |
| `banc-ids890.R` | Standalone root_890 backfill for the v888 → v890 migration. |
| `banc-updateids.R` | Thin wrapper around `bancr:::banctable_updateids()` — refreshes stale root_ids in SeaTable. |
| `banc-delete.R` | Sweep per-neuron file trees (OBJ / SWC / split / synapses / images) for files whose root_id no longer appears in SeaTable; reports coverage. |
| `banc-update-seatable.R` | Master orchestrator: updates root_ids, snapshots, merges `banc_ids.csv` + `banc_meta.csv` into SeaTable. |
| `banc-update-matches.R` | Push NBLAST + PNG match columns to SeaTable; detect conflicts. |
| `banc-update-celltypes.R` | Derive cell types from cross-dataset matches and push the dataset-specific + cascade columns to SeaTable. |
| `banc-update-status.R` | Compute + push status flags (`LR_TYPE_CONFLICT`, `SIDE_CONFLICT`, `UNROOTED`, `TOO_SMALL`, `TRACING_ISSUE_RESOLVED`, etc.) and side assignments. |
| `banc-update-metrics.R` | Join per-metric feathers and push metric columns to SeaTable. |
| `banc-update-ntpred.R` | Push neurotransmitter predictions to SeaTable (`neurotransmitter_predicted_v{2,3}`). |
| `banc-dcv-density.R` | Push per-neuron soma DCV count + density to SeaTable (joined on `nucleus_id`). |

## Executable for users?

No — needs SeaTable + CAVE credentials and HMS-O2 paths.

See [`banc/README.md`](../README.md), [`o2/production/o2_banc_update.sh`](../../o2/production/), top-level [`README.md`](../../README.md).
