# `banc/nblast/` — cross-dataset + mirror NBLAST

Computes pairwise NBLAST scores between BANC and each reference connectome (FAFB v783, MANC v1.2.1, FANC v1.116, hemibrain v1.2.1, maleCNS v0.9), plus left-right mirror NBLAST and within-BANC native NBLAST. Renders per-query review PNGs for manual inspection. Compiles all sources into per-dataset feathers and pushes the match columns to CAVE annotation tables.

Sits after [`banc/transforms/`](../transforms/) (which transforms reference meshes / skeletons into BANC space) and feeds [`banc/matching/`](../matching/), [`alignment/`](../../alignment/), and [`banc/update/banc-update-matches.R`](../update/).

## Scripts

| File | Purpose |
|---|---|
| `banc-make-proofread-ids.R` | Emit a `proofread.txt` of BANC root_ids; each downstream script honours `BANC_TEST_IDS_FILE`. |
| `banc-fafb-nblast.R` | Pairwise BANC ↔ FAFB v783 NBLAST. Resumable. |
| `banc-manc-nblast.R` | Pairwise BANC ↔ MANC v1.2.1 NBLAST. |
| `banc-fanc-nblast.R` | Pairwise BANC ↔ FANC v1.116 NBLAST (VNC + ANs + DNs + sensory + efferent). |
| `banc-hemibrain-nblast.R` | Pairwise BANC ↔ hemibrain v1.2.1 NBLAST (non-VNC only, native + mirrored). |
| `banc-malecns-nblast.R` | Pairwise BANC ↔ maleCNS v0.9 NBLAST (per super_class, motor + visceral first). |
| `banc-nblast-lr.R` | BANC self-NBLAST after thin-plate-spline left-right mirror. |
| `banc-nblast-native.R` | BANC self-NBLAST in native space (ipsi sister cells, partial duplicates). |
| `banc-{fafb,manc,fanc,hemibrain,malecns,lr}-nblast-images.R` | Render per-query NBLAST review PNGs (mesh comparisons). |
| `banc-nblast-search.R` | Manual: find BANC neurons matching a target cell_type via compiled NBLAST. |
| `banc-nblast-wrong-matches.R` | Process `*_PNG_MATCH_WRONG` SeaTable flags into all artefacts. Run before `banc-nblast-compile.R`. |
| `banc-nblast-compile.R` | Consolidate per-query NBLAST CSVs into per-dataset feathers. Largest script in the dir. |
| `banc-nblast-cave.R` | Sync compiled NBLAST match feathers to CAVE `cell_match` tables. |
| `banc-nblast-plot.R` | NBLAST score distribution + per-region matching inventory plots. |
| `banc-sort-folders.R` | Bidirectional rsync between O2 working tree + fileserve mirror. |

## Executable for users?

No — needs HMS-O2 paths, CAVE credentials, and the transformed-meshes tree from `banc/transforms/`. Compiled NBLAST feathers are published via GCS (`banc_<dataset>_<ver>_nblast.feather`) — fetch those rather than re-running.

See [`banc/README.md`](../README.md), [`docs/png-matching.md`](../../docs/png-matching.md), and top-level [`README.md`](../../README.md).
