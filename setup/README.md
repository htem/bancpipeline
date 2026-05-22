# `setup/` — one-time environment + reference assets

One-shot helpers for standing up the BANC pipeline on a fresh machine (HMS-O2 or local) and reference assets used by downstream scripts. Not part of any recurring chain.

## Files

| File | Purpose |
|---|---|
| `hms_setup.R` | Interactive recipe for setting up the HMS-O2 R / Python environment: module loads, library installs, CAVE / SeaTable / neuprint credential configuration, smoke tests. |
| `create_cave_table.py` | One-off helper to create the BANC `cell_match` CAVE annotation table (NBLAST + PNG-review matches push into this). |
| `paper_colours.csv` | Canonical hex codes for paper figures (super_class, region, neurotransmitter palettes). |
| `imported_meshes/` | Brain / VNC neuropil mesh assets used by visualisation scripts. |

## Executable for users?

Partial. `hms_setup.R` is HMS-O2-specific (lists `module load gcc/14.2.0` etc.) but the recipe is also a useful checklist for other clusters. `create_cave_table.py` is a one-off and was run when the BANC CAVE instance was first stood up.

See top-level [`README.md`](../README.md) and [`docs/data-layout.md`](../docs/data-layout.md).
