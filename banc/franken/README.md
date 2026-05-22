# `banc/franken/` — frankenbrain composite-dataset assembly

Builds the "franken-brain" meta CSV: a single cross-dataset table that merges FlyWire (FAFB) + MANC metadata, using BANC-mediated cell_type bridges. The franken CSV is what `franken_meta()` loads everywhere, so the rest of the pipeline can reason about brain + VNC cell types as a single vocabulary. Seed-label organisation (descending / ascending hierarchical clusters) is the second concern.

Sits upstream of `banc/annotations/`, `banc/matching/`, and `banc/meta/` — they all consume `franken_meta()`.

## Scripts

- `banc-frankenbrain.R` — Build the franken-brain meta CSV (FAFB + MANC merge). Canonical reference loaded by `franken_meta()`.
- `franken-seeds.R` — Organise descending / ascending neuron seed labels (`SD_*` / `ED_*` → `DN_*`, `SA_*` / `EA_*` → `AN_*`) and summarise cluster composition into seed worksheets consumed by `BANC-project`.

## Executable for users?

No — needs HMS-O2 paths and SeaTable credentials.

See [`banc/README.md`](../README.md).
