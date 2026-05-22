# `banc/meta/` — master metadata assembly

Builds the master `banc_meta.csv` from L2 skeleton metrics + NBLAST + the `banc_ids.csv` identity table. This CSV is then enriched by [`banc/share/banc-data.R`](../share/) into the released `banc_888_meta.feather`. The compile here covers the everything-from-disk pass; the share-side wraps it for publication.

Sits after [`banc/metrics/`](../metrics/) (writes the per-metric feathers) and [`banc/nblast/banc-nblast-compile.R`](../nblast/) (writes the compiled NBLAST feathers), and feeds [`banc/share/banc-data.R`](../share/).

## Scripts

- `banc-meta.R` — Compile `banc_meta.csv` from L2 metrics + NBLAST CSVs + IDs. Derives `auto:*` cell types above score thresholds.
- `banc-hemilineages.R` — One-off generator for the hemilineage mapping worksheet (Ito-Lee + Hartenstein).
- `banc-meta-fix.R` — Ad-hoc spot-fixer for franken / flywire / MANC metadata bugs. Not part of any automated run.

## Executable for users?

No — needs HMS-O2 paths and SeaTable credentials.

See [`banc/README.md`](../README.md), [`docs/data-products.md`](../../docs/data-products.md).
