# `banc/betweenness/` — all-to-all + sensory→effector betweenness centrality

Source-target betweenness centrality on the BANC edgelist via `networkx`. Produces the `banc_888_betweenness_*.csv` artefacts in the public release. Sits after `banc/metrics/banc-calculate-connectivity.R` (which writes the edgelist) and feeds `banc/share/banc-data.R` (which packages the CSV into the release tree).

## Scripts

- `banc-betweenness.py` — Compute source-target betweenness for one of three task presets (all-to-all, afferent → efferent, others). Reads the edgelist parquet + meta feather, writes a CSV.
- `banc-betweenness-run.R` — Dispatcher: invokes `banc-betweenness.py` from R, auto-detecting the data dir (O2 versioned tree vs local GCS mirror) and resolving `--source v2|v3`.

## Executable for users?

Partial. With local copies of `banc_888_meta.feather` + `banc_888_edgelist_simple_v{2,3}.feather` (downloadable from GCS), `banc-betweenness.py` runs anywhere with `networkx` + `pyarrow`. The `-run.R` wrapper is most useful on O2.

See [`banc/README.md`](../README.md), top-level [`README.md`](../../README.md), and [`docs/data-products.md`](../../docs/data-products.md).
