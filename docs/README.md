# `docs/` — companion documentation

Long-form companions to the top-level [`README.md`](../README.md). Kept as separate files so the main README stays a high-level overview.

| File | Covers |
|---|---|
| [`code-style.md`](code-style.md) | Comment / header / docstring conventions for `.R` and `.py` files. UK English in prose. |
| [`data-layout.md`](data-layout.md) | Where every file lives on O2 + GCS; the `banc.*.save.path` env-var to path mapping. |
| [`data-products.md`](data-products.md) | Producer → output map. Which script writes which `banc_888_*.{feather,parquet,csv}`, and the downstream `BANC-project` consumer. |
| [`png-matching.md`](png-matching.md) | Historical PNG-review matching workflow (superseded for cell-typing by the annealed alignment in [`alignment/`](../alignment/), but its outputs still feed the released metadata as `INVESTIGATE` / `TRACING_ISSUE` / `NO_*_MATCH` tags). |

Per-column release-file schemas live in [`BANC-project/manuscript/print/dataverse/documentation/`](https://github.com/htem/BANC-project/tree/main/manuscript/print/dataverse/documentation), not here.
