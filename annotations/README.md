# `annotations/` — reference annotation artefacts + run logs

Reference annotations and run logs consumed or produced by `banc/annotations/` and `banc/meta/`. Distinct from [`banc/annotations/`](../banc/annotations/), which holds the *scripts* — this directory is the *data*.

## Files

| File | Purpose |
|---|---|
| `annotations.md` | Controlled-vocabulary reference: every column in `banc_meta` with allowed values + provenance. Maintained by hand. |
| `cell_info_annotations.xlsx` | Cell-info reference spreadsheet (NT controlled vocab, neuropeptides, hemilineage table) — input to `banc/annotations/banc-cell-info.R`. |
| `banc_meta_from_cell_info_update_log.txt` | Log of cell_info → SeaTable propagation runs. |
| `franken_meta_issues.txt` | Running notes on franken-brain metadata issues (consumed by `banc/franken/banc-frankenbrain.R` + `banc/annotations/franken-annotations-fix.R`). |
| `cell_ids_tag.txt` | Empty placeholder tag file. |

## Executable for users?

No — these are reference artefacts, not scripts.

See [`banc/annotations/`](../banc/annotations/) for the scripts that read these and the top-level [`README.md`](../README.md).
