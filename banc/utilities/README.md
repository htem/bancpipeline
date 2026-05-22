# `banc/utilities/` — plotting + format-conversion helpers

Manual analysis / plotting helpers and one-off format converters. Not part of any production chain — `banc/share/banc-data.R` does not source any of these.

## Scripts

| File | Purpose |
|---|---|
| `banc-plot.R` | Density / beeswarm / heatmap diagnostics for cross-dataset matching + NT predictions. Manual exploration. |
| `banc-dimorphic-density.R` | Combined VNC + brain density plots of sexually-dimorphic presynapses across BANC, MANC, FAFB, maleCNS. |
| `banc-csv-feather-convert.R` | One-time migration: split the legacy monolithic `banc_root_positions.csv` into per-metric feathers. |
| `malecns-3dprint-stl.py` | Export BANC + maleCNS neuropil surfaces as watertight STL files for 3D printing. |

## Executable for users?

Partial. `malecns-3dprint-stl.py` runs anywhere with `navis` + the published neuropil meshes. The R scripts mostly need HMS-O2 paths.

See [`banc/README.md`](../README.md).
