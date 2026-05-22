# `malecns/` — maleCNS v0.9 reference deposit

Sister-dataset publish pipeline: assembles the maleCNS v0.9 metadata, NBLAST scores and skeleton bundle and rsyncs it to `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09/`. Parallels `banc/share/banc-sjcabs.R`.

## Scripts

| File | Purpose |
|---|---|
| `malecns-meta.R` | Pull maleCNS v0.9 per-neuron metadata via `malecns` package and assemble the master meta CSV. Large file — contains nearly all maleCNS data-wrangling logic for the BANC project. |
| `malecns-skels.R` | Download maleCNS SWC skeletons into the BANC matching tree. |
| `malecns-sjcabs.R` | Assemble + publish the maleCNS v0.9 SJCABS-style data bundle to GCS. |

## Executable for users?

No — these need HMS-O2 access and Janelia neuPrint credentials.

See top-level [`README.md`](../README.md), [`docs/data-layout.md`](../docs/data-layout.md), and [`docs/data-products.md`](../docs/data-products.md).
