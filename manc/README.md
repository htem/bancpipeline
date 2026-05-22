# `manc/` — MANC v1.2.1 reference deposit

Sister-dataset publish pipeline: assembles the MANC v1.2.1 metadata, flow-centrality splits, edgelist, NBLAST scores and skeleton bundle and rsyncs it to `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/manc_121/`. Parallels `banc/share/banc-sjcabs.R`.

## Scripts

| File | Purpose |
|---|---|
| `manc-meta.R` | Pull MANC v1.2.1 per-neuron metadata + assemble the master meta CSV. |
| `manc-split.R` | Compute MANC flow-centrality axon/dendrite splits per neuron. |
| `manc-data.R` | Compile MANC release artefacts (metrics, edgelists) from split outputs. |
| `manc-sjcabs.R` | Assemble + publish the MANC v1.2.1 SJCABS-style data bundle to GCS. |

## Executable for users?

No — these need HMS-O2 access (read `/n/data1/.../malevnc/`, write to versioned save tree).

See top-level [`README.md`](../README.md), [`docs/data-layout.md`](../docs/data-layout.md), and [`docs/data-products.md`](../docs/data-products.md).
