# `hemibrain/` — hemibrain v1.2.1 reference deposit

Sister-dataset publish pipeline: assembles the hemibrain v1.2.1 metadata, edgelist, NBLAST scores and skeleton bundle and rsyncs it to `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/hemibrain_121/`. Parallels `banc/share/banc-sjcabs.R`. Most data-wrangling logic was originally developed in [`fafbpipeline`](https://github.com/flyconnectome/fafbpipeline); this directory is the GCS-publish-side mirror copied into bancpipeline so the public release owns the publish step.

## Scripts

- `hemibrain-meta.R` — Pull hemibrain per-neuron metadata via `neuprintr` and assemble the master meta CSV.
- `hemibrain-sjcabs.R` — Assemble + publish the hemibrain v1.2.1 SJCABS-style data bundle (meta + edgelist + NBLAST + synapses + skeletons) to GCS.

## Executable for users?

Partial. `hemibrain-sjcabs.R` needs HMS-O2 access (reads `/n/data1/.../hemibrain/`). `hemibrain-meta.R` runs anywhere with neuprint credentials but is rarely useful outside the pipeline.

See top-level [`README.md`](../README.md), [`docs/data-layout.md`](../docs/data-layout.md), and [`docs/data-products.md`](../docs/data-products.md).
