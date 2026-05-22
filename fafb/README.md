# `fafb/` — FAFB v783 reference deposit

Sister-dataset publish pipeline: assembles the FAFB (FlyWire v783) metadata, edgelist, NBLAST scores and skeleton bundle and rsyncs it to `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/fafb_783/`. Parallels `banc/share/banc-sjcabs.R`. Most data-wrangling logic was originally developed in [`fafbpipeline`](https://github.com/flyconnectome/fafbpipeline); this directory is the GCS-publish-side mirror copied into bancpipeline so the public release owns the publish step.

## Scripts

| File | Purpose |
|---|---|
| `fafb-meta.R` | Pull FAFB v783 per-neuron metadata from `franken_meta()` + flywire CAVE; assemble the master meta CSV used by NBLAST + alignment. |
| `fafb-sjcabs.R` | Assemble + publish the FAFB v783 SJCABS-style data bundle (meta + edgelist + NBLAST + synapses + skeletons) to GCS. |
| `sort-connectome-data.R` | Organise connectome match files on disk for downstream consumption. |

## Executable for users?

Partial. `fafb-sjcabs.R` needs HMS-O2 access (reads `/n/data1/.../malevnc/`, `/n/data1/.../fafb/`). `fafb-meta.R` runs anywhere with `flywire` CAVE credentials, but is rarely useful outside the pipeline.

See top-level [`README.md`](../README.md) for the publish flow, [`docs/data-layout.md`](../docs/data-layout.md) for storage paths, and [`docs/data-products.md`](../docs/data-products.md) for the published file schema.
