# `banc/share/` — the publish boundary

**The single boundary between the internal pipeline and the public GCS bucket.** Every output that ships externally is written or rsynced by a script in here. Reads the per-metric / per-NBLAST / per-split intermediates from disk, assembles the versioned release artefacts (`banc_888_meta.feather`, `banc_888_metrics.feather`, `banc_888_synapses_v{2,3}_enriched.parquet`, `banc_888_edgelist_simple_v{2,3}.feather`, `banc_888_edgelist_split_v2.feather`, the NT predictions CSV, the per-region cutout neuron / synapse / connectivity feathers), and pushes them to `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/`.

Sits at the bottom of the pipeline — consumes outputs from every other `banc/<subdir>/` and feeds GCS, the Harvard Dataverse mirror, and the FlyWire Codex import.

## Scripts

| File | Purpose |
|---|---|
| `banc-data.R` | Assemble the versioned BANC release: meta, metrics, synapses, edgelists. The biggest script in the dir; one-stop release driver. |
| `banc-sjcabs.R` | Export per-neuron SWCs, region cut-outs, neuropil meshes for collaborator share. |
| `banc-sjcabs-upload.R` | Sync exported BANC / maleCNS / MANC connectomes to GCS. |
| `banc-nblast-share-gcs.R` | Publish per-dataset NBLAST feathers + reviewed CSVs to GCS. |
| `banc-publish-synapse-lookups.R` | Publish version-stable syn_id → neuropil/region/side parquets. |
| `banc-export-skeletons.R` | Package per-neuron SWCs into the versioned release tree (detailed preferred, L2 fallback). |

## Executable for users?

No — needs HMS-O2 paths + gsutil credentials with write access to the public bucket. End users should download the published artefacts directly from GCS / Dataverse rather than re-run these scripts.

See [`banc/README.md`](../README.md), top-level [`README.md`](../../README.md), and [`docs/data-products.md`](../../docs/data-products.md).
