# Data layout — O2, GCS, Dataverse

Where every file the pipeline reads or writes lives. Companion to [`docs/data-products.md`](data-products.md) (which goes file-by-file with schemas).

The pipeline is HMS-O2-rooted: scripts read/write under `/n/data1/hms/neurobio/wilson/banc/` for working data, and under `/n/data1/hms/neurobio/wilson/connectomes/banc/banc_<version>/` for versioned published outputs. The published outputs are mirrored to GCS, and a frozen snapshot is deposited on Harvard Dataverse.

## O2 — working data

Root: `/n/data1/hms/neurobio/wilson/banc/` (env var: `banc.save.path`).

```
/n/data1/hms/neurobio/wilson/banc/
├── meta/                     # banc_meta.csv + snapshots/ + reference metas (fanc, hemibrain, …)
│                             # (banc.meta.save.path)
├── connectivity/             # per-version synapse parquets, edgelists, lookups (intermediates)
│                             # (banc.connectivity.save.path)
├── synapses/                 # raw synapse downloads from CAVE (per version)
│                             # (banc.synapses.save.path)
├── synapses_v3/              # v3-specific (raw shards, processed batches)
├── obj/                      # per-neuron .obj meshes (banc.obj.save.path)
├── swc/                      # detailed SWC skeletons (banc.swc.save.path)
├── l2/                       # L2 SWC skeletons (banc.l2swc.save.path)
├── split/                    # detailed axon/dendrite split outputs (banc.split.save.path)
├── l2split/                  # L2 axon/dendrite split outputs (banc.l2split.save.path)
│                             # ├── synapses/{root_id}.csv   per-neuron labelled synapses
│                             # ├── swc/{root_id}.swc        labelled skeletons
│                             # └── metrics/banc_metrics_N.csv
├── metrics/                  # per-neuron metric tables (banc.metrics.save.path)
├── nt/                       # neurotransmitter CSVs (banc.nt.save.path)
├── influence/                # all-to-all influence intermediates (banc.influence.save.path)
├── matching/                 # cross-dataset PNG match images (see docs/png-matching.md)
│                             # mirror/, fafb/, manc/, fanc/, hemibrain/, malecns/
│                             # each with images/ + correct/ subtrees
├── nblast_<dataset>/         # per-target NBLAST result CSVs
└── ...
```

Notes:

- `banc.save.path` is set in `banc/banc-startup.R:248`. Off-O2 it does not resolve; off-O2 work uses `--source gcs` flags to fetch from GCS.
- The `obj/`, `swc/`, `l2/`, `split/`, `l2split/` trees can each hold ~150k per-neuron files. They are not synced anywhere as raw files — only the zipped exports below.

## O2 — versioned published outputs

Root: `/n/data1/hms/neurobio/wilson/connectomes/banc/banc_<version>/` (env var: `banc.versioned.save.path`; current version 888).

```
/n/data1/hms/neurobio/wilson/connectomes/banc/banc_888/
├── banc_888_meta.feather                       # master metadata (188k rows × 79 cols)
├── banc_888_metrics.feather                    # per-neuron metrics subset
├── banc_888_synapses_v2.parquet                # raw synapse table (v2)
├── banc_888_synapses_v2_enriched.parquet       # enriched (NT + neuropil + compartment)
├── banc_888_synapses_v3.parquet                # raw synapse table (v3)
├── banc_888_synapses_v3_enriched.parquet       # enriched (v3)
├── banc_888_edgelist_simple_v2.feather         # neuron-to-neuron, count >= 5
├── banc_888_edgelist_simple_v3.feather         # neuron-to-neuron, count >= 10
├── banc_888_edgelist_split_v2.feather          # compartment-resolved
├── banc_888_neurotransmitter_prediction_v2.csv # per-neuron NT
├── banc_<dataset>_<ver>_nblast.feather         # one per cross-dataset (fafb / manc / fanc / hemibrain / malecns)
├── banc_mirror_nblast.feather                  # left-right NBLAST
├── influence/all_to_all/chunk_NNNN.parquet     # all-to-all influence shards
├── influence_all_to_effector_subclass.parquet
├── influence_sensory_subclass_to_all.parquet
├── banc_banc_space_swc/{root_id}.swc           # L2 skeletons (zip on Dataverse)
└── (other versioned exports)
```

The versioned tree is what gets rsynced to GCS + frozen for the Dataverse deposit.

## GCS — public bucket

URL: `https://console.cloud.google.com/storage/browser/lee-lab_brain-and-nerve-cord-fly-connectome/` (also `gs://lee-lab_brain-and-nerve-cord-fly-connectome/`). Public, no auth required.

```
gs://lee-lab_brain-and-nerve-cord-fly-connectome/
├── compiled_data/banc_888/                     # versioned outputs (mirrors O2's banc_888/ tree)
├── compiled_data/banc_<ver>/                   # one tree per BANC version
├── compiled_data/fafb_783/                     # FAFB inputs the BANC alignment uses
├── compiled_data/manc_v1.2.1/                  # MANC inputs
├── compiled_data/<other-datasets>/             # hemibrain, malecns, fanc
├── synapses/v1.1/                              # per-version raw synapse + NT prediction parquets
├── synapses/v2.0/
├── synapses/v3.0/
├── neuron_connectivity/v888/                   # human-readable CAVE synapse exports
├── neuron_annotations/v888/                    # CAVE-derived annotation parquets
│                                               # (cell_info, codex_annotations, cell_representative_point,
│                                               #  somas_v1, backbone_proofread, peripheral_nerves, etc.)
├── nblast/                                     # shared NBLAST result feathers
└── matching/banc_*.zip                         # PNG matching archives (see docs/png-matching.md)
```

This is the "live" public mirror — it can update past the paper snapshot as the project evolves.

## Harvard Dataverse — frozen DOI snapshot

DOI: <https://doi.org/10.7910/DVN/7WTH1N>. CC BY 4.0.

The Dataverse deposit is a frozen snapshot of the GCS bucket at the paper-release point. Same filenames, same schemas. The canonical filename index is [`BANC-project/manuscript/print/banc_data_locations.md`](https://github.com/htem/BANC-project/blob/main/manuscript/print/banc_data_locations.md), and per-file column schemas live one directory deeper at `BANC-project/manuscript/print/dataverse/documentation/<filename>.md`.

The paper Methods cites each Dataverse file with the `[filename.feather]` convention.

## bancpipeline ↔ outputs map (which script writes which file)

The pipeline-stage to output mapping is fully detailed in [`docs/data-products.md`](data-products.md). Headline:

- `banc/share/banc-data.R` — assembles `banc_888_meta.feather`, `banc_888_metrics.feather`, all enriched-synapse parquets, both edgelists, and the split edgelist. Called by the master rebuild orchestrator.
- `banc/metrics/banc-calculate-connectivity.R` — `banc_888_synapses_v2.parquet` + `banc_888_edgelist_simple_v2.feather` (and v3 variants).
- `banc/metrics/banc-calculate-synapses.R` — per-neuron synapse-count summary; feeds banc-data.R.
- `banc/metrics/banc-calculate-split.R` — flow-centrality compartment labels per synapse; `l2split/` tree.
- `banc/nblast/banc-<target>-nblast.R` — per-target NBLAST feathers (`banc_fafb_783_nblast.feather`, etc.).
- `banc/nblast/banc-nblast-compile.R` — final consolidated NBLAST feathers.
- `banc/clustering/banc-spectral-clustering-run.R` — `banc_888_cns_network_spectral_clustering_v2.csv`.
- `banc/influence/banc-build-influence.R` + `banc-aggregate-influence.R` — `influence/all_to_all/*.parquet` shards then aggregated.
- `banc/share/banc-nblast-share-gcs.R` — publishes NBLAST feathers to GCS.
- `banc/share/banc-export-skeletons.R` — packages L2 SWCs into the per-neuron release tree.
- `banc/share/banc-publish-synapse-lookups.R` — publishes per-version synapse → neuropil/region/side lookup parquets to GCS.
- `banc/share/banc-sjcabs.R` + `banc-sjcabs-upload.R` — packages the sjcabs collaborator-share zips.

## Re-running notes

- Most stages depend on `banc/banc-startup.R` resolving its O2 paths. Off-O2 they error.
- A few stages take `--source gcs` and can fetch from GCS instead: `alignment/banc-alignment-prep.R`, `banc/clustering/banc-spectral-clustering*.R/py`, `banc/betweenness/banc-betweenness-run.R`.
- The full v888 rebuild is `o2/production/o2_banc_v888_rebuild.sh` (~250 GB / priority partition); see `o2/README.md` for the recurring-chain summary.
