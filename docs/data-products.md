# Data products — file-by-file

Every public artefact bancpipeline produces, with the producer script, its consumer schema doc, and where it lives. Companion to [`docs/data-layout.md`](data-layout.md) (which gives the overall storage map). For column-by-column schemas, follow the **Schema** column to the per-file docs in `BANC-project/manuscript/print/dataverse/documentation/`.

The paper Methods cites each file with the `[filename.feather]` square-bracket convention; every cited filename below appears both at Harvard Dataverse <https://doi.org/10.7910/DVN/7WTH1N> and at `gs://lee-lab_brain-and-nerve-cord-fly-connectome/`.

## Headline tables

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_888_meta.feather` | `banc/share/banc-data.R` (Section 1 META) | `banc_888_meta.md` | 188 162 rows × 79 cols. The headline per-neuron metadata. Soma + root position, region, side, neuromere, hemilineage, cell-type hierarchy (super_class > cell_class > cell_sub_class > cell_type), cross-dataset matches (FAFB v783, MANC v1.2.1, FANC v1.116, hemibrain v1.2.1, maleCNS v0.9), AN/DN behaviour cluster + super_cluster, CNS-network membership, NT prediction + verified, morphology metrics. |
| `banc_888_metrics.feather` | `banc/share/banc-data.R` (Section 2 METRICS) | `banc_888_metrics.md` | Per-neuron quantitative subset: cable length, volume, synapse counts, mitochondria count + volume, flow-centrality segregation index, primary-dendrite width. |

## Connectivity

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_888_edgelist_simple_v2.feather` | `banc/metrics/banc-calculate-connectivity.R --source v2` (then aggregated by `banc/share/banc-data.R`) | `banc_888_edgelist_simple_v2.md` | Neuron-to-neuron edgelist (paper version). One row per ordered (pre, post) pair with `count`, `norm` (fraction of post's total input), per-neuron totals. v2 = `size ≥ 5`. |
| `banc_888_edgelist_simple_v3.feather` | same, `--source v3` | `banc_888_edgelist_simple_v3.md` | v3 = `size ≥ 10`. NT predictions are still computed on v2. |
| `banc_888_edgelist_split_v2.feather` | `banc/share/banc-data.R` (Section 5 SPLIT) | `banc_888_edgelist_split_v2.md` | Compartment-resolved edgelist. Axon / dendrite / soma / primary-dendrite labels on both sides via flow-centrality. |
| `banc_888_synapses_v2_enriched.parquet` | `banc/share/banc-data.R` (Section 3 SYNAPSES) | `banc_888_synapses_v2_enriched.md` | 218 460 852-row per-synapse table. Pre/post root IDs, 3D coordinates, neuropil, region, per-synapse NT classifier output, compartment labels. |
| `banc_888_synapses_v3_enriched.parquet` | same, `--source v3` | `banc_888_synapses_v3_enriched.md` | v3 variant. 259 409 001 synaptic links. |
| `banc_888_synapses_v2_human_readable.csv.gz` | `banc/metrics/banc-calculate-connectivity.R` upstream | `banc_888_synapses_v2_human_readable.md` | Raw CAVE export with human-readable column headers. |
| `banc_888_synapses_v3_human_readable.csv.gz` | same | `banc_888_synapses_v3_human_readable.md` | v3 variant. |

## Neurotransmitter

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_888_neurotransmitter_prediction_v2.csv` | `banc/metrics/banc-calculate-ntpred.R` | `banc_888_neurotransmitter_prediction_v2.md` | Per-neuron NT summary. Per-NT presynapse counts, argmax class, confidence, cell-type-consensus check. v2 = post-synapse-size threshold ≥ 5. |
| `banc_nt_prediction_w_sizethresh_5_11102025.parquet` | external (synister_banc by K. Dasari) — bancpipeline consumes | `banc_nt_prediction_w_sizethresh_5_11102025.md` | Per-synapse NT classification output from the 3D-CNN classifier (eight-class softmax). v2. |
| `banc_nt_prediction_v3_w_sizethresh_10_05042026.parquet` | external | `banc_nt_prediction_v3_w_sizethresh_10_05042026.md` | v3 variant. |
| `banc_nt_ground_truth.csv` | hand-curated; ingested by `banc-calculate-ntpred.R` | `banc_nt_ground_truth.md` | 60 394 ground-truth-labelled neurons across 3 379 cell types from the literature + cross-matched FAFB/MANC/hemibrain. |

## NBLAST cross-dataset

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_fafb_783_nblast.feather` | `banc/nblast/banc-fafb-nblast.R` → `banc-nblast-compile.R` | `banc_fafb_783_nblast.md` | Pairwise NBLAST score table BANC ↔ FAFB v783. |
| `banc_manc_v1.2.1_nblast.feather` | `banc/nblast/banc-manc-nblast.R` → compile | `banc_manc_v1.2.1_nblast.md` | BANC ↔ MANC v1.2.1. |
| `banc_fanc_1116_nblast.feather` | `banc/nblast/banc-fanc-nblast.R` → compile | `banc_fanc_1116_nblast.md` | BANC ↔ FANC v1.116. |
| `banc_hemibrain_v1.2.1_nblast.feather` | `banc/nblast/banc-hemibrain-nblast.R` → compile | `banc_hemibrain_v1.2.1_nblast.md` | BANC ↔ hemibrain v1.2.1. |
| `banc_malecns_v0.9_nblast.feather` | `banc/nblast/banc-malecns-nblast.R` → compile | `banc_malecns_v0.9_nblast.md` | BANC ↔ maleCNS v0.9. |
| `banc_mirror_nblast.feather` | `banc/nblast/banc-nblast-lr.R` → compile | `banc_mirror_nblast.md` | BANC self-matched after thin-plate-spline mirror registration. |
| `banc_native_nblast.feather` | `banc/nblast/banc-nblast-native.R` → compile | `banc_native_nblast.md` | BANC self-NBLAST in native space (no template-space registration). |

NBLAST publication to GCS is handled by `banc/share/banc-nblast-share-gcs.R`.

## Influence

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `influence/all_to_all/chunk_NNNN.parquet` | `banc/influence/banc-build-influence.R` (SLURM-array sharded) | `influence_all_to_all.md` | All-to-all steady-state activation matrix, partitioned. Built via PETSc/SLEPc sparse solve. |
| `influence_all_to_effector_subclass.parquet` | `banc/influence/banc-aggregate-influence.R` | `influence_all_to_effector_subclass.md` | Pooled source / target aggregation. Influence onto effector cell sub-classes. |
| `influence_sensory_subclass_to_all.parquet` | same | `influence_sensory_subclass_to_all.md` | Pooled aggregation. Influence from sensory sub-classes. |

## Graph statistics

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_888_betweenness_all_to_all_v2.csv` | `banc/betweenness/banc-betweenness-run.R --source v2` | `banc_888_betweenness.md` | Unnormalised, unweighted, directed all-to-all betweenness centrality (Brandes algorithm via `igraph`). Paper Methods §"Betweenness centrality". |
| `banc_888_betweenness_all_to_all_v3.csv` | same, v3 | same | v3 variant. |
| `banc_888_betweenness_afferent_to_efferent_v2.csv` | same | same | Source–target betweenness restricted to sensory → effector pairs. |
| `banc_888_betweenness_afferent_to_efferent_v3.csv` | same, v3 | same | v3 variant. |

## CNS networks / clustering

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_888_cns_network_spectral_clustering_v2.csv` | `banc/clustering/banc-spectral-clustering-run.R --source v2 --cluster-count 13 --min-connection-strength 1` | `banc_888_cns_network_spectral_clustering_v2.md` | 13 graph-Laplacian clusters + UMAP coords (n_neighbors=100, cosine, min_dist=0, seed=3). Paper Methods §"Spectral clustering". |
| `banc_888_cns_network_spectral_clustering_v3.csv` | same, v3 | (forthcoming) | v3 variant. |

## Skeletons + meshes + MIPs

| File | Producer script | Schema doc | Brief |
|---|---|---|---|
| `banc_swc_skeletons.zip` | `banc/share/banc-export-skeletons.R` | `banc_swc_skeletons.md` | One `.swc` per neuron in BANC voxel space (L2 skeletons via pcg_skel). |
| `banc_neuron_meshes.zip` | `banc/metrics/banc-obj.R` → packaging | `banc_neuron_meshes.md` | One `.obj` mesh per proofread neuron. |
| `banc_neuropil_meshes.zip` | `banc/transforms/` neuropil-mesh generation | `banc_neuropil_meshes.md` | Closed CNS surface + per-neuropil sub-meshes (alpha-shape on synapse cloud). |
| `banc_color_mips.zip` | `banc/nblast/banc-fafb-nblast-images.R` etc. (template-space MIPs) | `banc_color_mips.md` | Color-depth MIPs registered to JRC2018_Unisex_20x_HR (brain) + JRC2018_VNC_Unisex_40x_DS (VNC). NeuronBridge-compatible. |
| `banc_microCT.zip` | external (Zeiss microCT) | `banc_microCT.md` | Pre-sectioning microCT volume of the resin-embedded sample (Methods §"Specimen"). |
| `banc_template_spaces.zip` | `banc/transforms/` | `banc_template_spaces.md` | BANC ↔ JRC2018F / JRC2018VNCF template registrations + the BANC space definition. |
| `registration_brain_jrc2018f.zip` | `banc/transforms/` (elastix) | `registration_brain_jrc2018f.md` | Brain elastix registration to JRC2018F. |
| `registration_vnc_jrc2018vncf.zip` | `banc/transforms/` (elastix) | `registration_vnc_jrc2018vncf.md` | VNC elastix registration to JRC2018VNCF. |
| `banc_mitochondria_v1.feather` | `banc/metrics/banc-calculate-volumes.R` (mitochondria branch) | `banc_mitochondria_v1.md` | Per-mitochondrion detections + volumes (~40 M mitochondria). |

## Annotation parquets (CAVE-derived)

| File | Producer | Schema doc | Brief |
|---|---|---|---|
| `backbone_proofread.parquet` | CAVE export (BANC project) | `backbone_proofread.md` | Per-neuron proofread-status snapshot at v888. |
| `cell_info.parquet` | CAVE community annotations | `cell_info.md` | Free-form community cell-info annotations. |
| `cell_representative_point.parquet` | CAVE | `cell_representative_point.md` | One representative xyz per neuron, keys to the master codex_annotations table. |
| `codex_annotations.parquet` | CAVE (curated by core BANC team) | `codex_annotations.md` | Authoritative annotation table federated to FlyWire Codex. |
| `somas_v1.parquet` | CAVE (somas_v1a + somas_v1b corrections) | `somas_v1.md` | Per-soma midpoint coordinates. |
| `peripheral_nerves.parquet` | CAVE | `peripheral_nerves.md` | Per-nerve cross-section neuron-profile annotations (47 seed planes). |
| `proofreading_notes.parquet` | CAVE | `proofreading_notes.md` | "Roughly proofread" classifications + tracing-issue notes. |
| `neck_connective_y92500.parquet`, `neck_connective_y121000.parquet` | CAVE | `neck_connective_y92500.md` / `neck_connective_y121000.md` | Two cross-section planes through the neck connective; basis for AN/DN enumeration. |
| `synapse_neuropil_lookup_v2.parquet` | `banc/share/banc-publish-synapse-lookups.R` | `synapse_neuropil_lookup_v2.md` | Per-synapse `syn_id → neuropil` lookup (syn_id is stable across segmentation versions). |
| `synapse_neuropil_lookup_v3.parquet` | same, v3 | `synapse_neuropil_lookup_v3.md` | v3 variant. |

## Supplemental + reference

| File | Producer | Schema doc | Brief |
|---|---|---|---|
| `2024-09-20_aelysia_synapse_sample_complete.csv` | Aelysia LTD manual review (Methods §"Synapse detection evaluation") | `aelysia_synapse_review_sample.md` | 4 648 synapse manual-review sample. Used to set the v2 size-threshold ≥ 5. |
| `banc_problem_regions.csv` | hand-curated | `banc_problem_regions.md` | Known dataset artefact regions (misalignment, data loss) flagged at xyz boxes. Supplemental Data 10 in the paper. |
| `behavior.zip` | pre-EM Y-maze trials of the BANC fly | `behavior.md` | Pre-fixation Y-maze handedness data (~582 choices). |

## Code archives bundled into the Dataverse deposit

| Zip | Source repo | Schema doc |
|---|---|---|
| `bancpipeline.zip` | `htem/bancpipeline` (this repo) | `bancpipeline_archive.md` |
| `banc_project_archive.zip` | `htem/BANC-project` | `banc_project_archive.md` |
| `bancr_archive.zip` | `natverse/bancr` | `bancr_archive.md` |
| `synister_banc.zip` | `htem/synister_banc` | `synister_banc_archive.md` |
| `connectome_influence_calculator.zip` | `jdrugo/connectome-influence-calculator` | `connectome_influence_calculator_archive.md` |
| `influencer.zip` | `natverse/influencer` | `influencer_archive.md` |
| `nat_ggplot_archive.zip` | `natverse/nat.ggplot` | `nat_ggplot_archive.md` |
| `fly_connectome_data_tutorial_archive.zip` | `sjcabs/fly_connectome_data_tutorial` | `fly_connectome_data_tutorial_archive.md` |
| `drosophila_neurotransmitters_archive.zip` | `funkelab/drosophila_neurotransmitters` | `drosophila_neurotransmitters_archive.md` |
| `drosophila_neuropeptides_archive.zip` | (companion to NT archive) | `drosophila_neuropeptides_archive.md` |
| `the_banc_fly_connectome_archive.zip` | `jasper-tms/the-BANC-fly-connectome` | `the_banc_fly_connectome_archive.md` |
| `banc_python_archive.zip` | `pypi.org/project/banc` | `banc_python_archive.md` |

## Versions to know

- **Materialization version**: `v888` (paper, snapshot 2026-04-17). `v626` was the preprint version — superseded; do not use.
- **Synapse version**: `v2` (paper, `size ≥ 5`). `v3` (`size ≥ 10`) is for future work. NT predictions are still computed on `v2`.
- **NBLAST reference versions** are listed in each per-target row above; they pin to the matching upstream connectome release.
