# `banc/` — main BANC pipeline tree

This is the working tree of the BANC pipeline. Every script under `banc/` operates on the BANC connectome itself: pulling raw outputs from CAVE + SeaTable, computing per-neuron metrics, building synaptic edgelists, matching neurons to reference connectomes, and assembling the public release artefacts in `banc/share/`.

The cross-dataset alignment algorithm (the paper's headline algorithmic contribution) lives at the top-level [`alignment/`](../alignment/) — it's been promoted out of `banc/` because it's intended to be reusable for any pair of connectomes, not BANC-specific.

For where each output lands on GCS / Dataverse, see [`docs/data-products.md`](../docs/data-products.md). For overall pipeline flow + GCS layout, see the top-level [`README.md`](../README.md).

## Layout

| Subdir | Role in the pipeline |
|---|---|
| [`update/`](update/) | SeaTable read/write — ID refresh from CAVE, push curated cell-type / status updates back to `banc_meta` |
| [`metrics/`](metrics/) | Per-neuron computation: L2 skeletons, regions, synapses, axon/dendrite split, cable length, volume, neurotransmitter prediction |
| [`nblast/`](nblast/) | Cross-dataset NBLAST against FAFB, MANC, FANC, hemibrain, maleCNS, plus left-right mirror; compile + publish |
| [`matching/`](matching/) | Cross-validation of NBLAST + connectivity-based matches; sexually dimorphic / abdominal neuromere typing |
| [`clustering/`](clustering/) | Spectral clustering of the CNS connectivity graph (k = 13, paper Methods §"Spectral clustering") |
| [`influence/`](influence/) | All-to-all influence scoring via PETSc/SLEPc, with SLURM-array sharding + aggregation |
| [`betweenness/`](betweenness/) | All-to-all + sensory→effector betweenness centrality (Brandes algorithm in `igraph`) |
| [`transforms/`](transforms/) | BANC ↔ JRC2018F / JRCVNC2018F registration; cross-dataset mesh + skeleton bridging; Neuroglancer upload |
| [`meta/`](meta/) | Master metadata assembly (`banc-meta.R`), hemilineage curation, meta fix-ups |
| [`annotations/`](annotations/) | Cell-type curation + one-shot annotation passes (community annotations, sensory classification, KC / PN by connectivity, banc-tracing-* family) |
| [`share/`](share/) | **The publish boundary.** Every output that ships to the public GCS bucket is written or rsynced by a script in here. |
| [`franken/`](franken/) | Frankenbrain composite-dataset assembly + seed labels |
| [`utilities/`](utilities/) | Plotting + format conversion helpers (not part of any production chain) |
| [`legacy/`](legacy/) | Archived superseded scripts — kept for paper-reproducibility cross-references; intentionally undocumented |

## Top-level files

- `banc-startup.R` — paths, libraries, helpers, `banc.keys` private-ID loader. **Sourced by every script** as `source("banc/banc-startup.R")`.
- `banc-functions.R` — shared utilities (`banc_filter_neurons()`, `banctable_query_cached()`, snapshot fallback, etc.).
- `banc-test.R` — ad-hoc smoke tests; not part of any pipeline run.
- `load-keys.R` — reads `data/private/keys.csv` (gitignored) into the global `banc.keys` list.

## Executing scripts

Scripts in `banc/` assume the working directory is the **repo root** (so `source("banc/banc-startup.R")` resolves). On HMS-O2 they are invoked via sbatch wrappers under [`o2/production/`](../o2/production/). Off-O2 most stages will hit "data missing" — only the `--source gcs` subset (alignment-prep, spectral-clustering, betweenness) runs without lab paths.

Per-file headers (in the `#'` roxygen-style block at the top of each script) state: what each script reads, what it writes, which sbatch wrapper invokes it, which BANC-project consumer uses its output, and which paper-Methods section it implements.
