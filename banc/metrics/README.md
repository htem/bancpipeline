# `banc/metrics/` — per-neuron computation: skeletons, regions, synapses, splits, NT

The compute-heavy core of the pipeline. For each proofread BANC neuron, computes (or fetches) L2 skeletons, OBJ meshes, nuclei + root positions, region of innervation, per-neuron volume, axon / dendrite flow-centrality splits, the versioned synapse parquet + simple edgelist, per-neuron NT predictions, and the cross-version synapse → neuropil lookup. Every output is consumed downstream by [`banc/update/`](../update/), [`banc/meta/`](../meta/), and ultimately [`banc/share/banc-data.R`](../share/).

Sits between [`banc/update/banc-ids.R`](../update/) (defines which neurons exist) and [`banc/meta/banc-meta.R`](../meta/) + [`banc/share/banc-data.R`](../share/) (which assemble the release artefacts).

## Scripts

| File | Purpose |
|---|---|
| `banc-l2.R` | Acquire + reroot per-neuron L2 skeletons via `pcg_skel`. |
| `banc-obj.R` | Download per-neuron OBJ meshes. |
| `banc-nuclei.R` | CAVE nuclei → nm coords + L/R side + containing neuropil. |
| `banc-calculate-skeletons.R` | Skeletonise OBJ meshes into radius-bearing SWCs via `fafbseg::skeletor()`. |
| `banc-calculate-root-positions.R` | Compute root / soma position per neuron (nucleus preferred). |
| `banc-calculate-regions.R` | Assign region of innervation (`ventral_nerve_cord` / `central_brain` / `optic_lobe` / `brain` / `rind`). |
| `banc-calculate-volumes.R` | Per-neuron volume via the CAVE L2 cache. |
| `banc-calculate-l2-metrics.R` | Merge per-neuron L2 + split metrics into a single feather. |
| `banc-calculate-connectivity.R` | Build the versioned BANC synapse parquet + simple edgelist (v2 / v3). |
| `banc-calculate-split.R` | Flow-centrality axon/dendrite splits per proofread neuron. |
| `banc-calculate-synapses.R` | Per-neuron synapse / mitochondria metrics, preferring split-CSV compartment labels. |
| `banc-calculate-neuropil-inclusion.R` | Assign neuropil / region / side to each synapse via stable `id` lookup + alpha-shape residuals. |
| `banc-calculate-ntpred.R` | Per-neuron NT predictions from the `synister_banc` per-synapse classifier. |
| `banc-calculate-completion.R` | Synapse-capture rates across v1, v2, v3 detection rounds. |
| `banc-extract-synapse-lookups.R` | Extract version-agnostic synapse → neuropil lookup (run before a rebuild to preserve assignments). |
| `banc-synapses-v3.R` | Process v3 synapse predictions: download → region/neuropil → svid→root_id → capture vs v2. |
| `banc-synapses-v3-optimised.R` | Bbox-accelerated v3 classifier; cuts point-in-surface tests by ~3 orders of magnitude per neuropil. |
| `banc-v3-synapse-sample.R` | Regenerate the v3 spatial sample CSV (standalone of completion script's v3 block). |

## Executable for users?

No — needs HMS-O2 paths and CAVE + SeaTable credentials. Output feathers / parquets are released via GCS / Dataverse — fetch those rather than re-running.

See [`banc/README.md`](../README.md) and [`docs/data-products.md`](../../docs/data-products.md).
