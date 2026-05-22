# `banc/transforms/` — cross-dataset mesh + skeleton registration

For each reference connectome (FAFB, MANC, FANC, hemibrain, maleCNS), transforms its neuron meshes + skeletons into BANC space via the relevant registration chain (e.g. `FlyWire → JRC2018F → BANC` for FAFB, `MANC → JRCVNC2018F → BANC` for MANC). Output is the per-version transformed-mesh / skel tree that feeds [`banc/nblast/`](../nblast/) (which needs all neurons in a common space for NBLAST).

Also handles the Neuroglancer side: uploads precomputed meshes + per-dataset `segment_properties` JSON to the public GCS Neuroglancer layer for browsing at [ng.banc.community](https://ng.banc.community/view).

Sits between the reference-dataset publish dirs ([`fafb/`](../../fafb/), [`manc/`](../../manc/), [`fanc/`](../../fanc/), [`hemibrain/`](../../hemibrain/), [`malecns/`](../../malecns/)) and [`banc/nblast/`](../nblast/).

## Scripts

| File | Purpose |
|---|---|
| `banc-fafb-mesh-transform.R` | Transform FAFB v783 meshes into BANC space (FlyWire → JRC2018F → BANC). |
| `banc-fafb-skel-transform.R` | Transform FAFB v783 skeletons into BANC space. |
| `banc-manc-mesh-transform.R` | Transform MANC meshes into BANC space (MANC → JRCVNC2018F → BANC). |
| `banc-manc-skel-transform.R` | Transform MANC skeletons into BANC space. |
| `banc-fanc-mesh-transform.R` | Transform FANC v1.116 meshes into BANC space (FANC → JRCVNC2018F → BANC). |
| `banc-fanc-skel-transform.R` | Transform FANC v1.116 skeletons into BANC space. |
| `banc-hemibrain-mesh-transform.R` | Transform hemibrain meshes into BANC space (JRCFIB2018F → JRC2018F → BANC). |
| `banc-hemibrain-skel-transform.R` | Transform hemibrain skeletons into BANC space. |
| `banc-malecns-mesh-transform.R` | Transform maleCNS meshes into BANC space via Python `navis` + `flybrains` TPS (no R-native registration yet). |
| `banc-ngl-upload.R` | Upload Neuroglancer-compatible meshes + per-dataset `segment_properties` JSON to GCS. |
| `banc-publish-segment-properties.R` | Refresh per-dataset `segment_properties` JSON only (lighter-weight than the full upload). |

## Executable for users?

No — needs HMS-O2 paths, mesh / skeleton fetch credentials, and gsutil write access for the Neuroglancer publish step.

See [`banc/README.md`](../README.md), top-level [`README.md`](../../README.md).
