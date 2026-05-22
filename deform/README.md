# `deform/` — Deformetrica shape-analysis runs

Cross-dataset shape-statistics experiments using [Deformetrica](https://www.deformetrica.org/). Builds the template + atlas inputs (mesh / skeleton point sets) for Deformetrica's geodesic shape modelling, runs the optimisation, and pulls the resulting VTK shapes back into R for inspection. Not part of the production publishing chain — outputs are research artefacts rather than released data products.

## Scripts

| File | Purpose |
|---|---|
| `deformetrica-fafb-brain.R` | Prepare FAFB brain mesh point sets and drive a Deformetrica brain-template fit. |
| `deformetrica-manc-vnc.R` | Prepare MANC VNC mesh point sets and drive a Deformetrica VNC-template fit. |
| `deformetrica-lr.R` | Build left ↔ right point correspondences for Deformetrica mirror-symmetry analysis. |
| `transfer_vtk.R` | Pull Deformetrica VTK output back into R, link to source neurons. |
| `deform_fafb.sh` | Run Deformetrica from the FAFB conda env (`deformetrica estimate model.xml ...`). |

The `example/`, `mesh/`, `skel/` subdirectories hold input templates / sample data used by the scripts.

## Executable for users?

No — needs HMS-O2 paths, a working Deformetrica conda environment, and lab-internal mesh inputs.

See top-level [`README.md`](../README.md).
