# `banc/annotations/` — cell-type curation + per-task annotation passes

One-shot and recurring curation scripts that derive or fix annotations on BANC neurons: community match imports, KC / PN naming by connectivity, sensory classification, cell-info publishing, plus the per-task tracing-* family split out from the legacy `banc-tracing.R` mega-script. Most write back to SeaTable `banc_meta`; some publish reference CSVs to GCS.

Sits between `banc/update/` (which pushes column updates to SeaTable) and `banc/meta/` (which compiles the final meta table) in the production chain.

## Scripts

| File | Purpose |
|---|---|
| `banc-cell-info.R` | Compile + publish brain hemilineage / NT / neuropeptide reference CSVs to GCS. |
| `banc-community.R` | Build the BANC ↔ MANC per-neuron community lookup (NBLAST + manual matches joined). |
| `banc-connectivity-review-groups.R` | Assemble cross-matched neuron groups (AL, CX, brain motor) into manual-review CSVs. |
| `banc-fix-major-cell-type-errors.R` | Integrate external tracing decisions (ORN + others). Manual spot-fixer. |
| `banc-kc-by-connectivity.R` | Name Kenyon cells by dominant DAN / MBON connectivity (`KCab` / `KCa'b'` / `KCg` / `KC`). |
| `banc-pn-by-connectivity.R` | Name projection neurons by glomerular upstream connectivity. |
| `banc-sensory-classify.R` | Classify proofread sensory neurons by NBLAST against named GT (sets `NBLAST:<type>`). |
| `banc-sensory-jo-orn.R` | Find missing Johnston's-organ sensory neurons + apply ORN cell-type corrections. |
| `banc-tracing-cns-network.R` | CNS-cluster / CNS-network curation + snapshot restore. |
| `banc-tracing-matches.R` | Cross-dataset match updates + Hampel-lab tracing integration. |
| `banc-tracing-png-tags.R` | Apply tracing-issue tags driven by PNG review folder contents. |
| `banc-tracing-status.R` | Status updates + roughly-proofread sensory marking. |
| `banc-vnc-retyping.R` | Cluster + retype female-specific efferent VNC neurons (anchored by typed neurons). |
| `franken-annotations-fix.R` | Standardise franken-brain annotations against published FAFB Supplemental table. |

## Executable for users?

No — every script reads / writes SeaTable + reads HMS-O2 paths. Invoke from the repo root via `Rscript banc/annotations/<file>.R` only with credentials configured.

See [`banc/README.md`](../README.md), [`banc/meta/`](../meta/), [`annotations/`](../../annotations/) (reference artefacts the cell-info scripts consume).
