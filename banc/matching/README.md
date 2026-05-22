# `banc/matching/` — cross-validation + correction of NBLAST / connectivity matches

Post-hoc validation and correction of cell-type matches: detects mismatches between Codex connectivity alignment scores and SeaTable `fafb_cell_type`, audits `fafb_match` quality via NBLAST top-10, propagates classification cascades, and handles VNC sexually-dimorphic / abdominal-neuromere typing. Produces review CSVs with neuroglancer links + comparison PNGs; on re-run with reviewed CSVs, pushes accepted changes back to SeaTable.

Sits between [`banc/nblast/`](../nblast/) (NBLAST results) + [`alignment/`](../../alignment/) (connectivity alignment scores) on the input side, and [`banc/update/`](../update/) (the SeaTable-push side) on the output side.

## Scripts

| File | Purpose |
|---|---|
| `banc-al-connectivity-matches.R` | Detect BANC↔FAFB type mismatches for AL neurons via Codex scores; emits review CSV + neuroglancer links. |
| `banc-al-review-images.R` | Render three-mesh comparison PNGs for AL mismatches (BANC blue / SeaTable-type red / connectivity-matched green). |
| `banc-al-type-changes.R` | AL local + projection-neuron type-change detection. |
| `banc-cx-connectivity-matches.R` | Same as AL but for central-complex neurons. |
| `banc-cx-fix.R` | Re-pick best FAFB + hemibrain matches per CX neuron from NBLAST top-10. |
| `banc-cx-review-images.R` | CX comparison PNGs. |
| `banc-vnc-type-changes.R` | Detect VNC cell_type changes vs the new VNC typology (intrinsic + effector). |
| `banc-vnc-type-review-images.R` | VNC type-change PNGs via MANC NBLAST + MANC→BANC registration. |
| `banc-vnc-sexually-dimorphic.R` | Propagate `sexually_dimorphic` labels across BANC, MANC, FAFB, maleCNS. |
| `banc-bm-asymmetry.R` | Diagnose BM sensory L/R asymmetry; propose retypes from downstream connectivity. |
| `banc-class-fixes.R` | Propagate downstream classification cascade for bates-sourced neurons (FAFB primary, maleCNS fallback). |
| `banc-fafb-cell-type-fixes.R` | Audit `fafb_match` quality via NBLAST top-10; auto-resolve same-type swaps; flag the rest. |
| `banc-update-abdominal-neuromere-types.R` | Push FAbG-project abdominal-neuromere types into SeaTable. |

## Executable for users?

No — needs SeaTable + CAVE credentials and HMS-O2 paths.

See [`banc/README.md`](../README.md), [`alignment/README.md`](../../alignment/README.md), and [`docs/png-matching.md`](../../docs/png-matching.md).
