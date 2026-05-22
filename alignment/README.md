# `alignment/` — cross-dataset cell-type alignment by annealed morphology + connectivity

Iterative algorithm that transfers cell-type labels from a fully-annotated **target** connectome into a partially-annotated **query** connectome by **annealing between morphology and connectivity**. In the BANC paper, the query is BANC and the target is FAFB/FlyWire (Methods §"Automated typing by morphology and connectivity"). The algorithm is target-dataset-agnostic — any densely-reconstructed reference with a published cell-type vocabulary can be slotted in (MANC, hemibrain, maleCNS, future connectomes).

This module is the flagship algorithmic contribution of `bancpipeline`. It is region-agnostic and intended to be reused for other cross-dataset alignment problems: different region cuts, different reference connectomes, or any pair of densely-reconstructed graphs that share a cell-type vocabulary.

> **Terminology used throughout this doc and the code:**
> - **query dataset** = the connectome whose neurons receive types (here: BANC).
> - **target dataset** = the reference connectome that provides the type vocabulary, the seed labels and the connectivity centroids (here: FAFB).
>
> Some legacy variable names in the code still use `banc.*` / `fafb.*` literals. Generalising the codebase so an arbitrary target dataset can be slotted in is tracked in `upgrade_plan.md` (Phase D.6).

---

## The problem

Cell-type identification across connectomes has historically pulled in two directions, neither of which finishes the job:

- **Morphology (NBLAST + 1:1 matching)** — fast, works well for distinctive types, but *fails systematically* where multiple cell types share arbor geometry and differ only in their synaptic partners. This is pervasive among columnar optic-lobe types (Mi, Tm, T4/T5 variants) and many central-brain interneurons.
- **Within-dataset label propagation (e.g. NTAC, Schwartzman et al. 2025)** — uses connectivity directly, but propagates labels *within a single connectome's graph*. It needs many seed examples per type and struggles when the reference graph is only partially labelled. It does not transfer across datasets.

Human annotators resolve the ambiguous cases by combining both — morphology first (narrows the candidate set), then connectivity (decides among look-alikes). No published algorithm did this end-to-end across datasets at the time of the BANC paper. This module is that algorithm.

---

## The approach: morphology + connectivity, annealed

The core idea is a **scheduled blend** of two scoring channels per iteration:

```
score = α · connectivity_score  +  (1 − α) · NBLAST_score
```

where α is annealed from balanced to connectivity-dominant over 80 iterations:

```
α:    0.05 ──────────────────────────────────► 0.95
      ┊ NBLAST anchors assignment              ┊ connectivity refines
      ┊ (early — high-confidence morphology)   ┊ (late — circuit context wins)
iter:  0                                      80
```

A second temperature parameter τ anneals the softmax sharpness in parallel (soft early, hard commitments late). The schedule is a simulated-annealing / label-propagation hybrid (cf. Zhu & Ghahramani 2002): early iterations explore under a strong morphological anchor, late iterations commit on connectivity to disambiguate look-alike types.

### What each channel does

- **Morphology (NBLAST)** — supplies initial type *hypotheses* per neuron (Tier-2 seeds from the best FAFB match) and constrains the search space (low-NBLAST candidates pruned). It anchors the algorithm against runaway connectivity-only drift.
- **Connectivity** — at each iteration, BANC neurons' input + output profiles are aggregated to FAFB cell-type centroids (mean input/output, 1-hop and 2-hop) and scored against fixed FAFB type centroids via a **cosine-seeded weighted-Jaccard ensemble** metric. The 2-hop term, with full weight, is key — it captures circuit context that single-hop misses.

### Capacity-constrained greedy assignment

Soft scores are converted to hard assignments via a two-pass greedy procedure that respects **per-type capacity limits derived from FAFB left-right count variability** (with a +10% slack). First pass: neurons with usable NBLAST data are matched to *individual* FAFB neurons (each FAFB neuron claimable by at most one BANC neuron). Second pass: remaining neurons are assigned to *types* by soft probability ranking, capped at capacity.

This constraint is what stops a popular type from absorbing everything it scores plausibly against — a failure mode for naive nearest-centroid methods. It also pre-bakes the biological prior that BANC and FAFB should have roughly the same number of each type.

### Curated guardrails

- **1,125 known-incorrect NBLAST matches** are excluded from the assignment pool (curated list maintained in `banc-alignment-false-positives.R`).
- A **soma-presence rule** prevents a neuron with a detected nucleus from being assigned to a type whose somata lie outside the CNS volume (e.g. photoreceptors, lamina neurons).
- **Early stopping**: best iteration by holdout accuracy is retained.

---

## Why this combination is unusual

The morphology+connectivity pairing is rare in connectomics tooling. Most published alternatives sit firmly in one camp:

| Approach | Channel | Cross-dataset? | Limitation this hits |
|---|---|---|---|
| NBLAST + manual review | morphology only | yes | misses connectivity-distinguished types |
| NTAC | connectivity only | no — within-graph | needs dense seeds; degrades on partial graphs |
| Brute-force connectivity centroid matching | connectivity only | yes in principle | no morphology prior → high collapse risk |
| Manual expert review | both, mentally | yes | doesn't scale to ~10⁵ neurons |
| **This algorithm** | **both, annealed** | **yes** | combines anchors and refinements; needs paired NBLAST + edgelist + seed set |

The schedule itself is the contribution: starting with morphology-balanced scoring gives the algorithm a sane anchor; ramping toward connectivity-dominant lets it refine the columnar / interneuron ambiguities that morphology can't resolve.

---

## Headline result (paper, BANC right optic lobe)

| Method | Holdout accuracy | Mi1 anchor type accuracy |
|---|---|---|
| NBLAST-only baseline | 43.9 % | — |
| NTAC (within-dataset, ≤9% seeds) | 65.0 % | — |
| **This algorithm** | **82.7 %** | **90.4 %** |

Spearman correlation between BANC and FAFB per-type neuron counts: **0.938**. Full numbers, parameter sweeps and ablations: [`experiment-log.md`](experiment-log.md).

---

## Pipeline

```
       prep                              run                   validate / plots
   ┌──────────┐                      ┌──────────┐               ┌─────────────┐
   │ presets/ │ ── inputs ──►        │ alignment │ ── CSV ──►   │  validate    │
   │ <region>/│                      │   run.py  │              │   plots      │
   │ prep.R   │                      │           │              │ discrepancies│
   └──────────┘                      └──────────┘               └─────────────┘
                                                                       │
                                                                       ▼
                                                                update-seatable
```

## Quickstart

```bash
# Optic lobe (right) — the paper's primary demo
Rscript alignment/banc-alignment-prep.R --region optic-lobe right --source gcs --syn-source v2
python alignment/banc-alignment-run.py \
  --side right --metric ensemble --max-iter 80 \
  --hop2-weight 1.0 --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --ensemble-blend 0.3 --manual-labels Mi1 --nt-weight 0
Rscript alignment/banc-alignment-validate.R right

# Whole brain (both hemispheres) — the production paper run
Rscript alignment/banc-alignment-prep.R --region whole-brain both --source gcs --syn-source v2
python alignment/banc-alignment-run.py \
  --side both --bilateral \
  --data-dir data/whole_brain_alignment --file-prefix brain \
  --metric ensemble --max-iter 80 \
  --hop2-weight 1.0 --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --ensemble-blend 0.3 --manual-labels Mi1 --nt-weight 0
Rscript alignment/banc-alignment-validate.R both
```

## Scripts

Top-level scripts in `alignment/` are region-agnostic. Region-specific configuration lives in `presets/<region>/`.

| Script | Purpose |
|---|---|
| `banc-alignment-prep.R --region <NAME>` | Dispatch to `presets/<NAME>/prep.R`. Loads BANC + FAFB metadata, filters to the region's neuron pool, builds the seed-label table, computes per-type capacity. |
| `banc-alignment-run.py` | The algorithm: iterative NBLAST + connectivity ensemble with α / τ annealing, capacity-constrained greedy assignment. Region-agnostic. |
| `banc-alignment-validate.R` | Holdout accuracy, NT consistency, mismatch detection against existing SeaTable annotations. |
| `banc-alignment-plots.R` | Validation plots (per-type accuracy, confusions, NT consistency). |
| `banc-alignment-discrepancies.{R,py}` | Super_class / region discrepancy CSV with Neuroglancer links for manual review. |
| `banc-alignment-update-seatable.R --region <NAME>` | Dispatch to `presets/<NAME>/update-seatable.R`. Pushes accepted alignment results to SeaTable. |
| `banc-alignment-sweep.py` | Two-stage hyperparameter sweep (metric, α schedule, NBLAST threshold). |
| `banc-alignment-sweep-eval.R` | Aggregate and plot sweep results. |
| `banc-alignment-ntac.py` | NTAC baseline (single-connectome label propagation). |
| `banc-alignment-diff.R` | Row-level diff between two alignment-result CSVs (e.g. dry-run vs. live). |
| `banc-alignment-fix-seeds.R` | One-shot: re-seed cell types after a SeaTable curation pass. |
| `banc-alignment-fill-cell-type.R` | Post-process: assign cell types where the algorithm produced no confident match. |
| `banc-alignment-fill-super-class.R` | Post-process: backfill super_class consistent with assigned cell_type. |
| `banc-alignment-false-positives{,-append}.R` | Maintain the curated blacklist of known-incorrect NBLAST matches. |
| `alignment-data-sources.R` | Shared helper: resolve BANC + FAFB data paths from `--source local|gcs`. |
| `alignment_splits.py` | Shared helper: side / data-split utilities. |
| `assessment/` | Quantitative benchmark harness comparing this algorithm to NBLAST-only and NTAC. |

## Presets

| Preset | Scope | Output dir | Notes |
|---|---|---|---|
| `optic-lobe` | BANC optic lobe + visual projection / centrifugal; FAFB equivalent. | `data/optic_lobe/` | Paper's primary demo. Intrinsic-typed neurons held out for validation. |
| `whole-brain` | BANC central brain + optic lobe + ANs/DNs (excludes glia, lamina). | `data/whole_brain_alignment/` (or `$BANC_WB_OUTPUT_DIR`). | Production paper run; wider FAFB-vocabulary check on BANC seed labels. |

### Adding a new preset

1. Pick a region name (lowercase-with-dashes), e.g. `central-complex`.
2. Create `presets/central-complex/prep.R` — copy an existing preset and adapt: filter spec (which `region` / `super_class` / `cell_type` values), output dir, filename infix, seed-label rules (which super_classes are anchors vs. holdout), per-type capacity formula.
3. Optionally create `presets/central-complex/update-seatable.R` for the SeaTable-push side.
4. Run: `Rscript alignment/banc-alignment-prep.R --region central-complex ...`.

Presets are sandboxed — they source `banc/banc-startup.R` and `alignment/alignment-data-sources.R` and otherwise stand alone.

## Algorithm parameters (paper run)

| Parameter | Value | Notes |
|---|---|---|
| Metric | `ensemble` | Cosine seeds → weighted-Jaccard re-rank of top-30 |
| `hop2_weight` | `1.0` | Equal weight on 1-hop and 2-hop type-connectivity |
| `nblast_threshold` | `0.15` | Below-threshold matches dropped from Tier-2 init |
| `alpha_start → alpha_end` | `0.05 → 0.95` | NBLAST ↔ connectivity anneal over 80 iter |
| `tau_start` | `4.0` | Softmax temperature anneal (sharp early → hard commitments late) |
| `ensemble_blend` | `0.3` | Cosine weight in the ensemble metric (Jaccard 0.7) |
| `manual_labels` | `Mi1` | Manually-anchored example types |
| `nt_weight` | `0` | NT ignored — connectivity encodes NT implicitly. NT-aware variants underperformed (`experiment-log.md`). |
| Per-type capacity | FAFB L-R count + 10 % slack | Greedy assignment respects these caps |

Full ablations and rationale: [`experiment-log.md`](experiment-log.md).

## Output schema

Each run writes to `data/<region>/`. Key outputs:

- `*_alignment_results.csv` — per-neuron: assigned cell_type, confidence, runner-up type / score, best FAFB NBLAST match.
- `*_holdout_accuracy.csv` — held-out neurons' predicted vs. true type, with per-type accuracy.
- `*_nt_mismatches.csv` — predicted NT vs. consensus NT of the assigned FAFB type.
- `*_mismatches.csv` — assignments that disagree with existing SeaTable `fafb_cell_type`.

`banc-alignment-update-seatable.R --region <NAME>` reads `*_alignment_results.csv` and pushes accepted assignments to SeaTable columns (`fafb_match`, `cell_type`, `fafb_cell_type`, `super_class`).

## Further reading

- Paper Methods §"Automated typing by morphology and connectivity" (Bates, Phelps, Kim, Yang et al., bioRxiv 10.1101/2025.07.31.667571).
- [`method.txt`](method.txt) — methods-section excerpt, with citations to Zhu & Ghahramani (2002), Costa et al. (2016), Schwartzman et al. (2025), Dorkenwald et al. (2024), Schlegel et al. (2024).
- [`experiment-log.md`](experiment-log.md) — full experiment history: every parameter sweep, every ablation, what worked and what didn't.
- [`optic-lobe-matching-plan.md`](optic-lobe-matching-plan.md) — original design doc.
