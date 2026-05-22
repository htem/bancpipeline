# Optic Lobe Alignment: Experiment Log

## Run 1: Baseline (cosine + soft assignments, no capacity constraint)

**Parameters**: tau 2.0->0.1, alpha 0.5->0.95, nblast_threshold 0.1, no capacity constraint
**Changes from original**: Replaced weighted Jaccard with cosine similarity. Added soft assignments with temperature annealing. Added NBLAST hard constraint + adaptive alpha regularizer.

| Metric | Init | Best iter (2) | Iter 15 | Iter 30 |
|--------|------|---------------|---------|---------|
| Holdout | 30.4% | **49.1%** | 35.7% | 3.0% |
| Mi1 | 0.0% | 38.1% | 54.2% | 0.0% |
| NT | 95.8% | 84.4% | 70.6% | 53.8% |
| Type rho | 0.801 | 0.807 | 0.636 | 0.228 |

**Problem**: Algorithm collapses — rare types (LPLC4: 54 FAFB, assigned 9022 BANC; cL16: 2 FAFB, assigned 8229 BANC) vacuum up thousands of neurons. Temperature anneals too fast, NBLAST threshold too permissive (669/669 types allowed per neuron).

---

## Run 2: Capacity constraint + slower annealing + stronger NBLAST

**Parameters**: tau 2.0->0.5 (slower anneal), alpha 0.3->0.8 (more NBLAST weight), nblast_threshold 0.25, capacity 1.5x FAFB count, 5 Sinkhorn rounds
**Changes**: Raised NBLAST threshold to 0.25 (664/669 types allowed vs 669/669). Slowed temperature annealing (tau_end 0.5 vs 0.1). Weighted NBLAST higher (alpha_start 0.3 = 70% NBLAST initially). Added Sinkhorn-like type capacity scaling.

| Metric | Init | Best iter (19) | Iter 25 | Iter 30 |
|--------|------|----------------|---------|---------|
| Holdout | 30.4% | **44.7%** | 43.3% | 35.3% |
| Mi1 | 0.0% | **80.4%** | 87.3% | 87.2% |
| NT | 95.8% | **88.3%** | 87.9% | 71.0% |
| Type rho | 0.801 | 0.208 | 0.213 | 0.212 |

**Mi1 target reached** (80.4% at iter 19, peaks at 88.9% at iter 28).

**Remaining problems**:
1. Type count correlation (rho) terrible (~0.21) — Sinkhorn capacity constraint not biting. TmY5a: 573 FAFB -> 6992 BANC at iter 20. The Sinkhorn modifies soft_probs but argmax still over-assigns. Need hard capacity enforcement.
2. Overall holdout (44.7%) lower than Run 1 (49.1%). The stronger NBLAST weighting helps Mi1 but slightly hurts overall.
3. Degradation after iter 22 — holdout starts dropping, changes increase, NT drops.

**Diagnosis**: The Sinkhorn operates on soft probability mass, but hard assignment (argmax) ignores it. Even if Sinkhorn perfectly balances probabilities, the argmax picks the highest-probability type regardless of how many neurons are already assigned there. Need greedy capacity-constrained assignment at the hard level.

---

## Run 3: Greedy capacity-constrained assignment (breakthrough)

**Parameters**: Same as Run 2 (tau 2.0->0.5, alpha 0.3->0.8, nblast_threshold 0.25, capacity 1.5x)
**Key change**: Replaced Sinkhorn soft scaling with greedy capacity-constrained hard assignment. After computing softmax probabilities, assign neurons to types in order of confidence (highest first), consuming type slots. Once a type's capacity (1.5x FAFB count) is exhausted, remaining neurons fall to their next-best type. Feedback: blend 50% hard-constrained + 50% original soft probs for next iteration's profiles.

| Metric | Init | Best iter (14) | Iter 20 | Iter 30 |
|--------|------|----------------|---------|---------|
| Holdout | 30.4% | **75.7%** | 74.9% | 71.3% |
| Mi1 | 0.0% | **98.5%** | 98.6% | 98.7% |
| NT | 95.8% | **94.8%** | 94.7% | 94.4% |
| Type rho | 0.801 | **0.923** | 0.930 | 0.953 |

**Transformative improvement** across all metrics:
- Mi1: 98.5% (far exceeds 90% target, stable throughout)
- Holdout: 75.7% (up from 44.7% in Run 2, 49.1% in Run 1)
- NT: 94.8% (near-initial levels, stable)
- Type count rho: 0.923 → 0.953 (up from 0.21 — distribution now matches FAFB)
- No degradation collapse — graceful decline from 75.7% to 71.3% over 30 iters
- Still ~9% of neurons changing per iteration (doesn't converge to <1%)

**Remaining issues**:
1. Large columnar types still under-assigned: T4a-d (~300-600 vs ~1400-1500 FAFB), T5a-d (~250-600 vs ~1400-1500), C2 (~140 vs 1477), C3 (~340 vs 1532). These represent the majority of optic lobe neurons and are systematically ~50-80% under-filled.
2. ~9% churn per iteration suggests the algorithm hasn't found a stable equilibrium — connectivity signal alone may not distinguish these highly similar columnar types.
3. The under-assignment of T/C types likely means many neurons that SHOULD be T4/T5/C2/C3 are being assigned to less common types that happen to have higher NBLAST or cosine scores.

**Why it works**: The capacity constraint prevents winner-take-all collapse. Without it, a type with a slightly higher-than-average centroid similarity to many neurons can "claim" thousands of neurons (see TmY5a in Run 2). With capacity, once TmY5a fills its ~860 slots, the remaining neurons must find their actual types. This cascading displacement is what makes the algorithm actually work — it forces specificity.

---

## Run 5: Type frequency prior (prior_weight=0.5)

**Parameters**: Same as Run 3 + prior_weight=0.5
**Change**: Added log-prior from FAFB type frequency, zero-mean so common types get bonus.

| Metric | Best iter |
|--------|-----------|
| Holdout | 30.4% (iter 0 — never improved) |
| Mi1 | 94.2% |

**Failed**: Prior range [-1.56, +5.17] with weight 0.5 = [-0.78, +2.59] added to scores. Much larger than cosine scores [0,1]. Common types completely dominate — all neurons get assigned to T3/T4c regardless of connectivity.

---

## Run 6: Small prior (prior_weight=0.03)

**Parameters**: Same as Run 3 + prior_weight=0.03

| Metric | Best iter (14) |
|--------|----------------|
| Holdout | 74.7% |
| Mi1 | 98.7% |
| NT | 94.7% |
| Type rho | 0.900 |

Negligible effect — essentially same as Run 3 (75.7%). The issue isn't scoring bias.

---

## Run 7: NT hard constraint + higher NBLAST weight

**Parameters**: alpha 0.2->0.6 (80% NBLAST initially, 40% at end), nblast_threshold 0.25, NT hard constraint (neurons with known NT can only be assigned to types with matching consensus NT)

| Metric | Best iter (23) |
|--------|----------------|
| Holdout | **75.9%** |
| Mi1 | **99.1%** |
| NT | **99.6%** |
| Type rho | 0.925 |

**Analysis**:
- NT consistency jumped to 99.6% (from 94.8%) — the constraint nearly eliminates NT mismatches
- Mi1 improved to 99.1% (from 98.5%)
- Holdout marginally better (75.9% vs 75.7%)
- Convergence improved: only 5% churn (vs 9% in Run 3), peaks later (iter 23 vs iter 14)
- But holdout is plateauing around 75-76% — this appears to be the ceiling for type-level connectivity matching

**Diagnosis**: The remaining ~24% errors are likely within-family confusions (T4a↔T4b, T5a↔T5b, etc.) where types have nearly identical type-level connectivity profiles. These can only be resolved by:
1. Neuron-level position information (which column the neuron is in)
2. Neuron-level graph matching (matching to individual FAFB neurons, not type centroids)
3. Coarse-to-fine matching (first assign type family, then refine)

---

## Run 4: Weighted Jaccard A/B comparison (all else identical to Run 3)

**Parameters**: Same as Run 3 but metric=jaccard instead of cosine
**Purpose**: Direct A/B test of similarity metric

| Metric | Cosine (Run 3, iter 14) | Jaccard (Run 4, iter 18, still improving) |
|--------|------------------------|-------------------------------------------|
| Holdout | 75.7% | **79.2%** |
| Mi1 | **98.5%** | 71.4% |
| NT | 94.8% | 93.9% |
| Type rho | 0.923 | 0.916 |
| Speed | 2s/iter | **70s/iter** |

**Key finding**: Jaccard gives better overall holdout (+3.5pp and still climbing) but worse Mi1. Jaccard better distinguishes cross-family types; cosine better resolves within-family subtypes. Jaccard is ~35x slower.

---

## Confusion Matrix Analysis (Run 3, best iteration)

Total holdout: 27,039. Correct: 19,653 (72.7%). Wrong: 7,386.

- **Within-family errors** (e.g., T4a↔T4c): 2,044 (27.7% of errors)
- **Cross-family errors** (e.g., C3→Mi4): 5,342 (72.3% of errors)

If ALL within-family fixed: 80.2%. Top attractor types: T4c (727 stolen from T4a/b/d), T5c (635 from T5a/b/d).
Top cross-family: C3→Mi4 (126), Tm4→Tm9 (88), Tm4→Tm2 (84), Mi15→Dm2 (72).

**Implication**: 90% requires fixing cross-family confusions — types with genuinely similar connectivity profiles. May need neuron-level matching or position features.

---

## NBLAST-only baseline

Best-match: 43.9%, Top-3: 66.1%, Mi1: 27.6%. Our algorithm adds +35pp from connectivity.

---

## Runs 8-12: Neuron-level matching experiments

| Run | Change | Holdout | Mi1 |
|-----|--------|---------|-----|
| 8-9 | Neuron-level NBLAST matching | 78.2% | 93.5% |
| 10 | Annealing feedback blend | 78.1% | 94.1% |
| 11 | 2-hop profiles | 78.6% | 93.4% |
| 12 | Individual FAFB profile matching | 76.6% | 92.2% |

Individual FAFB profiles were WORSE — noisier than centroids. Centroids smooth out variation.

---

## Run 13: Cosine-seeded Jaccard ensemble (new best)

**Approach**: Fast cosine for all 669 centroids, then re-rank top-30 with weighted Jaccard. Blend 50/50. Gets Jaccard's magnitude sensitivity without 70s/iter.

| Metric | Best (iter 40, still climbing) |
|--------|-------------------------------|
| Holdout | **81.3%** |
| Mi1 | **93.4%** |
| NT | 99.6% |
| Type rho | 0.929 |
| Speed | 20s/iter |

New best. Still improving at iter 40. The ensemble captures both angular (cosine) and magnitude (Jaccard) sensitivity.

---

## Parameter Sweep (fresh SeaTable data, per-type capacity)

Re-pulled metadata from SeaTable (2026-03-31). Added per-type capacity from FAFB left-right count variability (capacity = max(L,R) + max(|L-R|, 10% of max)). Mean capacity/max ratio: 1.51x but now correctly per-type.

### Stage 1: One-at-a-time screening (cosine, 15 iters)

| Parameter | Best value | Holdout | vs baseline |
|-----------|-----------|---------|-------------|
| `capacity_scale` | 1.0 | 81.5% | +0.0 (per-type capacity dominates) |
| **`hop2_weight`** | **1.0** | **82.4%** | **+0.9pp** |
| `nblast_threshold` | 0.15 | 81.8% | +0.3pp |
| `tau_start` | 1.0 | 81.6% | +0.1pp |
| `alpha_start` | 0.5 | 81.6% | +0.1pp |
| `alpha_end` | 0.8 | 81.5% | +0.0 |

### Stage 2: Ensemble validation (25 iters, best Stage 1 params)

| Config | Holdout | Mi1 | NT |
|--------|---------|-----|------|
| ensemble (best_s1, top_k=30, blend=0.5) | 84.1% | 93.0% | 99.7% |
| ensemble (top_k=50) | 84.1% | 93.0% | 99.7% |
| **ensemble (blend=0.3 = 70% Jaccard)** | **84.3%** | **93.2%** | **99.7%** |
| ensemble (defaults) | 83.7% | 93.3% | 99.7% |
| ensemble (blend=0.7 = 70% cosine) | 83.8% | 93.3% | 99.7% |
| ensemble (top_k=10) | 83.1% | 92.7% | 99.7% |

### Best config (right side, 60 iters)

Parameters: metric=ensemble, hop2_weight=1.0, nblast_threshold=0.15, tau_start=1.0, alpha_start=0.05, alpha_end=0.8, ensemble_blend=0.3, per-type capacity from FAFB L-R variability.

| Metric | Value |
|--------|-------|
| Holdout | **83.9%** |
| Mi1 | 93.0%+ |
| NT | 99.7% |
| Time | 1359s (23 min) |

Key insights:
- `hop2_weight=1.0` is the single biggest improvement (+0.9pp from screening, much more in ensemble)
- `ensemble_blend=0.3` (Jaccard-heavy) slightly better than 0.5
- Per-type capacity from FAFB L-R variability is more principled than fixed multiplier
- Fresh SeaTable data added ~400 neurons and more Tier 1 labels

---

## Pool correction: remove non-optic neurons, add L1-L5/R7/R8

**Critical bug fix**: The `root_region` filter was pulling in 1,727 non-optic neurons (descending, central_brain_intrinsic, sensory, motor) whose soma happened to be near an optic neuropil. Fix: root_region only includes neurons WITHOUT an existing super_class. Also removed exclusion of L1-L5/R7/R8 (these exist in BANC).

Impact: discrepancies dropped from 6.0% to 0.8% (right) and 30% to 1.1% (left).

### Final results (corrected pool, both sides)

| Metric | Right | Left |
|--------|-------|------|
| BANC pool | 41,558 | 36,909 |
| FAFB pool | 56,719 | 57,271 |
| Holdout accuracy | **82.7%** (25,941 eval) | 55.8% (859 eval, too small) |
| Mi1 accuracy | **90.4%** (855 neurons) | 71.4% (7 neurons) |
| NT consistency | 99.6% | 99.3% |
| Type count rho | 0.938 | 0.932 |
| FAFB match rate | 95.2% | 81.8% |
| Discrepancies | 235 (0.8%) | 42 (1.1%) |

Left side has similar alignment quality (rho 0.932 ≈ right 0.938) despite sparse holdout validation.

---

## Final configuration: re-sweep + manual labels + auto: fix (2026-04-01)

**Additional fixes applied:**
- GT cell_type: use cell_type, fall back to fafb_cell_type, drop "auto:" prefixed types
- R7/R8 FAFB neurons get histamine NT override
- NBLAST refreshed from GCS
- NGL links generated in R via fafbseg (base URL 5370752139264000)
- Discrepancies exclude "auto:" cell types

**Re-sweep results (Stage 1, cosine, 15 iters):**

| Parameter | Best | Holdout |
|-----------|------|---------|
| hop2_weight | 1.0 | 81.7% |
| nblast_threshold | 0.15 | 80.4% |
| tau_start | 4.0 | 80.2% |
| alpha_start | 0.5 | 80.2% |
| alpha_end | 0.95 | 80.1% |
| Combined (cosine) | — | 82.2% |

**Final runs: ensemble + Mi1 manual labels**

Parameters: metric=ensemble, hop2_weight=1.0, nblast_threshold=0.15, tau_start=4.0, alpha_start=0.05, alpha_end=0.95, ensemble_blend=0.3, --manual-labels Mi1, per-type capacity.

| Metric | Right | Left |
|--------|-------|------|
| Pool | 41,576 | 36,924 |
| Holdout (excl Mi1) | **82.7%** (25,086) | 55.4% (856) |
| Mi1 | 855/855 anchored | 7/7 anchored |
| Type count rho | **0.940** | **0.934** |
| FAFB match rate | 94.3% | 81.7% |
| Discrepancies | 129 (0.3%) | 37 (0.1%) |

Discrepancies dropped further after excluding "auto:" cell types from the comparison.

---

## NT weight optimization

Swept nt_weight (0, 0.25, 0.5, 0.75, 1.0) with cosine, 15 iters:
| nt_weight | Holdout |
|-----------|---------|
| **0.0** | **81.4%** |
| 0.25 | 81.3% |
| 0.5 | 81.0% |
| 0.75 | 81.0% |
| 1.0 | 80.8% |

**Finding**: Ignoring NT entirely gives the best accuracy. NT predictions are ~10% wrong, causing more harm than benefit. Connectivity implicitly encodes NT information.

---

## NTAC comparison (within-dataset baseline)

NTAC (Schwartzman et al., 2025, doi:10.1038/s41467-025-68044-1) run on BANC graph alone.

Initial attempt with ~9% seeds (matching our anchor setup): 0.0% — NTAC needs more seeds.

Re-run with 80% random seeds (fair setup for within-dataset propagation):
- **NTAC: 65.0% holdout**, 86.3% Mi1 (30 iterations, 411s)
- **Our method: 83.1% holdout**, Mi1 100% (with different holdout split)
Fair comparison with identical 80/20 random holdout (seed=42, 8,195-8,199 eval neurons):

| Method | Holdout | Mi1 | Uses FAFB? |
|--------|---------|-----|-----------|
| **Our method** | **84.6%** | 84.9% | Yes |
| **NTAC** | **65.0%** | 86.3% | No |

Our method's **+19.6pp advantage** comes from cross-dataset FAFB connectivity centroids. Mi1 is similar (~85%) — NTAC handles distinctive types well. The gap is on less distinctive columnar types where FAFB reference is essential.

---

## Individual FAFB profile scoring (2026-04-02)

Vectorized individual profile scoring: for each BANC neuron with NBLAST data, compare its connectivity profile against individual FAFB neurons' profiles (not just type centroids). Blends centroid and best individual scores: `(1-w)*centroid + w*best_individual`. Top-20 NBLAST matches per neuron to control memory.

### Sweep (no manual labels, 30 iters, right side)

| ind_weight | ind_metric | Holdout | Mi1 | NT | rho |
|-----------|-----------|---------|-----|------|-----|
| 0 (baseline) | — | 82.7% | 78.9% | 93.2% | 0.937 |
| 0.05 | ensemble | 82.9% | 80.7% | 93.4% | 0.939 |
| 0.10 | ensemble | 83.1% | 82.5% | 93.6% | 0.935 |
| 0.15 | ensemble | 83.2% | 82.8% | 93.6% | 0.941 |
| 0.20 | ensemble | 83.3% | 83.9% | 93.8% | 0.932 |
| **0.30** | **ensemble** | **83.4%** | **84.5%** | **93.8%** | **0.935** |
| 0.50 | ensemble | 83.4% | 85.6% | 93.9% | 0.938 |
| 0.10 | cosine | 82.7% | 90.4% | 94.3% | 0.940 |
| 0.20 | cosine | 82.6% | 90.7% | 94.2% | 0.941 |
| 0.30 | cosine | 82.6% | 90.6% | 94.2% | 0.941 |
| 0.10 | jaccard | 82.9% | — | — | — |
| 0.20 | jaccard | 82.5% | — | — | — |
| 0.30 | jaccard | 82.4% | — | — | — |

**Key findings:**
- Ensemble individual scoring best for holdout (+0.7pp), peaks at ind_weight=0.3
- Cosine-only individual: no holdout help but massive Mi1 improvement (+11.5pp)
- Jaccard-only individual: slightly worse than ensemble, degrades at higher weights
- The holdout improvement comes specifically from Jaccard component of ensemble
- Cosine component helps Mi1 by better discriminating individual morphological identity

---

## Final configuration (2026-04-02)

Parameters: metric=ensemble, hop2_weight=1.0, nblast_threshold=0.15, tau_start=4.0, alpha_start=0.05, alpha_end=0.95, ensemble_blend=0.3, nt_weight=0, **ind_weight=0.3**, --manual-labels Mi1, per-type capacity from FAFB L-R variability.

| Metric | Right |
|--------|-------|
| Pool | 41,576 |
| Holdout (excl Mi1) | **83.6%** (25,081) |
| Mi1 | 855/855 anchored |
| NT | 94.3% |
| Type count rho | 0.941 |
| Best iter | 54 (of 60) |

---

## Optic lobe full suite (2026-04-03, fresh SeaTable)

Re-pulled SeaTable data (+202 new annotations). ind_weight=0.5, ensemble, 60 iters.

| Config | Right | Left |
|--------|-------|------|
| No manual labels | **83.5%** (25,936 eval) | 54.4% (913 eval) |
| Mi1 manual labels | 83.2% (25,081 eval) | 54.6% (906 eval) |
| All GT labels | 98.6% typed, rho 0.940 | 82.5% typed |
| Stratified holdout 0.2 | 39.0% (8,413 eval) | 93.1% (6,324 eval) |

Left side accuracy limited by sparse holdout (913 neurons).

---

## Whole-brain bilateral alignment (2026-04-03)

127,281 BANC neurons (both sides), 126,895 FAFB neurons (deduplicated), 7,725 types.
Profile dimension: 30,900 (with hop2) or 15,450 (without).
Holdout: 59,821 neurons (central_brain_intrinsic + optic_lobe_intrinsic).

### Run 1: cosine, no hop2, no bilateral, no conf_gate

| Metric | Value |
|--------|-------|
| Best holdout | **69.3%** (iter ~50) |
| Final holdout (iter 100) | 67.3% |
| Mi1 | 76.4% |
| NT | 88.0% |
| Type rho | 0.974 |
| Time | 30,368s (8.4 hours), ~270s/iter |

Iter 1 breakdown: central_brain_intrinsic **69.0%** (23,705), optic_lobe_intrinsic 51.4% (22,111).
Central brain types are easier to match (more distinctive connectivity).
Algorithm oscillated from iter ~50 onward (21% churn, holdout declining).

### Run 2: cosine + hop2 + bilateral + conf_gate (partial, killed at iter 9)

67.9% at iter 9 (vs 67.7% for Run 1 at same iter). Slightly better trajectory from hop2/bilateral/conf_gate. ~207s/iter. Would need full 100 iters to compare properly.

### Scalability findings

- **Ensemble metric infeasible locally**: Jaccard re-ranking at 30.9k dims takes ~60 min/iter (vs ~4 min for cosine). Needs O2.
- **2-hop profiles OOM'd with non-chunked**: (127k × 30.9k) = 15 GB per array. Fixed by chunked computation (2-hop computed on-the-fly per 5k-neuron chunk).
- **Individual profiles (ind_weight>0) impractical for whole brain**: (127k × 30.9k) = 15 GB extra. Use ind_weight=0.
- **FAFB pool deduplication required**: bilateral mode has duplicate fafb_ids across hemispheres.
