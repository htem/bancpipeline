# Optic Lobe Cross-Dataset Cell Type Matching: BANC ↔ FAFB

## Context

BANC optic lobe intrinsic neurons (~60k per side) lack cell_type annotations. FAFB/FlyWire has established cell_type labels for the same neuron types. We need to transfer FAFB type labels to BANC neurons. NBLAST alone is weak for optic lobe (many fine, spatially repeated types) and cross-dataset graph alignment tools like NTAC only work within a single connectome. We'll build a custom iterative type-level connectivity matching algorithm that uses NBLAST for initialization and refines via connectivity profiles, processing one side at a time.

## Algorithm: Iterative Type-Connectivity Alignment

### Why not NTAC?
NTAC assigns types within a single connectome graph — it needs a single connected adjacency matrix. BANC and FAFB are separate graphs with no cross-dataset synapses. NTAC can't leverage FAFB's connectivity structure to type BANC neurons.

### Our approach
Cross-dataset alignment via iterative refinement of type-level connectivity profiles. The shared vocabulary is cell types: if a BANC neuron has the same pattern of connections to the same types as FAFB neurons of type X, it's probably type X.

### Algorithm

```
INPUTS:
  - BANC adjacency matrix (60k x 60k sparse, normalized weights)
  - FAFB adjacency matrix (60k x 60k sparse, normalized weights)
  - FAFB type labels (ground truth, fixed throughout)
  - BANC seeds: neurons with known cell_type (fixed anchors)
  - NBLAST scores: BANC neuron -> FAFB neuron similarity (medium sparse)

INITIALIZE BANC type assignments:
  Tier 1: BANC neurons with cell_type in SeaTable -> keep as fixed anchors
  Tier 2: BANC neurons with NBLAST scores -> assign type of best-scoring
           FAFB match (soft, fully overridable)
  Tier 3: BANC neurons with neither -> unassigned (typed by connectivity
           once neighbours are typed, after first iterations)

COMPUTE FAFB type centroids (done once, these are fixed):
  For each FAFB type T:
    input_centroid[T]  = mean over FAFB neurons of type T of their
                         input-by-type profile
    output_centroid[T] = mean over FAFB neurons of type T of their
                         output-by-type profile

ITERATE until convergence:
  1. Build type-connectivity profile for each BANC neuron:
     For neuron n with neighbours {m1, m2, ...}:
       input_profile[t]  = sum of normalized input weights from
                           neighbours assigned to type t
       output_profile[t] = sum of normalized output weights to
                           neighbours assigned to type t
     -> Each BANC neuron gets a sparse vector over ~800 FAFB types

  2. Score each BANC neuron against each FAFB type centroid:
     score(n, T) = weighted_jaccard(profile(n), centroid(T))
     where weighted_jaccard(A, B) = sum(min(a_i, b_i)) / sum(max(a_i, b_i))

  3. Reassign non-anchor BANC neurons to best-scoring FAFB type
     (Tier 1 anchors stay fixed)

  4. Check convergence: if < 1% of neurons changed type, stop

OUTPUT:
  Per BANC neuron: assigned_cell_type, best_score, runner_up_type,
  runner_up_score, best_fafb_match (specific FAFB neuron of assigned type
  with highest NBLAST score, where available)
```

### Scalability

The bottleneck is Step 2: scoring 60k BANC neurons against ~800 FAFB type centroids. Each neuron's profile and each centroid are sparse vectors over ~800 types. This is a single sparse matrix multiply:

```
profiles (60k x 800) @ centroids.T (800 x 800) -> scores (60k x 800)
```

This takes seconds in scipy. Step 1 (building profiles) is also a sparse matrix multiply:

```
assignment_matrix (60k x 800, one-hot) -> used to aggregate adjacency by type
profiles = adjacency (60k x 60k) @ assignment_matrix (60k x 800) -> (60k x 800)
```

Total per iteration: ~seconds. With ~10-20 iterations to convergence: under a minute. **No GPU needed.**

### How NBLAST is used

NBLAST is used **only for initialization** (Tier 2 neurons). It provides a starting point so that the first iteration has informative type-connectivity profiles. After that, connectivity drives all reassignments — a neuron with a strong NBLAST match to type X can be reassigned to type Y if connectivity evidence is stronger. NBLAST scores are medium sparse (not every BANC neuron has them), which is fine — Tier 3 neurons get typed purely from connectivity context after the first few iterations stabilize.

For the final `fafb_match` column (specific FAFB neuron ID), NBLAST is used post-hoc: among all FAFB neurons of the assigned type, pick the one with the highest NBLAST score to the BANC neuron (where available).

## Data Sources

### BANC
- **Metadata**: `banctable_query()` -- root_id, supervoxel_id, side, region, super_class, cell_type, top_nt
- **Edgelist**: `banc_850_edgelist_simple.feather`
  - O2: `/n/data1/hms/neurobio/wilson/connectomes/banc/banc_850/banc_850_edgelist_simple.feather`
  - GCS: `gs://brain-and-nerve-cord_exports/processed_data/banc/banc_850/banc_850_edgelist_simple.feather`
- **NBLAST scores**: `banc_fafb_783_nblast.feather` at `{banc.meta.save.path}/`

### FAFB
- **Metadata**: `franken_meta()` from SeaTable -- fafb_id, cell_type, side, super_class, top_nt
- **Edgelist**: `fafb_783_simple_edgelist.feather`
  - O2 (fafbpipeline output): `/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/`
  - O2 (banc connectivity dir): `/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_edgelist.feather`
  - GCS: `gs://brain-and-nerve-cord_exports/processed_data/fafb/fafb_783_simple_edgelist.feather`
- **Local fafbpipeline clone**: `/Users/abates/projects/flyconnectome/fafbpipeline/`

### Data loading pattern
Follow `banc-connectivity-review-groups.R` lines 11-91: use `gcs_cache()` helper when running locally (downloads from GCS to temp cache), direct O2 paths when on cluster.

```r
on_o2 <- nzchar(Sys.getenv("SLURM_JOB_ID")) || grepl("o2.rc.hms.harvard.edu", Sys.info()["nodename"])
if (on_o2) {
  banc.el <- arrow::read_feather("/n/data1/hms/neurobio/wilson/connectomes/banc/banc_850/banc_850_edgelist_simple.feather")
  fafb.el <- arrow::read_feather("/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_edgelist.feather")
} else {
  banc.el <- arrow::read_feather(gcs_cache("gs://brain-and-nerve-cord_exports/processed_data/banc/banc_850/banc_850_edgelist_simple.feather"))
  fafb.el <- arrow::read_feather(gcs_cache("gs://brain-and-nerve-cord_exports/processed_data/fafb/fafb_783_simple_edgelist.feather"))
}
```

## Neuron Pool Filters

**Include** (both datasets):
- `region == "optic_lobe"` OR `region == "optic_lobes"`
- `super_class %in% c("visual_centrifugal", "visual_projection", "optic_lobe_intrinsic")`

**Exclude** (both datasets):
- `super_class %in% c("glia", "not_a_neuron", "trachea")`
- Lamina types: `cell_type %in% c("R1-6", "R7", "R8", "L1", "L2", "L3", "L4", "L5", "Lai")` or `cell_class` containing "photoreceptor" or "lamina_intrinsic"

**Side filter**: Process right side first (`side == "right"`), then repeat for left.

## Pipeline Overview

### Files to create

| File | Language | Purpose |
|------|----------|---------|
| `alignment/banc-alignment-prep.R --region optic-lobe` | R | Data preparation (Phase 0) |
| `alignment/banc-alignment-run.py` | Python | Iterative connectivity alignment (Phase 1) |
| `alignment/banc-alignment-validate.R` | R | Validation and mismatch analysis (Phase 2) |
| `alignment/banc-optic-lobe-update.R` | R | SeaTable update (Phase 3) |
| `o2_banc_optic_match.sh` | Bash | SLURM orchestration |

---

## Phase 0: Data Preparation (R)

### Step 0.1: Build optic lobe neuron sets
- BANC: filter `banctable_query()` by region/super_class/side, exclude glia/lamina
- FAFB: filter `franken_meta()` identically
- Include `top_nt` / `neurotransmitter_predicted` columns for validation
- Output: `data/optic_lobe/banc_optic_{side}_meta.csv`, `data/optic_lobe/fafb_optic_{side}_meta.csv`

### Step 0.2: Extract edgelists
- Load full edgelists (GCS when local, O2 path on cluster)
- Filter to edges where **at least one** endpoint is in the optic neuron pool (keeps cross-region connectivity signal)
- Apply `count >= 2` threshold
- Output: `data/optic_lobe/banc_optic_{side}_edgelist.feather`, `data/optic_lobe/fafb_optic_{side}_edgelist.feather`

### Step 0.3: Extract NBLAST scores
- Filter `banc_fafb_783_nblast.feather` to rows where BANC neuron is in the optic pool
- Keep all FAFB matches (not just optic — a BANC optic neuron's best NBLAST hit might be an FAFB optic neuron not in our filtered set)
- Medium sparse: many BANC optic neurons will have no NBLAST scores
- Output: `data/optic_lobe/banc_fafb_optic_{side}_nblast.csv` (banc_id, fafb_id, nblast_score)

### Step 0.4: Identify seed labels
- **Tier 1 (anchors)**: BANC neurons with existing `cell_type` in SeaTable — split 80% anchor / 20% holdout
- **Tier 2 (NBLAST-initialized)**: BANC neurons with NBLAST scores but no cell_type — assign type of best FAFB match
- **Tier 3 (unassigned)**: remaining neurons — start with no type, typed purely from connectivity
- Output: `data/optic_lobe/banc_optic_{side}_seeds.csv` (root_850, cell_type, tier, is_holdout)

---

## Phase 1: Iterative Connectivity Alignment (Python)

### Step 1.1: Environment setup
```bash
pip install scipy numpy pandas pyarrow
```
No GPU needed. No external alignment package — this is our own algorithm.

### Step 1.2: Load data and build sparse matrices
- Load feather edgelists and metadata CSVs via pyarrow/pandas
- Build BANC and FAFB adjacency matrices as CSR (scipy.sparse)
- Build FAFB type assignment matrix (neuron x type, one-hot)
- Memory: 60k x 60k sparse ≈ 50-200MB

### Step 1.3: Compute FAFB type centroids
- `fafb_profiles = fafb_adj @ fafb_type_matrix` (60k x T sparse)
- For each type T: `centroid[T] = mean(profiles of FAFB neurons of type T)` (T x T)
- Compute separate input and output centroids
- These are fixed throughout — FAFB types are ground truth

### Step 1.4: Initialize BANC assignments
- Tier 1: from seeds CSV (fixed anchors)
- Tier 2: from NBLAST — for each neuron, type of the best-scoring FAFB match
- Tier 3: unassigned (null type)

### Step 1.5: Iterate
```python
for iteration in range(max_iter):
    # Build BANC type assignment matrix from current assignments
    banc_type_matrix = build_type_matrix(banc_assignments, type_index)

    # Type-connectivity profiles: adjacency @ type_matrix
    banc_input_profiles  = banc_adj.T @ banc_type_matrix   # (60k x T)
    banc_output_profiles = banc_adj   @ banc_type_matrix   # (60k x T)

    # Score each BANC neuron against each FAFB type centroid
    # Using weighted Jaccard on concatenated [input, output] profiles
    scores = weighted_jaccard_matrix(
        banc_profiles, fafb_centroids)  # (60k x T)

    # Reassign (skip Tier 1 anchors)
    new_assignments = types[scores.argmax(axis=1)]
    new_assignments[anchor_mask] = anchor_types  # keep anchors fixed

    # Convergence check
    changed = (new_assignments != banc_assignments).sum()
    pct = changed / len(banc_assignments) * 100
    print(f"Iteration {iteration}: {changed} changed ({pct:.1f}%)")
    banc_assignments = new_assignments
    if pct < 1.0:
        break
```

### Step 1.6: Output
- Per BANC neuron: root_850, assigned_cell_type, best_score, runner_up_type, runner_up_score, confidence (best - runner_up), tier
- For neurons with NBLAST data: best_fafb_match = FAFB neuron of assigned type with highest NBLAST score
- Output: `data/optic_lobe/banc_optic_{side}_alignment_results.csv`

---

## Phase 2: Validation (R)

### Step 2.1: Holdout accuracy
- Compare assignments against 20% holdout Tier 1 seeds (these were hidden from the algorithm)
- Report overall accuracy, per-type accuracy, confusion matrix
- Compare against NBLAST-only baseline (what type would NBLAST alone assign?)

### Step 2.2: Neurotransmitter consistency check
- For each BANC neuron with a predicted NT (`top_nt`), check that its assigned FAFB type has a compatible NT
- Override: any neuron with `cell_class == "photoreceptor"` or `cell_type %in% c("R1-6", "R7", "R8")` gets NT = "histamine" regardless of prediction (histamine is exclusive to photoreceptors in Drosophila)
- Mismatches are strong evidence of a wrong assignment (e.g., glutamatergic vs cholinergic is a hard constraint)
- Flag NT-inconsistent assignments as low-confidence
- Independent validation signal: NT prediction is from synapse ultrastructure, not connectivity

### Step 2.3: Mismatch detection
- Follow `banc-al-connectivity-matches.R` pattern
- Join results to SeaTable, identify conflicts with existing annotations
- Generate neuroglancer URLs for manual review
- Output: `data/codex/optic_lobe_connectivity_mismatches.csv`

---

## Phase 3: SeaTable Update (R)

- Follow `banc-al-connectivity-matches.R` update pattern
- High-confidence assignments (high score, NT-consistent) auto-accepted
- Mismatches and low-confidence flagged for review
- Updates: `fafb_match`, `fafb_cell_type`, `cell_type` columns

---

## SLURM Specifications

### Data prep (R)
```bash
#SBATCH -c 1 -t 0-02:00 -p short --mem=64G
```

### Alignment (Python) -- CPU only
```bash
#SBATCH -c 4 -t 0-01:00 -p short --mem-per-cpu=16G
```

### Validation (R)
```bash
#SBATCH -c 1 -t 0-04:00 -p short --mem=32G
```

No GPU needed. The iterative algorithm is sparse matrix multiplies — seconds per iteration, <1 minute total.

## Scalability Notes

- 60k x 60k sparse adjacency: ~50-200MB RAM
- 60k x 800 type profiles: ~10MB dense, much less sparse
- 800 x 800 centroids: trivial
- Per-iteration cost: two sparse matrix multiplies + weighted Jaccard scoring = seconds
- 10-20 iterations to convergence: under 1 minute total
- If needed, partition by neuropil (medulla ~30k, lobula ~15k, lobula plate ~15k) — but shouldn't be necessary

## Key References

### Code templates in this repo
- `banc/annotations/banc-connectivity-review-groups.R` -- GCS/O2 data loading, edgelist filtering, metadata joining
- `banc/matching/banc-al-connectivity-matches.R` -- mismatch detection, NGL URLs, SeaTable update
- `banc/clustering/banc-spectral-clustering.py` -- Python-on-O2 with sparse matrices, argparse
- `banc/banc-startup.R` -- all path definitions
- `fafb/fafb-meta.R` -- FAFB optic lobe metadata

### External
- FAFB data pipeline: `/Users/abates/projects/flyconnectome/fafbpipeline/`

## Verification

1. Run prep script, check pool sizes (~60k per side per dataset)
2. Run alignment on right side, check holdout accuracy (target: >80%)
3. Check NT consistency rate (target: >90% of typed neurons have matching NT)
4. Spot-check 20 assignments via neuroglancer links
5. Compare against NBLAST-only baseline
6. If satisfactory, run left side and push to SeaTable
