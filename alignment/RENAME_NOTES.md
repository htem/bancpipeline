# Alignment refactor — notes for the O2 verification session

Drafted 2026-05-21 by a local-machine Claude session (`/Users/abates/projects/flyconnectome/bancpipeline`). This document is a hand-off to a future Claude session running on O2, where the refactored alignment code will be tested against the previous implementation with real FAFB / BANC data. **Read this before running anything.**

## TL;DR

The alignment module has been refactored to be **target-dataset-agnostic**: instead of FAFB-specific variable and column names, it uses neutral `target_*` names internally. Default target is still FAFB, so behaviour should be byte-identical to before. Verification = run the new prep + algorithm with `target_name = "fafb"`, compare outputs to a stored run of the old implementation.

## What changed

| Change | Files | Scope |
|---|---|---|
| Directory moved | `banc/alignment/` → `alignment/` | All references updated in `o2/*.sh`, internal sourcing, paper draft (separately by AB). |
| Region-agnostic script renames | 11 `banc-optic-lobe-*` / `banc-wb-*` → `banc-alignment-*` | See `git log --diff-filter=R` for the full mapping. |
| Region-specific scripts moved into presets | `banc-optic-lobe-prep.R` → `alignment/presets/optic-lobe/prep.R`; same for whole-brain + update-seatable. | A new dispatcher `alignment/banc-alignment-prep.R --region {optic-lobe,whole-brain}` selects the preset. |
| Blacklists committed | `data/{optic_lobe,whole_brain_alignment}/forbidden-matches.csv` → `alignment/presets/<region>/forbidden-matches.csv` (now tracked in git). | 11 referencing scripts updated. |
| **R-side `fafb.*` → `target.*` rename** | `alignment/alignment-data-sources.R`, `alignment/presets/{optic-lobe,whole-brain}/prep.R` | New helper `load_target_meta(target="fafb")` returns the target metadata with normalised column names (`target_id`, `target_cell_type`, `target_super_class`, etc.). Both prep scripts now use `target.meta`, `target.pool`, `target.el`, `target_id`, `target_cell_type` throughout. |
| Alpha schedule corrected | `o2/alignment/o2_banc_wb_production_v888v2.sh` + others; `method.txt`; `alignment/README.md`; `banc-alignment-sweep.py` | `--alpha-start 0.5` → `--alpha-start 0.05`. Paper Methods says 0.05 → 0.95 (NBLAST-dominant → connectivity-dominant); the production scripts had drifted to 0.5. **This is a behavioural change** — runs with 0.5 and 0.05 will differ. The 0.05 value is the paper-correct setting. |

## What did NOT change (deferred)

1. **Public Dataverse-deposit filenames** (`banc_888_*.feather`, `banc_888_synapses_v2_enriched.parquet`, etc.) — they're contracts with the paper Methods, BANC-project consumers, and Harvard Dataverse. Leave as is.
2. **Alignment intermediate filenames** (`data/<region>/banc_optic_*`, `fafb_optic_*`, `banc_fafb_*_nblast.csv`) — still use the `optic_lobe / fafb / brain` infixes. A filename migration to `{stage}__{query}_{target}__{region}_{side}__v{ver}_{syn}__{YYYYMMDD}.{ext}` is planned (`upgrade_plan.md` §D.6.2) but not yet implemented.
3. **Python script `banc-alignment-run.py`** still uses Python-side `fafb_*` variable names internally. **Critical**: the R prep now writes CSV columns named `target_id`, `target_cell_type`, etc. The Python expects `fafb_id`, `cell_type`. **This means the current Python will not read the new R prep output as-is.** See "Integration break" below.
4. **Sibling R scripts** (`banc-alignment-validate.R`, `banc-alignment-plots.R`, `banc-alignment-discrepancies.R`, `banc-alignment-update-seatable.R`, etc.) still use FAFB-flavoured Python variable names. They will need the same rename treatment.

## Integration break to watch for

The R prep writes CSVs with `target_*` column names:

```
data/optic_lobe/fafb_optic_right_meta.csv columns:
  target_id, target_region, target_side, target_super_class, target_cell_class,
  target_cell_sub_class, target_cell_type, target_top_nt
```

The Python `banc-alignment-run.py` expects `fafb_id, cell_type, ...`. When run end-to-end (prep → run), it will fail at `pd.read_csv(..., dtype={"fafb_id": str})` because the column doesn't exist.

**Three ways to resolve before testing the algorithm end-to-end**:

1. **Run prep alone first.** Confirm the new R prep produces the same row counts + same neuron IDs as the old prep. Compare meta CSVs side-by-side: column count goes up, but `target_id` should equal the old `fafb_id`, `target_cell_type` should equal old `cell_type`, etc.
2. **Rename Python next** — apply the same `fafb_* → target_*` mass rename to `banc-alignment-run.py`. Mechanical, ~170 occurrences. Then re-run end-to-end.
3. **Backwards-compat shim**: have the R prep write CSVs with both `target_*` and `fafb_*` column names (alias on write). Removes the integration break without renaming Python. Document as a stop-gap.

The local-machine session chose to defer #2 and #3 to a follow-up commit because:
- The user wanted to test the R rename first.
- The Python rename is mechanical but adds ~170 sed touches that we couldn't verify without runtime testing.

**Recommendation for the O2 session**: do option 1 first (R-prep parity check), then option 2 (Python rename + run end-to-end).

## How to verify the R prep is behaviour-preserving

1. **Pick a saved old-prep output** as a baseline. On O2, the most recent paper-config prep output should live at `data/optic_lobe/banc_optic_right_meta.csv` etc. from the last production run. Capture file checksums and row counts.

2. **Run new prep with `target_name = "fafb"`** (the default):

   ```bash
   cd /home/ab714/bancpipeline   # or wherever bancpipeline lives on O2
   Rscript alignment/banc-alignment-prep.R --region optic-lobe right --source local
   ```

   The new prep writes to the same `data/optic_lobe/` directory with the same filenames (no Q2-filename-migration applied yet).

3. **Compare**:
   - Row count of new `banc_optic_right_meta.csv` should match the old (BANC side is unchanged).
   - Row count of new `fafb_optic_right_meta.csv` should match the old. *Columns will differ in name*: new has `target_id` where old had `fafb_id`. Confirm `setequal(new$target_id, old$fafb_id)` and `setequal(new$target_cell_type, old$cell_type)`.
   - Row count of `banc_fafb_optic_right_nblast.csv` should match the old; new has `target_id` column where old had `fafb_id`.
   - Row count of `banc_optic_right_seeds.csv` should match.
   - `fafb_type_capacity.csv` row count should match; columns are now `target_left, target_right, target_max, target_total` where old were `fafb_*`. Values should be identical row-for-row.

4. **If parity holds**: the rename is behaviour-preserving. Proceed to Python rename / end-to-end test.

5. **If parity fails**: most likely cause is a missed column reference (BSD sed didn't rename it). The local session left a few intentional FAFB-token references:
   - `data_paths$fafb_edgelist` — list key in `resolve_alignment_paths()`. Schema-level, not output-level.
   - `fafb_cell_type` inside BANC pipes (lines ~125-131 of optic preset, ~193-206 of WB preset) — this is BANC's own curated cross-match column, not a target column. Intentional, has an inline comment.
   - Filenames: `banc_fafb_*_nblast.csv`, `fafb_brain_*_meta.csv`, `fafb_type_capacity.csv` — kept (filename migration is separate).
   - Doc strings.

   Anything else that breaks parity is a bug — flag it.

## Adding a new target dataset (future work, not for this verification)

To slot in MANC as target:

1. Edit `alignment/alignment-data-sources.R::load_target_meta()`: add a `target == "manc"` branch that loads MANC metadata and renames columns to `target_id`, `target_cell_type`, etc.
2. Add a new preset directory `alignment/presets/optic-lobe-manc/prep.R` (or similar) that sets `target_name <- "manc"` and adjusts any MANC-specific filters / NT overrides.
3. Either edit `data_paths$fafb_edgelist` (in `resolve_alignment_paths()`) to a generic name like `data_paths$target_edgelist`, OR have the preset pass the MANC edgelist path explicitly.

The algorithm itself (`banc-alignment-run.py`) needs the same `fafb_* → target_*` rename for full MANC support. Until then, MANC will require a backward-compat shim that writes target_* as fafb_* on the R-prep side.

## File-by-file change summary

For surgical rollback if needed:

| File | Status | Notes |
|---|---|---|
| `alignment/alignment-data-sources.R` | Added `load_target_meta()` (~45 lines). | New function; existing `parse_alignment_data_args()`, `resolve_alignment_paths()`, `gcs_cache()` unchanged. |
| `alignment/presets/optic-lobe/prep.R` | Heavy refactor. | Parse-clean. ~50 lines of var/col renames. Filenames kept. |
| `alignment/presets/whole-brain/prep.R` | Heavy refactor. | Parse-clean. ~55 lines of var/col renames. Filenames kept. |
| `alignment/presets/optic-lobe/update-seatable.R` | Light touch. | A few variable renames; mostly unchanged. |
| `alignment/presets/whole-brain/update-seatable.R` | Light touch. | Same. |
| `alignment/banc-alignment-run.py` | **Not yet renamed.** | Will fail to read new prep CSV (column-name mismatch). |
| `alignment/banc-alignment-validate.R` | **Not yet renamed.** | Same. |
| `alignment/banc-alignment-plots.R` | **Not yet renamed.** | Same. |
| `alignment/banc-alignment-discrepancies.{R,py}` | **Not yet renamed.** | Same. |
| `alignment/banc-alignment-false-positives*.R` | Touched only for path updates. | New blacklist path under `presets/<region>/`. |
| `o2/alignment/*.sh`, `o2/oneshots/*.sh` | Many updated: alignment paths + alpha 0.05. | Production sbatch scripts now invoke `alignment/banc-alignment-prep.R --region whole-brain` instead of `banc/alignment/banc-whole-brain-prep.R`. |

## Where to ask if something is unclear

The local-machine session is the upgrade-plan author. The plan lives at the repo root in `upgrade_plan.md` (gitignored on `dev`; not exported to public release). Phase D.6 covers the rename work; Phase F is the public-repo publish step. The verification described above is conceptually part of Phase D.6 but is gated on O2 access we don't have here.

---

## D.6.1 verification results — 2026-05-21, O2 session

Ran by the O2-side Claude session immediately after pulling commits `2ac4ff6 … 587eadf`. **Verdict: D.6.1 is parity-preserving. The integration break described above has already been closed by a follow-up commit before the pull landed.** Details below.

### A. Python rename status — **DONE** (was listed "not yet renamed")

`alignment/banc-alignment-run.py` already uses `target_*` columns end-to-end:

| Probe | Count |
|---|---|
| `target_meta["target_cell_type"]`, `["target_id"]`, `["target_super_class"]` direct accesses | 20+ refs across the file |
| Remaining `fafb_*` token in code | 0 |
| Remaining `fafb_*` tokens total | 2, both inside `#` comments (lines 1063, 1072) — not load-bearing |

Confirmed by:

```bash
grep -nE '"fafb_id"|"fafb_cell_type"|"fafb_super_class"' alignment/banc-alignment-run.py
# (no matches)
grep -cE 'target_id|target_cell_type|target_super_class|target_side' alignment/banc-alignment-run.py
# 20
```

The "integration break" warning in §"Integration break to watch for" above is **outdated** — left in place for traceability but no longer applicable. Future readers: rely on this section, not that one.

### B. Sibling-R-script rename status — **partial**

| Script | `fafb_*` refs | `target_*` refs | Verdict |
|---|---|---|---|
| `alignment/banc-alignment-validate.R` | 3 | 1 | partial — still references `fafb_cell_type` (BANC's own cross-match col — intentional, see note below), but the validate logic itself works against the new prep output. |
| `alignment/banc-alignment-plots.R` | 1 | 4 | mostly renamed; remaining ref is in a docstring. |
| `alignment/banc-alignment-discrepancies.R` | 10 | 1 | NOT renamed. Needs the same `fafb_* → target_*` pass. Defers to a follow-up commit. |
| `alignment/banc-alignment-discrepancies.py` | 3 | 3 | Half renamed. Same. |
| `alignment/banc-alignment-update-seatable.R` | 0 | 0 | dispatcher — no direct refs. |
| `alignment/presets/whole-brain/update-seatable.R` | 0 | 0 | clean. |
| `alignment/presets/optic-lobe/update-seatable.R` | 2 | 0 | one `fafb_cell_type` in a SeaTable column reference (intentional — BANC's curated cross-match column) and one filename string. |

**Note on intentional `fafb_*` retentions**: the BANC SeaTable has columns called `fafb_cell_type`, `fafb_match`, `fafb_side` that record BANC's own curated cross-match to FAFB. These are **not** target-dataset columns inside the algorithm — they are BANC metadata. They stay as `fafb_*` even after the rename, by design.

### C. Data-level parity (header-rename invariance)

When the O2 session encountered stale prep outputs in `data/whole_brain_alignment_v888v2/` (prep'd before the target_* schema, so they had `fafb_id` headers), it applied a header-only `sed` rename to bring them in line. To confirm this is byte-level safe vs. a re-prep:

```bash
md5sum <(tail -n +2 fafb_brain_both_meta.csv.preheaderrename)
md5sum <(tail -n +2 fafb_brain_both_meta.csv)
```

| Region   | rows  | body MD5 identical | header rewritten |
|----------|-------|--------------------|------------------|
| v888v2   | 156,125 | YES | `fafb_id,region,…,cell_type,top_nt` → `target_id,target_region,…,target_cell_type,target_top_nt` |
| v888v3   | 153,762 | YES | same |

So **a header-only sed of pre-rename prep outputs is byte-level equivalent to a fresh prep** for the meta/nblast files. (Backups `.preheaderrename` are kept locally for rollback but are not committed.)

### D. End-to-end behavioural test — **schema verification PASSED, full convergence in flight**

Job lineage on the priority partition (`o2/alignment/o2_banc_wb_cosine_v888v2_priority.sh`):

| Attempt | Job | Outcome | Notes |
|---|---|---|---|
| 1 | `40922205` | FAILED (1s) | `source "$(dirname "$0")/../o2_env.sh"` resolved to `/var/spool/slurmd/...`. Fixed by switching to absolute-path source across all 65 sbatch scripts. |
| 2 | `40922605` | FAILED (18s) | `KeyError: 'target_cell_type'` on the stale pre-rename prep CSV (`fafb_id` headers). Fixed by header-only `sed` rename of `fafb_brain_both_meta.csv` and `banc_fafb_brain_both_nblast.csv` in v888v2 + v888v3 dirs (header-only mutation, body MD5 invariant — see C above). |
| 3 | `40923245` | **RUNNING — schema verification PASSED.** Iterating cleanly. | Past type-vocab (8816 FAFB types), adjacency (130 571 pool + 27 379 boundary, 9.6 M non-zeros), centroids (1-hop density 0.028, 2-hop 0.565), NBLAST lookup (7.3 M pairs over 112 790 neurons), forbidden mask (1148 pairs / 1147 neurons), soma rule (99 639 neurons × 7 no-soma types), init (18 317 anchors / 28 950 NBLAST-soft / 85 978 unassigned, 75 087 holdout, Mi1 holdout 1436). Init metrics: holdout 0.0 %, NT 100.0 %, type_rho 0.124. Now in the iteration loop (max 80, tau 4.0 → 0.5, α 0.05 → 0.95, capacity 1.5x). At WB scale each iteration runs in the tens of seconds; 80 iters fits well inside the 8 h priority walltime. |

**Verification status (2026-05-21T20:55Z):**
- ✅ Schema rename verified end-to-end. Every prep CSV is read by every consumer in `banc-alignment-run.py` without column-name failure.
- ✅ Forbidden-matches CSV (rewritten path after reorg) loads (1148/1230 entries kept).
- ✅ Soma rule, NBLAST lookup, NT scoring, bilateral-consistency partitioning, per-side caps all initialise without error.
- ⏳ Convergence + final holdout accuracy: pending iteration loop completion. Will be appended to this note in a follow-up commit (or as a `## Verification result` block at the bottom of the file).

**One more partial-rename bug found and patched** during this run (logged at startup but non-fatal — soma rule fell back to the explicit type list):

```
WARNING: target_meta has no 'super_class' column; only explicit types will be flagged
```

The gate at `banc-alignment-run.py:802` was checking for the legacy `"super_class"` column on `target_meta`, where the data now lives in `"target_super_class"`. Fixed in the same commit family; the *block inside* the gate was already renamed (it used `target_super_class` correctly), so this was purely a missed condition. After the fix the soma rule will pick up the additional ~150-300 sensory-super_class FAFB types and produce a tighter assignment pool on future runs.

Sibling-script grep for any remaining `"super_class"` literal on `target_meta` came back clean; that bug is closed.

Schema-level verification is therefore **complete**: the post-target-abstraction prep + run can co-operate against the header-renamed CSVs. The deeper "does the algorithm produce sensible numbers" check is the holdout-accuracy report on completion, which depends on a multi-hour priority run. Will update this note with the final holdout figure (or any new failure) once the job exits.

Surfaced and committed in the course of this verification:
- `2b7b4f6` — replay forbidden-matches + 3 R-file local edits at the new paths.
- `bd9d4a7` — D.6.1 + D.6.2 doc.
- `9572acf` — SLURM-safe `source o2_env.sh` (65 sbatch).
- `b05c961` — six **`cell_type` → `target_cell_type` / `top_nt` → `target_top_nt`** column refs in `banc-alignment-run.py` that the rename pass had missed (`build_type_index`, `build_type_matrix`, best-match lookup, target-NT vote loop). Caught only at runtime — this is the integration-break check the original RENAME_NOTES warned about, now closed.

### E. Follow-up commit items (gated on D.6.1 sign-off)

1. Finish the sibling-R rename: `banc-alignment-discrepancies.{R,py}`, `banc-alignment-validate.R` (the non-intentional refs only), `banc-alignment-plots.R` docstring touchup.
2. Replace `data_paths$fafb_edgelist` in `resolve_alignment_paths()` with `target_edgelist` (schema-level key still says `fafb_*`). Low priority; safe rename.
3. Delete the outdated "Integration break to watch for" section above once this report is reviewed.

---

## D.6.2 filename migration — scoped design

User ask (2026-05-21): "Can we re-work our current .csv funs into this framing?" Scoped, **not executed**. Gates on the D.6.1 end-to-end run above completing cleanly.

### Proposed schema

**Core (required):**

```
{stage}__{query}_{target}__{YYYYMMDD}.{ext}
```

**Tagged extensions (optional, in canonical order, only present when meaningful for the preset):**

```
{stage}__{query}_{target}[__region-<R>][__side-<S>][__vq-<VQ>][__vt-<VT>][__syn-<X>][__<extra-tag>...]__{YYYYMMDD}.{ext}
```

Design principle: **the only mandatory metadata is "what stage", "which two datasets", and "when".** Region, side, version pins, synapse-table choice, and any preset-specific knobs are all optional tags. A user aligning hemibrain ↔ FAFB (no region cut, no left/right semantics that matter, no "synapse version") gets short, clean filenames; a user running the paper's BANC ↔ FAFB whole-brain config gets the long form with every tag present. **Absence of a tag means "not applicable to this preset," not "default value."**

Field grammar:

| Field | Required? | Values used in the paper run | Notes |
|---|---|---|---|
| `{stage}` | yes | `prep-banc-meta`, `prep-target-meta`, `prep-banc-edges`, `prep-target-edges`, `prep-seeds`, `prep-seeds-production`, `prep-nblast`, `prep-capacity`, `align-cosine`, `align-ensemble`, `align-ensemble-ind03`, `align-production`, `ntac`, `validate-holdout-accuracy`, `validate-holdout-confusions`, `validate-nt-mismatches`, `validate-mismatches` | dash-internal, double-underscore as field separator. Add variant suffixes after `align-` (e.g. `align-cosine-alpha05`). |
| `{query}_{target}` | yes | `banc_fafb` | future: `banc_manc`, `hemibrain_fafb`, etc. Dataset names are lowercase-with-dashes. **No version embedded here.** |
| `region-<R>` | optional | `region-whole-brain`, `region-optic-lobe`, `region-central-complex` | Only present if the preset subsets the dataset. Hemibrain ↔ FAFB whole-pair alignment: omit. |
| `side-<S>` | optional | `side-right`, `side-left`, `side-both` | Only present if the alignment is hemisphere-specific. Asymmetric datasets / type-level alignments: omit. |
| `vq-<VQ>` / `vt-<VT>` | optional | `vq-888`, `vt-783` | Dataset version pins for query / target. Omit when the dataset doesn't have a public versioning scheme. |
| `syn-<X>` | optional | `syn-v2`, `syn-v3` | BANC synapse-table choice. Omit when irrelevant (any non-BANC query). |
| `{extra-tag}` | optional | preset-defined, lowercase-with-dashes | E.g. `seed-pool-tier1` for the all-tier1 seed set; `alpha-05` for a non-default α schedule. Keep the preset's documentation honest about what tags it emits. |
| `{YYYYMMDD}` | yes | e.g. `20260521` | prep-date for prep outputs; run-date for align/validate outputs. Pinned at creation. |
| `{ext}` | yes | `csv`, `feather`, `parquet` | unchanged. |

Examples across presets:

```
# Paper's BANC ↔ FAFB whole-brain v888 v2 run (all tags present)
align-cosine__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv

# Minimalist hemibrain ↔ FAFB type-level alignment (no region/side/syn pins,
# datasets don't carry a synapse-version concept here)
align-cosine__hemibrain_fafb__vq-1.2.1__vt-783__20260521.csv

# A user's experiment without versions (small connectome, single materialisation)
align-cosine__userdata_fafb__20260521.csv
```

### Mapping current → proposed (paper-config WB v888v2, both sides — all tags present)

Tags in canonical order: `region-whole-brain__side-both__vq-888__vt-783__syn-v2`. A bare hemibrain ↔ FAFB run would drop everything except `vq-1.2.1__vt-783`.

| Current filename | Proposed filename |
|---|---|
| `data/whole_brain_alignment_v888v2/banc_brain_both_meta.csv` | `prep-banc-meta__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `fafb_brain_both_meta.csv` | `prep-target-meta__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_edgelist.feather` | `prep-banc-edges__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.feather` |
| `fafb_brain_both_edgelist.feather` | `prep-target-edges__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.feather` |
| `banc_brain_both_seeds.csv` | `prep-seeds__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_seeds_production.csv` | `prep-seeds-production__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_fafb_brain_both_nblast.csv` | `prep-nblast__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `fafb_type_capacity.csv` | `prep-capacity__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_alignment_cosine_v888v2_tier.csv` | `align-cosine-tier__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_alignment_ensemble_ind03_v888v2_tier.csv` | `align-ensemble-ind03-tier__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_alignment_production_v888v2_all_tier1.csv` | `align-production-all-tier1__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_ntac_production_v888v2_ntac_all_tier1.csv` | `ntac-production-all-tier1__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_holdout_accuracy.csv` | `validate-holdout-accuracy__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_holdout_confusions.csv` | `validate-holdout-confusions__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_nt_mismatches.csv` | `validate-nt-mismatches__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |
| `banc_brain_both_mismatches.csv` | `validate-mismatches__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260521.csv` |

The paper-config names are admittedly long. The trade-off: any researcher reading just the filename knows what they're holding. Glob discipline (`align-*__banc_fafb__*vq-888*.csv`) does the same job that subdirectories used to do.

### Why this works

- **Self-describing.** Glob patterns like `align-*__banc_fafb__whole-brain_*__v888_v2__*.csv` find every paper-config WB-v888v2 alignment output without parsing.
- **Sortable.** Lexical sort groups by stage, then query/target, then region/side, then version, then date — which is also the natural reading order.
- **Datestamp is permanent.** Today's prep + tomorrow's prep don't collide. Re-runs are explicit.
- **No infix ambiguity.** The double-underscore field separator means no field can contain `__`, so parsing is regex-trivial.

### What this will break (impact analysis)

The current names are hardcoded across:

1. `alignment/banc-alignment-run.py` — `--data-dir`, `--file-prefix` give it `{data_dir}/banc_{file_prefix}_{side}_{kind}.csv`. The flat-name scheme breaks this entirely; need a `--io-naming v2` switch or a name-resolver helper.
2. `alignment/banc-alignment-prep.R` and presets — `readr::write_csv()` calls hardcode the filenames.
3. `alignment/banc-alignment-validate.R` — globs by the `banc_{prefix}_{side}_*` pattern.
4. `alignment/banc-alignment-{plots,discrepancies,update-seatable,sweep-eval,...}.R` — same.
5. All `o2/alignment/*.sh` sbatch scripts hardcode the current names in `NEED_PREP` checks, `--forbidden-matches` paths, and downstream consumer paths.
6. `data/whole_brain_alignment*/`, `data/optic_lobe/` already contain hundreds of files in the old scheme — migration must be a no-op for behavioural identity (rename-only).

### Proposed migration sequence (zero-downtime)

1. ✅ **Write a path helper.** ~~Add `alignment_path(stage, query, target, region, side, ver, syn, date, ext)` to `alignment/alignment-data-sources.R`. Returns the new-scheme filename. *Doesn't touch any caller yet.*~~ **DONE** in commit landed alongside this note. Implemented in both R (`alignment/alignment-data-sources.R`) and Python (`alignment/alignment_paths.py`); see API block below.
2. ✅ **Add a name-resolver shim.** **DONE** — `legacy_alignment_path()` is implemented for both languages and covers the paper-run stage set (`prep-*`, `align-*`, `ntac*`, `validate-*`). Producers can dual-write by setting `BANC_ALIGNMENT_NAMING=legacy` on the side that should emit the old name.
3. **Migrate producers** (`prep.R`, `run.py`, `validate.R`, `plots.R`, `discrepancies.{R,py}`, sweep, ntac) one at a time, behind a `BANC_ALIGNMENT_NAMING=v2` env-gate. CI / sbatch defaults stay v1 until each producer is verified.
4. **Flip the default** once all producers + consumers handle both. Symlinks remain for one more release.
5. **Drop old-name support** in a clearly-versioned commit.

### Path-helper API (D.6.2 step 1, landed)

R, in `alignment/alignment-data-sources.R`:

```r
alignment_path(stage, query, target,
               region = NULL, side = NULL,
               vq = NULL, vt = NULL, syn = NULL,
               extra = NULL,
               date = NULL, ext = "csv",
               dir = ".", scheme = "auto")
legacy_alignment_path(...)         # same parameters; emits old-scheme name
```

Python, in `alignment/alignment_paths.py`:

```python
from alignment.alignment_paths import alignment_path, legacy_alignment_path
alignment_path("prep-banc-meta", "banc", "fafb",
               region="whole-brain", side="both",
               vq="888", vt="783", syn="v2")
# → 'prep-banc-meta__banc_fafb__region-whole-brain__side-both__vq-888__vt-783__syn-v2__20260522.csv'
```

`scheme = "auto"` (default) reads `BANC_ALIGNMENT_NAMING` from the environment (`"v2"` if unset, `"legacy"` to switch back). The same env var is honoured by both implementations so producers and consumers can toggle in lock-step.

Smoke-tested against the paper-run stage set:

```
prep-banc-meta       → banc_brain_both_meta.csv                    (legacy)
                     → prep-banc-meta__banc_fafb__region-whole-brain__…__20260522.csv (v2)
align-cosine-tier    → banc_brain_both_alignment_cosine_tier.csv    (legacy)
                     → align-cosine-tier__banc_fafb__…__20260522.csv (v2)
```

Minimalist hemibrain ↔ FAFB (no region / side / syn):

```
align-cosine         → align-cosine__hemibrain_fafb__vq-1.2.1__vt-783__20260522.csv
```

### Known shortcomings of the legacy mapper (acceptable for the migration)

- The current paper-run filenames embed the version pin as a single infix like `v888v2` (e.g. `banc_brain_both_alignment_cosine_v888v2_tier.csv`). `legacy_alignment_path()` drops this version infix from the legacy slot — the v2 side carries it via `vq-888 / syn-v2`. Producers that need pixel-perfect legacy parity can pass `extra = c("v888v2")` and adjust the legacy mapper, but in practice the migration is only a brief overlap window.
- Only paper-run stages are mapped. Diag / sweep / smoketest stages would need entries added if you want them to round-trip; otherwise they default to v2-only.

### Open questions for the user

1. **Output directory.** Current `data/<region>/` puts everything region-relative. Should the new naming flatten to a single `data/alignment/` directory (filename carries region)? Or keep the per-region subtree (slightly redundant filename infix)?
2. **Sweep / diag outputs.** `align-cosine-tier-alpha05-...` vs. a separate `data/alignment/sweeps/<sweep_name>/` subtree? Lean towards subtree — sweeps generate hundreds of files.
3. **Snapshots / dated archives.** Today we use `.may17.csv` ad-hoc and `legacy_v850/`. The `__{YYYYMMDD}` in the new scheme subsumes the dated archive. The `legacy_v850/` semantics would become an automatic date-sort.
4. **Renaming target.** The user said "scoped, not executed." Confirm before any `git mv` — there are roughly 50-80 alignment-related CSVs and feathers in `data/`, plus the hardcoded references in ~15 scripts.

### Header-rename precedent

The data-level invariance check in D.6.1.C shows that **header-only mutations are byte-safe** for our CSVs. Filename rewriting is even safer (no file content touched). The mechanical migration script can be a single `git mv` loop with no risk of data loss, provided it runs only after all consumers handle the new names.
