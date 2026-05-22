# bancpipeline code-comment style guide

`bancpipeline` is a collection of scripts, not an R package — we do not export functions, build man pages, or run `R CMD check`. But we want the codebase to read like a documented package. We borrow the [roxygen2](https://roxygen2.r-lib.org/) tag conventions because they're familiar to R users from the natverse, and we apply the equivalent Google-style docstring conventions to Python files.

This document is the contract. New scripts and edits to retained scripts should follow it. Pre-existing scripts will be brought into compliance during the Phase E documentation pass (see `upgrade_plan.md`).

## Principles

1. **Explain why and where, not what.** A reader can read the code to learn what it does; comments exist for context they can't recover (paper-figure cross-reference, hidden constraint, the BANC-project consumer of this output, the reason a default is what it is).
2. **No narrating comments.** `# filter the data` adds nothing. Remove on sight.
3. **One file header at the top, function-level docs above each non-trivial function.** Skip docs for one-line helpers whose name + signature are self-evident.
4. **Cite the paper and the sister repos liberally.** Most files map to a Methods section in Bates, Phelps, Kim, Yang et al. and produce an output consumed by `BANC-project`. State both.
5. **Link, don't duplicate.** If the canonical schema for a data product lives in `BANC-project/manuscript/print/dataverse/documentation/<file>.md`, reference that file by path rather than copying its column list inline.
6. **British English in human-language prose.** Comments, docstrings, README sections, and commit messages use UK spelling — *colour, optimisation, behaviour, analyse, organise, neighbour, recognise, fibre, centre, modelling, normalised*. The lead author writes in UK English and the codebase should match. Code identifiers (variable names, function names, column names) are unaffected: existing API names like `color_mips`, `normalize()`, `cluster_center` keep their incoming spelling because renaming them would break callers and published schemas. So: `# Compute behavioural colour assignments` in a comment, but `assign_color_mips(...)` as the function call. When introducing *new* identifiers in new code, prefer the UK spelling where it doesn't break existing patterns. The paper itself uses American spelling in places (it's published in a US journal) — quoting paper Methods section names verbatim is fine; everything else defaults UK.

## File header — R

Top of every `.R` script that is invoked directly (i.e. via `Rscript` or `source()` from another script):

```r
#' <script-name> — <one-line purpose, imperative tense>
#'
#' <Optional paragraph if the one-liner isn't enough. Cite paper figure
#' or methods section if the script's existence is paper-motivated.>
#'
#' @section Reads:
#'   - `<path or table>` — <one-line description>
#'   - SeaTable `banc_meta`: cols `<col1>`, `<col2>`
#'   - CAVE table `<name>` (materialisation v888)
#'   - env var `<NAME>` — <how it gates behaviour>
#'
#' @section Writes:
#'   - `<path>` — <one-line description>
#'   - SeaTable `banc_meta`: cols `<col1>`
#'
#' @section CLI:
#'   --source {v2,v3}   default v3; which synapse version to read
#'   --dry-run          don't push to SeaTable
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_synapses.sh`; cron daily 01:00
#'
#' @section Used by:
#'   BANC-project/R/startup/banc-edgelist.R (loads via construct_path)
#'   BANC-project/R/figures/panel_an_dn_connectivity.R
#'
#' @section Schema:
#'   BANC-project/manuscript/print/dataverse/documentation/banc_888_edgelist_simple_v2.md
#'
#' @section Paper:
#'   Methods §"Synapse detection", §"Annotation taxonomy".
#'   File cited as `[banc_888_edgelist_simple_v2.feather]`.
#'
#' @section Notes:
#'   - <gotcha worth flagging; bug workarounds; non-obvious defaults>
```

Drop sections that don't apply (e.g. a pure analysis script with no SeaTable side effects can omit `Writes` and `Used by` if it produces nothing downstream).

## File header — Python

Use a module-level docstring with the same sections as a Google-style block. Tags become section headings.

```python
"""<one-line purpose>.

<Optional paragraph.>

Reads:
  - <path or table>: <description>
  - env var BANC_INFLUENCE_SHARD_IDX: SLURM array index

Writes:
  - <path>: <description>

CLI:
  --source {v2,v3}   default v3.

Invoked by:
  o2/production/o2_banc_influence.sh

Used by:
  BANC-project/R/figures/panel_all_to_all_influence.R

Schema:
  BANC-project/manuscript/print/dataverse/documentation/influence_all_to_all.md

Paper:
  Methods §"Influence". File cited as [influence/all_to_all/].

Notes:
  - <gotcha>
"""
```

## Function-level — R

Above each non-trivial function (anything longer than ~10 lines or with non-obvious signature):

```r
#' Short imperative title (one line).
#'
#' Optional longer description if the title isn't enough. Mention any
#' coupling to other scripts or to BANC-project consumers.
#'
#' @param x         <type> — <what it is>
#' @param threshold <type> — <what it controls; default rationale>
#' @return          <type and shape> — <what it means>
#' @seealso         `banc/metrics/banc-calculate-connectivity.R`
filter_valid_neurons <- function(meta, threshold = 5) {
  ...
}
```

Function-level docs **never use** `@export`, `@noRd`, `@rdname`, `@inheritParams`, `@import` / `@importFrom`, `@useDynLib`. Those are package-build mechanics and we are not a package.

`@examples` blocks are allowed but not required. Prefer pointing to a calling script over inlining an example.

## Function-level — Python

Google-style docstring. Same content as R, different syntax.

```python
def filter_valid_neurons(meta: pd.DataFrame, threshold: int = 5) -> pd.DataFrame:
    """Drop glia, debris, and duplicate neurons.

    Args:
        meta: BANC metadata table; must have columns `super_class`, `proofread`.
        threshold: minimum synapse count to keep. Default 5 matches paper.

    Returns:
        Filtered metadata, same columns. Rows reduced.

    See Also:
        `banc/metrics/banc-calculate-connectivity.R` for the consumer.
    """
```

## Inline comments

Default to none. Add an inline comment only when removing it would confuse a careful reader.

Keep:
- `# Paper Fig. 3e uses this threshold; see Methods §"Naming AN/DN clusters".`
- `# nolint start: long line below is a hardcoded CAVE column list — keep aligned`
- `# Wang et al. 2024 use 0.3 here; we use 0.5 because v888 synapses are denser.`
- `# v2 voxel = 16x16x45 nm; v3 = same units but smaller threshold (>=10).`

Delete:
- `# loop over neurons`
- `# read the file`
- `# return result`

## Paper cross-references

When the script implements something the paper describes, cite the Methods section by name. Don't quote line numbers — those drift.

> **Paper:** Methods §"Spectral clustering". File cited as `[banc_888_meta.feather]`.

When the paper cites a specific *bancpipeline* path (e.g. Methods cites `bancpipeline/alignment/`), that path is part of the publication record. Don't rename without thinking.

## BANC-project cross-references

The downstream consumer for most bancpipeline outputs lives in `/Users/papers/BANC-project` (paper analysis + figure code). When you can verify a specific consumer file references a specific bancpipeline output, cite it in the header's `@section Used by:` line. Paths are stable contracts:

```
@section Used by:
  BANC-project/R/startup/banc-meta.R (loads banc_888_meta.feather)
  BANC-project/R/figures/panel_an_dn_connectivity.R
```

Do not fabricate consumers. If you haven't read the consumer file, leave `Used by:` off — better silence than a wrong reference.

## Schema authority

For files that appear in the Dataverse release (`banc_888_*.{feather,parquet,csv}`, `banc_*_nblast.feather`, influence parquets, etc.), the canonical per-column schema lives in `BANC-project/manuscript/print/dataverse/documentation/<filename>.md`. Always link the relevant `.md` from the producer's header. If bancpipeline and that `.md` disagree, *either* can be the bug — flag the discrepancy in the script header `@section Notes:` and resolve it in a separate change.

## Linting / enforcement

Phase F (publish) does not block on this. Phase E (per-file documentation pass) brings the kept-script set into compliance. There is no lint rule for these tags; we rely on review.

A `grep` to find scripts that haven't been touched yet:

```bash
grep -L "^#' .* — " banc/**/*.R       # R scripts missing a header title line
grep -L '^"""' banc/**/*.py            # Python scripts missing a module docstring
```
