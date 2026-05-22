# `exploration/` — historical exploratory analyses

This directory holds **legacy code from exploratory analyses** carried out by members of the Wilson Lab during BANC development. It is a **historical record only**:

- The scripts here are **not part of the published bancpipeline release** and were not used to produce the paper's reported results.
- They are kept so that the lineage of decisions during the project is reproducible — e.g. why certain DN clusters were investigated, how an early connectivity-clustering attempt was structured before the alignment algorithm in `alignment/` took its place.
- They may not run on the current data layout. They reference older path conventions, older BANC versions (v626 / v850), or external GCS / Drive locations that may have moved.

For the paper-run pipeline, see the top-level [`README.md`](../README.md). For the canonical cell-type assignment algorithm, see [`alignment/README.md`](../alignment/README.md).

## Subdirectories

- `R/` — exploratory R scripts (DN / AN connectivity, neck-connective tracing, abdominal-neuromere influence, NBLAST cluster experiments, etc.).
- `python/` — exploratory Python (early cascade-model implementations, spectral-propagation tests, FAFB example analyses, Jupyter notebooks).
- `matlab/` — exploratory MATLAB (heatmap + histogram plots from early DN-influence work).

## Conventions

- Filenames use hyphens (matching the rest of the codebase) for `.R`, `.py`, `.sh`. MATLAB `.m` files keep their original names because MATLAB requires the filename to match the function name.
- Every code file begins with a four-line legacy-record header in the appropriate comment syntax (`###` for R, `#` for Python / shell, `%` for MATLAB).
- Jupyter notebooks (`.ipynb`) are not headered programmatically — they live under `exploration/python/` and inherit the legacy-record status from this README.
