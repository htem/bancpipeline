# `banc/legacy/` — archived superseded scripts

Scripts that have been replaced by current-pipeline equivalents but are retained for paper-reproducibility cross-references. Intentionally undocumented at the per-script level. Do not invoke for new work — every file here has a current-pipeline replacement.

Common replacements:

| Legacy file | Replaced by |
|---|---|
| `banc-build-influence.R`, `banc-influence.R`, `banc-all-by-all-influence.R` | [`banc/influence/banc-build-influence.R`](../influence/) |
| `banc-meta.R`, `banc-meta-integrate.R` | [`banc/meta/banc-meta.R`](../meta/), [`banc/share/banc-data.R`](../share/) |
| `banc-data.R`, `banc-build.R`, `banc-sjcabs.R` | [`banc/share/banc-data.R`](../share/), [`banc/share/banc-sjcabs.R`](../share/) |
| `banc-nblast-share.R` | [`banc/share/banc-nblast-share-gcs.R`](../share/), [`banc/nblast/banc-nblast-compile.R`](../nblast/) |
| `banc-split.R`, `banc-skel.R`, `banc-roots.R`, `banc-regions.R`, `banc-ntpred.R` | [`banc/metrics/banc-calculate-*.R`](../metrics/) |
| `banc-an-dn-connectivity.R`, `banc-synapse-proportion-plot.R`, `banc-assess-synapses.R` | exploratory; superseded by `banc/metrics/banc-calculate-completion.R` + analysis in [`BANC-project`](https://github.com/htem/BANC-project) |
| `fix-batch-003-duplicates.R` | one-off completed migration |

## Executable for users?

No — these reference old paths and BANC versions (v626 / v746 / v821 / v850) and may not run against the current data layout. They are kept solely so the pre-publication lineage of paper-era code remains accessible.

See [`banc/README.md`](../README.md). For working equivalents, follow the table above.
