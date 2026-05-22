# `o2/` — SLURM batch scripts

This directory holds every `sbatch` job script used by `bancpipeline` on HMS-O2. The pipeline is HMS-O2-specific; this layout is internal infrastructure rather than user-facing API. Paper readers do not need any of these — see the top-level [`README.md`](../README.md) instead.

## Layout

```
o2/
├── o2_env.sh              # module loads + env vars; sourced by every script
├── production/            # load-bearing recurring chain (cron / orchestrators / paper artefacts)
├── alignment/             # cross-dataset cell-type alignment runs (whole-brain + optic-lobe)
└── oneshots/              # diagnostics, sweeps, smoketests, one-time rebuilds
```

Counts (post-reorg 2026-05-21): production 40, alignment 11, oneshots 15.

## How scripts find the env

Every job sources the module environment with:

```bash
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
```

**Use the absolute path.** The earlier reorg used `source "$(dirname "$0")/../o2_env.sh"` on the theory that `$0` would resolve to the script's repo path. It does not — SLURM copies the submitted script to `/var/spool/slurmd/job<JID>/slurm_script` before executing, so `$(dirname "$0")` becomes `/var/spool/slurmd/job<JID>` and the relative `../o2_env.sh` misses by a mile (caught and fixed across all 65 sbatch scripts on 2026-05-21). The absolute form above is invariant under SLURM's script-spool behaviour.

`Rscript` calls inside each script still use repo-root-relative paths (`Rscript banc/...`), so the leading `cd /home/ab714/bancpipeline` is mandatory.

## Conventions

- Walltime: aim for ~125% of expected runtime.
- Memory: specify with `--mem=XG` (total) or `--mem-per-cpu=XG`. Default 1G is almost always too low.
- Partition: `short` (≤ 12h), `medium` (≤ 5d), `long` (≤ 30d), `priority` (≤ 30d, 2-job cap), `transfer` (for GCS rsync).
- stdout / stderr: `-o file_%j.out` `-e file_%j.err`.

See `o2_env.sh` for the module-load preamble (`gcc/14.2.0`, `R/4.4.2`, `cmake/3.31.2`, `java/jdk-23.0.1`, plus the UDUNITS2 env vars).

## Re-pointing cron entries (post-reorg)

`ab714` has no active crontab on O2 as of 2026-05-21 (verified: `crontab -l` reports "no crontab for ab714"). The recurring chain is launched manually or via the `o2_banc_update.sh` orchestrator rather than cron. If/when cron is reintroduced, the recurring set would be:

```cron
sbatch o2/production/o2_banc_updateids.sh
sbatch o2/production/o2_banc_nblast.sh
sbatch o2/production/o2_banc_metrics.sh
sbatch o2/production/o2_banc_synapses.sh
sbatch o2/production/o2_banc_data_push.sh
sbatch o2/production/o2_banc_update.sh
sbatch o2/production/o2_banc_results.sh
```

## Key recurring chains

- **Daily update**: `o2_banc_updateids.sh` → `o2_banc_nblast.sh` (refresh ID + NBLAST).
- **Full v888 rebuild**: `o2_banc_v888_rebuild.sh` (master orchestrator; ~250 GB / priority partition).
- **Influence**: `o2_banc_influence.sh` → `o2_banc_aggregate_influence.sh` (SLURM-array shard then aggregate).
- **Alignment**: `alignment/o2_banc_wb_prep_v888v2.sh` → `alignment/o2_banc_wb_production_v888v2.sh` → `alignment/o2_banc_wb_push.sh`.

## Deleted in Phase B (2026-05-21)

Six superseded variants removed to reduce ambiguity:

- `o2_banc_data_push_v2only.sh` — `data_push.sh` already orchestrates both v2 and v3.
- `o2_banc_ids888.sh`, `o2_banc_ids890.sh` — one-off ID migrations now subsumed by `o2_banc_updateids.sh`.
- `o2_banc_synapses_v3_optimised.sh` — `o2_banc_synapses_v3.sh` is the canonical path.
- `o2_banc_wb_cosine_v850.sh` — preprint-version sanity check, paper uses v888.
- `o2_banc_wb_production.sh` — older driver; `o2_banc_wb_production_v888v2.sh` is the paper-config production driver.

History for these files is in the pre-release `wilson-lab/bancpipeline` archive.
