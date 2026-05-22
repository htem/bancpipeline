#!/bin/bash
#SBATCH -J banc_ngl_upload
#SBATCH -c 4
#SBATCH -t 0-12:00                              # short partition cap
#SBATCH -p short
#SBATCH --mem=16G
#SBATCH --array=0-19                            # 20 disjoint shards
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_ngl_upload_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_ngl_upload_%A_%a.err

# Parallel upload of BANC v888 neuron meshes to
# gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/.
#
# Each array task takes a disjoint slice of the proofread+roughly_proofread
# v888 id list, keyed by (last 6 digits of root_888) %% SLURM_ARRAY_TASK_COUNT
# == SLURM_ARRAY_TASK_ID — see banc/transforms/banc-ngl-upload.R for the
# shard logic. Disjoint slices guarantee no two workers ever touch the same
# GCS object or local .obj scratch file, so no flock is needed.
#
# Re-run: just `sbatch o2_banc_ngl_upload_array.sh`. Each task's per-id
# `gsutil stat` pre-flight skips ids already on neuron_meshes/, so re-runs
# auto-pick-up where the prior batch left off.
#
# Resize: change `--array=0-(N-1)` and the script re-shards on next launch.

set -uo pipefail

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

mkdir -p jobs

# Export so the R script can read both task index AND total task count
# (sbatch sets SLURM_ARRAY_TASK_ID and SLURM_ARRAY_TASK_COUNT in the task env).
export SLURM_ARRAY_TASK_ID
export SLURM_ARRAY_TASK_COUNT

echo "=== banc-ngl-upload shard ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} @ $(date) ==="
echo "host: $(hostname)  job: ${SLURM_JOB_ID}  array: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

Rscript banc/transforms/banc-ngl-upload.R

echo "=== shard ${SLURM_ARRAY_TASK_ID} done @ $(date) ==="
