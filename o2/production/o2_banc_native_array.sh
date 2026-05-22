#!/bin/bash
#SBATCH -c 10                          # cores per array task
#SBATCH -t 0-12:00                     # 12h short (each shard fits)
#SBATCH -p short
#SBATCH --mem-per-cpu=8G               # 80G per task
#SBATCH --array=0-9                    # 10 disjoint shards
#SBATCH -o jobs/banc_native_arr_%A_%a.out
#SBATCH -e jobs/banc_native_arr_%A_%a.err

set -e  # surface Rscript non-zero so SLURM marks task FAILED instead of COMPLETED 0:0
start=$(date +%s)
cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

echo "Native BANC self-NBLAST — array shard ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}"
# trap so even if Rscript dies during shutdown (free(): invalid pointer) we still
# report a meaningful elapsed time before the trap propagates the error.
trap 'echo "Shard ${SLURM_ARRAY_TASK_ID} aborted after $(( $(date +%s) - start )) seconds"' ERR
Rscript banc/nblast/banc-nblast-native.R

echo "Shard ${SLURM_ARRAY_TASK_ID} done in $(( $(date +%s) - start )) seconds"
