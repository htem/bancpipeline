#!/bin/bash
#SBATCH -c 8
#SBATCH -t 5-00:01
#SBATCH -p long
#SBATCH --mem=249G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_ntac_rerun_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_ntac_rerun_%A_%a.err
#SBATCH --array=0-2
###############################################################################
# Rerun of NTAC shards 10-12 from o2_banc_wb_sweep.sh — these OOM-killed at
# 160G even with --mem-per-cpu=10G. Bumped to a flat --mem=300G on priority.
#
# Array index -> condition:
#   0 -> ntac 80% seeds  (sweep_ntac_a080_ho020)
#   1 -> ntac 90% seeds  (sweep_ntac_a090_ho010)
#   2 -> ntac 100% seeds (sweep_ntac_a100_ho000, no holdout)
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

NTAC_PY=/home/ab714/ntac-env/bin/python
module load gcc/14.2.0 python/3.13.1

WB_DIR=data/whole_brain_alignment
NTAC_COMMON=(
  --side both
  --data-dir "$WB_DIR" --file-prefix brain
  --max-iter 80
)

idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID unset}"
echo "=== WB NTAC rerun shard $idx (job ${SLURM_ARRAY_JOB_ID:-n/a}) ==="
start=$(date +%s)

case "$idx" in
  0) "$NTAC_PY" alignment/banc-alignment-ntac.py "${NTAC_COMMON[@]}" \
         --stratified-holdout 0.2 --output-suffix sweep_ntac_a080_ho020 ;;
  1) "$NTAC_PY" alignment/banc-alignment-ntac.py "${NTAC_COMMON[@]}" \
         --stratified-holdout 0.1 --output-suffix sweep_ntac_a090_ho010 ;;
  2) "$NTAC_PY" alignment/banc-alignment-ntac.py "${NTAC_COMMON[@]}" \
         --no-holdout            --output-suffix sweep_ntac_a100_ho000 ;;
  *) echo "FATAL: unknown shard idx $idx" >&2; exit 3 ;;
esac

echo "=== Shard $idx done ($(( $(date +%s) - start ))s) ==="
