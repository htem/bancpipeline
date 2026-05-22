#!/bin/bash
#SBATCH -c 8
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=96G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_optic_sweep_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_optic_sweep_%A_%a.err
#SBATCH --array=0-12
# Memory note: NTAC shards 9-12 (anchor frac >= 50%) OOM-killed at 32G on
# job 36588536 — NTAC's footprint grows with seed density. Bumped to 96G;
# wasteful for align shards but keeps the array uniform.
###############################################################################
# Optic lobe alignment sweep — cheap sanity-check mirror of o2_banc_wb_sweep.sh.
#
# Same 13-condition matrix as the WB sweep, but on optic_lobe data. Iteration
# cost is ~20s (vs 30+ min on WB), so this finishes in under an hour per shard
# and gives us a scaling curve we can compare the WB sweep against.
#
# Run this first:
#   Rscript alignment/banc-alignment-prep.R --region optic-lobe       both --source local
#   Rscript alignment/banc-alignment-false-positives.R
# Then: sbatch o2_banc_optic_sweep.sh
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python
NTAC_PY=/home/ab714/ntac-env/bin/python
module load gcc/14.2.0 python/3.13.1

OPTIC_DIR=data/optic_lobe
FORBIDDEN=$OPTIC_DIR/forbidden-matches.csv

for required in \
  "$OPTIC_DIR/banc_optic_both_meta.csv" \
  "$OPTIC_DIR/fafb_optic_both_meta.csv" \
  "$OPTIC_DIR/banc_optic_both_edgelist.feather" \
  "$OPTIC_DIR/fafb_optic_both_edgelist.feather" \
  "$OPTIC_DIR/banc_optic_both_seeds.csv" \
  "$OPTIC_DIR/fafb_type_capacity.csv" \
  "$FORBIDDEN"
do
  if [[ ! -f "$required" ]]; then
    echo "FATAL: missing prep file $required — run optic prep + false-positives first" >&2
    exit 2
  fi
done

ALIGN_COMMON=(
  --side both --bilateral
  --data-dir "$OPTIC_DIR" --file-prefix optic
  --metric ensemble
  --hop2-weight 1.0
  --nblast-threshold 0.15
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95
  --ensemble-blend 0.3
  --nt-weight 0 --ind-weight 0.5
  --max-iter 60
  --forbidden-matches "$FORBIDDEN"
)

NTAC_COMMON=(
  --side both
  --data-dir "$OPTIC_DIR" --file-prefix optic
  --max-iter 60
)

run_align_holdout() {
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" --stratified-holdout "$1" --output-suffix "$2"
}
run_align_holdout_anchored() {
  # --manual-labels Mi1: anchor only Mi1 type to its true cell_type. Bare
  # --manual-labels would anchor ALL holdout types and game the eval, masking
  # algorithm failures (caught 2026-04-25 via the round-of-diags showing the
  # iter-1 collapse was just anchor-echo, not a real algorithm regression).
  # Per alignment/experiment-log.md "Final configuration (2026-04-02)".
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" --stratified-holdout "$1" --manual-labels Mi1 --output-suffix "$2"
}
run_align_all_anchors() {
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" --manual-labels Mi1 --output-suffix "$1"
}
run_ntac_holdout() {
  "$NTAC_PY" alignment/banc-alignment-ntac.py \
    "${NTAC_COMMON[@]}" --stratified-holdout "$1" --output-suffix "$2"
}
run_ntac_no_holdout() {
  "$NTAC_PY" alignment/banc-alignment-ntac.py \
    "${NTAC_COMMON[@]}" --no-holdout --output-suffix "$1"
}

idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID unset}"
echo "=== Optic sweep shard $idx (job ${SLURM_ARRAY_JOB_ID:-n/a}) ==="
start=$(date +%s)

case "$idx" in
  0)  run_align_holdout          1.0 sweep_align_a000_ho100 ;;
  1)  run_align_holdout_anchored 0.9 sweep_align_a010_ho090 ;;
  2)  run_align_holdout_anchored 0.8 sweep_align_a020_ho080 ;;
  3)  run_align_holdout_anchored 0.5 sweep_align_a050_ho050 ;;
  4)  run_align_holdout_anchored 0.2 sweep_align_a080_ho020 ;;
  5)  run_align_holdout_anchored 0.1 sweep_align_a090_ho010 ;;
  6)  run_align_all_anchors          sweep_align_a100_ho000 ;;
  7)  run_ntac_holdout           0.9 sweep_ntac_a010_ho090  ;;
  8)  run_ntac_holdout           0.8 sweep_ntac_a020_ho080  ;;
  9)  run_ntac_holdout           0.5 sweep_ntac_a050_ho050  ;;
  10) run_ntac_holdout           0.2 sweep_ntac_a080_ho020  ;;
  11) run_ntac_holdout           0.1 sweep_ntac_a090_ho010  ;;
  12) run_ntac_no_holdout            sweep_ntac_a100_ho000  ;;
  *)  echo "FATAL: unknown shard idx $idx" >&2; exit 3 ;;
esac

echo "=== Shard $idx done ($(( $(date +%s) - start ))s) ==="
