#!/bin/bash
#SBATCH -c 16
#SBATCH -t 3-00:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=10G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_sweep_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_sweep_%A_%a.err
#SBATCH --array=0-12
###############################################################################
# Whole-brain alignment sweep — 13 conditions.
#
# Before launching this array, run:
#   sbatch o2_banc_wb_probe.sh          # (preps data once + measures ensemble cost)
#
# This script assumes data/whole_brain_alignment/ is already populated by the
# probe job's prep step. Each array task runs ONE condition and writes its
# output to a unique suffix. No SeaTable writeback.
#
# Conditions (array index -> method, anchor fraction):
#   0  align  0%   (all holdout, Tier 2 NBLAST seeds retained)
#   1  align  10%
#   2  align  20%
#   3  align  50%
#   4  align  80%
#   5  align  90%
#   6  align  100% (manual-labels, no holdout)
#   7  ntac   10%
#   8  ntac   20%
#   9  ntac   50%
#   10 ntac   80%
#   11 ntac   90%
#   12 ntac   100% (no-holdout)
#
# NOTE: NTAC still lacks --forbidden-matches, so its runs skip the
# forbidden-pair mask. Stratified holdout is now shared with align.py via
# alignment/alignment_splits.py, so NTAC and align at the same fraction
# hold out the exact same neurons. NTAC writes is_holdout to its output CSV
# so the eval script can compute holdout accuracy directly.
#
# Align params match experiment-log.md Run 13 + individual-profile blend:
#   ensemble metric, hop2_weight 1.0, nblast_threshold 0.15, tau_start 4.0,
#   alpha 0.05 -> 0.95, ensemble_blend 0.3, nt_weight 0, ind_weight 0.5.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Separate Python environments:
#   ALIGN_PY: reticulate conda env (Python 3.9) — used by align.py
#   NTAC_PY:  dedicated venv under O2 python/3.13.1 — needed because the
#             ntac package on PyPI requires Python >= 3.10
ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python
NTAC_PY=/home/ab714/ntac-env/bin/python
# NTAC venv's libpython3.13.so lives under the module tree, so we need the
# module loaded in this shell for LD_LIBRARY_PATH to resolve it.
module load gcc/14.2.0 python/3.13.1

WB_DIR=data/whole_brain_alignment
FORBIDDEN=$WB_DIR/forbidden-matches.csv

# Sanity check: prep outputs must already exist (run the probe first).
for required in \
  "$WB_DIR/banc_brain_both_meta.csv" \
  "$WB_DIR/fafb_brain_both_meta.csv" \
  "$WB_DIR/banc_brain_both_edgelist.feather" \
  "$WB_DIR/fafb_brain_both_edgelist.feather" \
  "$WB_DIR/banc_brain_both_seeds.csv" \
  "$WB_DIR/fafb_type_capacity.csv" \
  "$FORBIDDEN"
do
  if [[ ! -f "$required" ]]; then
    echo "FATAL: missing prep file $required — run o2_banc_wb_probe.sh first" >&2
    exit 2
  fi
done

ALIGN_COMMON=(
  --side both --bilateral
  --data-dir "$WB_DIR" --file-prefix brain
  --metric ensemble
  --hop2-weight 1.0
  --nblast-threshold 0.15
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95
  --ensemble-blend 0.3
  --nt-weight 0
  # ind_weight=0.0 for WB: individual-profile scoring OOMs at 127k×35k scale
  # (confirmed on probe 36588518). Ensemble metric is still on.
  --ind-weight 0.0
  --max-iter 50
  --checkpoint-every 5
  --forbidden-matches "$FORBIDDEN"
)

NTAC_COMMON=(
  --side both
  --data-dir "$WB_DIR" --file-prefix brain
  --max-iter 80
)

run_align_holdout() {
  local frac="$1"; local suffix="$2"
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" \
    --stratified-holdout "$frac" \
    --output-suffix "$suffix"
}

run_align_holdout_anchored() {
  # --manual-labels Mi1: anchor only Mi1 type. Bare --manual-labels would anchor
  # ALL holdouts and game the eval (caught 2026-04-25). Per experiment-log.md
  # "Final configuration (2026-04-02)".
  local frac="$1"; local suffix="$2"
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" \
    --stratified-holdout "$frac" --manual-labels Mi1 \
    --output-suffix "$suffix"
}

run_align_all_anchors() {
  local suffix="$1"
  "$ALIGN_PY" alignment/banc-alignment-run.py \
    "${ALIGN_COMMON[@]}" \
    --manual-labels Mi1 \
    --output-suffix "$suffix"
}

run_ntac_holdout() {
  local frac="$1"; local suffix="$2"
  "$NTAC_PY" alignment/banc-alignment-ntac.py \
    "${NTAC_COMMON[@]}" \
    --stratified-holdout "$frac" \
    --output-suffix "$suffix"
}

run_ntac_no_holdout() {
  local suffix="$1"
  "$NTAC_PY" alignment/banc-alignment-ntac.py \
    "${NTAC_COMMON[@]}" \
    --no-holdout \
    --output-suffix "$suffix"
}

idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID unset — submit via sbatch --array}"
echo "=== WB sweep shard $idx (job ${SLURM_ARRAY_JOB_ID:-n/a}) ==="
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
