#!/bin/bash
#SBATCH -c 12
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem-per-cpu=10G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v2_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v2_%j.err
###############################################################################
# Whole-brain cosine alignment test on v888 + v2 edgelists.
#
# Goal: reproduce the prior 72.6% holdout (run_no_manual_labels.log, Apr 7).
# Match that run's anchoring regime: NO --stratified-holdout flag — the
# default tier-based init in align.py keys on seeds.tier==1 & is_holdout==FALSE
# from the prep-time seeds file (~16,822 anchors at WB scale). Passing
# --stratified-holdout 1.0 (prior wrapper) zeros all anchors and capped the
# 38353394 run at 52.1%.
#
# Differences vs o2_banc_wb_sweep.sh:
#   - --metric cosine    (sweep uses ensemble + ensemble_blend 0.3)
#   - --max-iter 80      (sweep capped at 50; prior best peaked iter 30+)
#   - no --ensemble-blend
#   - no --manual-labels (clean baseline; matches the prior 72.6% run)
#   - no --stratified-holdout (use tier-based default; matches prior 72.6%)
#   - ind_weight stays 0
#
# Runs in a parallel data dir so it doesn't disturb the v850 sweep dir.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

WB_DIR=data/whole_brain_alignment_v888v2
FORBIDDEN_SRC=alignment/presets/whole-brain/forbidden-matches.csv

mkdir -p "$WB_DIR"

# Re-prep into the parallel directory using v888 + v2 edgelists, but only
# if the outputs aren't already present. The prep R script crashed with
# `free(): invalid pointer` on interpreter teardown in job 38292258 *after*
# writing every required output; rerunning prep just rolls the dice on
# whether the abort kills the wrapper before Phase 2 starts. The existence
# checks below are the source of truth — Phase 2 won't run if anything's
# missing.
NEED_PREP=0
for f in banc_brain_both_meta.csv fafb_brain_both_meta.csv \
         banc_brain_both_edgelist.feather fafb_brain_both_edgelist.feather \
         banc_brain_both_seeds.csv fafb_type_capacity.csv \
         banc_fafb_brain_both_nblast.csv; do
  [[ -f "$WB_DIR/$f" ]] || NEED_PREP=1
done

if (( NEED_PREP )); then
  echo "=== Phase 1: prep (banc_888 + v2 edgelists) ==="
  BANC_WB_OUTPUT_DIR="$WB_DIR" \
  BANC_SYN_SOURCE=v2 \
    Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || \
    echo "  prep R exited non-zero — checking outputs..."
else
  echo "=== Phase 1: prep outputs present, skipping ==="
fi

# Forbidden matches are dataset-version-agnostic (they're FAFB root_ids paired
# with cell types) — reuse the merged 1186-entry list from the main dir.
[[ -f "$WB_DIR/forbidden-matches.csv" ]] || cp "$FORBIDDEN_SRC" "$WB_DIR/forbidden-matches.csv"

for required in \
  "$WB_DIR/banc_brain_both_meta.csv" \
  "$WB_DIR/fafb_brain_both_meta.csv" \
  "$WB_DIR/banc_brain_both_edgelist.feather" \
  "$WB_DIR/fafb_brain_both_edgelist.feather" \
  "$WB_DIR/banc_brain_both_seeds.csv" \
  "$WB_DIR/fafb_type_capacity.csv" \
  "$WB_DIR/forbidden-matches.csv"
do
  if [[ ! -f "$required" ]]; then
    echo "FATAL: prep did not produce $required" >&2
    exit 2
  fi
done
echo "  prep outputs verified"

echo "=== Phase 2: cosine align (tier-anchored, 80 iters) ==="
"$ALIGN_PY" alignment/banc-alignment-run.py \
  --side both --bilateral \
  --data-dir "$WB_DIR" --file-prefix brain \
  --metric cosine \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --nt-weight 0 \
  --ind-weight 0.0 \
  --max-iter 80 \
  --checkpoint-every 5 \
  --forbidden-matches "$WB_DIR/forbidden-matches.csv" \
  --output-suffix cosine_v888v2_tier

echo "=== Done ==="
