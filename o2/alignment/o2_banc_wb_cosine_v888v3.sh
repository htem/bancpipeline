#!/bin/bash
#SBATCH -c 12
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem-per-cpu=10G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v3_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v3_%j.err
###############################################################################
# Whole-brain cosine alignment test on v888 + v3 edgelists.
#
# Apples-to-apples comparison against o2_banc_wb_cosine_v888v2.sh — same
# params, only the connectivity source differs (v3 = Zetta synapse export
# read directly from GCS, no MIN_SIZE filter at edgelist build time vs the
# v2 default size>=5).
#
# Tier-anchored regime (no --stratified-holdout); see v2 wrapper for the
# full param rationale.
#
# Runs in a parallel data dir so it doesn't disturb the v2 or v850 dirs.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

WB_DIR=data/whole_brain_alignment_v888v3
FORBIDDEN_SRC=alignment/presets/whole-brain/forbidden-matches.csv

mkdir -p "$WB_DIR"

NEED_PREP=0
for f in banc_brain_both_meta.csv fafb_brain_both_meta.csv \
         banc_brain_both_edgelist.feather fafb_brain_both_edgelist.feather \
         banc_brain_both_seeds.csv fafb_type_capacity.csv \
         banc_fafb_brain_both_nblast.csv; do
  [[ -f "$WB_DIR/$f" ]] || NEED_PREP=1
done

if (( NEED_PREP )); then
  echo "=== Phase 1: prep (banc_888 + v3 edgelists) ==="
  BANC_WB_OUTPUT_DIR="$WB_DIR" \
  BANC_SYN_SOURCE=v3 \
    Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || \
    echo "  prep R exited non-zero — checking outputs..."
else
  echo "=== Phase 1: prep outputs present, skipping ==="
fi

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
  --output-suffix cosine_v888v3_tier

echo "=== Done ==="
