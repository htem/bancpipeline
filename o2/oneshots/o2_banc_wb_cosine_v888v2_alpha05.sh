#!/bin/bash
#SBATCH -c 12
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem-per-cpu=10G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v2_alpha05_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_cosine_v888v2_alpha05_%j.err
###############################################################################
# A/B variant of o2_banc_wb_cosine_v888v2.sh: alpha-start lowered to 0.05.
#
# Hypothesis: starting NBLAST-heavy (alpha=0.05 → connectivity weight 5%,
# NBLAST 95%) lets the system stabilise on morphology while connectivity
# scores are still noisy from sparse early labels, then anneals to the same
# connectivity-dominant 0.95 endpoint. May reduce early cell-type churn at
# the cost of more total NBLAST exposure across the linear schedule.
#
# Reuses v888v2 prep cache produced by 38508224. Same suffix layout so output
# lands beside the baseline for direct comparison; suffix `_alpha05` distinguishes.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-12}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

WB_DIR=data/whole_brain_alignment_v888v2

for required in \
  "$WB_DIR/banc_brain_both_meta.csv" \
  "$WB_DIR/fafb_brain_both_meta.csv" \
  "$WB_DIR/banc_brain_both_edgelist.feather" \
  "$WB_DIR/fafb_brain_both_edgelist.feather" \
  "$WB_DIR/banc_brain_both_seeds.csv" \
  "$WB_DIR/fafb_type_capacity.csv" \
  "$WB_DIR/banc_fafb_brain_both_nblast.csv" \
  "$WB_DIR/forbidden-matches.csv"
do
  [[ -f "$required" ]] || { echo "FATAL: missing $required" >&2; exit 2; }
done
echo "  prep cache verified — skipping Phase 1"

echo "=== Phase 2: cosine align (tier-anchored, alpha 0.05→0.95, 80 iters) ==="
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
  --output-suffix cosine_v888v2_tier_alpha05

echo "=== Done ==="
