#!/bin/bash
#SBATCH -c 16
#SBATCH -t 2-06:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=15G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_ensemble_ind_v888v2_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_ensemble_ind_v888v2_%j.err
###############################################################################
# Whole-brain alignment with ensemble metric + ind_weight on v888 + v2.
#
# Validates the deferred-subset FAFB individual-profile rewrite (commit 4e14de3)
# at WB scale. ind_weight has been forced to 0 at WB until now because the
# full (n_fafb_pool, 4*n_types) materialisation OOMed; the rewrite densifies
# only the NBLAST-reachable FAFB subset (~14.5 GB at WB v888v2). Job 38353405
# was OOM-killed mid-iter on the 160 GB short node — bumped to priority
# 16c × 15G = 240 GB to clear working-memory headroom.
#
# Tier-anchored regime to match the prior 72.6% baseline; --stratified-holdout
# 1.0 was the bug that capped 38353394 (cosine) at 52.1%.
#
# Wall bumped 24h → 54h (2-06:00) after job 38508315 timed out at iter 44/80
# still climbing (best 64.1%). Iter time ~32min × 80 iter ≈ 43h; 54h ≈ 125%.
#
# Reuses the v888 + v2 prep produced by o2_banc_wb_cosine_v888v2.sh (Phase 1).
# Submit with --dependency=afterok:<cosine_job_id> if running back-to-back.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

WB_DIR=data/whole_brain_alignment_v888v2

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
    echo "FATAL: missing prep file $required — run o2_banc_wb_cosine_v888v2.sh first" >&2
    exit 2
  fi
done

echo "=== ensemble + ind_weight=0.3 align (tier-anchored, 80 iters) ==="
"$ALIGN_PY" alignment/banc-alignment-run.py \
  --side both --bilateral \
  --data-dir "$WB_DIR" --file-prefix brain \
  --metric ensemble \
  --ensemble-blend 0.3 \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --nt-weight 0 \
  --ind-weight 0.3 \
  --max-iter 80 \
  --checkpoint-every 5 \
  --forbidden-matches "$WB_DIR/forbidden-matches.csv" \
  --output-suffix ensemble_ind03_v888v2_tier

echo "=== Done ==="
