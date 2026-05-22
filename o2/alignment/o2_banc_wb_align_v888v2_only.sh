#!/bin/bash
#SBATCH -c 16
#SBATCH -t 2-12:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=15G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_align_v888v2_only_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_align_v888v2_only_%j.err
###############################################################################
# Alignment-only counterpart to o2_banc_wb_production_v888v2.sh. Skips Stage 1
# (prep) on the assumption that o2_banc_wb_prep_v888v2_priority.sh has just
# produced fresh banc_brain_both_* + seeds_production.csv. Runs the same
# production config (--manual-labels all + ensemble + ind_weight=0.3 + 80 iters)
# so align + NTAC can run in parallel after prep.
###############################################################################
set -euo pipefail
ulimit -c 0
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python
WB_DIR=data/whole_brain_alignment_v888v2

# Required prep outputs (from o2_banc_wb_prep_v888v2_priority.sh)
for required in \
  "$WB_DIR/banc_brain_both_meta.csv" \
  "$WB_DIR/fafb_brain_both_meta.csv" \
  "$WB_DIR/banc_brain_both_edgelist.feather" \
  "$WB_DIR/fafb_brain_both_edgelist.feather" \
  "$WB_DIR/banc_brain_both_seeds.csv" \
  "$WB_DIR/banc_brain_both_seeds_production.csv" \
  "$WB_DIR/banc_fafb_brain_both_nblast.csv" \
  "$WB_DIR/fafb_type_capacity.csv" \
  "$WB_DIR/forbidden-matches.csv"
do
  if [[ ! -f "$required" ]]; then
    echo "FATAL: missing prep file $required — run o2_banc_wb_prep_v888v2_priority.sh first" >&2
    exit 2
  fi
done

# Swap seeds.csv → seeds_production.csv for the align run; restore on exit
cp -p "$WB_DIR/banc_brain_both_seeds.csv" \
      "$WB_DIR/banc_brain_both_seeds.precapture.csv"
cp -p "$WB_DIR/banc_brain_both_seeds_production.csv" \
      "$WB_DIR/banc_brain_both_seeds.csv"
trap 'cp -p "$WB_DIR/banc_brain_both_seeds.precapture.csv" "$WB_DIR/banc_brain_both_seeds.csv" || true' EXIT

echo "=== Align: ensemble + ind_weight=0.3 + --manual-labels all (80 iters) ==="
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
  --manual-labels all \
  --max-iter 80 \
  --checkpoint-every 5 \
  --forbidden-matches "$WB_DIR/forbidden-matches.csv" \
  --output-suffix production_v888v2_all_tier1

echo "=== Done ==="
