#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2-12:00
#SBATCH -p priority
#SBATCH --mem=250G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_ntac_v888v2_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_ntac_v888v2_%j.err
###############################################################################
# WB NTAC run on v888v2 data — "100% anchors" / no-holdout (the NTAC analogue
# of the production --manual-labels all alignment in
# o2_banc_wb_align_v888v2_only.sh). Reads the prep outputs from
# o2_banc_wb_prep_v888v2_priority.sh.
#
# Memory: WB-scale NTAC at 100% anchors OOM'd at 160G even with
# --mem-per-cpu=10G (see o2_banc_wb_ntac_rerun.sh). Flat 250G on priority
# (MaxMemPerNode=257G — 300G was rejected as un-allocatable).
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

NTAC_PY=/home/ab714/ntac-env/bin/python
module load gcc/14.2.0 python/3.13.1

WB_DIR=data/whole_brain_alignment_v888v2

for required in \
  "$WB_DIR/banc_brain_both_meta.csv" \
  "$WB_DIR/fafb_brain_both_meta.csv" \
  "$WB_DIR/banc_brain_both_edgelist.feather" \
  "$WB_DIR/fafb_brain_both_edgelist.feather" \
  "$WB_DIR/banc_brain_both_seeds.csv" \
  "$WB_DIR/fafb_type_capacity.csv"
do
  if [[ ! -f "$required" ]]; then
    echo "FATAL: missing prep file $required — run o2_banc_wb_prep_v888v2_priority.sh first" >&2
    exit 2
  fi
done

start=$(date +%s)
echo "=== WB NTAC v888v2 --no-holdout (all GT as anchors) ==="
"$NTAC_PY" alignment/banc-alignment-ntac.py \
  --side both \
  --data-dir "$WB_DIR" --file-prefix brain \
  --max-iter 80 \
  --no-holdout \
  --output-suffix production_v888v2_ntac_all_tier1

echo "=== Done ($(( ($(date +%s) - start)/60 )) min) ==="
