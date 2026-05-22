#!/bin/bash
#SBATCH -c 8
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem-per-cpu=8G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_ol_ind_smoketest_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_ol_ind_smoketest_%j.err
###############################################################################
# OL ind_weight=0.5 smoke test — confirms commit 4e14de3's deferred-subset
# materialisation runs end-to-end on real data (OL bilateral prep already on
# disk). Compares against banc_optic_both_alignment_with_rules_holdout.csv
# baseline (82.4% holdout, Apr 10 from local run).
#
# 10 iterations is enough to see holdout climb and confirm no crash. Full
# rerun (60 iter) is a separate concern.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

"$ALIGN_PY" alignment/banc-alignment-run.py \
  --side both --bilateral \
  --metric ensemble \
  --ensemble-blend 0.3 \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --nt-weight 0 \
  --ind-weight 0.5 \
  --max-iter 10 \
  --forbidden-matches alignment/presets/optic-lobe/forbidden-matches.csv \
  --output-suffix ind_subset_smoke
