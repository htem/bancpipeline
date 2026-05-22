#!/bin/bash
#SBATCH -c 16
#SBATCH -t 3-00:00                           # generous — we do not yet know WB+ensemble wall-time
#SBATCH -p medium
#SBATCH --mem=160G                           # extra headroom for ensemble re-rank buffers
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_probe_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_probe_%j.err
###############################################################################
# WB ensemble probe run — de-risks the full sweep.
#
# Runs whole-brain prep (local O2 data) once, refreshes forbidden-matches, then
# executes a single whole-brain alignment at 100% anchors with the current best
# params INCLUDING the ensemble metric. The goal is to measure per-iter wall
# time and peak memory so we can size the sweep array correctly.
#
# Outputs:
#   data/whole_brain_alignment/banc_brain_both_alignment_probe_ensemble_g.csv
#
# Inspect .err log after completion:
#   grep -E 'iter .* took|MaxRSS' jobs/banc_wb_probe_*.err
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python

echo "=== WB ensemble probe (job ${SLURM_JOB_ID:-n/a}) ==="
start=$(date +%s)

# Some R/C library tears down with `free(): invalid pointer` on Rscript exit
# even after every output file has been written successfully. We tolerate the
# crash by gating on file existence rather than on exit code.
require_files() {
  for f in "$@"; do
    if [[ ! -f "$f" ]]; then
      echo "FATAL: missing prep output $f" >&2
      exit 4
    fi
  done
}

echo "[prep_wb]"
Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || true
require_files \
  data/whole_brain_alignment/banc_brain_both_meta.csv \
  data/whole_brain_alignment/fafb_brain_both_meta.csv \
  data/whole_brain_alignment/banc_brain_both_edgelist.feather \
  data/whole_brain_alignment/fafb_brain_both_edgelist.feather \
  data/whole_brain_alignment/banc_brain_both_seeds.csv \
  data/whole_brain_alignment/fafb_type_capacity.csv
echo "  prep_wb outputs verified"

echo "[refresh_forbidden]"
Rscript alignment/banc-alignment-false-positives.R || true
require_files \
  alignment/presets/whole-brain/forbidden-matches.csv
echo "  refresh_forbidden outputs verified"

echo "[align_wb_ensemble_g]"
# ind_weight=0.0 for WB: individual-profile scoring materialises BANC×types
# dense buffers that don't fit in memory at 127k neurons × 35k types even at
# 160G (OOM confirmed on job 36588518). Ensemble metric is still on.
"$ALIGN_PY" alignment/banc-alignment-run.py \
  --side both --bilateral \
  --data-dir data/whole_brain_alignment --file-prefix brain \
  --metric ensemble \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --ensemble-blend 0.3 \
  --nt-weight 0 --ind-weight 0.0 \
  --max-iter 80 \
  --manual-labels \
  --forbidden-matches alignment/presets/whole-brain/forbidden-matches.csv \
  --output-suffix probe_ensemble_g

echo "=== Done ($(( $(date +%s) - start ))s) ==="
echo "Probe output: data/whole_brain_alignment/banc_brain_both_alignment_probe_ensemble_g.csv"
