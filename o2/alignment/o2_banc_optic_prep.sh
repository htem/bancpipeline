#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=32G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_optic_prep_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_optic_prep_%j.err
###############################################################################
# Optic lobe alignment prep (local v850 data) + forbidden-matches refresh.
# Prerequisite for o2_banc_optic_sweep.sh.
#
#   sbatch o2_banc_optic_prep.sh
#   # when done:
#   sbatch o2_banc_optic_sweep.sh
#
# Or chain with dependency:
#   PREP=$(sbatch --parsable o2_banc_optic_prep.sh)
#   sbatch --dependency=afterok:$PREP o2_banc_optic_sweep.sh
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "=== Optic prep (job ${SLURM_JOB_ID:-n/a}) ==="
start=$(date +%s)

# Tolerate Rscript shutdown crashes (free(): invalid pointer in C library
# teardown after all files are written) — gate on file existence instead.
require_files() {
  for f in "$@"; do
    if [[ ! -f "$f" ]]; then
      echo "FATAL: missing prep output $f" >&2
      exit 4
    fi
  done
}

echo "[prep_optic]"
Rscript alignment/banc-alignment-prep.R --region optic-lobe both --source local || true
require_files \
  data/optic_lobe/banc_optic_both_meta.csv \
  data/optic_lobe/fafb_optic_both_meta.csv \
  data/optic_lobe/banc_optic_both_edgelist.feather \
  data/optic_lobe/fafb_optic_both_edgelist.feather \
  data/optic_lobe/banc_optic_both_seeds.csv \
  data/optic_lobe/fafb_type_capacity.csv
echo "  prep_optic outputs verified"

echo "[refresh_forbidden]"
Rscript alignment/banc-alignment-false-positives.R || true
require_files \
  alignment/presets/optic-lobe/forbidden-matches.csv \
  alignment/presets/whole-brain/forbidden-matches.csv
echo "  refresh_forbidden outputs verified"

echo "=== Done ($(( $(date +%s) - start ))s) ==="
