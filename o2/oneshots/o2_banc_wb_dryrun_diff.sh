#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J banc_wb_dryrun_diff
#SBATCH -o data/scheduled_runs/wb_dryrun_diff_%j.out
#SBATCH -e data/scheduled_runs/wb_dryrun_diff_%j.err
###############################################################################
# Dry-run of banc-alignment-update-seatable.R --region whole-brain AND a row-level diff of the new
# alignment results CSV against the snapshot of the previous run's results.
# Triggered after the running production v888v2 alignment job completes:
#   sbatch --dependency=afterok:<production_jobid> o2_banc_wb_dryrun_diff.sh
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
mkdir -p data/scheduled_runs

NEW_CSV=data/whole_brain_alignment_v888v2/banc_brain_both_alignment_production_v888v2_all_tier1.csv
OLD_CSV=data/whole_brain_alignment_v888v2/banc_brain_both_alignment_production_v888v2_all_tier1.may17.csv
DIFF_OUT=data/scheduled_runs/wb_dryrun_diff_${SLURM_JOB_ID}_rowdiff.csv

echo "=== Stage A: alignment-results diff (new vs previous tier1) ==="
if [[ -f "$NEW_CSV" && -f "$OLD_CSV" ]]; then
  Rscript alignment/banc-alignment-diff.R "$NEW_CSV" "$OLD_CSV" "$DIFF_OUT"
else
  echo "  Missing $NEW_CSV or $OLD_CSV — skipping result-level diff"
fi

echo ""
echo "=== Stage B: dry-run of banc-alignment-update-seatable.R --region whole-brain ==="
Rscript alignment/banc-alignment-update-seatable.R --region whole-brain --dry-run "$NEW_CSV"
