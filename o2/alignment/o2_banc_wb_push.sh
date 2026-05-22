#!/bin/bash
#SBATCH -p priority
#SBATCH -t 0-04:00
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J banc_wb_push
#SBATCH -o data/scheduled_runs/wb_push_%j.out
#SBATCH -e data/scheduled_runs/wb_push_%j.err
###############################################################################
# Live push of banc-alignment-update-seatable.R --region whole-brain for the v888v2 production alignment.
# Writes only fafb_alignment_* + fafb_ntac_cell_type columns; takes a
# full-table snapshot before any write (preserve-on-blank merge applied).
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
mkdir -p data/scheduled_runs

NEW_CSV=data/whole_brain_alignment_v888v2/banc_brain_both_alignment_production_v888v2_all_tier1.csv

echo "=== Live push of banc-alignment-update-seatable.R --region whole-brain ==="
echo "Input: $NEW_CSV"
Rscript alignment/banc-alignment-update-seatable.R --region whole-brain "$NEW_CSV"
