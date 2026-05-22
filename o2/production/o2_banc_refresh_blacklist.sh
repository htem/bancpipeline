#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:45
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -J banc_refresh_blacklist
#SBATCH -o data/scheduled_runs/refresh_blacklist_%j.out
#SBATCH -e data/scheduled_runs/refresh_blacklist_%j.err
###########################################################
### Pull fafb_alignment_decision == F rows from SeaTable
### and UNION them onto the existing forbidden-matches.csv
### blacklists in
###   data/optic_lobe/
###   data/whole_brain_alignment/
###########################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
mkdir -p data/scheduled_runs
Rscript alignment/banc-alignment-false-positives-append.R
