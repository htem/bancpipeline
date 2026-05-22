#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem=16G
#SBATCH -o jobs/test_l2_10ids_%j.out
#SBATCH -e jobs/test_l2_10ids_%j.err
# One-shot diagnostic: try fetching L2 skels for the 10 SeaTable neurons that
# came back without a region today (2026-05-15).
set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0
Rscript jobs/test_l2_10ids.R
