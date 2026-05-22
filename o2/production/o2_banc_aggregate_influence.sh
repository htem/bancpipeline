#!/bin/bash
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH --mem=64G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_aggregate_influence_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_aggregate_influence_%j.err

start=$(date +%s)

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

Rscript banc/influence/banc-aggregate-influence.R
Rscript banc/influence/banc-sync-influence.R

echo "script completed in $(($(date +%s) - start)) seconds"
