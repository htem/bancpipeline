#!/bin/bash
#SBATCH -c 10                          # cores
#SBATCH -t 2-00:00                     # 2 days (LR script took similar)
#SBATCH -p medium
#SBATCH --mem-per-cpu=8G               # 80G total — banc.dps for ~173k neurons is large
#SBATCH -o jobs/banc_native_%j.out
#SBATCH -e jobs/banc_native_%j.err

start=$(date +%s)
cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "Native BANC self-NBLAST"
Rscript banc/nblast/banc-nblast-native.R

echo "Done in $(( $(date +%s) - start )) seconds"
