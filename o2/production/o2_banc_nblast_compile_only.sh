#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH -o jobs/nblast_compile_%j.out
#SBATCH -e jobs/nblast_compile_%j.err

# Standalone re-run of banc-nblast-compile.R only (no NBLAST recompute, no
# CAVE push, no GCS share). The script wraps its body in local({...}); local
# does NOT swallow errors but the per-target tryCatch() blocks turn errors
# into warnings on stderr. After the run, scan the .err file for "Error",
# "Warning", "skipping", "failed".

set -e  # surface any non-zero from Rscript itself
start=$(date +%s)

cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "Compiling NBLAST results into feather files"
Rscript banc/nblast/banc-nblast-compile.R

end=$(date +%s)
echo "compile completed in $((end-start)) seconds"
