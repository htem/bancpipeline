#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-12:00
#SBATCH -p medium
#SBATCH --mem=128G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_betweenness_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_betweenness_%j.err
#
# Re-run betweenness centrality on the latest BANC v888 edgelists.
# Runs v3 + v2 sequentially (each ~hours; igraph is single-threaded for
# betweenness so running them concurrently wouldn't gain).
# Both an undated (canonical) and dated CSV are written by banc-betweenness.py.

set -uo pipefail
cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

start=$(date +%s)
echo "### banc betweenness re-run starting $(date) ###"

echo "--- v3 ---"
Rscript banc/betweenness/banc-betweenness-run.R both --source v3 \
  > jobs/banc_betweenness_v3_${SLURM_JOB_ID}.log 2>&1
RC_V3=$?
echo "v3 exit=$RC_V3"

echo "--- v2 ---"
Rscript banc/betweenness/banc-betweenness-run.R both --source v2 \
  > jobs/banc_betweenness_v2_${SLURM_JOB_ID}.log 2>&1
RC_V2=$?
echo "v2 exit=$RC_V2"

end=$(date +%s)
echo "Total runtime: $((end-start))s"

if [ $RC_V3 -ne 0 ] || [ $RC_V2 -ne 0 ]; then
  echo "One or both runs failed: v3=$RC_V3 v2=$RC_V2"
  exit 1
fi
