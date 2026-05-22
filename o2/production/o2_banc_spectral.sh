#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_spectral_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_spectral_%j.err
#
# Re-run spectral clustering on the latest BANC v888 + SeaTable inclusions.
# Runs v3 + v2 in parallel via background subshells (same pattern used in
# o2/o2_banc_v888_rebuild.sh GROUP B1).

set -uo pipefail
cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

start=$(date +%s)
echo "### banc spectral re-run starting $(date) ###"

# Run v3 and v2 concurrently. Each reads its own edgelist; independent memory.
Rscript banc/clustering/banc-spectral-clustering.R --source v3 \
  > jobs/banc_spectral_v3_${SLURM_JOB_ID}.log 2>&1 &
PID_V3=$!

Rscript banc/clustering/banc-spectral-clustering.R --source v2 \
  > jobs/banc_spectral_v2_${SLURM_JOB_ID}.log 2>&1 &
PID_V2=$!

# Wait, capture exit codes from both regardless of which fails first.
wait $PID_V3; RC_V3=$?
wait $PID_V2; RC_V2=$?

echo "v3 exit=$RC_V3, v2 exit=$RC_V2"
end=$(date +%s)
echo "Total runtime: $((end-start))s"

# Surface failures via job exit code so SLURM marks the job FAILED.
if [ $RC_V3 -ne 0 ] || [ $RC_V2 -ne 0 ]; then
  echo "One or both clusterings failed: v3=$RC_V3 v2=$RC_V2"
  exit 1
fi
