#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-04:00
#SBATCH -p priority
#SBATCH --mem=64G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_aggregate_influence_priority_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_aggregate_influence_priority_%j.err
###############################################################################
# Priority-partition copy of o2_banc_aggregate_influence.sh.
# Aggregates the 277 chunk_*.parquet shards in banc_888/influence/all_to_all/
# into the sensory → all + all → effector subclass tables. Both outputs
# are written to influence/ and pushed to GCS by banc-aggregate-influence.R.
#
# Chain: depends on o2_banc_data_push_priority.sh (needs fresh
# banc_888_meta.feather for the sensory/effector lookups).
###############################################################################
set -euo pipefail

start=$(date +%s)

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

Rscript banc/influence/banc-aggregate-influence.R
Rscript banc/influence/banc-sync-influence.R

echo "script completed in $(($(date +%s) - start)) seconds"
