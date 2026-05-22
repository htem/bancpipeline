#!/bin/bash
#SBATCH -c 4                               # 4 cores (calc is ~single-threaded, need RAM)
#SBATCH -t 2-00:00                         # 2 days (generous margin over ~28h expected per shard)
#SBATCH -p medium                          # Medium partition (>12h, up to 5 days)
#SBATCH --mem=64G                          # LU factorization of ~155K sparse matrix
#SBATCH --array=0-11                       # 12 shards covering the 230 remaining chunks
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_influence_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_influence_%A_%a.err

# -----------------------------------------------------------------------------
# BANC all-to-all influence — RESUBMISSION as a 12-way SLURM job array.
#
# Context:
#   - 70 of 300 chunks are already written under
#       /n/data1/hms/neurobio/wilson/connectomes/banc/banc_850/influence/all_to_all/
#   - The SQLite at banc_850_proofread_influence.sqlite (~601 MB) is complete,
#     so no shard needs to rebuild it. Shard 0 sees it exists and proceeds;
#     shards 1..11 fall through the "wait for SQLite" poll immediately.
#   - Each chunk takes ~90 minutes of wall time.
#
# Sharding strategy:
#   BANC_INFLUENCE_SHARD_TOTAL=12, BANC_INFLUENCE_SHARD_IDX=$SLURM_ARRAY_TASK_ID.
#   The R script assigns batch b to shard ((b - 1) %% 12). Each shard walks
#   the full 1..300 batch space but only computes its own stride, so no two
#   shards ever write the same chunk_NNNN.parquet — write-safe by construction.
#
# Resume behaviour:
#   overwrite_chunks = FALSE in the R script, so every shard simply skips any
#   chunk_NNNN.parquet that already exists on disk. The 70 finished chunks are
#   therefore skipped instantly, and each shard only does real work on its
#   remaining assigned batches. Across 12 shards the ~230 remaining chunks are
#   distributed as ~19 chunks/shard; at ~90 min/chunk that is ~28.5 h per shard,
#   well under the 2-day time limit.
#
# Submission:
#   sbatch o2_banc_influence_array.sh
#
# (The --array=0-11 and BANC_INFLUENCE_SHARD_TOTAL=12 are baked in above, so no
#  extra flags are needed. All 12 array tasks can run concurrently; O2 will
#  schedule them as resources permit.)
#
# After all shards finish:
#   Rscript banc/influence/banc-sync-influence.R
# -----------------------------------------------------------------------------

export BANC_INFLUENCE_SHARD_TOTAL=12
export BANC_INFLUENCE_SHARD_IDX="${SLURM_ARRAY_TASK_ID}"

echo "RUNNING BANC ALL-TO-ALL INFLUENCE CALCULATION (array resubmission)"
echo "  SLURM_JOB_ID=${SLURM_JOB_ID}  SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}  SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "  shard ${BANC_INFLUENCE_SHARD_IDX} of ${BANC_INFLUENCE_SHARD_TOTAL}"

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "calculating all-to-all influence scores (shard ${BANC_INFLUENCE_SHARD_IDX}/${BANC_INFLUENCE_SHARD_TOTAL})"
# Tolerate R's libc teardown crash (free(): invalid pointer on R_CleanUp)
# which fires after all chunk_NNNN.parquet files are saved.
Rscript banc/influence/banc-build-influence.R || true

# Report
end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime
