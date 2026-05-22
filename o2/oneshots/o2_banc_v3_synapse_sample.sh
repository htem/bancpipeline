#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem-per-cpu=48G  # 96G total — 24G/cpu OOM'd on 39969154 (one-pass 12-col read peaked >48G)
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_v3_synapse_sample_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_v3_synapse_sample_%j.err
#SBATCH -J banc_v3_sample

# -----------------------------------------------------------------------------
# Standalone v3 synapse sample regenerator. Replaces the broken Apr 24 sample
# (Coord 1 in nm, Coord 2 empty, Type=Point, ID empty) with the correct
# Coord 1 = postsyn / Coord 2 = presyn in BANC raw voxel, Type=Line, ID=syn_id.
#
# Runs banc/metrics/banc-v3-synapse-sample.R, which uses narrow parquet
# reads (12 cols, with a 2-pass fallback) to dodge the SLURM-only arrow C++
# URI bug that crashes the same step inside banc-calculate-completion.R.
# -----------------------------------------------------------------------------

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

echo "### banc-v3-synapse-sample ###"
echo "  cores=${SLURM_CPUS_PER_TASK:-?}  mem_per_cpu=${SLURM_MEM_PER_CPU:-?}"

Rscript banc/metrics/banc-v3-synapse-sample.R

end=`date +%s`
echo "v3 synapse sample finished in: $((end-start)) seconds"
