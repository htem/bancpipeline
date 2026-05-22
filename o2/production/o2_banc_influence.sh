#!/bin/bash
#SBATCH -c 4                               # Request cores (influence calc is single-threaded, but need RAM)
#SBATCH -t 1-00:00                         # Runtime: 24 hours (155K neurons * 0.4s + factorization)
#SBATCH -p medium                          # Medium partition for >12h jobs
#SBATCH --mem=64G                          # LU factorization of ~155K sparse matrix needs substantial RAM
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_influence_%A_%a.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_influence_%A_%a.err
# Job 38760185 (2026-05-04) all 4 shards failed in <2min on a-16-{108,109,139}
# with OMPI's PMI2_Init failing under SLURM. Apr 13 successful runs landed on
# different nodes (a-16-65, a-17-1XX). Excluding the known-failing trio.
#SBATCH --exclude=compute-a-16-108,compute-a-16-109,compute-a-16-139

# -----------------------------------------------------------------------------
# Sharded run via SLURM job arrays:
#
#   sbatch --array=0-3 --export=ALL,BANC_INFLUENCE_SHARD_TOTAL=4 o2_banc_influence.sh
#
# Each array task processes a disjoint stride of chunk indices, so chunks
# never collide. Shard 0 builds the SQLite if needed; shards 1..N-1 wait
# (up to 4h) for it to appear, so submit them all at once.
#
# Single-job (un-sharded) run is still supported — just `sbatch` with no
# --array, and BANC_INFLUENCE_SHARD_TOTAL defaults to 1.
#
# After all shards finish, run the sync script (separate sbatch / login):
#
#   Rscript banc/influence/banc-sync-influence.R
#
# -----------------------------------------------------------------------------

export BANC_INFLUENCE_SHARD_IDX="${SLURM_ARRAY_TASK_ID:-0}"
export BANC_INFLUENCE_SHARD_TOTAL="${BANC_INFLUENCE_SHARD_TOTAL:-1}"

echo "RUNNING BANC ALL-TO-ALL INFLUENCE CALCULATION"
echo "  shard ${BANC_INFLUENCE_SHARD_IDX} of ${BANC_INFLUENCE_SHARD_TOTAL}"

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

# OpenMPI PMI workaround: PETSc (used by influence_calculator_py) calls
# MPI_Init unconditionally; OMPI auto-detects SLURM and tries to attach to
# SLURM's PMI server. The conda-env-bundled OMPI was apparently rebuilt
# without SLURM PMI support, so MPI_Init fails fast with "PMI2_Init failed"
# and crashes Python before any chunks are written (job 38760185 / 38788333
# both died in <2min on different nodes — cluster-wide, not node-specific).
# Run Rscript in a subshell with PMI env vars stripped so OMPI can't see
# the SLURM allocation and falls back to single-process mode. The MCA
# settings are belt+braces in case OMPI rediscovers SLURM via other paths.
export OMPI_MCA_pmix=^slurm
export OMPI_MCA_psec=native
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_plm=^slurm
export PMIX_MCA_psec=native

echo "calculating all-to-all influence scores"
# Subshell with PMI handles stripped — see comment above. Parent shell keeps
# its full SLURM env so things like sacct, scontrol, etc. still work.
(
  unset SLURM_PMI_FD SLURM_PMI_RANK SLURM_PMI_SIZE \
        PMI_FD PMI_RANK PMI_SIZE \
        PMIX_RANK PMIX_NAMESPACE PMIX_SERVER_URI PMIX_SERVER_TMPDIR \
        SLURM_TASK_PID SLURM_LOCALID SLURM_PROCID SLURM_STEPID \
        SLURM_STEP_ID SLURM_STEP_NUM_NODES SLURM_STEP_NODELIST
  Rscript banc/influence/banc-build-influence.R
)

# Report
end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime
