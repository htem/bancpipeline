#!/bin/bash
#SBATCH -c 4                                  # 4 cores (mostly I/O + classify_status)
#SBATCH -t 0-04:00                            # 4 hours (load 3 parquets + 12 CSVs)
#SBATCH -p short                              # short partition
#SBATCH --mem-per-cpu=48G                     # 4 * 48 = 192G total (was 80G; v3 spatial OOM'd
                                              # after adding mean_score+X/Y/Z cols for sample)
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_completion_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_completion_%j.err

# -----------------------------------------------------------------------------
# BANC synapse completion metrics — v1 vs v2 vs v3 capture rates.
#
# Pulls v888 v1 + v2 synapse parquets from GCS, reads v3 from local
# (banc_<ver>_synapses_v3.parquet, produced by o2_banc_synapses_v3_optimised.sh
# Job 2), classifies each synapse's pre/post root_id as
# proofread / identified / fragment, and writes 12 capture-rate CSVs +
# a 3-way comparison summary so we can judge which detection captures the
# most identified connectivity.
#
# Resource notes:
#   - The biggest in-memory object is the v3 parquet (~6-8G when read). With
#     the v1+v2 parquets + the SeaTable banc.meta + identified.ids vector,
#     80G total is comfortable. Bump if you see Killed in the .err.
#   - GCS pull is idempotent (cached under synapses_v3/cache/v<ver>_inputs/).
#     First run downloads ~12G (v1=3.7G, v2=3.4G, v3 already local).
#
# Submission:
#   # Default: target = banc.version (currently 888)
#   sbatch o2_banc_completion.sh
#
#   # Override (e.g. compute v850 metrics):
#   BANC_V3_TARGET_VERSION=850 sbatch o2_banc_completion.sh
#
# Outputs in /n/data1/.../banc/synapses_v3/capture_rates/:
#   banc_<ver>_v1_<gross|inout|region|neuropil>_capture_rates.csv
#   banc_<ver>_v2_<gross|inout|region|neuropil>_capture_rates.csv
#   banc_<ver>_v3_<gross|inout|region|neuropil>_capture_rates.csv
#   banc_<ver>_v1_v2_v3_summary.csv
#   banc_<ver>_v1_v2_v3_region_summary.csv
# -----------------------------------------------------------------------------

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

echo "### banc completion: v1 vs v2 vs v3 capture rates ###"
echo "  cores=${SLURM_CPUS_PER_TASK:-?}  mem_per_cpu=${SLURM_MEM_PER_CPU:-?}"
echo "  BANC_V3_TARGET_VERSION=${BANC_V3_TARGET_VERSION:-<default>}"

Rscript banc/metrics/banc-calculate-completion.R

end=`date +%s`
runtime=$((end-start))
echo "completion job finished in: $runtime seconds"
