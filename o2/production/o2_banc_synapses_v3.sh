#!/bin/bash
#SBATCH -c 8                                 # cores (region/neuropil classification)
#SBATCH -t 02-00:00                          # 2 days
#SBATCH -p priority
#SBATCH --mem-per-cpu=16G                    # 8 * 16 = 128 GB total
#SBATCH -o jobs/banc_synapses_v3_%j.out
#SBATCH -e jobs/banc_synapses_v3_%j.err

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "### banc-synapses-v3: process v3 predictions and compute capture rates ###"
Rscript banc/metrics/banc-synapses-v3.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: $runtime seconds"
