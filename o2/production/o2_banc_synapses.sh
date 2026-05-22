#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 00-600:00                      # Runtime in D-HH:MM format
#SBATCH -p priority                       # Partition to run in
#SBATCH --mem-per-cpu=250G                # Memory per core
#SBATCH -o jobs/banc_synapses_%j.out      # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/banc_synapses_%j.err      # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "calculating per-neuron synapse metrics and compartment labels"
Rscript banc/metrics/banc-calculate-synapses.R

echo "downloading synapse table and building connectivity"
Rscript banc/metrics/banc-calculate-connectivity.R

echo "calculating neuropil inclusion and completion metrics"
Rscript banc/metrics/banc-calculate-neuropil-inclusion.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

