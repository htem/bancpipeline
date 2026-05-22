#!/bin/bash
#SBATCH -c 10                              # Request cores
#SBATCH -t 0-96:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem-per-cpu=10G                  # Memory per core
#SBATCH -o jobs/banc_meshes_%j.out                # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/banc_meshes_%j.err                # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "saving banc .obj files"
Rscript banc/metrics/banc-obj.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

