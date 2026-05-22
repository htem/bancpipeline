#!/bin/bash
#SBATCH -c 10                              # Request cores
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=8G                  # Memory per core
#SBATCH -o jobs/banc_split_%j.out          # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/banc_split_%j.err          # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "getting BANC meshes"
Rscript banc/metrics/banc-obj.R

echo "skeletonising BANC neurons"
Rscript banc/metrics/banc-calculate-skeletons.R

echo "splitting BANC neurons"
Rscript banc/metrics/banc-calculate-split.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

