#!/bin/bash
#SBATCH -c 1                          # Request cores
#SBATCH -t 00-12:00                  # Runtime in D-HH:MM format
#SBATCH -p short                     # Partition to run in
#SBATCH --mem-per-cpu=250G            # Memory per core
#SBATCH -o jobs/manc_skel_%j.out           # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/manc_skel_%j.err           # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "do MANC split"
Rscript manc/manc-split.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

