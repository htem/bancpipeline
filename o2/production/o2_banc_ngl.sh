#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 0-96:00                        # Runtime in D-HH:MM format
#SBATCH -p medium                         # Partition to run in
#SBATCH --mem-per-cpu=10G                # Memory per core
#SBATCH -o jobs/banc_ngl_%j.out            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/banc_ngl_%j.err            # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "transforming maleCNS meshes to BANC space"
Rscript banc/transforms/banc-malecns-mesh-transform.R

echo "transforming FAFB-783 meshes to BANC space"
Rscript banc/transforms/banc-fafb-mesh-transform.R

echo "transforming MANC meshes into BANC space"
Rscript banc/transforms/banc-manc-mesh-transform.R

echo "transforming hemibrain meshes to BANC space"
Rscript banc/transforms/banc-hemibrain-mesh-transform.R

echo "uploading data for neuroglancer scenes"
Rscript banc/transforms/banc-ngl-upload.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

