#!/bin/bash
# Per-dataset NBLAST + image generation
# Usage: sbatch [SBATCH_OPTIONS] o2/o2_banc_nblast_dataset.sh <dataset>
#
# Recommended SBATCH settings per dataset:
#   fafb:      -c 1  -t 0-96:00  -p medium  --mem-per-cpu=25G
#   fanc:      -c 1  -t 0-12:00  -p short   --mem-per-cpu=10G
#   hemibrain: -c 1  -t 0-96:00  -p medium  --mem-per-cpu=25G
#   manc:      -c 10 -t 0-96:00  -p medium  --mem-per-cpu=5G
#
#SBATCH -o jobs/banc_nblast_%j.out
#SBATCH -e jobs/banc_nblast_%j.err

DATASET="${1:-${DATASET}}"

if [ -z "$DATASET" ]; then
    echo "Usage: sbatch [OPTIONS] o2/o2_banc_nblast_dataset.sh <dataset>"
    echo "Datasets: fafb, fanc, hemibrain, manc"
    exit 1
fi

start=$(date +%s)
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

case "$DATASET" in
    fafb)
        echo "NBLASTing BANC vs FAFB"
        Rscript banc/nblast/banc-fafb-nblast.R
        echo "Creating 2D matching images for BANC-FAFB"
        Rscript banc/nblast/banc-fafb-nblast-images.R
        ;;
    fanc)
        echo "Transforming FANC meshes into BANC space"
        Rscript banc/transforms/banc-fanc-mesh-transform.R
        echo "Getting FANC metadata"
        Rscript fanc/fanc-meta.R
        echo "Transforming FANC skeletons into BANC space"
        Rscript banc/transforms/banc-fanc-skel-transform.R
        echo "NBLASTing BANC vs FANC"
        Rscript banc/nblast/banc-fanc-nblast.R
        echo "Creating 2D matching images for BANC-FANC"
        Rscript banc/nblast/banc-fanc-nblast-images.R
        ;;
    hemibrain)
        echo "NBLASTing BANC vs Hemibrain"
        Rscript banc/nblast/banc-hemibrain-nblast.R
        echo "Creating 2D matching images for BANC-Hemibrain"
        Rscript banc/nblast/banc-hemibrain-nblast-images.R
        ;;
    manc)
        echo "NBLASTing BANC vs MANC"
        Rscript banc/nblast/banc-manc-nblast.R
        echo "Creating 2D matching images for BANC-MANC"
        Rscript banc/nblast/banc-manc-nblast-images.R
        ;;
    malecns)
        echo "NBLASTing BANC vs maleCNS"
        Rscript banc/nblast/banc-malecns-nblast.R
        echo "Creating 2D matching images for BANC-maleCNS"
        Rscript banc/nblast/banc-malecns-nblast-images.R
        ;;
    lr)
        echo "Left-right mirror NBLAST (intra-BANC)"
        Rscript banc/nblast/banc-nblast-lr.R
        ;;
    *)
        echo "Unknown dataset: $DATASET"
        echo "Valid datasets: fafb, fanc, hemibrain, manc, malecns, lr"
        exit 1
        ;;
esac

end=$(date +%s)
echo "script completed in: $((end - start)) seconds"
