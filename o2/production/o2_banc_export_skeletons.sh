#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=24G
#SBATCH -o jobs/banc_export_skeletons_%j.out
#SBATCH -e jobs/banc_export_skeletons_%j.err

# Exports per-neuron SWCs for compiled_data/banc_<ver>/.
# Detailed skeleton (banc/swc/<id>.swc) preferred; falls back to L2
# (banc/l2/<id>.swc); SWC name suffix records which (_skeleton|_l2).
# Filter to valid banc_<ver> ids from banc_<ver>_meta.feather.
# Then rsyncs to gs://lee-lab_brain-and-nerve-cord-fly-connectome/.
#
# Designed to chain after o2_banc_data_push.sh so it sees the freshest meta.

set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

start=$(date +%s)
Rscript banc/share/banc-export-skeletons.R
end=$(date +%s)
echo "Total: $(( (end-start)/60 )) min"
