#!/bin/bash
#SBATCH -J banc_synapse_prop_plot
#SBATCH -c 4
#SBATCH -t 0-04:00
#SBATCH -p priority
#SBATCH --mem=64G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_synapse_proportion_plot_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_synapse_proportion_plot_%j.err

# Standalone plot regeneration for banc_synapse_proportion_on_cell.{png,pdf}
# in BANC-project/figures/figure_1/links/supplement. Backup for when the
# foreground Rscript run is delayed or fails.

set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

start=$(date +%s)
echo "=== banc-synapse-proportion-plot @ $(date) on $(hostname) ==="
Rscript banc/legacy/banc-synapse-proportion-plot.R
end=$(date +%s)
echo "completed in $((end-start))s"
