#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-18:00
#SBATCH -p medium
#SBATCH --mem=150G
#SBATCH -o jobs/banc_assess_synapses_plot_%j.out
#SBATCH -e jobs/banc_assess_synapses_plot_%j.err

# Regenerate EDF 1d (banc_synapse_proportion_on_cell.png) with the corrected
# Y axis (proportion 0-1 instead of percent 0-100). Just the first section of
# banc/legacy/banc-assess-synapses.R — line 138 (ggsave) is the last line we
# need; we stop with quit() before the downstream analyses begin.

set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

Rscript -e '
  options(warn = 1)
  # Source the script but quit after the first plot is saved.
  # Use a sentinel to halt execution at the right line.
  lines <- readLines("banc/legacy/banc-assess-synapses.R")
  stop_idx <- grep("Examine manually reviewed synapse sample", lines)[1] - 2L
  if (is.na(stop_idx) || stop_idx <= 0L)
    stop("Could not find section divider; aborting.")
  tmp <- tempfile(fileext = ".R")
  writeLines(lines[seq_len(stop_idx)], tmp)
  cat(sprintf("Running first %d lines of banc-assess-synapses.R...\n", stop_idx))
  source(tmp, echo = FALSE)
  cat("\nPlot saved.\n")
'
