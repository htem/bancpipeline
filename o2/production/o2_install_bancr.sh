#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o /home/ab714/bancpipeline/jobs/install_bancr_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/install_bancr_%j.err

set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

echo "=== Installing natverse/bancr ==="
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
old <- packageVersion("bancr")
cat("Before:", as.character(old), "\n")
remotes::install_github("natverse/bancr", upgrade = "never", quiet = FALSE)
new <- packageVersion("bancr")
cat("After:", as.character(new), "\n")
# Smoke-test banc_nuclei
library(bancr)
nuc <- tryCatch(banc_nuclei(rawcoords = TRUE, table = "both"),
                error = function(e) { cat("banc_nuclei error:", conditionMessage(e), "\n"); NULL })
if (!is.null(nuc)) cat("banc_nuclei OK, rows:", nrow(nuc), "\n")
'
echo "=== Done ==="
