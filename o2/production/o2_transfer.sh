#!/bin/bash

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh

# Connect to transfer node
# ssh $(whoami)@transfer.rc.hms.harvard.edu 'bash -s; exit' < o2_banc_results.sh
# ssh $(whoami)@transfer.rc.hms.harvard.edu 'Rscript -e "$HOME/bancpipeline/banc/meta/banc-sort-folders.R" && exit'
Rscript -e 'banc/meta/banc-sort-folders.R'