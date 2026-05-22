# -----------------------------------------------------------------
# Legacy exploratory analysis by a member of the Wilson Lab.
# Historical record only — not part of the bancpipeline release.
# -----------------------------------------------------------------

#!/bin/bash
#SBATCH -p short
#SBATCH -c 1                         # 1 CPU core
#SBATCH --mem=15G                    # 15GB RAM
#SBATCH -t 02:00:00                   # 2hr
#SBATCH --job-name="vscodetunnel"    
#SBATCH --mail-user=jfan@g.harvard.edu  # Email to which notifications will be sent

# Load required modules
module load gcc/9.2.0
module load miniconda3/23.1.0

# Activate your conda environment
source activate banc
#python batch_cascade.py
#python ../real_data_cascade_main.py --mask_A "MBON09" --mask_B "descending" --max_timestep 10
#python ../real_data_cascade_main.py --mask_A "mechanosensory" --mask_B "MBON" --mask_key_A "modality" --mask_key_B "cell_type" --max_timestep 10
# Keep the job running
sleep 2h