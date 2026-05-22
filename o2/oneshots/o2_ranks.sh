#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 0-04:00                        # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem-per-cpu=25G                 # Memory per core
#SBATCH -o /home/zaa827/bancpipeline/jobs/o2_ranks_%j.out # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/zaa827/bancpipeline/jobs/o2_ranks_%j.err # File to which STDERR will be written, including job ID (%j)

# data lives at:
# /n/data1/hms/neurobio/wilson/banc/connectivity/frankenbrain_v1.1_data.sqlite

echo "RUNNING NEURON RANKING ALGORITHM"
echo "usage: sbatch o2_ranks.sh -s olfactory"

start=`date +%s`

# Load O2 modules
echo "loading modules ..."
module purge
module load gcc/14.2.0
module load python/3.10.11

# Initialize SEED variable
SEED=""

# Parse command-line options
while getopts ":s:" opt; do
  case $opt in
    s)
      SEED="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if SEED was provided
if [ -z "$SEED" ]; then
  echo "Please provide a SEED using -s"
  exit 1
fi

# Run code
echo "running algorithm ..." # script runs algorithm on given seeds, saves in bancpipeline/data/ranks/...[seed].csv
virtualenv jupyter_connectome_venv
python3 bancpipeline/analysis/python/neuron_ranking.py "$SEED"

# Push to github
cd bancpipeline

# Configure Git to use a non-interactive merge strategy
git config pull.rebase false
git config merge.conflictstyle diff3

# Pull changes from the remote repository
git pull origin main

# Check if there are conflicts
if git diff --name-only --diff-filter=U | grep -q .; then
    # Auto-resolve conflicts by favoring the remote changes
    git checkout --theirs .
    git add .
    git commit -m "Auto-resolved conflicts favoring remote changes"
fi

# Stage all changes
git add .

# Commit changes (if any)
if git diff --staged --quiet; then
    echo "No changes to commit"
else
    git commit -m "Auto-commit from O2: $(date)"
fi

# Push changes to the remote repository
git push origin main

# Report
end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime
