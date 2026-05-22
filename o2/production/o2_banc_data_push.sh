#!/bin/bash
#SBATCH -c 4                              # Cores
#SBATCH -t 0-18:00                        # Runtime — bumped from 4h after job 38760183
                                          # timed out reading split CSVs in v2 Section 5
                                          # (v3 + v2 each need ~2-3h for Section 5 alone).
                                          # O2 policy: <12h → short, >12h → medium.
#SBATCH -p medium                         # medium needed for >12h walltime band
#SBATCH --mem=150G                        # Total memory
#SBATCH -o jobs/banc_data_push_%j.out
#SBATCH -e jobs/banc_data_push_%j.err

# Runs banc/share/banc-data.R with both --source v3 and --source v2,
# producing all v888 versioned outputs in
# /n/data1/hms/neurobio/wilson/connectomes/banc/banc_888/
# and rsyncing to gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/
#
# Resilient to R cleanup aborts: banc-data.R can SIGABRT after writing all
# outputs (reticulate/arrow destructor race). Each run's success is judged
# by the presence of expected output files, not Rscript's exit code.

set -uo pipefail
trap 'echo "[$(date +%H:%M:%S)] caught signal — continuing if outputs present"' ERR

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

OUT_DIR=/n/data1/hms/neurobio/wilson/connectomes/banc/banc_888

verify_outputs() {
  local source="$1"
  local missing=0
  for f in banc_888_meta.feather \
           banc_888_metrics.feather \
           "banc_888_synapses_${source}_enriched.parquet" \
           "banc_888_edgelist_simple_${source}.feather"; do
    if [[ ! -s "$OUT_DIR/$f" ]]; then
      echo "  MISSING: $OUT_DIR/$f" >&2
      missing=$((missing+1))
    else
      local mtime
      mtime=$(stat -c %Y "$OUT_DIR/$f")
      if (( mtime < start )); then
        echo "  STALE (mtime < job start): $OUT_DIR/$f" >&2
        missing=$((missing+1))
      fi
    fi
  done
  return $missing
}

start=$(date +%s)

echo "=== banc-data.R --source v3 ==="
Rscript banc/share/banc-data.R --source v3 || echo "  v3 Rscript exit $? — checking outputs..."
if ! verify_outputs v3; then
  echo "FATAL: v3 outputs missing or stale" >&2
  exit 2
fi
echo "  v3 outputs verified"

echo "=== banc-data.R --source v2 ==="
Rscript banc/share/banc-data.R --source v2 || echo "  v2 Rscript exit $? — checking outputs..."
if ! verify_outputs v2; then
  echo "FATAL: v2 outputs missing or stale" >&2
  exit 2
fi
echo "  v2 outputs verified"

end=$(date +%s)
echo "Total: $(( (end-start)/60 )) min"
