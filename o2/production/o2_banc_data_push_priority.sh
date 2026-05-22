#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-08:00                        # 8h cap on priority (was 18h on medium)
#SBATCH -p priority
#SBATCH --mem=220G                        # bumped from 150G: job 40119383 hit
                                          # 153GB peak then OOM-killed in
                                          # Section 3 v3 NT-join (198M rows).
                                          # MaxMemPerNode=257G on priority.
#SBATCH -o jobs/banc_data_push_priority_%j.out
#SBATCH -e jobs/banc_data_push_priority_%j.err
###############################################################################
# Priority-partition copy of o2_banc_data_push.sh. Refreshes banc-data.R outputs
# (v3 + v2) after a SeaTable banc_meta update, so downstream priority-chained
# jobs (aggregate-influence, WB prep, WB align/NTAC) read fresh inputs.
#
# Same defensive output verification as the medium version — Rscript may exit
# non-zero from a reticulate/arrow destructor race after writing all outputs.
###############################################################################
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
