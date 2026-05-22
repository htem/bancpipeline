#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-18:00
#SBATCH -p medium
#SBATCH --mem=200G
#SBATCH -o jobs/banc_data_push_v3only_%j.out
#SBATCH -e jobs/banc_data_push_v3only_%j.err

# One-shot v3 rebuild after banc-data.R was patched to:
#  - persist v3 spatial-NN pre_label/post_label into the enriched parquet
#    (Section 3 v3 branch)
#  - always build v3 simple edgelist from banc.syns (Section 4: skip stale
#    source feather for v3)
#  - simplify Section 5 v3 to read persisted labels (no NN re-compute)
# v2 outputs are untouched (separate sbatch handles those on the regular cron).

set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

OUT_DIR=/n/data1/hms/neurobio/wilson/connectomes/banc/banc_888
start=$(date +%s)

echo "=== banc-data.R --source v3 (post-patch) ==="
Rscript banc/share/banc-data.R --source v3 || echo "  Rscript exit $? — checking outputs..."

missing=0
for f in "banc_888_synapses_v3_enriched.parquet" \
         "banc_888_edgelist_simple_v3.feather" \
         "banc_888_edgelist_split_v3.feather"; do
  if [[ ! -s "$OUT_DIR/$f" ]]; then
    echo "  MISSING: $OUT_DIR/$f" >&2; missing=$((missing+1))
  else
    mtime=$(stat -c %Y "$OUT_DIR/$f")
    if (( mtime < start )); then
      echo "  STALE: $OUT_DIR/$f (mtime $mtime < start $start)" >&2; missing=$((missing+1))
    fi
  fi
done

end=$(date +%s)
echo "Total: $(( (end-start)/60 )) min"
[[ $missing -eq 0 ]] || exit 2
echo "v3 outputs verified"
