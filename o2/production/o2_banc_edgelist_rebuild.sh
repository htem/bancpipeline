#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-15:00                        # 15h walltime on medium (each stage
                                          # ~30m-4h, 4 stages, leaves headroom)
#SBATCH -p medium
#SBATCH --mem=250G                        # data-push step OOMed at 150G last run
#SBATCH -o jobs/banc_edgelist_rebuild_%j.out
#SBATCH -e jobs/banc_edgelist_rebuild_%j.err
###############################################################################
# 2026-05-15 — rebuild simple + split + enriched outputs for v2 and v3 after
# tasks.md §18 filter overhaul:
#   - banc.size.threshold = 5 everywhere
#   - glia + trachea kept in root.ids
#   - v3 CSV reader now uses 17-col schema (fixes empty v3 simple edgelist)
#
# Stages:
#   1. banc-calculate-connectivity.R --source v2  → banc_888_synapses_v2.parquet
#                                                   + banc_888_edgelist_simple_v2.feather
#   2. banc-calculate-connectivity.R --source v3  → banc_888_edgelist_simple_v3.feather
#                                                   (was 802 B / 0 rows before the
#                                                   17-col schema fix)
#   3. banc-data.R --source v3                    → banc_888_synapses_v3_enriched.parquet
#                                                   + banc_888_edgelist_split_v3.feather
#                                                   + copies simple_v3
#   4. banc-data.R --source v2                    → banc_888_synapses_v2_enriched.parquet
#                                                   + banc_888_edgelist_split_v2.feather
#                                                   + copies simple_v2
#
# Does NOT re-run banc-calculate-neuropil-inclusion.R — that's parked in
# tasks.md §18 "Deferred" (would re-classify ~213M synapses to fix v2's
# missing canonical 'outside' labels; ~hours of alpha-shape work).
###############################################################################
set -uo pipefail
trap 'echo "[$(date +%H:%M:%S)] caught signal — continuing if outputs present"' ERR

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

CONN_DIR=/n/data1/hms/neurobio/wilson/banc/connectivity
OUT_DIR=/n/data1/hms/neurobio/wilson/connectomes/banc/banc_888

start=$(date +%s)

check_nonempty() {
  local f="$1"
  local min_bytes="${2:-1000}"  # default: must be at least 1 KB
  if [[ ! -s "$f" ]]; then
    echo "  MISSING: $f" >&2
    return 1
  fi
  local sz
  sz=$(stat -c %s "$f")
  if (( sz < min_bytes )); then
    echo "  TOO SMALL: $f ($sz bytes < $min_bytes)" >&2
    return 1
  fi
  local mtime
  mtime=$(stat -c %Y "$f")
  if (( mtime < start )); then
    echo "  STALE: $f (mtime < job start)" >&2
    return 1
  fi
  echo "  OK: $f ($sz bytes)"
  return 0
}

echo "=== Stage 1: banc-calculate-connectivity.R --source v2 ==="
t1=$(date +%s)
Rscript banc/metrics/banc-calculate-connectivity.R --source v2 \
  || echo "  v2 connectivity Rscript exit $? — checking outputs..."
check_nonempty "$CONN_DIR/banc_888_synapses_v2.parquet" 1000000 || exit 2
check_nonempty "$CONN_DIR/banc_888_edgelist_simple_v2.feather" 1000000 || exit 2
echo "  Stage 1 done in $(( ($(date +%s) - t1) / 60 )) min"

echo "=== Stage 2: banc-calculate-connectivity.R --source v3 ==="
t2=$(date +%s)
Rscript banc/metrics/banc-calculate-connectivity.R --source v3 \
  || echo "  v3 connectivity Rscript exit $? — checking outputs..."
# v3 simple edgelist must be substantially populated (was 802 B before fix);
# expect tens of MB minimum after the schema fix.
check_nonempty "$CONN_DIR/banc_888_edgelist_simple_v3.feather" 10000000 || exit 2
echo "  Stage 2 done in $(( ($(date +%s) - t2) / 60 )) min"

echo "=== Stage 3: banc-data.R --source v3 ==="
t3=$(date +%s)
Rscript banc/share/banc-data.R --source v3 \
  || echo "  v3 banc-data.R exit $? — checking outputs..."
check_nonempty "$OUT_DIR/banc_888_synapses_v3_enriched.parquet" 1000000000 || exit 2
check_nonempty "$OUT_DIR/banc_888_edgelist_simple_v3.feather" 10000000 || exit 2
check_nonempty "$OUT_DIR/banc_888_edgelist_split_v3.feather" 100000000 || exit 2
echo "  Stage 3 done in $(( ($(date +%s) - t3) / 60 )) min"

echo "=== Stage 4: banc-data.R --source v2 ==="
t4=$(date +%s)
Rscript banc/share/banc-data.R --source v2 \
  || echo "  v2 banc-data.R exit $? — checking outputs..."
check_nonempty "$OUT_DIR/banc_888_synapses_v2_enriched.parquet" 1000000000 || exit 2
check_nonempty "$OUT_DIR/banc_888_edgelist_simple_v2.feather" 1000000 || exit 2
check_nonempty "$OUT_DIR/banc_888_edgelist_split_v2.feather" 100000000 || exit 2
echo "  Stage 4 done in $(( ($(date +%s) - t4) / 60 )) min"

end=$(date +%s)
echo "Total: $(( (end-start)/60 )) min"
