#!/bin/zsh
# Showers-off SRC on/off samples at reduced primordial kT, plus the default-kT pair with the
# same seeds (identical events to run_fig3_patched.sh, now with kT / diquark-merge columns).
# Usage: run_primkt_study.sh [primkt sigma GeV] [rows per chunk] [chunks per sample]
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
PRIMKT=${1:-0.13}
NCHUNK_ROWS=${2:-1600000}
NCHUNKS=${3:-4}
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production_primkt
TAG=primkt${PRIMKT/./p}
mkdir -p $OUT/chunks

SAMPLES=(
  "showers_off_unchanged_defaultkt|sample=unchanged src=0 showers=0"
  "showers_off_altered_defaultkt|sample=altered src=1 showers=0"
  "showers_off_unchanged_$TAG|sample=unchanged src=0 showers=0 primkt=$PRIMKT"
  "showers_off_altered_$TAG|sample=altered src=1 showers=0 primkt=$PRIMKT"
)

for entry in $SAMPLES; do
  name=${entry%%|*}
  args=(${=entry#*|})
  for c in $(seq 0 $((NCHUNKS - 1))); do
    $BIN out=$OUT/chunks/${name}_c$c.csv $args nacc=$NCHUNK_ROWS seed=$((20261028 + c)) \
      eidoffset=$((c * 10000000)) > $OUT/chunks/${name}_c$c.log 2>&1 &
  done
done
wait

for entry in $SAMPLES; do
  name=${entry%%|*}
  head -1 $OUT/chunks/${name}_c0.csv > $OUT/$name.csv
  for c in $(seq 0 $((NCHUNKS - 1))); do tail -n +2 $OUT/chunks/${name}_c$c.csv >> $OUT/$name.csv; done
  echo "$name: $(($(wc -l < $OUT/$name.csv) - 1)) rows"
done
