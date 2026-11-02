#!/bin/zsh
# SRC kick of 3.74 GeV (twice the realized default primordial kT on the jet side), default
# primordial kT, junction remnant, showers off and on. Used for the shower-survival numbers and the
# absolute-vs-ratio figure (compare_primkt.py, plot_absolute_vs_ratio.py).
# Usage: run_kick3p74.sh [rows per chunk] [chunks per sample]
# Writes $STUDY_ROOT/production_primkt/<sample>.csv (+ per-chunk logs and summaries).
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
NCHUNK_ROWS=${1:-1600000}
NCHUNKS=${2:-4}
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production_primkt
mkdir -p $OUT/chunks

# name | arguments. The SRC-off arm does not depend on kt; it is generated with it for symmetry.
SAMPLES=(
  "showers_off_unchanged_kick3p74|sample=unchanged src=0 showers=0 kt=3.74"
  "showers_off_altered_kick3p74|sample=altered src=1 showers=0 kt=3.74"
  "showers_on_unchanged_kick3p74|sample=unchanged src=0 showers=1 kt=3.74"
  "showers_on_altered_kick3p74|sample=altered src=1 showers=1 kt=3.74"
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
