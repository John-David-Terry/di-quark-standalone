#!/bin/zsh
# Produce the Figure 3 shower study samples with fig3_patched (build.sh first).
# Usage: run_fig3_patched.sh [rows per chunk] [chunks per sample]
# Writes $STUDY_ROOT/production/<sample>.csv (+ per-chunk summaries) outside the repo.
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
NCHUNK_ROWS=${1:-1600000}
NCHUNKS=${2:-4}
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production
mkdir -p $OUT/chunks

# name | arguments. unchanged/altered pairs differ only in src; off/on pairs only in showers.
SAMPLES=(
  "showers_off_unchanged|sample=unchanged src=0 showers=0"
  "showers_off_altered|sample=altered src=1 showers=0"
  "showers_on_unchanged|sample=unchanged src=0 showers=1"
  "showers_on_altered|sample=altered src=1 showers=1"
  "showers_off_diquark_control|sample=unchanged src=0 showers=0 remnant=diquark"
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
