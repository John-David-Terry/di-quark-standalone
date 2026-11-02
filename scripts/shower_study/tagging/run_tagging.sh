#!/bin/zsh
# Samples for the leading-hadron flavour tag (Track D): struck u (struckid=2) and struck d
# (struckid=1), Figure 3 DIS window, SRC off, junction remnant, default primordial kT, showers off
# and on. Writes $STUDY_ROOT/production_leadpi/showers_{off,on}_struck{u,d}.csv for leadpi_tag.py.
# Usage: run_tagging.sh [events per chunk] [chunks per sample]   (defaults 250000, 4)
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
NROWS=${1:-250000}
NCHUNKS=${2:-4}
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production_leadpi
mkdir -p $OUT/chunks

for sh in 0 1; do
  for f in u d; do
    name=showers_$([ $sh = 1 ] && echo on || echo off)_struck$f
    id=$([ $f = u ] && echo 2 || echo 1)
    for c in $(seq 0 $((NCHUNKS - 1))); do
      $BIN out=$OUT/chunks/${name}_c$c.csv sample=unchanged src=0 showers=$sh struck=any struckid=$id \
        nacc=$NROWS seed=$((20261030 + c)) eidoffset=$((c * 10000000)) > $OUT/chunks/${name}_c$c.log 2>&1 &
    done
  done
done
wait

for sh in off on; do
  for f in u d; do
    name=showers_${sh}_struck$f
    head -1 $OUT/chunks/${name}_c0.csv > $OUT/$name.csv
    for c in $(seq 0 $((NCHUNKS - 1))); do tail -n +2 $OUT/chunks/${name}_c$c.csv >> $OUT/$name.csv; done
    echo "production_leadpi/$name: $(($(wc -l < $OUT/$name.csv) - 1)) events"
  done
done
