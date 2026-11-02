#!/bin/zsh
# Samples for the charge-charge correlator (QQC) analyses: struck valence u, SRC kick 3.74 GeV,
# default primordial kT, junction remnant, SRC off/on x showers off/on, with one row per final-state
# charged hadron (qqc=).
#   production_qqc_utag/: only events passing the experimental u-tag (leading current-hemisphere
#     charged hadron is pi+ or K+ with z_h > 0.5); the chunk summaries count all selected DIS events.
#     Used by qqc_utag.py, qqc_chi.py, recoil_charge.py.
#   production_qqc/: all selected events (untagged comparison in qqc_utag.py).
# Usage: run_qqc.sh [tagged-sample rows per chunk] [untagged-sample rows per chunk] [chunks per sample]
# ("rows" = selected DIS events per chunk).
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
NROWS_TAG=${1:-2000000}
NROWS_UNTAG=${2:-500000}
NCHUNKS=${3:-4}
BIN=$STUDY_ROOT/bin/fig3_patched

for dir in production_qqc_utag production_qqc; do
  OUT=$STUDY_ROOT/$dir
  mkdir -p $OUT/chunks
  if [[ $dir == production_qqc_utag ]]; then extra=(utagonly=1); nrows=$NROWS_TAG; else extra=(); nrows=$NROWS_UNTAG; fi
  for sh in 0 1; do
    for src in 0 1; do
      name=showers_$([ $sh = 1 ] && echo on || echo off)_src$src
      arm=$([ $src = 1 ] && echo altered || echo unchanged)
      for c in $(seq 0 $((NCHUNKS - 1))); do
        $BIN out=$OUT/chunks/${name}_c$c.csv qqc=$OUT/chunks/${name}_qqc_c$c.csv sample=$arm src=$src \
          showers=$sh struck=u kt=3.74 $extra nacc=$nrows seed=$((20261028 + c)) \
          eidoffset=$((c * 10000000)) > $OUT/chunks/${name}_c$c.log 2>&1 &
      done
    done
  done
done
wait

for dir in production_qqc_utag production_qqc; do
  OUT=$STUDY_ROOT/$dir
  for f in $OUT/chunks/*_c0.csv; do
    [[ $f == *_qqc_c0.csv ]] && continue
    name=$(basename $f _c0.csv)
    for kind in "" _qqc; do
      head -1 $OUT/chunks/${name}${kind}_c0.csv > $OUT/${name}${kind}.csv
      for c in $(seq 0 $((NCHUNKS - 1))); do tail -n +2 $OUT/chunks/${name}${kind}_c$c.csv >> $OUT/${name}${kind}.csv; done
    done
    echo "$dir/$name: $(($(wc -l < $OUT/$name.csv) - 1)) events, $(($(wc -l < $OUT/${name}_qqc.csv) - 1)) hadrons"
  done
done
