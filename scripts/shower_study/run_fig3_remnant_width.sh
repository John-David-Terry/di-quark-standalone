#!/bin/zsh
# Figure 3 rerun with the SRC pair kick tied to the remnant intrinsic kT: in SRC events the partner
# d's intrinsic-kT width is doubled, BeamRemnants:primordialKTremnant 0.4 -> 0.8 GeV, by adding a
# correlated Gaussian kick of per-component width sqrt(0.8^2 - 0.4^2) = 0.693 GeV, +k to the
# struck u and -k to the d (transverse momentum conserved; SRC:kTdist = 0). All other intrinsic-kT
# parameters stay at the Pythia defaults. SpaceShower:pTmaxMatch is set for every sample (it only
# matters with showers on; 0 is the Pythia default, 2 is Pythia's DIS example).
# Usage: run_fig3_remnant_width.sh [rows per chunk] [chunks per sample] [kick width GeV] [pTmaxMatch]
# Writes $STUDY_ROOT/production_primkt/showers_{off,on}_{unchanged,altered}_<tag>.csv with
# tag ktw<width>_m<pTmaxMatch> (default ktw0p693_m0), for compare_primkt.py and
# plot_src_radiation.py.
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
NCHUNK_ROWS=${1:-1600000}
NCHUNKS=${2:-4}
KW=${3:-0.693}
PTM=${4:-0}
TAG=ktw${KW/./p}_m$PTM
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production_primkt
mkdir -p $OUT/chunks

COMMON="kt=$KW ktdist=0 ptmaxmatch=$PTM"
SAMPLES=(
  "showers_off_unchanged_$TAG|sample=unchanged src=0 showers=0 $COMMON"
  "showers_off_altered_$TAG|sample=altered src=1 showers=0 $COMMON"
  "showers_on_unchanged_$TAG|sample=unchanged src=0 showers=1 $COMMON"
  "showers_on_altered_$TAG|sample=altered src=1 showers=1 $COMMON"
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
