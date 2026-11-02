#!/bin/zsh
# Small event sample for leading_hadron_notebook.ipynb: Figure 3 DIS window, every struck quark flavour
# (struck=any), SRC off, with the Breit-frame four-momenta of all visible final-state particles
# (the jet-clustering inputs). Writes $STUDY_ROOT/notebook/events_showers_<on|off>.csv and
# particles_showers_<on|off>.csv. Needs the generator from build.sh.
# Usage: make_notebook_sample.sh [events] [showers 0|1]   (defaults 300000, 1; about 30 seconds)
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
N=${1:-300000}
SH=${2:-1}
NAME=showers_$([ $SH = 1 ] && echo on || echo off)
OUT=$STUDY_ROOT/notebook
mkdir -p $OUT
$STUDY_ROOT/bin/fig3_patched out=$OUT/events_$NAME.csv particles=$OUT/particles_$NAME.csv \
  sample=unchanged src=0 showers=$SH struck=any nacc=$N seed=20261101 > $OUT/generate_$NAME.log 2>&1
echo "wrote $OUT/events_$NAME.csv ($(($(wc -l < $OUT/events_$NAME.csv) - 1)) events)" \
  "and particles_$NAME.csv ($(($(wc -l < $OUT/particles_$NAME.csv) - 1)) particles)"
