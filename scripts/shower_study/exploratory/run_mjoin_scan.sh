#!/bin/zsh
# Scan FragmentationSystems:mJoinJunction (junction-quark -> diquark merging threshold) at fixed
# primordial kT, showers off, SRC on/off. Everything else as in fig3_patched defaults.
# Usage: run_mjoin_scan.sh [primkt sigma] [rows per chunk] [chunks per sample] [thresholds...]
set -e
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
PRIMKT=${1:-0.13}
NCHUNK_ROWS=${2:-800000}
NCHUNKS=${3:-2}
shift 3 2>/dev/null || true
THRESHOLDS=(${@:-1.0 0.5 0.2 0.0})
BIN=$STUDY_ROOT/bin/fig3_patched
OUT=$STUDY_ROOT/production_mjoin
mkdir -p $OUT/chunks

for mj in $THRESHOLDS; do
  for arm in unchanged altered; do
    src=$([ $arm = altered ] && echo 1 || echo 0)
    name=showers_off_${arm}_primkt${PRIMKT/./p}_mjoin${mj/./p}
    for c in $(seq 0 $((NCHUNKS - 1))); do
      $BIN out=$OUT/chunks/${name}_c$c.csv sample=$arm src=$src showers=0 primkt=$PRIMKT \
        mjoinjunction=$mj nacc=$NCHUNK_ROWS seed=$((20261028 + c)) \
        eidoffset=$((c * 10000000)) > $OUT/chunks/${name}_c$c.log 2>&1 &
    done
  done
done
wait

for f in $OUT/chunks/*_c0.csv; do
  name=$(basename $f _c0.csv)
  head -1 $f > $OUT/$name.csv
  for c in $(seq 0 $((NCHUNKS - 1))); do tail -n +2 $OUT/chunks/${name}_c$c.csv >> $OUT/$name.csv; done
  echo "$name: $(($(wc -l < $OUT/$name.csv) - 1)) rows"
done
