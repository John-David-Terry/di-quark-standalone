#!/bin/zsh
# Compile the generator (fig3_patched.cc) and the kick validation (pythia/validate_kick.cc) against
# the SRC-patched Pythia 8.311 (pythia/setup_pythia.sh). Binaries go to the data area, never into
# the repo.
set -e
PYTHIA8_DIR=${PYTHIA8_DIR:-$HOME/Data/pythia8311-SRCpatch}
STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
HERE=${0:A:h}
mkdir -p $STUDY_ROOT/bin
for src in $HERE/fig3_patched.cc $HERE/pythia/validate_kick.cc; do
  exe=$STUDY_ROOT/bin/${src:t:r}
  g++ -O2 -std=c++11 -W -Wall -I$PYTHIA8_DIR/include \
    -DXMLDIR="\"$PYTHIA8_DIR/share/Pythia8/xmldoc\"" \
    $src -o $exe \
    -L$PYTHIA8_DIR/lib -Wl,-rpath,$PYTHIA8_DIR/lib -lpythia8 -ldl
  echo "built $exe"
done
