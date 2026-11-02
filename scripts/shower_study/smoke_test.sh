#!/bin/zsh
# End-to-end check with tiny samples (a few minutes): build the driver, generate every sample at
# toy size, run every analysis. Checks that the chain runs; the numbers are not meaningful.
# Usage: smoke_test.sh [scratch dir]   (default: a new temporary directory; needs PYTHIA8_DIR built)
set -e
HERE=${0:A:h}
export STUDY_ROOT=${1:-$(mktemp -d)/fig3_shower_study_smoke}
mkdir -p $STUDY_ROOT
echo "smoke test in $STUDY_ROOT"
$HERE/build.sh
$HERE/run_fig3_patched.sh 4000 2
$HERE/run_kick3p74.sh 4000 2
$HERE/run_qqc.sh 20000 4000 2
$HERE/run_fig3_remnant_width.sh 10000 2 0.693 0
$HERE/run_analysis.sh
$HERE/run_analysis_remnant_width.sh ktw0p693_m0
echo "smoke test passed; outputs in $STUDY_ROOT"
