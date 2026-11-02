#!/bin/zsh
# Leading-hadron flavour tag (Track D), from the samples of run_tagging.sh: the u-tag "leading
# charged hadron of the Breit current hemisphere is a pi+ with SIDIS z_h > cut" at cut 0.5
# (and its jet version), and purity / efficiency scans in the z_h cut for pi+, pi+ or K+, and any
# positive hadron. Results go to $STUDY_ROOT/results_tagging/.
# Usage: run_analysis_tagging.sh            (Python: $PYTHON, default python3)
set -e
HERE=${0:A:h}
PY=${PYTHON:-python3}
export STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
$PY $HERE/leadpi_tag.py 0.5
$PY $HERE/leadpi_tag.py scan
$PY $HERE/leadpi_tag.py scan_positive
echo "tagging analysis done; results in $STUDY_ROOT/results_tagging"
