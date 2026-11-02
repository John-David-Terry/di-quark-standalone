#!/bin/zsh
# All analyses of the doubled remnant-width SRC kick samples (run_fig3_remnant_width.sh):
#   Figure 3 mixtures, flavour-tagged shifts and shower survival (compare_primkt.py), with the
#     Figure 3-layout plots;
#   the two-panel 1/sigma dsigma/dphi_hJ and R figure (plot_absolute_vs_ratio.py);
#   pi- distributions and SRC / no-SRC ratios without and with radiation (plot_src_radiation.py);
#   the charge ratio R_pi = pi-/pi+ and R_pi^SRC / R_pi^noSRC, overall and in bins of the
#     Breit-frame Delta eta_hJ: < 4, 4-5.5, >= 5.5 and 1-2.5, 2.5-4, 4-5.5 (charge_ratio.py).
# Results go to $STUDY_ROOT/results*_<tag>/ (and results_absolute_vs_ratio/).
# Usage: run_analysis_remnant_width.sh [tag]   (default ktw0p693_m0; Python: $PYTHON, default python3)
set -e
HERE=${0:A:h}
REPO=${HERE:h:h}
PY=${PYTHON:-python3}
export STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
TAG=${1:-ktw0p693_m0}
CAPTION=$'$x > 0.05$\n$4~\\mathrm{GeV} < Q < 5~\\mathrm{GeV}$\n$2.0 < \\eta_h^{\\mathrm{lab}} < 3.5$\n$z_{LC} > 0.1$\n$E_e \\times E_p = 5 \\times 41~\\mathrm{GeV}^2$'

# Figure 3 layout (Track A plot script, same arguments as AGENT_INSTRUCTIONS.md A5).
paper_figure() {  # <ok-row csv> <figure dir> <filename suffix>
  mkdir -p $2
  $PY $REPO/scripts/analysis/plot_jet_hadron_pi_pm_observables.py --csv $1 --figure-dir $2 \
    --layout combined --only-ok-rows --combined-only angle \
    --compare-pooled-altered-fractions 0,0.05,0.10 --pooled-balance-seed 42 \
    --combined-filename-suffix $3 --no-title --caption "$CAPTION" > $2/plot$3.log
}

ARGS=(--tag $TAG --label "d width 0.4 -> 0.8 GeV" --showers on
      --baseline $TAG --baseline-label "d width 0.4 -> 0.8 GeV" --baseline-showers off)
$PY $HERE/compare_primkt.py $ARGS --out-dir $STUDY_ROOT/results_remnant_width_$TAG
$PY $HERE/compare_primkt.py $ARGS --write-plot-inputs
for prefix in "" showers_on_; do
  paper_figure $STUDY_ROOT/production_primkt/fig3_plot_input_${prefix}${TAG}_antikt.csv \
    $STUDY_ROOT/results_remnant_width_$TAG/figures _patched_${prefix}${TAG}_antikt
done
$PY $HERE/plot_absolute_vs_ratio.py --tag $TAG
$PY $HERE/plot_src_radiation.py --tag $TAG
$PY $HERE/charge_ratio.py --tag $TAG
$PY $HERE/charge_ratio.py --tag $TAG --edges 1,2.5,4,5.5
echo "remnant-width analysis done; results in $STUDY_ROOT"
