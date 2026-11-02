#!/bin/zsh
# All analyses and figures of the shower study, from the samples of run_fig3_patched.sh,
# run_kick3p74.sh and run_qqc.sh. Results go to $STUDY_ROOT/results*/ (summary.txt + figures).
# Usage: run_analysis.sh            (Python: $PYTHON, default python3)
set -e
HERE=${0:A:h}
REPO=${HERE:h:h}
PY=${PYTHON:-python3}
export STUDY_ROOT=${STUDY_ROOT:-$HOME/Data/fig3_shower_study}
CAPTION=$'$x > 0.05$\n$4~\\mathrm{GeV} < Q < 5~\\mathrm{GeV}$\n$2.0 < \\eta_h^{\\mathrm{lab}} < 3.5$\n$z_{LC} > 0.1$\n$E_e \\times E_p = 5 \\times 41~\\mathrm{GeV}^2$'

# Figure 3 layout (Track A plot script, same arguments as AGENT_INSTRUCTIONS.md A5).
paper_figure() {  # <ok-row csv> <figure dir> <filename suffix>
  mkdir -p $2
  $PY $REPO/scripts/analysis/plot_jet_hadron_pi_pm_observables.py --csv $1 --figure-dir $2 \
    --layout combined --only-ok-rows --combined-only angle \
    --compare-pooled-altered-fractions 0,0.05,0.10 --pooled-balance-seed 42 \
    --combined-filename-suffix $3 --no-title --caption "$CAPTION" > $2/plot$3.log
}

echo "== Figure 3 rerun, paper kick (0.4 GeV), showers off vs on"
$PY $HERE/compare_fig3.py --write-plot-inputs
for sh in off on; do
  paper_figure $STUDY_ROOT/production/fig3_plot_input_showers_$sh.csv \
    $STUDY_ROOT/results/figures_showers_$sh _patched_showers_${sh}_antikt
done

echo "== SRC kick 3.74 GeV, showers off vs on"
KICK=(--tag kick3p74 --label "kick 3.74 GeV, mJoinJunction 1.0" --showers on
      --baseline kick3p74 --baseline-label "kick 3.74 GeV, mJoinJunction 1.0" --baseline-showers off)
$PY $HERE/compare_primkt.py $KICK --out-dir $STUDY_ROOT/results_showers_kick3p74
$PY $HERE/compare_primkt.py $KICK --write-plot-inputs
for prefix in "" showers_on_; do
  paper_figure $STUDY_ROOT/production_primkt/fig3_plot_input_${prefix}kick3p74_antikt.csv \
    $STUDY_ROOT/results_showers_kick3p74/figures _patched_${prefix}kick3p74_antikt
done
$PY $HERE/plot_absolute_vs_ratio.py

echo "== Charge-charge correlators (u-tagged)"
$PY $HERE/qqc_utag.py
$PY $HERE/qqc_chi.py
$PY $HERE/recoil_charge.py
echo "analysis done; results in $STUDY_ROOT"
