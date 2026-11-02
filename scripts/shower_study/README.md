# Track C — Figure 3 and charge correlators with parton showers

This track reruns the Figure 3 observable R_{π⁻/π⁺}(φ_hJ) with the SRC kick applied **inside**
Pythia's beam-remnant step, so that initial- and final-state radiation (ISR/FSR) act after the kick
and before hadronization. It answers the question of how much perturbative radiation washes out the
SRC signal, and adds the u-tagged charge-charge correlators (QQC) built on the same samples.

The Track A pipeline (`scripts/analysis/`) is not modified. Track C only imports its mixture and
ratio functions (read-only), so the 0/5/10% mixtures and R(φ_hJ) are computed with the original
code, and the Figure 3-style plots use the original plot script.

## How this differs from Track A

| | Track A (Figure 3) | Track C (this folder) |
|---|---|---|
| SRC kick | added by hand to the struck u and a split-off d, then partons re-hadronized; no showers possible | inside Pythia's beam-remnant step (`pythia/SRC_patch_v2.diff`); one `pythia.next()` does kick → shower → hadronize |
| Diquark breakup | surgery v2: (ud) → d + (us) + s̄ | stock junction remnant: `BeamRemnants:remnantMode = 1`, `beamJunction = on` |
| Kick | 0.4·sin θ along Breit x, θ ~ U[0, π] | same magnitude law, random azimuth (`SRC:kTdist = 2`), added to the primordial kT (`SRC:kickMode = 1`) |
| Jet | struck outgoing parton | Breit-frame e⁺e⁻ anti-kT (`ee_genkt`, p = −1, R = 1.0); the struck-parton φ_hJ is also written |
| Pythia | 8.31x Python bindings | patched Pythia 8.311 (C++) |

Unchanged from Figure 3: 5 × 41 GeV, 16 < Q² < 25 GeV², 0.05 < x ≤ 0.495, struck valence u,
leading charged target pion (Breit p_z > 0, proton-forward lab 2 < η < 3.5, z_LC ≥ 0.1), 40 bins
on [0, π], 0/5/10% mixtures with π⁺ from the SRC-off arm. Showers on/off differ only in
`PartonLevel:ISR/FSR` (QED showers off, `SpaceShower:dipoleRecoil = on` in both).

The patch kicks only events whose initiator is a valence u and whose remnant contains a separate d
(the partner). Every event of an SRC-on ("altered") sample is kicked; SRC fractions are set by the
0/5/10% mixtures (Figure 3) or compared as 0% vs 100% (QQC).

## Requirements

- macOS or Linux, a C++11 compiler, `zsh`, `curl`, `make`.
- Python 3 with `numpy`, `pandas`, `matplotlib`, `pyarrow` (`pip install -r scripts/shower_study/requirements.txt`).
  Tested with Python 3.12, numpy 2.4, pandas 3.0, matplotlib 3.10, pyarrow 24.
- Optional: a LaTeX installation (`latex`, `dvipng`) for LaTeX-rendered labels in the QQC figures;
  without it matplotlib's own math rendering is used.
- About 25 GB of disk for the full production (half of it is per-chunk CSVs that can be deleted
  after merging; keep the `*.summary.json` files), and about 5 CPU-hours (≈ 30–45 min on 16 cores).

Keep all data and the Pythia build **outside the repository** and outside cloud-synced folders,
as for Tracks A and B.

## Running

From the repository root:

```bash
export PYTHIA8_DIR=$HOME/Data/pythia8311-SRCpatch   # patched Pythia (outside the repo)
export STUDY_ROOT=$HOME/Data/fig3_shower_study      # all samples and results (outside the repo)
export PYTHON=python3                               # interpreter with the packages above

scripts/shower_study/pythia/setup_pythia.sh         # download 8.311 (29 MB), patch, build (~5 CPU-minutes)
scripts/shower_study/smoke_test.sh                  # optional: toy-size end-to-end check (a few minutes)

scripts/shower_study/build.sh                       # compile the generator and the kick validation
scripts/shower_study/run_fig3_patched.sh            # Figure 3 rerun, 0.4 GeV kick: 5 samples x 6.4M events
scripts/shower_study/run_kick3p74.sh                # 3.74 GeV kick: 4 samples x 6.4M events
scripts/shower_study/run_qqc.sh                     # QQC: 4 x 8M DIS events (u-tagged) + 4 x 2M (untagged)
scripts/shower_study/run_analysis.sh                # all summaries and figures of the samples above
```

Each production script takes the events per chunk and the number of chunks as arguments (defaults
reproduce the numbers below) and runs all chunks in parallel. Seeds are `20261028 + chunk` and
event ids are offset by `chunk × 10⁷`, so the same arguments give the same events on the same
platform. Every chunk writes `<chunk>.csv.summary.json` with all settings and counters.

SRC kick from the doubled remnant intrinsic kT (radiation-survival study, `AGENT_INSTRUCTIONS.md`
§C6): in SRC events the d's intrinsic-kT width goes from 0.4 to 0.8 GeV through a correlated
Gaussian kick of per-component width 0.693 GeV, +k on the struck u and −k on the d, with
`SpaceShower:pTmaxMatch = 0` everywhere. It has its own production and analysis steps:

```bash
scripts/shower_study/run_fig3_remnant_width.sh 4000000 8 0.693 0   # 4 samples x 32M events (≈ 4 CPU-hours)
scripts/shower_study/run_analysis_remnant_width.sh ktw0p693_m0     # all its analyses (≈ 10 minutes)
```

Optional comparison with the original Figure 3: after producing the Track A observables CSV
(`AGENT_INSTRUCTIONS.md` §A4), run
`$PYTHON scripts/shower_study/compact_original_fig3.py --csv <Track A jet_hadron_pi_pm_observables.csv>`
before `run_analysis.sh`; `compare_fig3.py` then adds the original curves.

## Outputs and reference results

All under `$STUDY_ROOT`. Each results folder has a `summary.txt` with every number below.

| Folder | Content |
|---|---|
| `results/` | Figure 3 rerun at the paper's 0.4 GeV kick, showers off vs on: R(φ_hJ) mixtures and cos-modulation fits (`compare_R_phi_hJ`), Figure 3-layout plots in `figures_showers_{off,on}/` |
| `results_showers_kick3p74/` | 3.74 GeV kick, showers off vs on: SRC shifts of ⟨cos φ_hJ⟩, R modulation, shower survival; Figure 3-layout plots in `figures/` |
| `results_absolute_vs_ratio/` | Two-panel figure: radiation changes the absolute φ_hJ shape, the ratio keeps the SRC signal |
| `results_qqc_utag/` | u-tagged azimuthal QQC: raw vs per-pair, tag-hadron inclusion variants, untagged comparison |
| `results_qqc_chi/` | Opening-angle QQC vs z = (1 + cos χ)/2 (Breit frame; γ*p c.m. frame as a diagnostic) |
| `results_recoil_charge/` | p_T-weighted recoil charge Q_rec vs tagging-hadron p_T |
| `results_remnant_width_<tag>/` | Doubled remnant-width SRC kick: Figure 3 mixtures, flavour-tagged shifts, shower survival; Figure 3-layout plots in `figures/` |
| `results_src_radiation_<tag>/` | π⁻ 1/σ dσ/dφ_hJ for no SRC / SRC × no radiation / radiation, and the SRC / no-SRC ratios with their radiation survival |
| `results_charge_ratio_<tag>/` | Charge ratio R_π = π⁻/π⁺ for the same four cases, R_π^SRC / R_π^noSRC, the Breit-frame Δη_hJ distribution, and the same ratios in three Δη_hJ bins |

Reference values from the full production (statistical errors):

| Quantity | Showers off | Showers on |
|---|---|---|
| Figure 3 rerun, 0.4 GeV kick, 10% mixture: R cos modulation (anti-kT) | −0.002 ± 0.003 | −0.002 ± 0.003 |
| 3.74 GeV kick, 10% mixture: R cos modulation, anti-kT / struck-parton jet | −0.031 ± 0.003 / −0.069 ± 0.004 | −0.027 / −0.061 |
| 3.74 GeV kick: flavour-tagged shift of ⟨cos φ_hJ⟩ (anti-kT) | −0.087 ± 0.001 | −0.081 ± 0.001 (survival 0.93 ± 0.02) |
| u-tagged QQC, tagging hadron × forward target, Δφ > 3π/4: SRC change per pair | −0.124 ± 0.002 | −0.121 ± 0.002 |
| Opening angle, 25% most back-to-back pairs: SRC change per pair, Breit / γ*p c.m. | −0.040 ± 0.004 / −0.132 ± 0.003 | −0.033 / −0.127 |
| Recoil charge Q_rec, full target hemisphere, SRC off / SRC on | +0.337 / +0.154 (± 0.001) | +0.328 / +0.159 |

For comparison, the original Figure 3 (Track A, no showers) gives a 10% R modulation of
−0.120 ± 0.004.

## Files

| File | Purpose |
|---|---|
| `pythia/setup_pythia.sh` | Download Pythia 8.311, verify the three files the patch touches (`pristine_8311.sha256`), apply `SRC_patch_v2.diff`, build |
| `pythia/SRC_patch_v2.diff` | The SRC patch (`patch -p1` in the Pythia source tree). Settings `SRC:on, splitDiquarkDIS, partnerId, kTdist, kT, kbarWidth, kickMode` are registered by the driver |
| `pythia/validate_kick.cc` | Kick validation: `validate_kick <primordial kT 0=zero,1=default> <showers 0/1> <src 0/1> <nEvents>` |
| `fig3_patched.cc` | Generator, event selection, Breit anti-kT jets, observables; options are listed at the top of the file (including `ktdist=` for the kick law and `ptmaxmatch=` for `SpaceShower:pTmaxMatch`) |
| `build.sh` | Compile `fig3_patched` and `validate_kick` into `$STUDY_ROOT/bin` |
| `run_fig3_patched.sh`, `run_kick3p74.sh`, `run_qqc.sh` | Production |
| `run_fig3_remnant_width.sh` | Production for the doubled remnant-width SRC kick: `[events per chunk] [chunks] [kick width] [pTmaxMatch]`, tag `ktw<width>_m<pTmaxMatch>` |
| `run_analysis.sh`, `smoke_test.sh` | All analyses; toy-size end-to-end check (covers every production and analysis script) |
| `run_analysis_remnant_width.sh` | All analyses of the doubled remnant-width SRC samples (`[tag]`, default `ktw0p693_m0`) |
| `compare_fig3.py` | Figure 3 rerun vs showers (and vs the original, if present) |
| `compare_primkt.py` | Compare two sample pairs in `production_primkt/` (used for the 3.74 GeV shower survival) |
| `plot_absolute_vs_ratio.py` | Absolute vs ratio figure (`--tag` picks the sample pair, default `kick3p74`; `--full` needs the exploratory samples) |
| `plot_src_radiation.py` | π⁻ distributions and SRC / no-SRC ratios with and without radiation (`--tag`, `--fraction`) |
| `charge_ratio.py` | Charge ratio R_π and its SRC / no-SRC double ratio, overall and in bins of the Breit-frame separation Δη_hJ = \|η_π − η_J\| (`--tag`; `--edges` all bin boundaries, default 0,4,5.5,inf, e.g. 1,2.5,4,5.5; `--sep dy` for rapidity) |
| `qqc_utag.py`, `qqc_chi.py`, `recoil_charge.py` | Charge-correlator analyses |
| `compact_original_fig3.py` | Optional: compact the original Figure 3 CSV for `compare_fig3.py` |
| `tagging/` | Leading-hadron flavour tag (Track D, section below): `run_tagging.sh`, `run_analysis_tagging.sh`, `leadpi_tag.py`, and the student notebook `leading_hadron_notebook.ipynb` with `make_notebook_sample.sh` |
| `exploratory/` | Studies that led here but are not part of the results above: reduced primordial kT and `mJoinJunction` scans, untagged QQC, back-to-back shape checks, and jet-charge / charge-flow flavour taggers (`flavor_tagging.py`, `charge_profile.py`, `charge_tag.py`). Each script documents the samples it needs |

## Track D — struck-quark flavour from the leading hadron

Every Track C observable that selects "u-quark" events relies on telling a struck u from a struck
d. Track D measures how well the **leading hadron** does that: when the struck quark fragments, the
hadron with the largest share of its momentum most often contains it, so a struck u tends to lead
with a π⁺ (u d̄) or K⁺ (u s̄) and a struck d with a π⁻. Track D uses the same patched Pythia and
generator as Track C, with SRC off and the Figure 3 DIS window, showers off and on. The
generator-level struck flavour is used only to measure efficiency and purity.

The tag: the **leading charged hadron of the Breit current hemisphere**, the one with the largest
SIDIS z_h = P·P_h / P·q among charged hadrons with Breit p_z < 0, is a π⁺ (or π⁺/K⁺, or any of
π⁺, K⁺, p) with z_h above a cut. No jet is needed. With π⁺ or K⁺ and z_h > 0.5 this is the u-tag of
the charge-charge correlators (`qqc_utag.py`, `qqc_chi.py`, `recoil_charge.py`). `leadpi_tag.py`
also reports a jet version (the highest-energy constituent of the current anti-kT jet is a π⁺
carrying more than the cut of the jet energy).

Efficiency is the fraction of all struck-u (struck-d) events that pass the tag; purity is the
fraction of tagged events with a struck u rather than a d, at the Figure 3 mix u : d = 0.742 : 0.116.

```bash
scripts/shower_study/tagging/run_tagging.sh            # struck-u and struck-d samples, showers off/on (4 x 1M events)
scripts/shower_study/tagging/run_analysis_tagging.sh   # efficiencies, purities, z_h-cut scans (seconds)
```

Results go to `results_tagging/`: `leadpi_tag_z0p5.txt` (the z_h > 0.5 tag, including the jet
version and the species of the leading hadron), `leadpi_tag_zscan.*` and `leadpos_tag_zscan.*`
(purity and efficiency vs the z_h cut for π⁺, π⁺ or K⁺, and any positive hadron).

Reference values (showers off / on):

| Leading charged hadron of the current hemisphere | u efficiency | d efficiency | u purity |
|---|---|---|---|
| no tag | 1 | 1 | 0.865 |
| π⁺ (no z_h cut) | 0.500 / 0.494 | 0.254 / 0.260 | 0.926 / 0.924 |
| π⁺, z_h > 0.5 | 0.113 / 0.100 | 0.031 / 0.028 | 0.959 / 0.958 |
| π⁺ or K⁺, z_h > 0.5 (the QQC u-tag) | 0.143 / 0.127 | 0.040 / 0.036 | 0.958 / 0.957 |
| π⁺, K⁺ or p, z_h > 0.5 | 0.183 / 0.162 | 0.062 / 0.056 | 0.950 / 0.949 |

**Student notebook.** `tagging/leading_hadron_notebook.ipynb` does the same with jet clustering, on
a small sample of already-generated events: it implements e⁺e⁻ anti-kT in a few lines of numpy,
checks that it reproduces the generator's jets (it does, event by event), takes the jet's leading
charged hadron by z_h, and reproduces the purity and efficiency plots above. With a z_h > 0.5 cut the
jet and whole-hemisphere tags agree (π⁺: 0.099 / 0.028, purity 0.957, showers on); without a cut the
jet version tags fewer events because in about a quarter of events the hemisphere's leading charged
hadron is not in the jet. Make its events once (about 20 seconds for 300k) and open it in Jupyter,
VS Code or Cursor (with `ipykernel` installed):

```bash
scripts/shower_study/tagging/make_notebook_sample.sh            # 300k events, showers on
scripts/shower_study/tagging/make_notebook_sample.sh 300000 0   # optional: showers off
```

The sample (`$STUDY_ROOT/notebook/`) has the usual event file plus one row per visible final-state
particle with its Breit-frame four-momentum (generator option `particles=`).
