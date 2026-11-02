#!/usr/bin/env python3
"""
Compare the SRC signal between two sample pairs in production_primkt/ (e.g. default vs reduced
primordial kT, or showers off vs on). Each pair is ``showers_{off,on}_{unchanged,altered}_{tag}``.

Uses the same R(phi_hJ) mixtures and signal summaries as compare_fig3.py (original Figure 3
mixture/ratio code). Adds, per arm: realized primordial kT (scattered system and remnant d, in
Pythia's hadronic rest frame), how often hadronization recombines the remnant u and d into a
diquark (status 74), and the SRC shift of <cos phi_hJ> split by recombined vs separate d.
Ends with the ratio second/first of the flavour-tagged shift and of the 10% R modulation
(the shower survival fraction when the two pairs differ only in showers).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fig3 import CENTERS, FRACS, STUDY, plt, r_curves, signal  # noqa: E402

COLS = ["event_id", "sample", "pion", "ok", "phi_hJ", "phi_hJ_parton", "kt_scat_GeV", "kt_d_GeV",
        "remnant_merged"]
JETS = (("anti-kT jet", "phi_hJ", "antikt"), ("struck-parton jet", "phi_hJ_parton", "parton"))


def load(name: str) -> pd.DataFrame:
    df = pd.read_csv(STUDY / "production_primkt" / f"{name}.csv", usecols=COLS)
    df["ok"] = df["ok"].astype(bool)
    return df


def cos_stats(v: np.ndarray) -> tuple:
    c = np.cos(v[np.isfinite(v)])
    return c.mean(), c.std() / math.sqrt(len(c)), len(c)


def ratio(a: tuple, b: tuple) -> tuple:
    r = a[0] / b[0]
    return r, abs(r) * math.hypot(a[1] / a[0], b[1] / b[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="primkt0p13", help="Sample tag compared with the baseline.")
    ap.add_argument("--label", default="reduced primordial kT", help="Description of --tag.")
    ap.add_argument("--showers", choices=("off", "on"), default="off", help="Showers for --tag.")
    ap.add_argument("--baseline", default="defaultkt", help="Sample tag compared against.")
    ap.add_argument("--baseline-label", default="default primordial kT", help="Description of --baseline.")
    ap.add_argument("--baseline-showers", choices=("off", "on"), default="off")
    ap.add_argument("--out-dir", type=Path, default=STUDY / "results_primkt")
    ap.add_argument("--write-plot-inputs", action="store_true",
                    help="Only write ok-row CSVs (anti-kT and parton phi_hJ) for "
                         "plot_jet_hadron_pi_pm_observables.py, then exit.")
    args = ap.parse_args()
    pairs = [(args.baseline, args.baseline_showers, args.baseline_label),
             (args.tag, args.showers, args.label)]

    if args.write_plot_inputs:
        for tag, sh, _ in pairs:
            both = pd.concat([load(f"showers_{sh}_unchanged_{tag}"), load(f"showers_{sh}_altered_{tag}")])
            both = both[both["ok"]]
            prefix = "" if sh == "off" else "showers_on_"
            for _, col, jet in JETS:
                path = STUDY / "production_primkt" / f"fig3_plot_input_{prefix}{tag}_{jet}.csv"
                out = both[["event_id", "sample", "pion", "ok"]].copy()
                out["phi_hJ"] = both[col].to_numpy()
                out.to_csv(path, index=False)
                print(f"wrote {path}")
        return
    args.out_dir.mkdir(parents=True, exist_ok=True)

    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    curves_by = {}
    results = []

    for tag, sh, desc in pairs:
        label = f"{desc} ({tag}, showers {sh})"
        unc = load(f"showers_{sh}_unchanged_{tag}")
        alt = load(f"showers_{sh}_altered_{tag}")
        res = {}
        say(f"== {label} ==")
        for arm, d in (("unchanged", unc), ("altered", alt)):
            ok = d[d["ok"]]
            say(f"  {arm:9s} eligible={len(d):,}  ok={len(ok):,} ({len(ok) / len(d):.3f})  "
                f"pi- frac={np.mean(ok['pion'] == 'piminus'):.3f}  "
                f"<|kT scat|>={d['kt_scat_GeV'].mean():.4f}  <|kT d|>={d['kt_d_GeV'].mean():.4f} GeV  "
                f"d recombined: all={d['remnant_merged'].mean():.3f}  ok pi-={ok.loc[ok['pion'] == 'piminus', 'remnant_merged'].mean():.3f}  "
                f"ok pi+={ok.loc[ok['pion'] == 'piplus', 'remnant_merged'].mean():.3f}")

        for jet, col, _ in JETS:
            say(f"  <cos phi_hJ>, {jet}:")
            for subset, sel in (("all", None), ("d separate", 0), ("d recombined", 1)):
                parts = []
                shifts = {}
                for pion in ("piminus", "piplus"):
                    stats = []
                    for d in (unc, alt):
                        m = d["ok"] & (d["pion"] == pion)
                        if sel is not None:
                            m &= d["remnant_merged"] == sel
                        if not m.any():
                            break
                        stats.append(cos_stats(d.loc[m, col].to_numpy(dtype=float)))
                    if len(stats) < 2:
                        break
                    (cu, eu, nu), (ca, ea, na) = stats
                    shifts[pion] = (ca - cu, math.hypot(ea, eu))
                    parts.append(f"{pion}: unch {cu:+.4f} alt {ca:+.4f} shift {ca - cu:+.4f}±{math.hypot(ea, eu):.4f} (n={nu:,}/{na:,})")
                if len(shifts) < 2:
                    say(f"    {subset:13s} (no events)")
                    continue
                tag_shift = shifts["piminus"][0] - shifts["piplus"][0]
                tag_err = math.hypot(shifts["piminus"][1], shifts["piplus"][1])
                if sel is None:
                    res[(jet, "flavour-tagged shift")] = (tag_shift, tag_err)
                    res[(jet, "pi- shift")] = shifts["piminus"]
                say(f"    {subset:13s} " + "   ".join(parts) + f"   flavour-tagged (pi- minus pi+) {tag_shift:+.4f}±{tag_err:.4f}")

        both = pd.concat([unc[unc["ok"]], alt[alt["ok"]]])
        for jet, col, _ in JETS:
            curves = r_curves(both, col)
            curves_by[(label, jet)] = curves
            say(f"  R(phi_hJ) mixtures, {jet}:")
            for f in FRACS[1:]:
                s = signal(curves, f)
                say(f"    f={f:.2f}: mod {s['mod']:+.4f}±{s['mod_err']:.4f}   b2b excess (phi>3pi/4) "
                    f"{s['b2b']:+.4f}±{s['b2b_err']:.4f}   n(pi+)={curves[f]['n_pip']:,} n(pi-)={curves[f]['n_pim']:,}")
                if f == 0.10:
                    res[(jet, "R mod at 10%")] = (s["mod"], s["mod_err"])
        results.append(res)
        say()

    say("== Ratio second / first (shower survival when only showers differ) ==")
    for key in results[0]:
        r, e = ratio(results[1][key], results[0][key])
        say(f"  {key[0]:18s} {key[1]:22s} {results[0][key][0]:+.4f} -> {results[1][key][0]:+.4f}   ratio {r:.3f}±{e:.3f}")

    (args.out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, len(curves_by), figsize=(3.6 * len(curves_by), 3.4), sharey=True)
    colors = {0.0: "gray", 0.05: "purple", 0.10: "red"}
    for ax, ((label, jet), curves) in zip(np.atleast_1d(axes), curves_by.items()):
        for f in FRACS:
            c = curves[f]
            ax.errorbar(CENTERS[c["valid"]], c["R"][c["valid"]], yerr=c["err"][c["valid"]], fmt="o",
                        ms=2.5, color=colors[f], label=f"$\\pi^-$ ({int(round(100 * f))}%)")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_title(f"{label}\n{jet}", fontsize=7)
        ax.set_xlabel(r"$\phi_{hJ}$")
    np.atleast_1d(axes)[0].set_ylabel(r"$R_{\pi^-/\pi^+}$")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out_dir / "R_phi_hJ_compare.pdf")
    fig.savefig(args.out_dir / "R_phi_hJ_compare.png", dpi=130)
    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
