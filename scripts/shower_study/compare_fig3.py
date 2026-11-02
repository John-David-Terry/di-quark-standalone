#!/usr/bin/env python3
"""
Compare R_{pi-/pi+}(phi_hJ) of the original Figure 3 with the patched-Pythia shower study.

The 0% / 5% / 10% mixtures and the ratio use the original code unchanged:
``subsample_altered_unchanged_for_pool`` (plot_jet_hadron_pi_pm_observables.py) and
``_ratio_series`` (jet_hadron_pi_pm_figures.py), 40 bins on [0, pi], seeds 42 + index, pi+
from unchanged rows only, each curve normalized to its own count.

Signal summaries (all from the per-bin ratios):
  * mod(f): weighted fit R(f)/R(0%) = a (1 + mod cos phi). mod < 0 means a back-to-back
    enhancement of pi- relative to the 0% baseline; it is the SRC signal.
  * b2b(f): mean of R(f)/R(0%) - 1 over phi > 3 pi/4.
Errors are per-bin Poisson errors propagated as if R(f) and R(0%) were independent (they share
the pi+ curve and part of the pi- rows, so these errors are conservative).
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
for p in (REPO / "src", REPO / "scripts" / "analysis"):
    sys.path.insert(0, str(p))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from jet_hadron_pi_pm_figures import JH, _ratio_series  # noqa: E402
from plot_jet_hadron_pi_pm_observables import subsample_altered_unchanged_for_pool  # noqa: E402

BINS = np.linspace(JH.ANGLE_RANGE[0], JH.ANGLE_RANGE[1], JH.ANGLE_BINS + 1)
CENTERS = 0.5 * (BINS[:-1] + BINS[1:])
FRACS = (0.0, 0.05, 0.10)
SEED = 42
STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
NEW_COLS = ["event_id", "sample", "pion", "ok", "phi_hJ", "phi_hJ_parton", "remnant_code"]


def load_new(name: str) -> pd.DataFrame:
    df = pd.read_csv(STUDY / "production" / f"{name}.csv", usecols=NEW_COLS)
    df = df[df["ok"].astype(bool)].copy()
    df["ok"] = True
    return df


def load_original() -> Optional[pd.DataFrame]:
    """Original Figure 3 rows (compact_original_fig3.py); None when that file is absent."""
    path = STUDY / "original/original_fig3_compact.parquet"
    if not path.is_file():
        return None
    df = pd.read_parquet(path)
    df = df[df["ok"]].copy()
    for c in ("sample", "pion"):
        df[c] = df[c].astype(str)
    return df


def r_curves(df: pd.DataFrame, phi_col: str) -> Dict[float, Dict[str, np.ndarray]]:
    """R(phi_hJ) for the Figure 3 mixtures, exactly as plot_jet_hadron_pi_pm_observables.py."""
    d = df[["event_id", "sample", "pion", "ok"]].copy()
    d["phi_hJ"] = df[phi_col].to_numpy()
    unc = d[d["sample"] == "unchanged"]
    vp = unc.loc[unc["pion"] == "piplus", "phi_hJ"].to_numpy(dtype=float)
    vp = vp[np.isfinite(vp)]
    out = {}
    for ix, f in enumerate(FRACS):
        if f == 0.0:
            mix, achieved = unc, 0.0
        else:
            mix, _, _, achieved = subsample_altered_unchanged_for_pool(d, f, seed=SEED + ix)
        vm = mix.loc[mix["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float)
        vm = vm[np.isfinite(vm)]
        _, ratio, err, valid = _ratio_series(vp, vm, BINS)
        out[f] = {"R": ratio, "err": err, "valid": valid, "n_pip": len(vp), "n_pim": len(vm),
                  "achieved": achieved}
    return out


def cos_fit(y: np.ndarray, e: np.ndarray) -> tuple:
    """Weighted fit y = a (1 + m cos phi); returns (m, sigma_m)."""
    ok = np.isfinite(y) & np.isfinite(e) & (e > 0)
    X = np.vstack([np.ones(ok.sum()), np.cos(CENTERS[ok])]).T
    W = 1.0 / e[ok] ** 2
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    a, c = cov @ (X.T @ (W * y[ok]))
    m = c / a
    rel2 = cov[1, 1] / c**2 + cov[0, 0] / a**2 - 2 * cov[0, 1] / (a * c)
    # A handful of populated bins can leave rel2 numerically negative: report the error as undefined.
    sm = abs(m) * math.sqrt(rel2) if rel2 >= 0 else float("nan")
    return m, sm


def signal(curves: Dict[float, Dict[str, np.ndarray]], f: float) -> Dict[str, float]:
    r0, rf = curves[0.0], curves[f]
    rel = rf["R"] / r0["R"]
    rel_err = rel * np.sqrt((rf["err"] / rf["R"]) ** 2 + (r0["err"] / r0["R"]) ** 2)
    m, sm = cos_fit(rel, rel_err)
    sel = CENTERS > 0.75 * math.pi
    w = 1.0 / rel_err[sel] ** 2
    b2b = float(np.sum(w * (rel[sel] - 1.0)) / np.sum(w))
    return {"mod": m, "mod_err": sm, "b2b": b2b, "b2b_err": float(1.0 / math.sqrt(np.sum(w)))}


def mean_cos(df: pd.DataFrame, phi_col: str, sample: str, pion: str) -> str:
    v = df.loc[(df["sample"] == sample) & (df["pion"] == pion), phi_col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    c = np.cos(v)
    return f"{c.mean():+.4f}±{c.std() / math.sqrt(len(c)):.4f} (n={len(c):,})"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nokick-csv", type=Path,
                    default=STUDY / "original_surgery_nokick/jet_hadron_pi_pm_observables_nokick.csv",
                    help="Original pipeline, surgery v2 with --delta-px 0 (altered rows only).")
    ap.add_argument("--out-dir", type=Path, default=STUDY / "results")
    ap.add_argument("--write-plot-inputs", action="store_true",
                    help="Also write ok-row CSVs (anti-kT phi_hJ) for plot_jet_hadron_pi_pm_observables.py.")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    orig = load_original()
    off = pd.concat([load_new("showers_off_unchanged"), load_new("showers_off_altered")])
    on = pd.concat([load_new("showers_on_unchanged"), load_new("showers_on_altered")])
    ctrl = load_new("showers_off_diquark_control")
    if args.write_plot_inputs:
        for tag, df in (("showers_off", off), ("showers_on", on)):
            path = STUDY / "production" / f"fig3_plot_input_{tag}.csv"
            df[["event_id", "sample", "pion", "ok", "phi_hJ"]].to_csv(path, index=False)
            print(f"wrote {path}")

    configs: Dict[str, tuple] = {}
    if orig is not None:
        configs["original Fig. 3 (surgery v2, parton jet)"] = (orig, "phi_hJ")
    else:
        print("original Figure 3 data not found; comparing the new samples only")
    configs.update({
        "new, showers off, parton jet": (off, "phi_hJ_parton"),
        "new, showers off, anti-kT jet": (off, "phi_hJ"),
        "new, showers on, anti-kT jet": (on, "phi_hJ"),
    })
    nokick: Optional[pd.DataFrame] = None
    if orig is not None and args.nokick_csv.is_file():
        # The zero-kick run covers only the first 800k original events (event_id < n_sub), so
        # both it and its reference use the unchanged rows of that same subset.
        nk = pd.read_csv(args.nokick_csv, usecols=["event_id", "sample", "pion", "ok", "phi_hJ"])
        n_sub = int(nk["event_id"].max()) + 1
        nk = nk[nk["ok"].astype(str).str.lower() == "true"].copy()
        nk["ok"] = True
        sub = orig[orig["event_id"] < n_sub][nk.columns]
        configs[f"original Fig. 3, first {n_sub // 1000}k events"] = (sub, "phi_hJ")
        nokick = pd.concat([sub[sub["sample"] == "unchanged"], nk])
        configs[f"original pipeline, zero kick, first {n_sub // 1000}k events"] = (nokick, "phi_hJ")

    lines: List[str] = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731

    say("== Sample composition (ok rows) ==")
    for label, df in (("original", orig), ("new showers off", off), ("new showers on", on),
                      ("new diquark control (off)", ctrl)):
        if df is None:
            continue
        for s in ("unchanged", "altered"):
            sub = df[df["sample"] == s]
            if len(sub):
                say(f"  {label:26s} {s:9s} ok rows={len(sub):>9,}  pi- fraction={np.mean(sub['pion'] == 'piminus'):.3f}")

    say("\n== <cos phi_hJ> per arm and pion ==")
    if orig is not None:
        say(f"  {'original unchanged (parton jet)':44s} pi+ {mean_cos(orig, 'phi_hJ', 'unchanged', 'piplus')}  pi- {mean_cos(orig, 'phi_hJ', 'unchanged', 'piminus')}")
        say(f"  {'original altered (parton jet)':44s} pi+ {mean_cos(orig, 'phi_hJ', 'altered', 'piplus')}  pi- {mean_cos(orig, 'phi_hJ', 'altered', 'piminus')}")
    c2101 = ctrl[ctrl["remnant_code"] == 2101]
    for label, df, col in (
        ("diquark control, ud0 only (parton jet)", c2101, "phi_hJ_parton"),
        ("diquark control, ud0+ud1 (parton jet)", ctrl, "phi_hJ_parton"),
        ("diquark control, ud0 only (anti-kT)", c2101, "phi_hJ"),
    ):
        say(f"  {label:44s} pi+ {mean_cos(df, col, 'unchanged', 'piplus')}  pi- {mean_cos(df, col, 'unchanged', 'piminus')}")
    for label, df, col in (("new off, parton jet", off, "phi_hJ_parton"), ("new off, anti-kT", off, "phi_hJ"),
                           ("new on, anti-kT", on, "phi_hJ"), ("new on, parton (post-FSR quark)", on, "phi_hJ_parton")):
        for s in ("unchanged", "altered"):
            say(f"  {label + ' ' + s:44s} pi+ {mean_cos(df, col, s, 'piplus')}  pi- {mean_cos(df, col, s, 'piminus')}")
    if nokick is not None:
        say(f"  {'original surgery v2, zero kick, altered':44s} pi+ {mean_cos(nokick, 'phi_hJ', 'altered', 'piplus')}  pi- {mean_cos(nokick, 'phi_hJ', 'altered', 'piminus')}")

    say("\n== R(phi_hJ) mixtures: SRC signal relative to the 0% baseline ==")
    say(f"  {'configuration':48s} {'f':>5s} {'n(pi+)':>9s} {'n(pi-)':>9s} {'achieved f':>10s} "
        f"{'mod (cos fit)':>17s} {'b2b excess, phi>3pi/4':>24s}")
    rows = []
    all_curves = {}
    for label, (df, col) in configs.items():
        curves = r_curves(df, col)
        all_curves[label] = curves
        for f in FRACS:
            c = curves[f]
            if f > 0:
                s = signal(curves, f)
                sig = f"{s['mod']:+.4f}±{s['mod_err']:.4f}   {s['b2b']:+.4f}±{s['b2b_err']:.4f}"
            else:
                m, sm = cos_fit(c["R"], c["err"])
                sig = f"(baseline R(0%) mod {m:+.4f}±{sm:.4f})"
            say(f"  {label:48s} {f:5.2f} {c['n_pip']:>9,} {c['n_pim']:>9,} {c['achieved']:>10.4f}   {sig}")
            for i in range(len(CENTERS)):
                rows.append({"configuration": label, "f": f, "phi_center": CENTERS[i],
                             "R": c["R"][i], "R_err": c["err"][i]})
    pd.DataFrame(rows).to_csv(args.out_dir / "r_phi_hJ_values.csv", index=False)

    off_lab, on_lab = "new, showers off, anti-kT jet", "new, showers on, anti-kT jet"
    say("\n== Shower effect on the SRC signal (anti-kT jets) ==")
    for f in FRACS[1:]:
        so, sn = signal(all_curves[off_lab], f), signal(all_curves[on_lab], f)
        ratio = sn["mod"] / so["mod"]
        rerr = abs(ratio) * math.sqrt((sn["mod_err"] / sn["mod"]) ** 2 + (so["mod_err"] / so["mod"]) ** 2)
        say(f"  f={f:.2f}: mod off={so['mod']:+.4f}±{so['mod_err']:.4f}  on={sn['mod']:+.4f}±{sn['mod_err']:.4f}"
            f"  on/off={ratio:.2f}±{rerr:.2f}")

    (args.out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Figure: R(phi_hJ) for each configuration, same colours as Figure 3.
    n = len(all_curves)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.4), sharey=True)
    colors = {0.0: "gray", 0.05: "purple", 0.10: "red"}
    for ax, (label, curves) in zip(np.atleast_1d(axes), all_curves.items()):
        for f in FRACS:
            c = curves[f]
            ax.errorbar(CENTERS[c["valid"]], c["R"][c["valid"]], yerr=c["err"][c["valid"]], fmt="o",
                        ms=2.5, color=colors[f], label=f"$\\pi^-$ ({int(round(100 * f))}%)")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_title(label, fontsize=8)
        ax.set_xlabel(r"$\phi_{hJ}$")
    np.atleast_1d(axes)[0].set_ylabel(r"$R_{\pi^-/\pi^+}$")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out_dir / "compare_R_phi_hJ.pdf")
    fig.savefig(args.out_dir / "compare_R_phi_hJ.png", dpi=130)
    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
