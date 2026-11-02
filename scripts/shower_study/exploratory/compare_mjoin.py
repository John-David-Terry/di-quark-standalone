#!/usr/bin/env python3
"""
Summarize the mJoinJunction scan (run_mjoin_scan.sh): per threshold, how often hadronization
fuses the remnant u and d into a diquark, Pythia failures/errors, pion yield, the SRC shift of
<cos phi_hJ> (pi-, pi+, flavour-tagged difference), and the Figure 3 R(phi_hJ) modulation.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_fig3 import STUDY, r_curves, signal  # noqa: E402

COLS = ["event_id", "sample", "pion", "ok", "phi_hJ", "phi_hJ_parton", "remnant_merged"]
DIR = STUDY / "production_mjoin"


def shift(unc: pd.DataFrame, alt: pd.DataFrame, col: str, pion: str) -> tuple:
    vals = []
    for d in (unc, alt):
        v = d.loc[d["ok"] & (d["pion"] == pion), col].to_numpy(dtype=float)
        c = np.cos(v[np.isfinite(v)])
        vals.append((c.mean(), c.std() / math.sqrt(len(c))))
    return vals[1][0] - vals[0][0], math.hypot(vals[0][1], vals[1][1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--primkt-tag", default="primkt0p13")
    ap.add_argument("--write-plot-inputs", action="store_true",
                    help="Only write ok-row CSVs (anti-kT and parton phi_hJ) per threshold for "
                         "plot_jet_hadron_pi_pm_observables.py, then exit.")
    args = ap.parse_args()

    tags = sorted({re.search(r"_mjoin([0-9p]+)\.csv$", p.name).group(1)
                   for p in DIR.glob(f"showers_off_unchanged_{args.primkt_tag}_mjoin*.csv")},
                  key=lambda t: -float(t.replace("p", ".")))
    if args.write_plot_inputs:
        for t in tags:
            both = pd.concat([pd.read_csv(DIR / f"showers_off_{arm}_{args.primkt_tag}_mjoin{t}.csv", usecols=COLS)
                              for arm in ("unchanged", "altered")])
            both = both[both["ok"].astype(bool)]
            for jet, col in (("antikt", "phi_hJ"), ("parton", "phi_hJ_parton")):
                path = DIR / f"fig3_plot_input_{args.primkt_tag}_mjoin{t}_{jet}.csv"
                out = both[["event_id", "sample", "pion", "ok"]].copy()
                out["phi_hJ"] = both[col].to_numpy()
                out.to_csv(path, index=False)
                print(f"wrote {path}")
        return
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say(f"mJoinJunction scan, showers off, {args.primkt_tag}")
    for t in tags:
        unc = pd.read_csv(DIR / f"showers_off_unchanged_{args.primkt_tag}_mjoin{t}.csv", usecols=COLS)
        alt = pd.read_csv(DIR / f"showers_off_altered_{args.primkt_tag}_mjoin{t}.csv", usecols=COLS)
        for d in (unc, alt):
            d["ok"] = d["ok"].astype(bool)
        fails, errors = 0, Counter()
        for arm in ("unchanged", "altered"):
            for js in DIR.glob(f"chunks/showers_off_{arm}_{args.primkt_tag}_mjoin{t}_c*.csv.summary.json"):
                fails += json.loads(js.read_text())["n_next_failed"]
            for lg in DIR.glob(f"chunks/showers_off_{arm}_{args.primkt_tag}_mjoin{t}_c*.log"):
                for line in lg.read_text(errors="ignore").splitlines():
                    if "PYTHIA Error" in line or "PYTHIA Abort" in line:
                        errors[line.split(":", 1)[0].replace(" PYTHIA Error in ", "").strip()] += 1
        n_ev = len(unc) + len(alt)
        say(f"\n== mJoinJunction = {t.replace('p', '.')} GeV ==")
        for arm, d in (("unchanged", unc), ("altered", alt)):
            ok = d[d["ok"]]
            say(f"  {arm:9s} eligible={len(d):,} ok={len(ok) / len(d):.3f} pi- frac={np.mean(ok['pion'] == 'piminus'):.3f}  "
                f"d fused into diquark: all={d['remnant_merged'].mean():.3f} ok pi-={ok.loc[ok['pion'] == 'piminus', 'remnant_merged'].mean():.3f}")
        say(f"  failed events: {fails} of {n_ev:,};  PYTHIA error lines: {sum(errors.values())} "
            + (str(dict(errors.most_common(4))) if errors else ""))
        both = pd.concat([unc[unc["ok"]], alt[alt["ok"]]])
        for jet, col in (("anti-kT", "phi_hJ"), ("parton", "phi_hJ_parton")):
            sm, em = shift(unc, alt, col, "piminus")
            sp, ep = shift(unc, alt, col, "piplus")
            s10 = signal(r_curves(both, col), 0.10)
            say(f"  {jet:8s} SRC shift <cos>: pi- {sm:+.4f}±{em:.4f}  pi+ {sp:+.4f}±{ep:.4f}  "
                f"flavour-tagged {sm - sp:+.4f}±{math.hypot(em, ep):.4f}   R mod(10%) {s10['mod']:+.4f}±{s10['mod_err']:.4f}")
    out = STUDY / "results_mjoin"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'summary.txt'}")


if __name__ == "__main__":
    main()
