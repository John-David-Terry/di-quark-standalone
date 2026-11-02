#!/usr/bin/env python3
"""
p_T-weighted recoil charge in u-tagged events (production_qqc_utag/, no regeneration).

  Q_rec = sum_b Q_b w_b / sum_b w_b,   w_b = -p_Tb . n_tag,

the charge of the target-side hadrons weighted by how much of their Breit transverse momentum
points opposite the tagging hadron. Parton-level expectation: a remnant recoiling as one ud unit
gives +1/3; the SRC partner d alone gives -1/3. Computed over the full target hemisphere and over
the original forward acceptance (2 < eta < 3.5, z_LC > 0.1), in bins of the tagging hadron's p_T.
Q_rec is a ratio of sample sums with event-level delta-method errors. Also reports the net charge
of the current and target hemispheres per tagged event (charge conservation check: sum = +1).
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
TAG = STUDY / "production_qqc_utag"
OUT = STUDY / "results_recoil_charge"
PT_BINS = ((0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, math.inf))


def ratio(a: np.ndarray, b: np.ndarray) -> tuple:
    R = a.sum() / b.sum()
    v = (a.var() - 2 * R * np.cov(a, b, bias=True)[0, 1] + R ** 2 * b.var()) / (len(a) * b.mean() ** 2)
    return R, math.sqrt(v)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say("Recoil charge Q_rec = sum Q_b w_b / sum w_b, w_b = -p_Tb . n_tag (Breit), u-tagged events")
    for sh in ("off", "on"):
        for src in (0, 1):
            had = pd.read_csv(TAG / f"showers_{sh}_src{src}_qqc.csv")
            ev_ids, ev = np.unique(had["event_id"].to_numpy(), return_inverse=True)
            n = len(ev_ids)
            had["ev"] = ev
            lead = had[(had["current"] == 1) & (had["is_lead"] == 1)].set_index("ev")
            cur, tgt = had[had["current"] == 1], had[had["current"] == 0]
            qcur = np.bincount(cur["ev"], weights=cur["charge"], minlength=n)
            qtgt = np.bincount(tgt["ev"], weights=tgt["charge"], minlength=n)
            say(f"\nshowers {sh}, SRC {'on' if src else 'off'}: {n:,} tagged events; <Q_current> {qcur.mean():.3f}, "
                f"<Q_target hemisphere> {qtgt.mean():.3f}, sum {np.mean(qcur + qtgt):.3f}")
            pt_tag = lead["pt"].reindex(np.arange(n)).to_numpy()
            acc = (tgt["eta_eic"] > 2) & (tgt["eta_eic"] < 3.5) & (tgt["zlc"] > 0.1)
            for lab, t in (("full target hemisphere", tgt), ("forward acceptance", tgt[acc])):
                w = -t["pt"].to_numpy() * np.cos(t["phi"].to_numpy() - lead.loc[t["ev"], "phi"].to_numpy())
                a = np.bincount(t["ev"], weights=t["charge"].to_numpy() * w, minlength=n)
                b = np.bincount(t["ev"], weights=w, minlength=n)
                R, e = ratio(a, b)
                parts = []
                for lo, hi in PT_BINS:
                    s = (pt_tag >= lo) & (pt_tag < hi)
                    Rb, eb = ratio(a[s], b[s])
                    parts.append(f"[{lo:g},{hi:g}) {Rb:+.3f}±{eb:.3f} ({s.mean():.2f})")
                say(f"  Q_rec {lab:23s} all {R:+.3f}±{e:.3f} | by tag p_T [GeV] (event fraction): " + ", ".join(parts))
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
