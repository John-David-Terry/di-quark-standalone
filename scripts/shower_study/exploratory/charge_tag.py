#!/usr/bin/env python3
"""
Event-by-event struck-quark charge tagging from the angular charge profile (charge_profile.py
inputs). Per event, discriminants built from the net charge of charged hadrons vs theta_hJ:
  * net charge in the current hemisphere;
  * net charge in cones theta_hJ < 0.5 and < 1.0 rad (all hadrons);
  * Fisher-optimal angular weights: w = Cov^-1 (mu_u - mu_d) on the 36-bin current-hemisphere
    net-charge vector, trained on even events and evaluated on odd events.
Reported: AUC (u vs d), and u/d efficiencies and purities at the Figure 3 struck-flavour mix
(u : d = 0.742 : 0.116 from the flavour-tagging sample; other flavours ignored).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from charge_profile import CENTERS, STUDY, load  # noqa: E402

F_U, F_D = 0.742, 0.116


def auc(s_u: np.ndarray, s_d: np.ndarray) -> float:
    r = pd.Series(np.concatenate([s_u, s_d])).rank(method="average").to_numpy()
    return float((r[: len(s_u)].sum() - len(s_u) * (len(s_u) + 1) / 2) / (len(s_u) * len(s_d)))


def working_points(s_u: np.ndarray, s_d: np.ndarray) -> list:
    """u-tag: s > cut; d-tag: s < cut. Report a few cuts with efficiencies and purities."""
    out = []
    for cut in np.unique(np.quantile(np.concatenate([s_u, s_d]), [0.2, 0.35, 0.5, 0.65, 0.8])):
        eu, mu = (s_u > cut).mean(), (s_d > cut).mean()
        ed, md = (s_d < cut).mean(), (s_u < cut).mean()
        pu = F_U * eu / (F_U * eu + F_D * mu) if eu + mu > 0 else np.nan
        pdd = F_D * ed / (F_D * ed + F_U * md) if ed + md > 0 else np.nan
        out.append((cut, eu, mu, pu, ed, md, pdd))
    return out


def main() -> None:
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say(f"Struck-quark charge tagging from the angular charge profile; purities use u:d = {F_U}:{F_D} "
        f"(no-tag purities u {F_U / (F_U + F_D):.3f}, d {F_D / (F_U + F_D):.3f})")
    for sh in ("off", "on"):
        data = {f: load(sh, f) for f in ("u", "d")}
        net = {f: data[f]["mats"][("current", "pos")] - data[f]["mats"][("current", "neg")] for f in data}
        full = {f: data[f]["mats"][("full", "pos")] - data[f]["mats"][("full", "neg")] for f in data}
        disc = {
            "current-hemisphere net charge": {f: net[f].sum(1) for f in data},
            "net charge, theta < 0.5 rad": {f: full[f][:, CENTERS < 0.5].sum(1) for f in data},
            "net charge, theta < 1.0 rad": {f: full[f][:, CENTERS < 1.0].sum(1) for f in data},
        }
        # Fisher weights from even events, applied to odd events.
        tr = {f: net[f][0::2].astype(np.float64) for f in data}
        te = {f: net[f][1::2].astype(np.float64) for f in data}
        cov = 0.5 * (np.cov(tr["u"], rowvar=False) + np.cov(tr["d"], rowvar=False))
        keep = np.diag(cov) > 1e-9
        w = np.zeros(len(CENTERS))
        w[keep] = np.linalg.solve(cov[np.ix_(keep, keep)], tr["u"].mean(0)[keep] - tr["d"].mean(0)[keep])
        disc["Fisher angular weights (test half)"] = {f: te[f] @ w for f in data}

        say(f"\n== showers {sh} ==")
        for name, s in disc.items():
            su, sd = s["u"], s["d"]
            say(f"  {name:36s} AUC {auc(su, sd):.3f}   <s>_u {su.mean():+.3f}  <s>_d {sd.mean():+.3f}")
        su, sd = disc["current-hemisphere net charge"]["u"], disc["current-hemisphere net charge"]["d"]
        say("  current-hemisphere net charge, distribution (fraction of events):")
        vals = np.arange(-3, 5)
        say("    Q    " + " ".join(f"{v:+6d}" for v in vals))
        say("    u    " + " ".join(f"{np.mean(su == v):6.3f}" for v in vals))
        say("    d    " + " ".join(f"{np.mean(sd == v):6.3f}" for v in vals))
        for name in ("current-hemisphere net charge", "Fisher angular weights (test half)"):
            say(f"  working points, {name}:")
            for cut, eu, mu, pu, ed, md, pdd in working_points(disc[name]["u"], disc[name]["d"]):
                say(f"    cut {cut:+.3f}: u-tag eff {eu:.3f} (d mistag {mu:.3f}) purity {pu:.3f} | "
                    f"d-tag eff {ed:.3f} (u mistag {md:.3f}) purity {pdd:.3f}")
    out = STUDY / "results_charge_profile" / "charge_tag_summary.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
