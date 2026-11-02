#!/usr/bin/env python3
"""
How well can the current jet tag the struck-quark flavour (u vs d)?

Input: production_tagging/showers_{off,on}_anystruck.csv from ``fig3_patched struck=any``
(Figure 3 DIS window, every struck quark/antiquark, SRC off). Two taggers:
  * leading pion: the highest-energy constituent of the current anti-kT jet (Breit frame) is a
    pi+ (u-tag) or a pi- (d-tag); any other leading particle leaves the event untagged;
  * jet charge Q_kappa = sum_i (E_i/E_J)^kappa Q_i, cut from above (u-tag) or below (d-tag).
Reported for all events with a current jet and for the Figure 3 analysis subset (accepted
target pion, ok rows). Efficiency = P(tag | true flavour); purity = P(true flavour | tag),
with all struck flavours (including antiquarks, s, c) in the denominator.
"""

from __future__ import annotations

import os
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
DIR = STUDY / "production_tagging"
OUT = STUDY / "results_tagging"
KAPPAS = (("0", "jetQ_k0"), ("0.3", "jetQ_k0p3"), ("0.5", "jetQ_k0p5"), ("1", "jetQ_k1"))
COLS = ["ok", "struck_id", "tagjet_E", "lead_id", "leadch_id"] + [c for _, c in KAPPAS]


def load(sh: str) -> pd.DataFrame:
    d = pd.read_csv(DIR / f"showers_{sh}_anystruck.csv", usecols=COLS)
    d["ok"] = d["ok"].astype(bool)
    return d[d["tagjet_E"].notna()].copy()


def auc(q_u: np.ndarray, q_d: np.ndarray) -> float:
    """P(Q of a random u-jet > Q of a random d-jet), ties counted half."""
    r = pd.Series(np.concatenate([q_u, q_d])).rank(method="average").to_numpy()
    return float((r[: len(q_u)].sum() - len(q_u) * (len(q_u) + 1) / 2) / (len(q_u) * len(q_d)))


def tag_stats(d: pd.DataFrame, u_tag: np.ndarray, d_tag: np.ndarray) -> dict:
    is_u, is_d = d["struck_id"].to_numpy() == 2, d["struck_id"].to_numpy() == 1
    return {
        "eff_u": u_tag[is_u].mean(), "eff_d": d_tag[is_d].mean(),
        "mis_u_as_d": d_tag[is_u].mean(), "mis_d_as_u": u_tag[is_d].mean(),
        "pur_u": is_u[u_tag].mean() if u_tag.any() else math.nan,
        "pur_d": is_d[d_tag].mean() if d_tag.any() else math.nan,
        "frac_u_tag": u_tag.mean(), "frac_d_tag": d_tag.mean(),
    }


def matched_cut(q: np.ndarray, is_flav: np.ndarray, target_eff: float, below: bool) -> float:
    """Jet-charge cut giving (as close as possible) the target efficiency for this flavour."""
    vals = np.sort(q[is_flav])
    k = int(round(target_eff * len(vals)))
    k = min(max(k, 1), len(vals) - 1)
    return vals[k - 1] + 1e-12 if below else vals[len(vals) - k] - 1e-12


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    data = {sh: load(sh) for sh in ("off", "on")}
    roc = {}

    for sh, dall in data.items():
        for subset, d in (("all events with a current jet", dall), ("Fig. 3 subset (accepted target pion)", dall[dall["ok"]])):
            sid = d["struck_id"].to_numpy()
            is_u, is_d = sid == 2, sid == 1
            say(f"\n== showers {sh}, {subset}: n={len(d):,} ==")
            fr = d["struck_id"].value_counts(normalize=True)
            say("  struck flavour fractions: " + ", ".join(f"{k}:{fr.get(v, 0):.3f}" for k, v in
                (("u", 2), ("d", 1), ("ubar", -2), ("dbar", -1), ("s", 3), ("sbar", -3), ("c", 4), ("cbar", -4))))
            say(f"  <E_J> = {d['tagjet_E'].mean():.2f} GeV")

            lead = d["lead_id"].to_numpy()
            s = tag_stats(d, lead == 211, lead == -211)
            say(f"  leading pi tag : u-tag eff {s['eff_u']:.3f} (d mistag {s['mis_d_as_u']:.3f}) purity {s['pur_u']:.3f} | "
                f"d-tag eff {s['eff_d']:.3f} (u mistag {s['mis_u_as_d']:.3f}) purity {s['pur_d']:.3f}")
            lch = d["leadch_id"].to_numpy()
            pos = np.isin(lch, [211, 321, 2212, -11, -13, 3222, -3112, -3312, 413, 411, 431])
            neg = np.isin(lch, [-211, -321, -2212, 11, 13, -3222, 3112, 3312, -413, -411, -431])
            s2 = tag_stats(d, pos, neg)
            say(f"  leading charged: u-tag eff {s2['eff_u']:.3f} (d mistag {s2['mis_d_as_u']:.3f}) purity {s2['pur_u']:.3f} | "
                f"d-tag eff {s2['eff_d']:.3f} (u mistag {s2['mis_u_as_d']:.3f}) purity {s2['pur_d']:.3f}")

            for kap, col in KAPPAS:
                qv = d[col].to_numpy()
                a = auc(qv[is_u], qv[is_d])
                cu = matched_cut(qv, is_u, s["eff_u"], below=False)
                cd = matched_cut(qv, is_d, s["eff_d"], below=True)
                su = tag_stats(d, qv > cu, np.zeros(len(d), bool))
                sd = tag_stats(d, np.zeros(len(d), bool), qv < cd)
                say(f"  jet charge k={kap:3s}: <Q>_u {qv[is_u].mean():+.3f}  <Q>_d {qv[is_d].mean():+.3f}  AUC(u vs d) {a:.3f} | "
                    f"at leading-pi efficiencies: u-tag eff {su['eff_u']:.3f} purity {su['pur_u']:.3f}, "
                    f"d-tag eff {sd['eff_d']:.3f} purity {sd['pur_d']:.3f}")
                if subset.startswith("all") and kap == "0.5":
                    thr = np.unique(np.quantile(qv, np.linspace(0, 1, 401)))
                    roc[sh] = (np.array([(qv[is_d] < t).mean() for t in thr]),
                               np.array([(qv[is_u] < t).mean() for t in thr]),
                               s["eff_d"], s["mis_u_as_d"], qv[is_u], qv[is_d])

    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Figure: jet-charge distributions (showers on) and d-tag ROC (showers off/on).
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    qu, qd = roc["on"][4], roc["on"][5]
    bins = np.linspace(-1.5, 1.5, 61)
    ax[0].hist(qu, bins=bins, density=True, histtype="step", lw=2, color="tab:red", label=r"struck $u$")
    ax[0].hist(qd, bins=bins, density=True, histtype="step", lw=2, color="tab:blue", label=r"struck $d$")
    ax[0].set_xlabel(r"jet charge $Q_{\kappa=0.5}$", fontsize=fs)
    ax[0].set_ylabel("normalized", fontsize=fs)
    ax[0].legend(fontsize=fs, frameon=False, loc="upper left")
    for sh, ls in (("off", "--"), ("on", "-")):
        ed, eu, pd_, pu_, _, _ = roc[sh]
        ax[1].plot(eu, ed, ls=ls, lw=2, color="k", label=rf"$Q_{{\kappa=0.5}}$, showers {sh}")
        ax[1].plot([pu_], [pd_], marker="o" if sh == "on" else "s", ms=10, color="tab:purple", ls="none",
                   label=rf"leading $\pi^-$, showers {sh}")
    ax[1].plot([0, 1], [0, 1], color="0.6", lw=1)
    ax[1].set_xlabel(r"$u$ mistagged as $d$", fontsize=fs)
    ax[1].set_ylabel(r"$d$-tag efficiency", fontsize=fs)
    ax[1].set_xlim(0, 1)
    ax[1].set_ylim(0, 1)
    ax[1].legend(fontsize=fs * 0.7, frameon=False, loc="lower right")
    for a in ax:
        a.tick_params(which="both", direction="in", labelsize=fs)
        a.tick_params(which="major", length=10, width=1.2)
        a.tick_params(which="minor", length=5, width=1.0)
        a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    fig.tight_layout()
    fig.savefig(OUT / "flavor_tagging.pdf")
    fig.savefig(OUT / "flavor_tagging.png", dpi=150)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
