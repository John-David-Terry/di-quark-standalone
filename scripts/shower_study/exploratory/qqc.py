#!/usr/bin/env python3
"""
Azimuthal current-target charge-charge correlation in the Breit frame:

  C_QQ(dphi) = (1/N_evt) sum_events sum_{a in current} sum_{b in target} Q_a Q_b
               delta(dphi - |phi_a - phi_b|_folded),

current = charged final-state hadrons with Breit pz < 0, target = Breit pz > 0. Normalized per
selected DIS event (events with no pairs included) and per radian; errors from event-to-event
fluctuations. Input: production_qqc/ (fig3_patched Figure 3 struck-u selection, qqc=...).

Target selections: "acceptance" = original target cuts (proton-forward lab 2 < eta < 3.5,
z_LC > 0.1) applied to all charged hadrons; "full" = whole target hemisphere.
Diagnostics: split by target-hadron charge sign; charge-blind pair counts; charges of the target
hadrons shuffled within each event (keeps all momenta, removes the direction-charge link);
p_T threshold scan; charge conservation.
"""

from __future__ import annotations

import os
import math
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
DIR = STUDY / "production_qqc"
OUT = STUDY / "results_qqc"
NB = 36
EDGES = np.linspace(0.0, np.pi, NB + 1)
WIDTH = EDGES[1] - EDGES[0]
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])
B2B = CENTERS > 0.75 * np.pi
NEAR = CENTERS < 0.25 * np.pi
SAMPLES = [(sh, src) for sh in ("off", "on") for src in (0, 1)]


def load(sh: str, src: int) -> tuple:
    ev = pd.read_csv(DIR / f"showers_{sh}_src{src}.csv", usecols=["event_id"])["event_id"].to_numpy()
    had = pd.read_csv(DIR / f"showers_{sh}_src{src}_qqc.csv")
    idx = pd.Series(np.arange(len(ev)), index=ev)
    had["ev"] = idx.loc[had["event_id"].to_numpy()].to_numpy()
    return len(ev), had.drop(columns="event_id")


def target_mask(had: pd.DataFrame, sel: str) -> np.ndarray:
    m = had["current"].to_numpy() == 0
    if sel == "acceptance":
        eta, z = had["eta_eic"].to_numpy(), had["zlc"].to_numpy()
        m &= (eta > 2.0) & (eta < 3.5) & (z > 0.1)
    return m


def pairs(had: pd.DataFrame, sel: str, ptmin: float, seed: int = 1) -> pd.DataFrame:
    pt = had["pt"].to_numpy()
    cur = had.loc[(had["current"].to_numpy() == 1) & (pt > ptmin), ["ev", "charge", "phi"]]
    tgt = had.loc[target_mask(had, sel) & (pt > ptmin), ["ev", "charge", "phi"]].sort_values("ev", kind="stable")
    # Within-event permutation of the target charges (momenta untouched).
    key = np.random.default_rng(seed).random(len(tgt))
    order = np.lexsort((key, tgt["ev"].to_numpy()))
    tgt = tgt.assign(charge_shuf=tgt["charge"].to_numpy()[order])
    m = cur.merge(tgt, on="ev", suffixes=("_a", "_b"), sort=True)
    d = np.abs(m["phi_a"].to_numpy() - m["phi_b"].to_numpy()) % (2 * np.pi)
    d = np.where(d > np.pi, 2 * np.pi - d, d)
    qa, qb = m["charge_a"].to_numpy(), m["charge_b"].to_numpy()
    return pd.DataFrame({
        "ev": m["ev"].to_numpy(), "bin": np.clip(np.digitize(d, EDGES) - 1, 0, NB - 1),
        "w": qa * qb, "w_tpos": np.where(qb > 0, qa * qb, 0), "w_tneg": np.where(qb < 0, qa * qb, 0),
        "w_count": np.ones(len(m)), "w_shuf": qa * m["charge_shuf"].to_numpy(),
    })


def hist_stats(n_evt: int, p: pd.DataFrame, col: str) -> tuple:
    """Per-bin mean per event per rad and its event-level error (events processed in blocks)."""
    ev, b, w = p["ev"].to_numpy(), p["bin"].to_numpy(), p[col].to_numpy()
    s1, s2 = np.zeros(NB), np.zeros(NB)
    step = 250000
    for lo in range(0, n_evt, step):
        i0, i1 = np.searchsorted(ev, [lo, lo + step])
        k = min(step, n_evt - lo)
        flat = np.bincount((ev[i0:i1] - lo) * NB + b[i0:i1], weights=w[i0:i1], minlength=k * NB).reshape(k, NB)
        s1 += flat.sum(0)
        s2 += (flat ** 2).sum(0)
    mean = s1 / n_evt
    err = np.sqrt(np.maximum(s2 / n_evt - mean ** 2, 0) / n_evt)
    return mean / WIDTH, err / WIDTH


def region_stats(n_evt: int, p: pd.DataFrame, col: str, region: np.ndarray) -> tuple:
    """Per-event integral of the correlation over a set of bins: mean and error."""
    sel = region[p["bin"].to_numpy()]
    s = np.bincount(p["ev"].to_numpy()[sel], weights=p[col].to_numpy()[sel], minlength=n_evt)
    return s.mean(), s.std() / math.sqrt(n_evt)


def diff(a: tuple, b: tuple) -> tuple:
    return a[0] - b[0], math.hypot(a[1], b[1])


def fmt(v: tuple) -> str:
    return f"{v[0]:+.5f}±{v[1]:.5f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    cols = ("w", "w_tpos", "w_tneg", "w_count", "w_shuf")
    hist = {}
    reg = {}
    scan = {}

    for sh, src in SAMPLES:
        n_evt, had = load(sh, src)
        q = had["charge"].to_numpy()
        ev = had["ev"].to_numpy()
        cur = had["current"].to_numpy() == 1
        qc = np.bincount(ev[cur], weights=q[cur], minlength=n_evt)
        qt = np.bincount(ev[~cur], weights=q[~cur], minlength=n_evt)
        say(f"\n== showers {sh}, SRC {'on' if src else 'off'}: {n_evt:,} events, {len(had) / n_evt:.2f} charged hadrons/event ==")
        say(f"  charge: <Q_current> {qc.mean():+.4f}  <Q_target> {qt.mean():+.4f}  per-event total == +1 in "
            f"{np.mean(qc + qt == 1):.4f} of events  (hemispheres: current = Breit pz<0, target = pz>0)")
        for sel in ("acceptance", "full"):
            p = pairs(had, sel, 0.0)
            say(f"  [{sel}] pairs/event {len(p) / n_evt:.3f}")
            for c in cols:
                hist[(sh, src, sel, c)] = hist_stats(n_evt, p, c)
                for rname, rmask in (("b2b", B2B), ("near", NEAR), ("all", np.ones(NB, bool))):
                    reg[(sh, src, sel, c, rname)] = region_stats(n_evt, p, c, rmask)
            for ptmin in (0.1, 0.2, 0.3, 0.5):
                pp = pairs(had, sel, ptmin)
                for c in ("w", "w_shuf"):
                    scan[(sh, src, sel, ptmin, c)] = region_stats(n_evt, pp, c, B2B)
            for c in ("w", "w_shuf"):
                scan[(sh, src, sel, 0.0, c)] = reg[(sh, src, sel, c, "b2b")]
            del p
        del had

    say("\n== SRC-induced change, Delta = (SRC on) - (SRC off), per event (integrals over dphi regions) ==")
    say("   b2b = dphi > 3pi/4, near = dphi < pi/4, all = [0, pi]")
    for sel in ("acceptance", "full"):
        for sh in ("off", "on"):
            say(f"  [{sel}, showers {sh}]")
            for c, lab in (("w", "C_QQ"), ("w_tpos", "target +"), ("w_tneg", "target -"),
                           ("w_count", "pair count"), ("w_shuf", "C_QQ, target charges shuffled")):
                base = {r: reg[(sh, 0, sel, c, r)] for r in ("b2b", "near", "all")}
                dd = {r: diff(reg[(sh, 1, sel, c, r)], reg[(sh, 0, sel, c, r)]) for r in ("b2b", "near", "all")}
                say(f"    {lab:30s} baseline b2b {fmt(base['b2b'])}   Delta b2b {fmt(dd['b2b'])}  "
                    f"near {fmt(dd['near'])}  all {fmt(dd['all'])}")
            dq = diff(reg[(sh, 1, sel, "w", "b2b")], reg[(sh, 0, sel, "w", "b2b")])
            ds = diff(reg[(sh, 1, sel, "w_shuf", "b2b")], reg[(sh, 0, sel, "w_shuf", "b2b")])
            say(f"    charge-sensitive part of Delta b2b (C_QQ minus shuffled): {fmt(diff(dq, ds))}")
    say("\n== Showers on / off for the SRC-induced back-to-back change ==")
    for sel in ("acceptance", "full"):
        for c, lab in (("w", "C_QQ"), ("w_tneg", "target -")):
            d_off = diff(reg[("off", 1, sel, c, "b2b")], reg[("off", 0, sel, c, "b2b")])
            d_on = diff(reg[("on", 1, sel, c, "b2b")], reg[("on", 0, sel, c, "b2b")])
            r = d_on[0] / d_off[0]
            say(f"  [{sel}] {lab:9s} Delta b2b off {fmt(d_off)}  on {fmt(d_on)}  ratio {r:.3f}±"
                f"{abs(r) * math.hypot(d_off[1] / d_off[0], d_on[1] / d_on[0]):.3f}")
    say("\n== p_T threshold scan (both hadrons, Breit p_T > cut): Delta b2b of C_QQ and of the shuffled version ==")
    for sel in ("acceptance", "full"):
        for sh in ("off", "on"):
            say(f"  [{sel}, showers {sh}] " + "  ".join(
                f"pT>{pt}: {fmt(diff(scan[(sh, 1, sel, pt, 'w')], scan[(sh, 0, sel, pt, 'w')]))} "
                f"(shuf {diff(scan[(sh, 1, sel, pt, 'w_shuf')], scan[(sh, 0, sel, pt, 'w_shuf')])[0]:+.5f})"
                for pt in (0.0, 0.1, 0.2, 0.3, 0.5)))
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fs = 16
    for sel in ("acceptance", "full"):
        fig, ax = plt.subplots(2, 3, figsize=(17, 9.5), sharex=True)
        for r, sh in enumerate(("off", "on")):
            a = ax[r, 0]
            for src, col, mk, lab in ((0, "0.35", "o", "SRC off"), (1, "tab:red", "s", "SRC on")):
                m, e = hist[(sh, src, sel, "w")]
                a.errorbar(CENTERS, m, yerr=e, fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            a = ax[r, 1]
            for c, col, mk, lab in (("w", "k", "o", r"$\Delta C_{QQ}$"),
                                    ("w_shuf", "tab:green", "^", "target charges shuffled")):
                m1, e1 = hist[(sh, 1, sel, c)]
                m0, e0 = hist[(sh, 0, sel, c)]
                a.errorbar(CENTERS, m1 - m0, yerr=np.hypot(e1, e0), fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            a = ax[r, 2]
            for c, col, mk, lab in (("w_tpos", "tab:orange", "o", r"positive target hadrons"),
                                    ("w_tneg", "tab:blue", "s", r"negative target hadrons")):
                m1, e1 = hist[(sh, 1, sel, c)]
                m0, e0 = hist[(sh, 0, sel, c)]
                a.errorbar(CENTERS, m1 - m0, yerr=np.hypot(e1, e0), fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            for c_ in (1, 2):
                ax[r, c_].axhline(0, color="k", lw=0.8)
            ax[r, 0].set_ylabel(rf"showers {sh}" + "\n" + "per event per rad", fontsize=fs)
        ax[0, 0].set_title(r"$C_{QQ}(\Delta\phi)$", fontsize=fs)
        ax[0, 1].set_title(r"$\Delta C_{QQ}$ = SRC on $-$ SRC off", fontsize=fs)
        ax[0, 2].set_title(r"$\Delta C_{QQ}$ by target-hadron charge", fontsize=fs)
        for a in ax.flat:
            a.set_xlim(0, np.pi)
            a.set_xticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
            a.set_xticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
            a.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(np.pi / 16))
            a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            a.tick_params(which="both", direction="in", labelsize=fs)
            a.tick_params(which="major", length=8, width=1.1)
            a.tick_params(which="minor", length=4, width=0.9)
        for c_ in range(3):
            ax[1, c_].set_xlabel(r"$\Delta\phi$ [rad]", fontsize=fs)
            ax[0, c_].legend(fontsize=fs * 0.8, frameon=False, loc="best")
        fig.tight_layout()
        fig.savefig(OUT / f"qqc_{sel}.pdf")
        fig.savefig(OUT / f"qqc_{sel}.png", dpi=130)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
