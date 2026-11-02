#!/usr/bin/env python3
"""
Angular charge profile around the current jet axis for struck u vs struck d (SRC off).

  C_f(theta) = < sum_h Q_h delta(theta - theta_hJ) >_f ,

with theta_hJ the 3D opening angle between a final-state charged hadron and the anti-kT jet axis
in the Breit frame, f the true struck-quark flavour. Densities are per selected event (events
with a current jet) and per radian, in identical bins for both flavours; positive and negative
hadron densities are shown separately and C_f = positive - negative. Errors come from the
event-to-event spread of each bin. Input: production_charge/ (fig3_patched struck=any
struckid=2|1 hadrons=...), same DIS window for both flavours.

Reported per shower setting:
  * integrated net charge in the current hemisphere and in cones theta < 0.5, 1.0 rad;
  * full-covariance chi2 between the u and d net-charge profiles (current hemisphere, full range);
  * the same after reweighting the d sample to the u sample's (x, Q^2) distribution.
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
DIR = STUDY / "production_charge"
OUT = STUDY / "results_charge_profile"
NB = 36
EDGES = np.linspace(0.0, np.pi, NB + 1)
WIDTH = EDGES[1] - EDGES[0]
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])
X_EDGES = np.geomspace(0.05, 0.495, 11)
Q2_EDGES = np.linspace(16.0, 25.0, 6)


def load(sh: str, f: str) -> dict:
    ev = pd.read_csv(DIR / f"showers_{sh}_struck{f}.csv", usecols=["event_id", "xB", "Q2", "tagjet_E"])
    ev = ev[ev["tagjet_E"].notna()].reset_index(drop=True)
    had = pd.read_csv(DIR / f"showers_{sh}_struck{f}_hadrons.csv")
    idx = pd.Series(np.arange(len(ev)), index=ev["event_id"].to_numpy())
    had = had[had["event_id"].isin(idx.index)]
    row = idx.loc[had["event_id"].to_numpy()].to_numpy()
    b = np.clip(np.digitize(had["theta_hJ"].to_numpy(), EDGES) - 1, 0, NB - 1)
    mats = {}
    for region, sel in (("current", had["current"].to_numpy() == 1), ("full", np.ones(len(had), bool))):
        for sign, s in (("pos", had["charge"].to_numpy() > 0), ("neg", had["charge"].to_numpy() < 0)):
            m = sel & s
            flat = np.bincount(row[m] * NB + b[m], minlength=len(ev) * NB).astype(np.float32)
            mats[(region, sign)] = flat.reshape(len(ev), NB)
    return {"ev": ev, "mats": mats}


def wstats(X: np.ndarray, w: np.ndarray) -> tuple:
    """Weighted per-bin mean, standard error of the mean, and covariance of the mean."""
    sw = w.sum()
    mu = (w[:, None] * X).sum(0) / sw
    n_eff = sw**2 / (w**2).sum()
    cov = np.zeros((X.shape[1], X.shape[1]))
    for i in range(0, len(X), 200000):
        D = X[i:i + 200000] - mu
        cov += (D * w[i:i + 200000, None]).T @ D
    cov /= sw
    cov_mean = cov / n_eff
    return mu, np.sqrt(np.diag(cov_mean)), cov_mean


def chi2(mu_a, cov_a, mu_b, cov_b) -> tuple:
    d = mu_a - mu_b
    c = cov_a + cov_b
    keep = np.diag(c) > 0
    v = float(d[keep] @ np.linalg.solve(c[np.ix_(keep, keep)], d[keep]))
    return v, int(keep.sum())


def kin_weights(ev_d: pd.DataFrame, ev_u: pd.DataFrame) -> np.ndarray:
    """Weights making the d sample's (x, Q^2) distribution match the u sample's."""
    hu, _, _ = np.histogram2d(ev_u["xB"], ev_u["Q2"], bins=[X_EDGES, Q2_EDGES])
    hd, _, _ = np.histogram2d(ev_d["xB"], ev_d["Q2"], bins=[X_EDGES, Q2_EDGES])
    ratio = np.divide(hu / hu.sum(), hd / hd.sum(), out=np.zeros_like(hu), where=hd > 0)
    ix = np.clip(np.digitize(ev_d["xB"], X_EDGES) - 1, 0, len(X_EDGES) - 2)
    iq = np.clip(np.digitize(ev_d["Q2"], Q2_EDGES) - 1, 0, len(Q2_EDGES) - 2)
    return ratio[ix, iq]


def style(a, fs):
    a.set_xlim(0, np.pi)
    a.set_xticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    a.set_xticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
    a.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(np.pi / 16))
    a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.tick_params(which="both", direction="in", labelsize=fs)
    a.tick_params(which="major", length=8, width=1.1)
    a.tick_params(which="minor", length=4, width=0.9)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    res = {}
    for sh in ("off", "on"):
        data = {f: load(sh, f) for f in ("u", "d")}
        wd = kin_weights(data["d"]["ev"], data["u"]["ev"])
        say(f"\n== showers {sh} ==")
        for f in ("u", "d"):
            ev = data[f]["ev"]
            say(f"  struck {f}: events with a current jet {len(ev):,}  <x>={ev['xB'].mean():.4f}  <Q2>={ev['Q2'].mean():.3f} GeV^2")
        say(f"  d -> u (x, Q2) reweighting: effective d sample {wd.sum() ** 2 / (wd ** 2).sum():,.0f} events, "
            f"reweighted <x>_d = {np.average(data['d']['ev']['xB'], weights=wd):.4f}")
        for region in ("current", "full"):
            stats = {}
            for f, w in (("u", None), ("d", None), ("d_rw", wd)):
                ff = "d" if f == "d_rw" else f
                m = data[ff]["mats"]
                ww = np.ones(len(data[ff]["ev"])) if w is None else w
                stats[f] = {k: wstats(m[(region, k)] if k != "net" else m[(region, "pos")] - m[(region, "neg")], ww)
                            for k in ("pos", "neg", "net")}
            res[(sh, region)] = stats
            say(f"  [{region} {'hemisphere' if region == 'current' else 'angular range'}]")
            for f in ("u", "d", "d_rw"):
                mu, err, _ = stats[f]["net"]
                pos, neg = stats[f]["pos"][0].sum(), stats[f]["neg"][0].sum()
                tot_err = math.sqrt(stats[f]["net"][2].sum())
                cones = []
                for tc in (0.5, 1.0):
                    k = CENTERS < tc
                    cones.append(f"Q(theta<{tc})={mu[k].sum():+.4f}±{math.sqrt(stats[f]['net'][2][np.ix_(k, k)].sum()):.4f}")
                say(f"    {f:5s} <n+>={pos:.3f} <n->={neg:.3f}  net charge {mu.sum():+.4f}±{tot_err:.4f}   " + "  ".join(cones))
            for a, b in (("u", "d"), ("u", "d_rw")):
                c2, ndf = chi2(stats[a]["net"][0], stats[a]["net"][2], stats[b]["net"][0], stats[b]["net"][2])
                say(f"    net-charge profile {a} vs {b}: chi2 = {c2:,.0f} / {ndf} bins")
            # Per-event separation of the integrated current-hemisphere (or full) net charge.
            qu = (data["u"]["mats"][(region, "pos")] - data["u"]["mats"][(region, "neg")]).sum(1)
            qd = (data["d"]["mats"][(region, "pos")] - data["d"]["mats"][(region, "neg")]).sum(1)
            pooled = math.sqrt(0.5 * (qu.var() + qd.var()))
            dprime = (qu.mean() - qd.mean()) / pooled if pooled > 0 else 0.0
            if dprime > 0:
                n5 = 2 * (5.0 / dprime) ** 2
                say(f"    per-event: <Q>_u - <Q>_d = {qu.mean() - qd.mean():+.4f}, per-event spread {pooled:.3f}, "
                    f"effect size {dprime:.3f} -> 5 sigma with ~{n5:,.0f} events per flavour")
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Figures: rows = showers off / on, columns = positive / negative / net; u vs d overlaid.
    fs = 16
    for region in ("current", "full"):
        fig, ax = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
        for r, sh in enumerate(("off", "on")):
            for c, (k, lab) in enumerate((("pos", r"positive hadrons"), ("neg", r"negative hadrons"),
                                          ("net", r"net charge $C_f(\theta)$"))):
                a = ax[r, c]
                for f, col, mk in (("u", "tab:red", "o"), ("d", "tab:blue", "s")):
                    mu, err, _ = res[(sh, region)][f][k]
                    a.errorbar(CENTERS, mu / WIDTH, yerr=err / WIDTH, fmt=mk, ms=3.5, color=col, ls="none",
                               label=rf"struck ${f}$")
                if k == "net":
                    a.axhline(0.0, color="k", lw=0.8)
                style(a, fs)
                if r == 0:
                    a.set_title(lab, fontsize=fs)
                if c == 0:
                    a.set_ylabel(rf"showers {sh}" + "\n" + r"per event per rad", fontsize=fs)
                if r == 1:
                    a.set_xlabel(r"$\theta_{hJ}$ [rad]", fontsize=fs)
        handles, labels = ax[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, fontsize=fs, frameon=False, loc="upper center", ncol=2)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(OUT / f"charge_profile_{region}.pdf")
        fig.savefig(OUT / f"charge_profile_{region}.png", dpi=130)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
