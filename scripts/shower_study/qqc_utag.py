#!/usr/bin/env python3
"""
u-tagged azimuthal QQC (current x forward-target charged-hadron pairs, Breit frame).

Tag (experimental): leading current-hemisphere charged hadron (largest z_h = P.P_h/P.q) is a
pi+ or K+ with z_h > 0.5; no generator-level flavour is used. Input: production_qqc_utag/
(fig3_patched utagonly=1, 8M selected DIS events per sample) and, for the untagged comparison,
production_qqc/ (2M DIS events per sample, all events).

Target region: original forward acceptance (proton-forward lab 2 < eta < 3.5, z_LC > 0.1).
Current-side sums for tagged events: all current charged hadrons ("incl. tag"), all except the
tagging hadron ("excl. tag"), and the tagging hadron alone ("tag only").

Per bin: raw C_QQ per event per rad, pair density N_pairs per event per rad, and the per-pair
average charge product C_QQ/N_pairs (event-level delta-method errors). Back-to-back = dphi > 3pi/4.
Significances are scaled to 10^6 selected DIS events so tagged and untagged compare fairly.
"""

from __future__ import annotations

import os
import json
import math
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
TAG = STUDY / "production_qqc_utag"
UNTAG = STUDY / "production_qqc"
OUT = STUDY / "results_qqc_utag"
NB = 36
EDGES = np.linspace(0.0, np.pi, NB + 1)
WIDTH = EDGES[1] - EDGES[0]
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])
B2B = CENTERS > 0.75 * np.pi


def n_dis(directory: Path, name: str) -> int:
    tot = 0
    for js in directory.glob(f"chunks/{name}_c*.csv.summary.json"):
        tot += int(json.loads(js.read_text())["n_rows"])
    return tot


def load(directory: Path, name: str) -> tuple:
    ev = pd.read_csv(directory / f"{name}.csv", usecols=["event_id"])["event_id"].to_numpy()
    had = pd.read_csv(directory / f"{name}_qqc.csv")
    idx = pd.Series(np.arange(len(ev)), index=ev)
    had["ev"] = idx.loc[had["event_id"].to_numpy()].to_numpy()
    if "is_lead" not in had:
        had["is_lead"] = 0
    return len(ev), had


def pairs(had: pd.DataFrame, cur_variant: str) -> pd.DataFrame:
    cur = had["current"].to_numpy() == 1
    lead = had["is_lead"].to_numpy() == 1
    cm = {"incl": cur, "excl": cur & ~lead, "only": cur & lead}[cur_variant]
    eta, z = had["eta_eic"].to_numpy(), had["zlc"].to_numpy()
    tm = (had["current"].to_numpy() == 0) & (eta > 2.0) & (eta < 3.5) & (z > 0.1)
    a = had.loc[cm, ["ev", "charge", "phi"]]
    b = had.loc[tm, ["ev", "charge", "phi"]]
    m = a.merge(b, on="ev", suffixes=("_a", "_b"), sort=True)
    d = np.abs(m["phi_a"].to_numpy() - m["phi_b"].to_numpy()) % (2 * np.pi)
    d = np.where(d > np.pi, 2 * np.pi - d, d)
    qb = m["charge_b"].to_numpy()
    w = m["charge_a"].to_numpy() * qb
    return pd.DataFrame({"ev": m["ev"].to_numpy(), "bin": np.clip(np.digitize(d, EDGES) - 1, 0, NB - 1),
                         "w": w, "w_tneg": np.where(qb < 0, w, 0.0), "n": np.ones(len(m))})


def stats(n_evt: int, p: pd.DataFrame) -> dict:
    """Per-bin and back-to-back means with event-level errors, incl. the per-pair ratio."""
    ev, b = p["ev"].to_numpy(), p["bin"].to_numpy()
    out = {}
    acc = {k: np.zeros(NB) for k in ("a", "n", "t", "aa", "nn", "an", "tt")}
    step = 250000
    for lo in range(0, n_evt, step):
        i0, i1 = np.searchsorted(ev, [lo, lo + step])
        k = min(step, n_evt - lo)
        idx = (ev[i0:i1] - lo) * NB + b[i0:i1]
        fa = np.bincount(idx, weights=p["w"].to_numpy()[i0:i1], minlength=k * NB).reshape(k, NB)
        fn = np.bincount(idx, weights=p["n"].to_numpy()[i0:i1], minlength=k * NB).reshape(k, NB)
        ft = np.bincount(idx, weights=p["w_tneg"].to_numpy()[i0:i1], minlength=k * NB).reshape(k, NB)
        acc["a"] += fa.sum(0); acc["n"] += fn.sum(0); acc["t"] += ft.sum(0)
        acc["aa"] += (fa ** 2).sum(0); acc["nn"] += (fn ** 2).sum(0); acc["an"] += (fa * fn).sum(0)
        acc["tt"] += (ft ** 2).sum(0)
    N = n_evt
    ma, mn, mt = acc["a"] / N, acc["n"] / N, acc["t"] / N
    va, vn, vt = acc["aa"] / N - ma ** 2, acc["nn"] / N - mn ** 2, acc["tt"] / N - mt ** 2
    cov = acc["an"] / N - ma * mn
    R = np.divide(ma, mn, out=np.full(NB, np.nan), where=mn > 0)
    vR = np.divide(va - 2 * R * cov + R ** 2 * vn, N * mn ** 2, out=np.full(NB, np.nan), where=mn > 0)
    out["C"] = (ma / WIDTH, np.sqrt(va / N) / WIDTH)
    out["N"] = (mn / WIDTH, np.sqrt(vn / N) / WIDTH)
    out["R"] = (R, np.sqrt(np.maximum(vR, 0)))
    # Back-to-back integrals from per-event region sums.
    sel = B2B[b]
    sa = np.bincount(ev[sel], weights=p["w"].to_numpy()[sel], minlength=N)
    sn = np.bincount(ev[sel], weights=p["n"].to_numpy()[sel], minlength=N)
    st = np.bincount(ev[sel], weights=p["w_tneg"].to_numpy()[sel], minlength=N)
    Rb = sa.mean() / sn.mean()
    vRb = (sa.var() - 2 * Rb * np.cov(sa, sn, bias=True)[0, 1] + Rb ** 2 * sn.var()) / (N * sn.mean() ** 2)
    out["C_b2b"] = (sa.mean(), sa.std() / math.sqrt(N))
    out["N_b2b"] = (sn.mean(), sn.std() / math.sqrt(N))
    out["T_b2b"] = (st.mean(), st.std() / math.sqrt(N))
    out["R_b2b"] = (Rb, math.sqrt(vRb))
    return out


def d(a: tuple, b: tuple) -> tuple:
    return a[0] - b[0], math.hypot(a[1], b[1])


def f(v: tuple, n: int = 5) -> str:
    return f"{v[0]:+.{n}f}±{v[1]:.{n}f}"


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
    ndis = {}
    for sh in ("off", "on"):
        for src in (0, 1):
            name = f"showers_{sh}_src{src}"
            ndis[("tag", sh, src)] = n_dis(TAG, name)
            ndis[("untag", sh, src)] = n_dis(UNTAG, name)
            n_t, had_t = load(TAG, name)
            q, cur, lead = had_t["charge"].to_numpy(), had_t["current"].to_numpy() == 1, had_t["is_lead"].to_numpy() == 1
            say(f"\n== showers {sh}, SRC {'on' if src else 'off'}: tagged {n_t:,} of {ndis[('tag', sh, src)]:,} DIS events "
                f"(efficiency {n_t / ndis[('tag', sh, src)]:.4f}); tagging hadron charge +1 in "
                f"{np.mean(q[lead & cur] == 1):.4f}; <Q_current> incl. tag {np.bincount(had_t['ev'][cur], weights=q[cur], minlength=n_t).mean():+.3f}, "
                f"excl. tag {np.bincount(had_t['ev'][cur & ~lead], weights=q[cur & ~lead], minlength=n_t).mean():+.3f}")
            for var in ("incl", "excl", "only"):
                res[("tag", var, sh, src)] = stats(n_t, pairs(had_t, var))
            del had_t
            n_u, had_u = load(UNTAG, name)
            res[("untag", "incl", sh, src)] = stats(n_u, pairs(had_u, "incl"))
            del had_u

    say("\n== Back-to-back (dphi > 3pi/4), SRC on - SRC off ==")
    say("   raw C_QQ per selected event; per-pair = average charge product; Z = significance per 10^6 selected DIS events")
    summary = {}
    for sample, var, lab in (("untag", "incl", "untagged, all current hadrons"), ("tag", "incl", "u-tagged, current incl. tag"),
                             ("tag", "excl", "u-tagged, current excl. tag"), ("tag", "only", "u-tagged, tagging hadron only")):
        say(f"  [{lab}]")
        for sh in ("off", "on"):
            r0, r1 = res[(sample, var, sh, 0)], res[(sample, var, sh, 1)]
            dC, dR, dN, dT = d(r1["C_b2b"], r0["C_b2b"]), d(r1["R_b2b"], r0["R_b2b"]), d(r1["N_b2b"], r0["N_b2b"]), d(r1["T_b2b"], r0["T_b2b"])
            scale = math.sqrt(1e6 / ndis[(sample, sh, 0)])
            zc, zr = dC[0] / dC[1] * scale, dR[0] / dR[1] * scale
            summary[(sample, var, sh)] = (dC, dR, zc, zr)
            say(f"    showers {sh}: baseline C {f(r0['C_b2b'])}  dC {f(dC)} (Z/1e6 DIS {zc:+.1f})  "
                f"baseline per-pair {f(r0['R_b2b'], 4)}  d(per-pair) {f(dR, 4)} (Z/1e6 DIS {zr:+.1f})  "
                f"dN_pairs {f(dN, 4)}  dC from negative target hadrons {f(dT)}")
        for k, (lab2, i) in enumerate((("raw", 0), ("per-pair", 1))):
            off, on = summary[(sample, var, "off")][i], summary[(sample, var, "on")][i]
            r = on[0] / off[0]
            say(f"    survival (on/off), {lab2}: {r:.3f}±{abs(r) * math.hypot(off[1] / off[0], on[1] / on[0]):.3f}")
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fs = 15
    # Main figure: u-tagged (current excl. tag), rows showers off / on.
    for var, tagname in (("excl", "excl_tag"), ("incl", "incl_tag")):
        fig, ax = plt.subplots(2, 4, figsize=(20, 9), sharex=True)
        for r, sh in enumerate(("off", "on")):
            for src, col, mk, lab in ((0, "0.35", "o", "SRC off"), (1, "tab:red", "s", "SRC on")):
                m, e = res[("tag", var, sh, src)]["C"]
                ax[r, 0].errorbar(CENTERS, m, yerr=e, fmt=mk, ms=3.5, color=col, ls="none", label=lab)
                m, e = res[("tag", var, sh, src)]["R"]
                ax[r, 2].errorbar(CENTERS, m, yerr=e, fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            for sample, v, col, mk, lab in (("tag", var, "tab:purple", "s", r"$u$-tagged"),
                                            ("untag", "incl", "0.5", "o", "untagged")):
                for c_, key in ((1, "C"), (3, "R")):
                    m1, e1 = res[(sample, v, sh, 1)][key]
                    m0, e0 = res[(sample, v, sh, 0)][key]
                    ax[r, c_].errorbar(CENTERS, m1 - m0, yerr=np.hypot(e1, e0), fmt=mk, ms=3.5, color=col,
                                       ls="none", label=lab)
            for c_ in range(4):
                ax[r, c_].axhline(0, color="k", lw=0.7)
            ax[r, 0].set_ylabel(f"showers {sh}", fontsize=fs)
        titles = (r"$C_{QQ}$ per tagged event per rad", r"$\Delta C_{QQ}$ = SRC on $-$ off",
                  r"$C_{QQ}/N_{\rm pairs}$ (charge product per pair)", r"$\Delta(C_{QQ}/N_{\rm pairs})$")
        for c_, t in enumerate(titles):
            ax[0, c_].set_title(t, fontsize=fs)
            ax[1, c_].set_xlabel(r"$\Delta\phi$ [rad]", fontsize=fs)
            ax[0, c_].legend(fontsize=fs * 0.8, frameon=False, loc="lower left")
        for a in ax.flat:
            style(a, fs)
        fig.tight_layout()
        fig.savefig(OUT / f"qqc_utag_{tagname}.pdf")
        fig.savefig(OUT / f"qqc_utag_{tagname}.png", dpi=120)
    # Tag-hadron bias figure: SRC-off baselines and SRC changes for the three current-side sums.
    fig, ax = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    for r, sh in enumerate(("off", "on")):
        for var, col, mk, lab in (("incl", "k", "o", "current incl. tag"), ("excl", "tab:blue", "s", "current excl. tag"),
                                  ("only", "tab:orange", "^", "tagging hadron only")):
            m0, e0 = res[("tag", var, sh, 0)]["R"]
            ax[r, 0].errorbar(CENTERS, m0, yerr=e0, fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            m1, e1 = res[("tag", var, sh, 1)]["C"]
            m0c, e0c = res[("tag", var, sh, 0)]["C"]
            ax[r, 1].errorbar(CENTERS, m1 - m0c, yerr=np.hypot(e1, e0c), fmt=mk, ms=3.5, color=col, ls="none", label=lab)
            m1, e1 = res[("tag", var, sh, 1)]["R"]
            ax[r, 2].errorbar(CENTERS, m1 - m0, yerr=np.hypot(e1, e0), fmt=mk, ms=3.5, color=col, ls="none", label=lab)
        for c_ in range(3):
            ax[r, c_].axhline(0, color="k", lw=0.7)
        ax[r, 0].set_ylabel(f"showers {sh}", fontsize=fs)
    for c_, t in enumerate((r"SRC off: $C_{QQ}/N_{\rm pairs}$", r"$\Delta C_{QQ}$ per tagged event per rad",
                            r"$\Delta(C_{QQ}/N_{\rm pairs})$")):
        ax[0, c_].set_title(t, fontsize=fs)
        ax[1, c_].set_xlabel(r"$\Delta\phi$ [rad]", fontsize=fs)
        ax[0, c_].legend(fontsize=fs * 0.8, frameon=False, loc="lower left")
    for a in ax.flat:
        style(a, fs)
    fig.tight_layout()
    fig.savefig(OUT / "qqc_utag_tag_bias.pdf")
    fig.savefig(OUT / "qqc_utag_tag_bias.png", dpi=120)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
