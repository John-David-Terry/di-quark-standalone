#!/usr/bin/env python3
"""
u-tagged QQC with the 3D EEC opening angle chi instead of the azimuthal separation dphi.

Same events and pairs as qqc_utag.py "tagging hadron only": current side = the tagging hadron
(leading current-hemisphere charged hadron by z_h, a pi+ or K+ with z_h > 0.5), target side =
charged hadrons in the original forward acceptance (lab 2 < eta < 3.5, z_LC > 0.1).
Input: production_qqc_utag/ (no regeneration). Breit 4-momenta are rebuilt exactly from the
stored quantities: p_T, phi, p_z = z_LC Q/(2x) (q.P_h = Q p_z, q.P = Q^2/(2x)) and
E = (z_h Q^2/(2x) + P_z p_z)/E_P.

  cos chi_ab = p_a . p_b / (|p_a| |p_b|),   back-to-back: chi -> pi, i.e. (1 + cos chi)/2 -> 0.

The Breit-frame chi is the observable. As a diagnostic of the longitudinal mismatch between the
current and target hadrons, chi is also computed in the gamma*p c.m. frame (a z-boost of the
Breit frame), where a back-to-back transverse kick of struck quark and remnant is back-to-back
in 3D. Raw C_QQ is per tagged event per unit ln(1 + cos chi); C_QQ/N_pairs is binning-independent.
For a fair chi vs dphi comparison, "back-to-back" is also defined as the 25% of SRC-off pairs
closest to back-to-back in each variable.
"""

from __future__ import annotations

import json
import math
import os
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

# Real LaTeX text when a TeX installation is available, matplotlib mathtext otherwise.
USETEX = shutil.which("latex") is not None
if USETEX:
    plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                         "text.latex.preamble": r"\usepackage{amsmath}"})

M_P = 0.938272
STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
TAG = STUDY / "production_qqc_utag"
OUT = STUDY / "results_qqc_chi"
U_EDGES = np.geomspace(1e-4, 10 ** -0.2, 20)   # (1 + cos chi)/2, up to chi ~ 75 deg (few pairs beyond)
U_CENT = np.sqrt(U_EDGES[:-1] * U_EDGES[1:])
LN_WIDTH = np.diff(np.log(U_EDGES))
A_EDGES = np.linspace(0.0, np.pi, 37)          # dphi and chi on a common linear axis
A_CENT = 0.5 * (A_EDGES[:-1] + A_EDGES[1:])
MIN_PAIR_FRAC = 0.002                          # linear-angle per-pair plot: drop bins with < 0.2% of the pairs
PI_TICKS = ([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi], [r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])


def n_dis(name: str) -> int:
    """Selected DIS events behind a utagonly sample (its events CSV holds the tagged ones only)."""
    return sum(int(json.loads(js.read_text())["n_rows"]) for js in TAG.glob(f"chunks/{name}_c*.csv.summary.json"))


def load_pairs(sh: str, src: int) -> tuple:
    name = f"showers_{sh}_src{src}"
    ev = pd.read_csv(TAG / f"{name}.csv", usecols=["event_id", "xB", "Q2"])
    had = pd.read_csv(TAG / f"{name}_qqc.csv")
    idx = pd.Series(np.arange(len(ev)), index=ev["event_id"].to_numpy())
    i = idx.loc[had["event_id"].to_numpy()].to_numpy()
    q, x = np.sqrt(ev["Q2"].to_numpy())[i], ev["xB"].to_numpy()[i]
    pz_p = q / (2 * x)
    e_p = np.hypot(pz_p, M_P)
    beta = (pz_p - q) / e_p                    # Breit -> gamma*p c.m. (P + q has no transverse momentum)
    pz = had["zlc"].to_numpy() * pz_p
    e = (had["zh"].to_numpy() * q * pz_p + pz_p * pz) / e_p
    had["ev"] = i
    had["pz"] = pz
    had["pzcm"] = (pz - beta * e) / np.sqrt(1 - beta ** 2)
    keep = ["ev", "charge", "phi", "pt", "pz", "pzcm"]
    lead = had.loc[(had["current"] == 1) & (had["is_lead"] == 1), keep]
    tgt = had.loc[(had["current"] == 0) & (had["eta_eic"] > 2.0) & (had["eta_eic"] < 3.5) & (had["zlc"] > 0.1), keep]
    m = lead.merge(tgt, on="ev", suffixes=("_a", "_b"), sort=True)
    dphi = np.abs(m["phi_a"].to_numpy() - m["phi_b"].to_numpy()) % (2 * np.pi)
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    pta, ptb = m["pt_a"].to_numpy(), m["pt_b"].to_numpy()
    out = {"ev": m["ev"].to_numpy(), "dphi": dphi, "w": m["charge_a"].to_numpy() * m["charge_b"].to_numpy()}
    for col, frame in (("pz", ""), ("pzcm", "cm")):
        pza, pzb = m[f"{col}_a"].to_numpy(), m[f"{col}_b"].to_numpy()
        cos = (pta * ptb * np.cos(dphi) + pza * pzb) / (np.hypot(pta, pza) * np.hypot(ptb, pzb))
        out["u" + frame] = np.clip((1 + cos) / 2, 1e-12, None)
        out["chi" + frame] = np.arccos(np.clip(cos, -1, 1))
        out["tha" + frame] = np.arctan2(pta, -pza)       # tagging hadron, polar angle from -z
        out["thb" + frame] = np.arctan2(ptb, pzb)        # target hadron, polar angle from +z
    xev = ev["xB"].to_numpy()
    return len(ev), pd.DataFrame(out), float(np.median(xev / (1 - xev)))


def binned(n_evt: int, ev: np.ndarray, b: np.ndarray, w: np.ndarray, nb: int) -> dict:
    """Per-bin raw mean per event, pair density, per-pair ratio; event-level (delta-method) errors."""
    ok = (b >= 0) & (b < nb)
    ev, b, w = ev[ok], b[ok], w[ok]
    acc = {k: np.zeros(nb) for k in ("a", "n", "aa", "nn", "an")}
    step = 250000
    for lo in range(0, n_evt, step):
        i0, i1 = np.searchsorted(ev, [lo, lo + step])
        k = min(step, n_evt - lo)
        idx = (ev[i0:i1] - lo) * nb + b[i0:i1]
        fa = np.bincount(idx, weights=w[i0:i1], minlength=k * nb).reshape(k, nb)
        fn = np.bincount(idx, minlength=k * nb).reshape(k, nb).astype(float)
        acc["a"] += fa.sum(0); acc["n"] += fn.sum(0)
        acc["aa"] += (fa ** 2).sum(0); acc["nn"] += (fn ** 2).sum(0); acc["an"] += (fa * fn).sum(0)
    N = n_evt
    ma, mn = acc["a"] / N, acc["n"] / N
    va, vn, cov = acc["aa"] / N - ma ** 2, acc["nn"] / N - mn ** 2, acc["an"] / N - ma * mn
    R = np.divide(ma, mn, out=np.full(nb, np.nan), where=mn > 0)
    vR = np.divide(va - 2 * R * cov + R ** 2 * vn, N * mn ** 2, out=np.full(nb, np.nan), where=mn > 0)
    return {"C": (ma, np.sqrt(va / N)), "N": (mn, np.sqrt(vn / N)), "R": (R, np.sqrt(np.maximum(vR, 0)))}


def region(n_evt: int, p: pd.DataFrame, sel: np.ndarray) -> dict:
    ev, w = p["ev"].to_numpy()[sel], p["w"].to_numpy()[sel]
    sa = np.bincount(ev, weights=w, minlength=n_evt)
    sn = np.bincount(ev, minlength=n_evt).astype(float)
    R = sa.mean() / sn.mean()
    vR = (sa.var() - 2 * R * np.cov(sa, sn, bias=True)[0, 1] + R ** 2 * sn.var()) / (n_evt * sn.mean() ** 2)
    return {"C": (sa.mean(), sa.std() / math.sqrt(n_evt)), "N": (sn.mean(), sn.std() / math.sqrt(n_evt)),
            "R": (R, math.sqrt(vR))}


def regions_for(p0: pd.DataFrame) -> list:
    """(key, label, variable, cut, below) with quantile cuts taken from the SRC-off pairs."""
    deg = lambda u: math.degrees(math.acos(2 * u - 1))  # noqa: E731
    u170, u175 = ((1 + math.cos(math.radians(a))) / 2 for a in (170, 175))
    qb, qc = np.quantile(p0["u"], 0.25), np.quantile(p0["ucm"], 0.25)
    phi75 = np.quantile(p0["dphi"], 0.75)
    return [
        ("phi_3pi4", "dphi > 3pi/4 (previous azimuthal analysis)", "dphi", 0.75 * math.pi, False),
        ("phi_q25", f"dphi, 25% most b2b (dphi > {math.degrees(phi75):.1f} deg)", "dphi", phi75, False),
        ("chi_q25", f"chi Breit, 25% most b2b (chi > {deg(qb):.1f} deg)", "u", qb, True),
        ("chi170", "chi Breit > 170 deg", "u", u170, True),
        ("chi175", "chi Breit > 175 deg", "u", u175, True),
        ("chicm_q25", f"chi g*p c.m., 25% most b2b (chi > {deg(qc):.1f} deg)", "ucm", qc, True),
        ("chicm170", "chi g*p c.m. > 170 deg", "ucm", u170, True),
        ("chicm175", "chi g*p c.m. > 175 deg", "ucm", u175, True),
    ]


def d(a, b):
    return a[0] - b[0], math.hypot(a[1], b[1])


def style(a, fs, minor_x=True):
    if minor_x:
        a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.tick_params(which="both", direction="in", labelsize=fs)
    a.tick_params(which="major", length=10, width=1.2)
    a.tick_params(which="minor", length=5, width=1.0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    data, hist, regs = {}, {}, {}
    for sh in ("off", "on"):
        for src in (0, 1):
            n, p, xr = load_pairs(sh, src)
            data[(sh, src)] = (n, p, xr)
            ev, w = p["ev"].to_numpy(), p["w"].to_numpy()
            for key in ("u", "ucm"):
                hist[(sh, src, key)] = binned(n, ev, np.digitize(p[key].to_numpy(), U_EDGES) - 1, w, len(U_CENT))
            for key in ("dphi", "chi", "chicm"):
                b = np.clip(np.digitize(p[key].to_numpy(), A_EDGES) - 1, 0, len(A_CENT) - 1)
                hist[(sh, src, key)] = binned(n, ev, b, w, len(A_CENT))

    say("Tagging hadron x forward-acceptance target hadrons, u-tagged events (production_qqc_utag)")
    say("Columns: pair fraction (SRC off) | per-pair C/N SRC off | d(C/N) | Z per 1e6 DIS | raw dC per tagged event | dN/N")
    for sh in ("off", "on"):
        n0, p0, xr = data[(sh, 0)]
        ndis = n_dis(f"showers_{sh}_src0")
        say(f"\n== showers {sh}: tagged events SRC off {n0:,}, SRC on {data[(sh, 1)][0]:,}; "
            f"pairs per tagged event {len(p0) / n0:.3f} / {len(data[(sh, 1)][1]) / data[(sh, 1)][0]:.3f} ==")
        for key, lab, var, cut, below in regions_for(p0):
            r = {}
            for src in (0, 1):
                n, p, _ = data[(sh, src)]
                x = p[var].to_numpy()
                r[src] = region(n, p, x < cut if below else x > cut)
            frac = np.mean(p0[var].to_numpy() < cut if below else p0[var].to_numpy() > cut)
            dR, dC, dN = d(r[1]["R"], r[0]["R"]), d(r[1]["C"], r[0]["C"]), d(r[1]["N"], r[0]["N"])
            regs[(sh, key)] = (lab.split(" (")[0] if "25%" in lab else lab, dR)
            say(f"  {lab:50s} {frac:5.3f} | {r[0]['R'][0]:+.3f} | {dR[0]:+.4f}±{dR[1]:.4f} | "
                f"{dR[0] / dR[1] * math.sqrt(1e6 / ndis):+5.1f} | {dC[0]:+.5f}±{dC[1]:.5f} | {dN[0] / r[0]['N'][0]:+.2f}")
        say("  longitudinal mismatch (medians; tagging hadron polar angle from -z, target hadron from +z):")
        for frame, nm in (("", "Breit"), ("cm", "g*p c.m.")):
            ta = {s: np.degrees(np.median(data[(sh, s)][1]["tha" + frame])) for s in (0, 1)}
            tb = np.degrees(np.median(p0["thb" + frame]))
            ratio = np.median(p0["thb" + frame] / p0["tha" + frame])
            reach = np.mean(np.pi - np.abs(p0["tha" + frame] - p0["thb" + frame]) > math.radians(170))
            say(f"    {nm:9s} theta_a {ta[0]:5.1f} (SRC on {ta[1]:5.1f}) deg, theta_b {tb:5.1f} deg, "
                f"theta_b/theta_a {ratio:.2f}; pairs able to reach chi > 170 deg at dphi = pi: {reach:.3f}")
        say(f"    parton-level expectation for a pure kT kick in the Breit frame: theta_b/theta_a = x/(1-x), median {xr:.2f}")
        b_phi = p0["dphi"].to_numpy() > np.quantile(p0["dphi"], 0.75)
        for key, nm in (("u", "Breit"), ("ucm", "g*p c.m.")):
            b_chi = p0[key].to_numpy() < np.quantile(p0[key], 0.25)
            say(f"    of the 25% most dphi-b2b pairs, {np.mean(b_chi[b_phi]):.3f} are in the 25% most chi-b2b ({nm})")

    say("\n== Survival: d(C/N) showers on / showers off ==")
    for key, *_ in regions_for(data[("off", 0)][1]):
        (lab, off), (_, on) = regs[("off", key)], regs[("on", key)]
        r = on[0] / off[0]
        say(f"  {lab:50s} {r:.2f}±{abs(r) * math.hypot(off[1] / off[0], on[1] / on[0]):.2f}")
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 1) raw C_QQ and C_QQ/N_pairs vs (1 + cos chi)/2, SRC on/off; rows: showers off/on.
    fs = 18
    fig, ax = plt.subplots(2, 2, figsize=(14, 10.5), sharex=True)
    for r_, sh in enumerate(("off", "on")):
        for src, col, mk, lab in ((0, "0.3", "o", "SRC off"), (1, "tab:red", "s", "SRC on")):
            m, e = hist[(sh, src, "u")]["C"]
            ax[r_, 0].errorbar(U_CENT, m / LN_WIDTH, yerr=e / LN_WIDTH, fmt=mk, ms=5, color=col, ls="none", label=lab)
            m, e = hist[(sh, src, "u")]["R"]
            ax[r_, 1].errorbar(U_CENT, m, yerr=e, fmt=mk, ms=5, color=col, ls="none", label=lab)
        ax[r_, 0].set_ylabel(r"$dC_{QQ}/d\ln(1+\cos\chi)$", fontsize=fs)
        ax[r_, 1].set_ylabel(r"$C_{QQ}/N_{\rm pairs}$", fontsize=fs)
        ax[r_, 0].text(0.05, 0.9, f"showers {sh}", transform=ax[r_, 0].transAxes, fontsize=fs)
    ax[0, 1].legend(fontsize=fs, frameon=False, loc="lower right")
    for a in ax.flat:
        a.set_xscale("log")
        a.set_xlim(1e-4, U_EDGES[-1])
        style(a, fs, minor_x=False)
    for a in ax[1]:
        a.set_xlabel(r"$(1+\cos\chi)/2$", fontsize=fs)
    fig.tight_layout()
    fig.savefig(OUT / "qqc_chi_raw_and_per_pair.pdf")
    fig.savefig(OUT / "qqc_chi_raw_and_per_pair.png", dpi=110)

    # 1b) vs z = (1 + cos chi)/2, SRC on/off; rows: showers off/on.
    dz = np.diff(U_EDGES)

    def z_panels(values, ylabel, fname, legend_loc="upper right", ylim=None):
        fs = 20
        fig, ax = plt.subplots(2, 1, figsize=(9, 11), sharex=True)
        for r_, sh in enumerate(("off", "on")):
            for src, col, mk, lab in ((0, "0.3", "o", "SRC off"), (1, "tab:red", "s", "SRC on")):
                m, e = values(sh, src)
                ax[r_].errorbar(U_CENT, m, yerr=e, fmt=mk, ms=6, color=col, ls="none", label=lab)
            ax[r_].set_ylabel(ylabel, fontsize=fs)
            ax[r_].text(0.95, 0.88, f"showers {sh}", transform=ax[r_].transAxes, fontsize=fs, ha="right")
            ax[r_].set_xscale("log")
            ax[r_].set_xlim(1e-4, U_EDGES[-1])
            if ylim:
                ax[r_].set_ylim(*ylim)
            style(ax[r_], fs, minor_x=False)
        ax[0].legend(fontsize=fs, frameon=False, loc=legend_loc, bbox_to_anchor=(0.0, 0.0, 1.0, 0.85))
        ax[1].set_xlabel(r"$z$", fontsize=fs)
        fig.tight_layout()
        fig.savefig(OUT / f"{fname}.pdf")
        fig.savefig(OUT / f"{fname}.png", dpi=110)

    def per_dz(sh, src, total=False):
        m, e = hist[(sh, src, "u")]["C"]
        n, p, _ = data[(sh, src)]
        norm = dz * (len(p) / n if total else 1.0)   # total: divide by all tag-target pairs per tagged event
        return m / norm, e / norm

    z_panels(per_dz, r"$dC_{QQ}/dz$", "qqc_chi_raw_z")
    z_panels(lambda sh, src: per_dz(sh, src, total=True), (r"$\textrm{QQC}(z)$" if USETEX else r"$\mathrm{QQC}(z)$"), "qqc_chi_raw_z_pairnorm")
    z_panels(lambda sh, src: hist[(sh, src, "u")]["R"], r"$C_{QQ}(z)/N_{\rm pair}(z)$", "qqc_chi_per_pair_z",
             legend_loc="upper center", ylim=(0, 0.45))

    # 2) SRC-induced change, showers off vs on.
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(16, 6.2))
    for sh, col, mk in (("off", "tab:red", "o"), ("on", "tab:blue", "s")):
        for c_, key, scale in ((0, "C", LN_WIDTH), (1, "R", 1.0)):
            m1, e1 = hist[(sh, 1, "u")][key]
            m0, e0 = hist[(sh, 0, "u")][key]
            ax[c_].errorbar(U_CENT, (m1 - m0) / scale, yerr=np.hypot(e1, e0) / scale, fmt=mk, ms=5,
                            color=col, ls="none", label=f"showers {sh}")
    ax[0].set_ylabel(r"$\Delta\, dC_{QQ}/d\ln(1+\cos\chi)$", fontsize=fs)
    ax[1].set_ylabel(r"$\Delta(C_{QQ}/N_{\rm pairs})$", fontsize=fs)
    ax[1].legend(fontsize=fs, frameon=False, loc="upper center")
    for a in ax:
        a.axhline(0, color="k", lw=0.8)
        a.set_xscale("log")
        a.set_xlim(1e-4, U_EDGES[-1])
        a.set_xlabel(r"$(1+\cos\chi)/2$", fontsize=fs)
        style(a, fs, minor_x=False)
    fig.tight_layout()
    fig.savefig(OUT / "qqc_chi_src_change.pdf")
    fig.savefig(OUT / "qqc_chi_src_change.png", dpi=110)

    # 3) Same events and pairs: SRC change per pair vs dphi, chi (Breit), chi (gamma*p c.m.).
    fs = 18
    fig, ax = plt.subplots(1, 3, figsize=(20, 6.2), sharey=True)
    for sh, col, mk in (("off", "tab:red", "o"), ("on", "tab:blue", "s")):
        for c_, key in enumerate(("dphi", "chi", "chicm")):
            m1, e1 = hist[(sh, 1, key)]["R"]
            m0, e0 = hist[(sh, 0, key)]["R"]
            ok = np.ones(len(A_CENT), bool)
            for src in (0, 1):
                nb = hist[(sh, src, key)]["N"][0]
                ok &= nb >= MIN_PAIR_FRAC * nb.sum()
            m1, e1, m0, e0 = (np.where(ok, v, np.nan) for v in (m1, e1, m0, e0))
            ax[c_].errorbar(A_CENT, m1 - m0, yerr=np.hypot(e1, e0), fmt=mk, ms=5, color=col, ls="none",
                            label=f"showers {sh}")
    for a, xl in zip(ax, (r"$\Delta\phi$", r"$\chi$ (Breit)", r"$\chi$ ($\gamma^*p$ c.m.)")):
        a.axhline(0, color="k", lw=0.8)
        a.set_xlim(0, np.pi)
        a.set_xticks(*PI_TICKS)
        a.set_xlabel(xl, fontsize=fs)
        style(a, fs)
    ax[0].set_ylabel(r"$\Delta(C_{QQ}/N_{\rm pairs})$", fontsize=fs)
    ax[2].legend(fontsize=fs, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "qqc_chi_vs_dphi.pdf")
    fig.savefig(OUT / "qqc_chi_vs_dphi.png", dpi=110)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
