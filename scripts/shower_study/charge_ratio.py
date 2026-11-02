#!/usr/bin/env python3
"""
SRC signal in the experimentally measurable charge ratio, and its rapidity-separation dependence.

  R_pi(phi_hJ) = [(1/N_pi-) dN_pi-/dphi_hJ] / [(1/N_pi+) dN_pi+/dphi_hJ]

with pi- and pi+ taken from the same sample (what an experiment measures), the manuscript's
leading-target-pion selection and kinematic cuts, the Breit anti-kT current jet, 40 bins on
[0, pi]; the ratio and its Poisson errors come from the manuscript code (_ratio_series).
Configurations: no SRC / SRC (the pure SRC-on sample) x no radiation / radiation.

Test 1: R_pi for the four configurations, and R_pi^SRC / R_pi^noSRC without and with radiation.
Test 2: the same in three bins of the Breit-frame separation between the pion and the jet:
        pseudorapidity Delta eta_hJ = |eta_pi - eta_J| (default, bins < 4, 4-5.5, >= 5.5) or
        rapidity |y_pi - y_J| (--sep dy). eta follows exactly from the stored Breit rapidity,
        p_T and energy: sinh(eta) = (m_T / p_T) sinh(y), m_T = E / cosh(y) (jet) or
        sqrt(m_pi^2 + p_T^2) (pion). The lab-frame acceptance is unchanged; the separation only
        sorts accepted events. --edges takes all bin boundaries (events outside them are not used);
        --edges tertiles uses the no-SRC, no-radiation tertiles instead.
The per-bin figure has three rows: pi- 1/sigma dsigma/dphi_hJ (unit area in each bin) for the four
configurations, R_pi, and R_pi^SRC / R_pi^noSRC.
Summary per double ratio: back-to-back excess (mean - 1 over phi_hJ > 3 pi/4), cos-modulation m of
a (1 + m cos phi_hJ), and their radiation survival (with / without radiation).
Input: production_primkt/showers_{off,on}_{unchanged,altered}_<tag>.csv (run_fig3_remnant_width.sh).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fig3 import BINS, CENTERS, STUDY, _ratio_series, cos_fit, plt  # noqa: E402

import matplotlib.ticker  # noqa: E402

PK = STUDY / "production_primkt"
COLS = ["pion", "ok", "phi_hJ", "jet_y_breit", "pion_y_breit", "jet_pjt_GeV", "jet_E_breit", "pion_pht_GeV"]
M_PI = 0.13957
STYLE = {"off": ("tab:red", "o", "no radiation"), "on": ("tab:blue", "s", "with radiation")}
PI_TICKS = ([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi], [r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])


def load(sh: str, src: int, tag: str) -> pd.DataFrame:
    arm = "altered" if src else "unchanged"
    d = pd.read_csv(PK / f"showers_{sh}_{arm}_{tag}.csv", usecols=COLS, engine="pyarrow")
    d = d[d["ok"].astype(bool)]
    yp, yj = d["pion_y_breit"].to_numpy(), d["jet_y_breit"].to_numpy()
    ptp, ptj = d["pion_pht_GeV"].to_numpy(), d["jet_pjt_GeV"].to_numpy()
    etap = np.arcsinh(np.sqrt(M_PI ** 2 + ptp ** 2) / ptp * np.sinh(yp))
    etaj = np.arcsinh(d["jet_E_breit"].to_numpy() / np.cosh(yj) / ptj * np.sinh(yj))
    return pd.DataFrame({"pion": d["pion"].to_numpy(), "phi": d["phi_hJ"].to_numpy(dtype=float),
                         "dy": np.abs(yp - yj), "deta": np.abs(etap - etaj), "eta_pi": etap, "eta_j": etaj})


def r_pi(d: pd.DataFrame) -> tuple:
    vp = d.loc[d["pion"] == "piplus", "phi"].to_numpy()
    vm = d.loc[d["pion"] == "piminus", "phi"].to_numpy()
    _, r, e, _ = _ratio_series(vp[np.isfinite(vp)], vm[np.isfinite(vm)], BINS)
    return r, e, len(vp), len(vm)


def pim_density(d: pd.DataFrame) -> tuple:
    """pi- 1/sigma dsigma/dphi_hJ (unit area on [0, pi]) with Poisson errors."""
    v = d.loc[d["pion"] == "piminus", "phi"].to_numpy()
    v = v[np.isfinite(v)]
    h, _ = np.histogram(v, bins=BINS)
    norm = len(v) * np.diff(BINS)
    return h / norm, np.sqrt(h) / norm


def b2b(r: np.ndarray, e: np.ndarray) -> tuple:
    sel = (CENTERS > 0.75 * math.pi) & np.isfinite(r)
    w = 1.0 / e[sel] ** 2
    return float(np.sum(w * (r[sel] - 1.0)) / np.sum(w)), float(1.0 / math.sqrt(np.sum(w)))


def ratio_of(a: tuple, b: tuple) -> tuple:
    r = a[0] / b[0]
    return r, r * np.sqrt((a[1] / a[0]) ** 2 + (b[1] / b[0]) ** 2)


def surv(on: tuple, off: tuple) -> tuple:
    s = on[0] / off[0]
    return s, abs(s) * math.hypot(on[1] / on[0], off[1] / off[0])


def style_axis(a, fs, xlabel=True):
    a.set_xlim(0.0, np.pi)
    a.set_xticks(*PI_TICKS)
    a.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(np.pi / 16))
    a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    if xlabel:
        a.set_xlabel(r"$\phi_{hJ}$ [rad]", fontsize=fs)
    a.tick_params(which="both", direction="in", labelsize=fs)
    a.tick_params(which="major", length=10, width=1.2)
    a.tick_params(which="minor", length=5, width=1.0)


def draw_abs(ax, res, fs, legend=True):
    """pi- 1/sigma dsigma/dphi_hJ for the four configurations."""
    for src in (0, 1):
        for sh, (col, mk, rad) in STYLE.items():
            h, e = res[("dens", sh, src)]
            ax.errorbar(CENTERS, h, yerr=e, fmt=mk, ms=4, color=col, mfc=col if src else "white", ls="none",
                        label=f"{'SRC' if src else 'no SRC'}, {rad}")
    if legend:
        ax.legend(fontsize=fs * 0.75, loc="center left", frameon=False)


def draw(ax_r, ax_d, res, fs, legend=True, legend_r=True):
    """Four R_pi curves on ax_r, the two SRC / no-SRC double ratios on ax_d."""
    for src in (0, 1):
        for sh, (col, mk, rad) in STYLE.items():
            r, e = res[(sh, src)][:2]
            ax_r.errorbar(CENTERS, r, yerr=e, fmt=mk, ms=4, color=col, mfc=col if src else "white", ls="none",
                          label=f"{'SRC' if src else 'no SRC'}, {rad}")
    for sh, (col, mk, rad) in STYLE.items():
        r, e = res[("dr", sh)]
        ax_d.errorbar(CENTERS, r, yerr=e, fmt=mk, ms=4, color=col, ls="none", label=rad)
    for a in (ax_r, ax_d):
        a.axhline(1.0, color="k", lw=0.8)
    if legend and legend_r:
        ax_r.legend(fontsize=fs * 0.75, loc="lower right", frameon=False)
    if legend:
        ax_d.legend(fontsize=fs * 0.75, loc="upper left", frameon=False)


def analyse(data: dict, sel: dict | None, say, label: str) -> dict:
    res = {}
    for sh in ("off", "on"):
        for src in (0, 1):
            d = data[(sh, src)] if sel is None else data[(sh, src)][sel[(sh, src)]]
            res[(sh, src)] = r_pi(d)
            res[("dens", sh, src)] = pim_density(d)
        res[("dr", sh)] = ratio_of(res[(sh, 1)], res[(sh, 0)])
    say(f"\n== {label} ==")
    for sh in ("off", "on"):
        for src in (0, 1):
            r, e, npl, nmi = res[(sh, src)]
            say(f"  R_pi {('SRC' if src else 'no SRC'):6s} {STYLE[sh][2]:15s} n(pi+) {npl:>10,} n(pi-) {nmi:>10,}  "
                f"b2b R-1 {b2b(r, e)[0]:+.4f}±{b2b(r, e)[1]:.4f}  cos mod {cos_fit(r, e)[0]:+.4f}±{cos_fit(r, e)[1]:.4f}")
    for sh in ("off", "on"):
        r, e = res[("dr", sh)]
        res[("b2b", sh)], res[("mod", sh)] = b2b(r, e), cos_fit(r, e)
        say(f"  R_SRC/R_noSRC {STYLE[sh][2]:15s} b2b excess {res[('b2b', sh)][0]:+.4f}±{res[('b2b', sh)][1]:.4f}  "
            f"cos mod {res[('mod', sh)][0]:+.4f}±{res[('mod', sh)][1]:.4f}")
    for key, lab in (("b2b", "b2b excess"), ("mod", "cos modulation")):
        s = surv(res[(key, "on")], res[(key, "off")])
        res[("surv", key)] = s
        say(f"  survival (with / without radiation), {lab:15s} {s[0]:.3f} ± {s[1]:.3f}")
    return res


def bin_name(lo: float, hi: float, var: str) -> str:
    if lo <= 0.0:
        return f"{var} < {hi:g}"
    if not np.isfinite(hi):
        return f"{var} >= {lo:g}"
    return f"{lo:g} <= {var} < {hi:g}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="ktw0p693_m0", help="production_primkt sample tag")
    ap.add_argument("--sep", choices=("deta", "dy"), default="deta",
                    help="Breit-frame separation: pseudorapidity (deta) or rapidity (dy)")
    ap.add_argument("--edges", default="0,4,5.5,inf", help="bin boundaries (e.g. 1,2.5,4,5.5), or 'tertiles'")
    args = ap.parse_args()
    out = STUDY / f"results_charge_ratio_{args.tag}"
    out.mkdir(parents=True, exist_ok=True)
    sym = {"deta": (r"\Delta\eta_{hJ}", "Delta eta_hJ"), "dy": (r"|y_\pi - y_J|", "|y_pi - y_J|")}[args.sep]
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say(f"Charge ratio R_pi = shape(pi-)/shape(pi+), anti-kT jet, tag {args.tag}; SRC = pure SRC-on sample")
    data = {(sh, src): load(sh, src, args.tag) for sh in ("off", "on") for src in (0, 1)}

    # Test 1.
    res = analyse(data, None, say, "Test 1: all events")
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    draw(ax[0], ax[1], res, fs)
    ax[0].set_ylabel(r"$R_\pi(\phi_{hJ})$", fontsize=fs)
    ax[1].set_ylabel(r"$R_\pi^{\rm SRC}/R_\pi^{\rm noSRC}$", fontsize=fs)
    for a in ax:
        style_axis(a, fs)
    fig.tight_layout()
    fig.savefig(out / "charge_ratio.pdf")
    fig.savefig(out / "charge_ratio.png", dpi=150)

    # Test 2: bins of the Breit-frame separation.
    if args.edges == "tertiles":
        edges = [0.0] + list(np.quantile(data[("off", 0)][args.sep], [1 / 3, 2 / 3])) + [np.inf]
    else:
        edges = [float(x) for x in args.edges.split(",")]
    suffix = "" if args.edges == ap.get_default("edges") else "_" + args.edges.replace(",", "-")
    names = [bin_name(lo, hi, sym[1]) for lo, hi in zip(edges[:-1], edges[1:])]
    say(f"\nTest 2: bins in {sym[1]} (Breit frame): " + ", ".join(names))
    say("  per configuration: mean, median, 5-95% range, fraction of events (and pi+ / pi- counts) per bin")
    for (sh, src), d in data.items():
        v = d[args.sep].to_numpy()
        frac, cnt = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (v >= lo) & (v < hi)
            frac.append(m.mean())
            cnt.append(f"{int((m & (d['pion'] == 'piplus')).sum()):,}/{int((m & (d['pion'] == 'piminus')).sum()):,}")
        q5, q50, q95 = np.quantile(v, [0.05, 0.5, 0.95])
        say(f"  {'SRC' if src else 'no SRC':6s} {STYLE[sh][2]:15s} mean {v.mean():.2f} median {q50:.2f} "
            f"[{q5:.2f}, {q95:.2f}]  fractions " + " ".join(f"{f:.3f}" for f in frac)
            + f" (outside {1 - sum(frac):.3f})  pi+/pi- " + "  ".join(cnt))
    m = data[("off", 0)]
    say(f"  no SRC, no radiation: median eta_pi {np.median(m['eta_pi']):+.2f}, median eta_J {np.median(m['eta_j']):+.2f} (Breit)")

    # Distribution of the separation, with the bin edges.
    fs = 20
    fig, a = plt.subplots(figsize=(8.5, 5.6))
    hb = np.linspace(0.0, 8.0, 81)
    for src in (0, 1):
        for sh, (col, mk, rad) in STYLE.items():
            v = data[(sh, src)][args.sep].to_numpy()
            h, _ = np.histogram(v, bins=hb)
            a.step(0.5 * (hb[1:] + hb[:-1]), h / (len(v) * np.diff(hb)), where="mid", color=col,
                   ls="-" if src else "--", lw=1.6, label=f"{'SRC' if src else 'no SRC'}, {rad}")
    for x in edges:
        if 0.0 < x < np.inf:
            a.axvline(x, color="k", lw=0.8, ls=":")
    a.set_xlabel(rf"${sym[0]}$ (Breit frame)", fontsize=fs)
    a.set_ylabel("normalized distribution", fontsize=fs)
    a.set_xlim(hb[0], hb[-1])
    a.legend(fontsize=fs * 0.7, frameon=False, loc="upper right")
    a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.tick_params(which="both", direction="in", labelsize=fs)
    a.tick_params(which="major", length=10, width=1.2)
    a.tick_params(which="minor", length=5, width=1.0)
    fig.tight_layout()
    fig.savefig(out / f"{args.sep}_distribution{suffix}.pdf")
    fig.savefig(out / f"{args.sep}_distribution{suffix}.png", dpi=150)

    bins = []
    for (lo, hi), name in zip(zip(edges[:-1], edges[1:]), names):
        sel = {k: (d[args.sep] >= lo) & (d[args.sep] < hi) for k, d in data.items()}
        bins.append((lo, hi, name, analyse(data, sel, say, f"Test 2: {name}")))
    say(f"\n== Test 2 summary: SRC enhancement near phi = pi (b2b excess of R_SRC/R_noSRC), bins of {sym[1]} ==")
    for lo, hi, name, r in bins:
        say(f"  {name:22s} no rad {r[('b2b', 'off')][0]:+.4f}±{r[('b2b', 'off')][1]:.4f}  "
            f"rad {r[('b2b', 'on')][0]:+.4f}±{r[('b2b', 'on')][1]:.4f}  "
            f"survival {r[('surv', 'b2b')][0]:.3f}±{r[('surv', 'b2b')][1]:.3f}  "
            f"(cos-mod survival {r[('surv', 'mod')][0]:.3f}±{r[('surv', 'mod')][1]:.3f})")
    (out / f"summary_{args.sep}{suffix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # No gaps between panels: y tick labels only in the left column, x tick labels only in the bottom
    # row, bin ranges above the top row, and the outermost ticks of each row pruned so labels of
    # neighbouring rows do not collide.
    fs = 32
    fig, ax = plt.subplots(3, len(bins), figsize=(8 * len(bins), 21), sharex=True, sharey="row")
    for j, (lo, hi, name, r) in enumerate(bins):
        draw_abs(ax[0, j], r, fs, legend=(j == 0))
        draw(ax[1, j], ax[2, j], r, fs, legend=(j == 0), legend_r=False)
        txt = "$" + bin_name(lo, hi, sym[0]).replace("<=", r"\leq").replace(">=", r"\geq") + "$"
        ax[0, j].set_title(txt, fontsize=fs, pad=16)
        for i in (0, 1, 2):
            style_axis(ax[i, j], fs, xlabel=(i == 2))
            ax[i, j].tick_params(which="both", top=True, right=True)
    ax[0, 0].set_ylim(0.0, None)
    for i in (0, 1, 2):
        ax[i, 0].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, prune="both"))
    ax[0, 0].set_ylabel(r"$\dfrac{1}{\sigma}\dfrac{d\sigma}{d\phi_{hJ}}$ ($\pi^-$)", fontsize=fs)
    ax[1, 0].set_ylabel(r"$R_\pi(\phi_{hJ})$", fontsize=fs)
    ax[2, 0].set_ylabel(r"$R_\pi^{\rm SRC}/R_\pi^{\rm noSRC}$", fontsize=fs)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.07, top=0.95, wspace=0.0, hspace=0.0)
    fig.savefig(out / f"charge_ratio_{args.sep}_bins{suffix}.pdf")
    fig.savefig(out / f"charge_ratio_{args.sep}_bins{suffix}.png", dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
