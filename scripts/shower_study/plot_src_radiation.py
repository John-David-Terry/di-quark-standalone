#!/usr/bin/env python3
"""
Isolate the SRC effect on phi_hJ and how much of it survives parton showers (anti-kT jet).

Left:  1/sigma dsigma/dphi_hJ of the leading target pi- (unit area on [0, pi], Figure 3 binning)
       for no SRC / SRC x no radiation / radiation.
Right: SRC / no SRC ratio of those normalized distributions, without and with radiation.

"SRC" is the SRC-on (altered) sample by default (--fraction 1); --fraction f uses the Figure 3
mixture with altered fraction f instead. Ratios are fitted to a (1 + m cos phi_hJ); the summary
gives m, the back-to-back excess (phi_hJ > 3 pi/4) and their radiation survival (with / without).
Errors are Poisson per bin; for f = 1 the two samples are independent, for f < 1 the mixture
shares the SRC-off events with the reference, so the quoted errors are conservative.
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
from compare_fig3 import CENTERS, STUDY, cos_fit, plt  # noqa: E402
from plot_absolute_vs_ratio import PK, density, load, pim  # noqa: E402

import matplotlib.ticker  # noqa: E402


def b2b(r: np.ndarray, e: np.ndarray) -> tuple:
    sel = CENTERS > 0.75 * math.pi
    w = 1.0 / e[sel] ** 2
    return float(np.sum(w * (r[sel] - 1.0)) / np.sum(w)), float(1.0 / math.sqrt(np.sum(w)))


def survival(on: tuple, off: tuple) -> tuple:
    r = on[0] / off[0]
    return r, abs(r) * math.hypot(on[1] / on[0], off[1] / off[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="ktw0p693_m0", help="production_primkt sample tag")
    ap.add_argument("--fraction", type=float, default=1.0, help="SRC fraction (1 = pure SRC sample)")
    args = ap.parse_args()
    out = STUDY / f"results_src_radiation_{args.tag}"
    out.mkdir(parents=True, exist_ok=True)
    fpct = f"{int(round(100 * args.fraction))}pct"

    dens, ratio, lines = {}, {}, []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say(f"pi- phi_hJ (anti-kT jet), sample tag {args.tag}, SRC fraction {args.fraction:g}")
    for sh in ("off", "on"):
        unc = load(PK / f"showers_{sh}_unchanged_{args.tag}.csv")
        alt = load(PK / f"showers_{sh}_altered_{args.tag}.csv")
        v0 = unc.loc[unc["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float)
        if args.fraction >= 1.0:
            v1 = alt.loc[alt["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float)
        else:
            v1 = pim(pd.concat([unc, alt]), args.fraction, 2)
        dens[(sh, 0)], dens[(sh, 1)] = density(v0), density(v1)
        (h0, e0), (h1, e1) = dens[(sh, 0)], dens[(sh, 1)]
        r = h1 / h0
        re = r * np.sqrt((e1 / h1) ** 2 + (e0 / h0) ** 2)
        ratio[sh] = (r, re)
        m = cos_fit(r, re)
        bb = b2b(r, re)
        c0, c1 = np.cos(v0[np.isfinite(v0)]), np.cos(v1[np.isfinite(v1)])
        dcos = (c1.mean() - c0.mean(), math.hypot(c1.std() / math.sqrt(len(c1)), c0.std() / math.sqrt(len(c0))))
        ratio[sh] += (m, bb, dcos)
        say(f"\nradiation {sh}: n(pi-) no SRC {len(v0):,}, SRC {len(v1):,}")
        say(f"  SRC / no SRC: cos modulation m = {m[0]:+.4f} ± {m[1]:.4f}, "
            f"back-to-back excess (phi > 3pi/4) {bb[0]:+.4f} ± {bb[1]:.4f}")
        say(f"  <cos phi_hJ> shift (SRC - no SRC) {dcos[0]:+.4f} ± {dcos[1]:.4f}")
    say("\nsurvival with radiation / without:")
    for i, lab in ((2, "cos modulation"), (3, "back-to-back excess"), (4, "<cos phi_hJ> shift")):
        s, se = survival(ratio["on"][i], ratio["off"][i])
        say(f"  {lab:22s} {s:.3f} ± {se:.3f}")
    (out / f"summary_{fpct}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    src_lab = "SRC" if args.fraction >= 1.0 else f"{int(round(100 * args.fraction))}% SRC"
    style = {"off": ("tab:red", "o", "no radiation"), "on": ("tab:blue", "s", "with radiation")}
    for src in (0, 1):
        for sh, (col, mk, rad) in style.items():
            h, e = dens[(sh, src)]
            lab = f"{src_lab if src else 'no SRC'}, {rad}"
            ax[0].errorbar(CENTERS, h, yerr=e, fmt=mk, ms=5, color=col, mfc=col if src else "white",
                           ls="none", label=lab)
    for sh, (col, mk, rad) in style.items():
        r, re = ratio[sh][:2]
        ax[1].errorbar(CENTERS, r, yerr=re, fmt=mk, ms=4, color=col, ls="none", label=rad)
    ax[0].set_ylabel(r"$\dfrac{1}{\sigma}\dfrac{d\sigma}{d\phi_{hJ}}$", fontsize=fs)
    ax[1].set_ylabel(f"{src_lab} / no SRC", fontsize=fs)
    ax[0].legend(fontsize=fs * 0.75, loc="upper left", frameon=False)
    ax[1].legend(fontsize=fs * 0.75, loc="upper left", frameon=False)
    ax[0].set_ylim(0.0, 1.1)
    ax[0].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    ax[0].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.05))
    ax[1].axhline(1.0, color="k", lw=0.8)
    ax[1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    for a in ax:
        a.set_xlabel(r"$\phi_{hJ}$ [rad]", fontsize=fs)
        a.set_xlim(0.0, np.pi)
        a.set_xticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
        a.set_xticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
        a.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(np.pi / 16))
        a.tick_params(which="both", direction="in", labelsize=fs)
        a.tick_params(which="major", length=10, width=1.2)
        a.tick_params(which="minor", length=5, width=1.0)
    fig.tight_layout()
    fig.savefig(out / f"src_radiation_{fpct}.pdf")
    fig.savefig(out / f"src_radiation_{fpct}.png", dpi=150)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
