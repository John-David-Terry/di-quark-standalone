#!/usr/bin/env python3
"""
Absolute phi_hJ distributions vs the ratio R_{pi-/pi+}(phi_hJ), from existing samples (anti-kT jet).

(a) normalized dsigma/dphi_hJ of the leading pi- without SRC, for several radiation and
    non-perturbative settings: the absolute shape moves a lot.
(b) R(phi_hJ) without SRC (0% mixture) for the same settings: stays flat at 1.
(c) SRC 3.74 GeV kick, default kT: change of the absolute pi- (and pi+) shape from showers
    vs from the 10% SRC mixture, each relative to the showers-off 0% curve.
(d) R(phi_hJ) at 0% and 10% for the same sample, showers off vs on: the SRC modulation survives.

Mixtures and ratios use the original Figure 3 code (via compare_fig3).

Default: a two-panel version with the 10% SRC mixture only, showers off vs on — (left) the
pi- distribution 1/sigma dsigma/dphi_hJ (unit area on [0, pi], as in Figure 3), (right) R(phi_hJ).
``--tag <production_primkt tag>`` picks another sample pair for it (default kick3p74);
``--full`` makes the four-panel figure above.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.ticker
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fig3 import (  # noqa: E402
    BINS, CENTERS, SEED, STUDY, cos_fit, plt, r_curves, signal, subsample_altered_unchanged_for_pool,
)

COLS = ["event_id", "sample", "pion", "ok", "phi_hJ"]
PK = STUDY / "production_primkt"
BASELINES = [  # (label, unchanged-arm file), all SRC off, anti-kT jet
    ("default $k_T$, showers off", PK / "showers_off_unchanged_defaultkt.csv"),
    ("default $k_T$, showers on", PK / "showers_on_unchanged_kick3p74.csv"),
    ("default $k_T$, mJoinJunction 0, showers off", PK / "showers_off_unchanged_kick3p74_mjoin0p0.csv"),
    ("default $k_T$, mJoinJunction 0, showers on", PK / "showers_on_unchanged_kick3p74_mjoin0p0.csv"),
    (r"primordial $k_T$ $\sigma=0.13$ GeV, showers off", PK / "showers_off_unchanged_primkt0p13.csv"),
    (r"$\sigma=0.13$ GeV, mJoinJunction 0, showers off",
     STUDY / "production_mjoin/showers_off_unchanged_primkt0p13_mjoin0p0.csv"),
]
def signal_files(tag: str) -> dict:
    """SRC-off / SRC-on sample pair per shower setting for a production_primkt tag."""
    return {sh: (PK / f"showers_{sh}_unchanged_{tag}.csv", PK / f"showers_{sh}_altered_{tag}.csv")
            for sh in ("off", "on")}


SIGNAL = signal_files("kick3p74")  # SRC kick 3.74 GeV, default kT, default hadronization
BINW = BINS[1] - BINS[0]


def load(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, usecols=COLS)
    d = d[d["ok"].astype(bool)].copy()
    d["ok"] = True
    return d


def density(v: np.ndarray) -> tuple:
    v = v[np.isfinite(v)]
    h, _ = np.histogram(v, bins=BINS)
    return h / (len(v) * BINW), np.sqrt(h) / (len(v) * BINW)


def pim(df: pd.DataFrame, f: float, ix: int) -> np.ndarray:
    """pi- phi_hJ values of the Figure 3 mixture with altered fraction f (f = 0: unchanged only)."""
    if f == 0.0:
        mix = df[df["sample"] == "unchanged"]
    else:
        mix = subsample_altered_unchanged_for_pool(df, f, seed=SEED + ix)[0]
    return mix.loc[mix["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float)


def two_panel(out: Path, tag: str = "kick3p74") -> None:
    """10% SRC, showers off vs on: absolute pi- shape and R (default: kick 3.74 GeV, default kT)."""
    sig = {sh: pd.concat([load(u), load(a)]) for sh, (u, a) in signal_files(tag).items()}
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    for sh, col, mk in (("off", "tab:red", "o"), ("on", "tab:blue", "s")):
        h, e = density(pim(sig[sh], 0.10, 2))
        ax[0].errorbar(CENTERS, h, yerr=e, fmt=mk, ms=3, color=col, ls="none", label=f"showers {sh}")
        c = r_curves(sig[sh], "phi_hJ")[0.10]
        ax[1].errorbar(CENTERS, c["R"], yerr=c["err"], fmt=mk, ms=3, color=col, ls="none", label=f"showers {sh}")
    ax[0].set_ylabel(r"$\dfrac{1}{\sigma}\dfrac{d\sigma}{d\phi_{hJ}}$", fontsize=fs)
    ax[1].set_ylabel(r"$R_{\pi^-/\pi^+}$ (10% SRC)", fontsize=fs)
    ax[1].legend(fontsize=fs, loc="upper left", frameon=False)
    for a in ax:
        a.set_xlabel(r"$\phi_{hJ}$ [rad]", fontsize=fs)
        a.set_xlim(0.0, np.pi)
        a.set_xticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
        a.set_xticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
        a.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(np.pi / 16))
        a.tick_params(which="both", direction="in", labelsize=fs, labelleft=True)
        a.tick_params(which="major", length=10, width=1.2)
        a.tick_params(which="minor", length=5, width=1.0)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    ax[0].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.05))
    ax[1].set_ylim(0.93, 1.09)
    ax[1].set_yticks([0.95, 1.00, 1.05])
    ax[1].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.01))
    fig.tight_layout()
    name = "absolute_vs_ratio_10pct" + ("" if tag == "kick3p74" else f"_{tag}")
    fig.savefig(out / f"{name}.pdf")
    fig.savefig(out / f"{name}.png", dpi=150)
    print(f"wrote {out / name}.png")


def main() -> None:
    out = STUDY / "results_absolute_vs_ratio"
    out.mkdir(parents=True, exist_ok=True)
    args = sys.argv[1:]
    if "--full" not in args:
        two_panel(out, args[args.index("--tag") + 1] if "--tag" in args else "kick3p74")
        return
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.2))
    colors = plt.cm.viridis(np.linspace(0.0, 0.9, len(BASELINES)))

    say("SRC off (0% mixture), anti-kT jet: absolute pi- shape vs R(0%) flatness")
    say(f"  {'configuration':52s} {'<cos phi_hJ> pi-':>18s} {'R(0%) cos modulation':>22s}")
    for (label, path), col in zip(BASELINES, colors):
        d = load(path)
        h, e = density(d.loc[d["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float))
        ax[0, 0].step(CENTERS, h, where="mid", color=col, lw=1.4, label=label)
        c0 = r_curves(d, "phi_hJ")[0.0]
        ax[0, 1].errorbar(CENTERS, c0["R"], yerr=c0["err"], fmt="o", ms=2.5, color=col, label=label)
        m, sm = cos_fit(c0["R"], c0["err"])
        v = d.loc[d["pion"] == "piminus", "phi_hJ"].to_numpy(dtype=float)
        cv = np.cos(v[np.isfinite(v)])
        say(f"  {label.replace('$', ''):52s} {cv.mean():+.4f}±{cv.std() / math.sqrt(len(cv)):.4f}   {m:+.4f}±{sm:.4f}")

    sig = {sh: pd.concat([load(u), load(a)]) for sh, (u, a) in SIGNAL.items()}
    ref, ref_e = density(pim(sig["off"], 0.0, 0))
    curves = {}
    for sh, f, ix, col, ls, lab in (
        ("on", 0.0, 0, "tab:blue", "-", r"$\pi^-$, showers on, no SRC (radiation only)"),
        ("off", 0.10, 2, "tab:red", "-", r"$\pi^-$, showers off, 10% SRC (SRC only)"),
        ("on", 0.10, 2, "tab:purple", "--", r"$\pi^-$, showers on, 10% SRC (both)"),
    ):
        h, e = density(pim(sig[sh], f, ix))
        r = h / ref
        re = r * np.sqrt((e / h) ** 2 + (ref_e / ref) ** 2)
        ax[1, 0].errorbar(CENTERS, r, yerr=re, fmt="o", ms=2.5, color=col, label=lab, ls="none")
        curves[(sh, f)] = r
    unc_off = sig["off"][sig["off"]["sample"] == "unchanged"]
    unc_on = sig["on"][sig["on"]["sample"] == "unchanged"]
    hp_off, ep_off = density(unc_off.loc[unc_off["pion"] == "piplus", "phi_hJ"].to_numpy(dtype=float))
    hp_on, ep_on = density(unc_on.loc[unc_on["pion"] == "piplus", "phi_hJ"].to_numpy(dtype=float))
    rp = hp_on / hp_off
    ax[1, 0].errorbar(CENTERS, rp, yerr=rp * np.sqrt((ep_on / hp_on) ** 2 + (ep_off / hp_off) ** 2),
                      fmt="s", ms=2.5, color="tab:cyan", label=r"$\pi^+$, showers on, no SRC (radiation only)", ls="none")

    say("\nAbsolute pi- shape, relative to showers off with no SRC (kick 3.74 GeV, default kT):")
    sel = CENTERS > 0.75 * math.pi
    for key, lab in (((("on", 0.0)), "radiation only (showers on, 0%)"),
                     ((("off", 0.10)), "SRC only (showers off, 10%)"),
                     ((("on", 0.10)), "both (showers on, 10%)")):
        r = curves[key]
        say(f"  {lab:34s} mean change at phi > 3pi/4: {np.mean(r[sel]) - 1:+.4f}   at phi < pi/4: {np.mean(r[CENTERS < 0.25 * math.pi]) - 1:+.4f}")

    say("\nR(phi_hJ), kick 3.74 GeV, default kT:")
    for (sh, col), mk in zip((("off", "tab:red"), ("on", "tab:purple")), ("o", "s")):
        c = r_curves(sig[sh], "phi_hJ")
        for f, alpha, lab in ((0.0, 0.45, "0%"), (0.10, 1.0, "10%")):
            ax[1, 1].errorbar(CENTERS, c[f]["R"], yerr=c[f]["err"], fmt=mk, ms=2.5, color=col, alpha=alpha,
                              ls="none", label=f"showers {sh}, {lab} SRC")
        s = signal(c, 0.10)
        m0, s0 = cos_fit(c[0.0]["R"], c[0.0]["err"])
        say(f"  showers {sh:3s}: R(0%) modulation {m0:+.4f}±{s0:.4f}   SRC modulation at 10% {s['mod']:+.4f}±{s['mod_err']:.4f}")

    titles = [r"(a) absolute: $\frac{1}{\sigma}\frac{d\sigma}{d\phi_{hJ}}$, leading $\pi^-$, no SRC",
              r"(b) ratio $R_{\pi^-/\pi^+}$, no SRC (0%)",
              r"(c) absolute $\pi^\pm$ shape / (showers off, no SRC); SRC kick 3.74 GeV",
              r"(d) ratio $R_{\pi^-/\pi^+}$, SRC kick 3.74 GeV: showers off vs on"]
    for a, t in zip(ax.flat, titles):
        a.set_title(t, fontsize=9)
        a.set_xlabel(r"$\phi_{hJ}$ [rad]")
        a.legend(fontsize=6.5, loc="upper left")
    ax[0, 1].axhline(1.0, color="k", lw=0.6)
    ax[0, 1].set_ylim(0.9, 1.1)
    ax[1, 0].axhline(1.0, color="k", lw=0.6)
    ax[1, 1].axhline(1.0, color="k", lw=0.6)
    ax[0, 0].set_ylabel(r"$\frac{1}{\sigma}\frac{d\sigma}{d\phi_{hJ}}$")
    ax[0, 1].set_ylabel(r"$R_{\pi^-/\pi^+}$")
    ax[1, 0].set_ylabel("ratio to showers off, no SRC")
    ax[1, 1].set_ylabel(r"$R_{\pi^-/\pi^+}$")
    fig.suptitle(r"Patched Pythia 8.311, Fig. 3 selection ($5\times41$ GeV, $4<Q<5$ GeV), anti-$k_T$ jet", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "absolute_vs_ratio.pdf")
    fig.savefig(out / "absolute_vs_ratio.png", dpi=140)
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
