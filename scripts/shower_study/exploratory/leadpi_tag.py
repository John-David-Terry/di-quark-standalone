#!/usr/bin/env python3
"""
u-jet tag "leading pi+ with z > cut", tested on separate struck-u and struck-d samples
(production_leadpi/: fig3_patched struck=any struckid=2|1, SRC off, showers off/on).

Two definitions of the leading pi+:
  * SIDIS: the charged hadron with the largest z_h = P.P_h / P.q in the Breit current hemisphere
    is a pi+ with z_h > cut (no jet needed);
  * jet: the highest-energy constituent of the current anti-kT jet is a pi+ carrying more than
    `cut` of the jet energy.
Efficiency = fraction of all events tagged; purity uses the Figure 3 struck mix u:d = 0.742:0.116.
Errors are binomial.
"""

from __future__ import annotations

import os
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STUDY = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study"))
DIR = STUDY / "production_leadpi"
F_U, F_D = 0.742, 0.116
COLS = ["ok", "xB", "lead_id", "lead_z", "leadcur_id", "leadcur_z"]


def rate(mask: np.ndarray) -> tuple:
    p = float(mask.mean())
    return p, math.sqrt(p * (1 - p) / len(mask))


def main(cut: float = 0.5) -> None:
    lines = []
    say = lambda s="": (print(s), lines.append(s))  # noqa: E731
    say(f"u-tag = leading pi+ with z > {cut}; purity at u:d = {F_U}:{F_D} (untagged u purity {F_U / (F_U + F_D):.3f})")
    for sh in ("off", "on"):
        d = {f: pd.read_csv(DIR / f"showers_{sh}_struck{f}.csv", usecols=COLS) for f in ("u", "d")}
        for f in d:
            d[f]["ok"] = d[f]["ok"].astype(bool)
        say(f"\n== showers {sh} ==")
        for subset in ("all events", "Fig. 3 subset (accepted target pion)"):
            s = {f: (d[f] if subset == "all events" else d[f][d[f]["ok"]]) for f in d}
            say(f"  [{subset}: {len(s['u']):,} u, {len(s['d']):,} d events]")
            for name, tag in (
                ("SIDIS z_h, current hemisphere", lambda t: (t["leadcur_id"] == 211) & (t["leadcur_z"] > cut)),
                ("jet energy fraction", lambda t: (t["lead_id"] == 211) & (t["lead_z"] > cut)),
            ):
                (eu, seu), (ed, sed) = rate(tag(s["u"]).to_numpy()), rate(tag(s["d"]).to_numpy())
                pur = F_U * eu / (F_U * eu + F_D * ed)
                say(f"    {name:30s} tags {eu:.4f}±{seu:.4f} of struck u, {ed:.4f}±{sed:.4f} of struck d   "
                    f"ratio {eu / ed:.2f}±{eu / ed * math.hypot(seu / eu, sed / ed):.2f}   u purity {pur:.3f}")
            lc = {f: s[f]["leadcur_id"] for f in s}
            say("    leading current-hemisphere charged hadron (no z cut), fraction of events:  "
                + "  ".join(f"{lab} u {np.mean(lc['u'] == pid):.3f} / d {np.mean(lc['d'] == pid):.3f}"
                            for lab, pid in (("pi+", 211), ("pi-", -211), ("K+", 321), ("p", 2212))))
        # x dependence of the SIDIS tag (all events, showers as above)
        say("    SIDIS tag vs x (all events):  " + "  ".join(
            f"x {lo:.2f}-{hi:.2f}: u {((d['u']['leadcur_id'] == 211) & (d['u']['leadcur_z'] > cut))[(d['u']['xB'] > lo) & (d['u']['xB'] <= hi)].mean():.3f}"
            f" d {((d['d']['leadcur_id'] == 211) & (d['d']['leadcur_z'] > cut))[(d['d']['xB'] > lo) & (d['d']['xB'] <= hi)].mean():.3f}"
            for lo, hi in ((0.05, 0.1), (0.1, 0.2), (0.2, 0.495))))
    out = STUDY / "results_tagging" / f"leadpi_tag_z{str(cut).replace('.', 'p')}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")


def scan() -> None:
    """Purity and tagged fraction vs the z_h cut (SIDIS definition, all events)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    cuts = np.round(np.arange(0.0, 0.901, 0.025), 3)
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    rows = []
    for sh, ls in (("off", "--"), ("on", "-")):
        d = {f: pd.read_csv(DIR / f"showers_{sh}_struck{f}.csv", usecols=["leadcur_id", "leadcur_z"]) for f in ("u", "d")}
        is_pip = {f: (d[f]["leadcur_id"] == 211).to_numpy() for f in d}
        z = {f: d[f]["leadcur_z"].to_numpy() for f in d}
        eff = {f: np.array([rate(is_pip[f] & (z[f] > c)) for c in cuts]) for f in d}
        eu, su, ed, sd = eff["u"][:, 0], eff["u"][:, 1], eff["d"][:, 0], eff["d"][:, 1]
        pur = F_U * eu / (F_U * eu + F_D * ed)
        # d(purity) = purity (1 - purity) * sqrt((s_u/e_u)^2 + (s_d/e_d)^2)
        perr = pur * (1 - pur) * np.sqrt((su / eu) ** 2 + (sd / ed) ** 2)
        ax[0].plot(cuts, pur, ls=ls, lw=2, color="k", label=f"showers {sh}")
        ax[0].fill_between(cuts, pur - perr, pur + perr, color="k", alpha=0.2, lw=0)
        for f, col in (("u", "tab:red"), ("d", "tab:blue")):
            e, s = eff[f][:, 0] * 100, eff[f][:, 1] * 100
            ax[1].plot(cuts, e, ls=ls, lw=2, color=col, label=rf"struck ${f}$, showers {sh}")
            ax[1].fill_between(cuts, e - s, e + s, color=col, alpha=0.2, lw=0)
        rows += [(sh, c, a, b, p) for c, a, b, p in zip(cuts, eu, ed, pur)]
    ax[0].axhline(F_U / (F_U + F_D), color="0.5", ls=":", lw=1.5, label="no tag")
    ax[0].set_ylabel(r"$u$ purity of leading-$\pi^+$ tag", fontsize=fs)
    ax[1].set_ylabel(r"events tagged [%]", fontsize=fs)
    ax[0].legend(fontsize=fs * 0.8, frameon=False, loc="lower right")
    ax[1].legend(fontsize=fs * 0.8, frameon=False, loc="upper right")
    for a in ax:
        a.set_xlabel(r"cut on $z_h$ of the leading $\pi^+$", fontsize=fs)
        a.set_xlim(0, 0.9)
        a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        a.tick_params(which="both", direction="in", labelsize=fs)
        a.tick_params(which="major", length=10, width=1.2)
        a.tick_params(which="minor", length=5, width=1.0)
    ax[1].set_ylim(bottom=0)
    fig.tight_layout()
    out = STUDY / "results_tagging"
    fig.savefig(out / "leadpi_tag_zscan.pdf")
    fig.savefig(out / "leadpi_tag_zscan.png", dpi=150)
    pd.DataFrame(rows, columns=["showers", "z_cut", "eff_u", "eff_d", "u_purity"]).to_csv(
        out / "leadpi_tag_zscan.csv", index=False)
    print(f"wrote {out / 'leadpi_tag_zscan.png'}")


def scan_positive() -> None:
    """Leading positive charged hadron (pi+, K+, p) vs leading pi+ vs leading pi+/K+, showers on."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    taggers = (
        (r"leading $\pi^+$", (211,), "-"),
        (r"leading $\pi^+$ or $K^+$", (211, 321), ":"),
        (r"leading positive ($\pi^+$, $K^+$, $p$)", (211, 321, 2212), "--"),
    )
    cuts = np.round(np.arange(0.0, 0.901, 0.025), 3)
    fs = 20
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    rows = []
    for sh in ("off", "on"):
        d = {f: pd.read_csv(DIR / f"showers_{sh}_struck{f}.csv", usecols=["leadcur_id", "leadcur_z"]) for f in ("u", "d")}
        for lab, ids, ls in taggers:
            sel = {f: d[f]["leadcur_id"].isin(ids).to_numpy() for f in d}
            z = {f: d[f]["leadcur_z"].to_numpy() for f in d}
            eff = {f: np.array([rate(sel[f] & (z[f] > c)) for c in cuts]) for f in d}
            eu, su, ed, sd = eff["u"][:, 0], eff["u"][:, 1], eff["d"][:, 0], eff["d"][:, 1]
            pur = F_U * eu / (F_U * eu + F_D * ed)
            rows += [(sh, lab.replace("$", "").replace("\\", ""), c, a, b, p) for c, a, b, p in zip(cuts, eu, ed, pur)]
            if sh != "on":
                continue
            perr = pur * (1 - pur) * np.sqrt((su / eu) ** 2 + (sd / ed) ** 2)
            ax[0].plot(cuts, pur, ls=ls, lw=2, color="k", label=lab)
            ax[0].fill_between(cuts, pur - perr, pur + perr, color="k", alpha=0.15, lw=0)
            for f, col in (("u", "tab:red"), ("d", "tab:blue")):
                ax[1].plot(cuts, eff[f][:, 0] * 100, ls=ls, lw=2, color=col)
    ax[0].axhline(F_U / (F_U + F_D), color="0.5", ls=(0, (1, 3)), lw=1.5, label="no tag")
    ax[1].plot([], [], color="tab:red", lw=2, label=r"struck $u$")
    ax[1].plot([], [], color="tab:blue", lw=2, label=r"struck $d$")
    ax[0].set_ylabel(r"$u$ purity of the tag", fontsize=fs)
    ax[1].set_ylabel(r"events tagged [%]", fontsize=fs)
    ax[0].legend(fontsize=fs * 0.75, frameon=False, loc="lower right")
    ax[1].legend(fontsize=fs * 0.8, frameon=False, loc="upper right")
    for a in ax:
        a.set_xlabel(r"cut on $z_h$ of the leading charged hadron", fontsize=fs)
        a.set_xlim(0, 0.9)
        a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        a.tick_params(which="both", direction="in", labelsize=fs)
        a.tick_params(which="major", length=10, width=1.2)
        a.tick_params(which="minor", length=5, width=1.0)
    ax[1].set_ylim(bottom=0)
    fig.tight_layout()
    out = STUDY / "results_tagging"
    fig.savefig(out / "leadpos_tag_zscan.pdf")
    fig.savefig(out / "leadpos_tag_zscan.png", dpi=150)
    pd.DataFrame(rows, columns=["showers", "tagger", "z_cut", "eff_u", "eff_d", "u_purity"]).to_csv(
        out / "leadpos_tag_zscan.csv", index=False)
    print(f"wrote {out / 'leadpos_tag_zscan.png'}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "scan":
        scan()
    elif len(sys.argv) > 1 and sys.argv[1] == "scan_positive":
        scan_positive()
    else:
        main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.5)
