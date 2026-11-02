"""Back-to-back shape diagnostics for the opening-angle QQC: which choice suppresses z -> 0?"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

M_P = 0.938272
TAG = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study")) / "production_qqc_utag"
EDGES = np.geomspace(1e-4, 1.0, 13)
CENT = np.sqrt(EDGES[:-1] * EDGES[1:])


def load(name):
    ev = pd.read_csv(TAG / f"{name}.csv", usecols=["event_id", "xB", "Q2"])
    had = pd.read_csv(TAG / f"{name}_qqc.csv")
    idx = pd.Series(np.arange(len(ev)), index=ev["event_id"].to_numpy())
    i = idx.loc[had["event_id"].to_numpy()].to_numpy()
    q, x = np.sqrt(ev["Q2"].to_numpy())[i], ev["xB"].to_numpy()[i]
    had["ev"] = i
    had["px"] = had["pt"] * np.cos(had["phi"])
    had["py"] = had["pt"] * np.sin(had["phi"])
    had["pz"] = had["zlc"] * q / (2 * x)
    return len(ev), had


def zvals(a, b):
    m = a.merge(b, on="ev", suffixes=("_a", "_b"))
    pa = m[["px_a", "py_a", "pz_a"]].to_numpy()
    pb = m[["px_b", "py_b", "pz_b"]].to_numpy()
    cos = (pa * pb).sum(1) / (np.linalg.norm(pa, axis=1) * np.linalg.norm(pb, axis=1))
    return np.clip((1 + cos) / 2, 1e-12, None), (m["charge_a"] * m["charge_b"]).to_numpy()


def report(label, n, z, w):
    cnt = np.histogram(z, EDGES)[0] / n
    c = np.histogram(z, EDGES, weights=w)[0] / n
    dz, dlnz = np.diff(EDGES), np.diff(np.log(EDGES))
    print(f"  {label}")
    print("    z          " + " ".join(f"{v:8.1e}" for v in CENT))
    print("    dN/dz      " + " ".join(f"{v:8.3f}" for v in cnt / dz))
    print("    dN/dlnz    " + " ".join(f"{v:8.4f}" for v in cnt / dlnz))
    print("    dC/dz      " + " ".join(f"{v:8.3f}" for v in c / dz))
    print(f"    peak of dN/dz at z = {CENT[np.argmax(cnt / dz)]:.1e}, of dN/dlnz at z = {CENT[np.argmax(cnt / dlnz)]:.1e}")


for name in sys.argv[1:]:
    n, had = load(name)
    cur, tgt_all = had[had["current"] == 1], had[had["current"] == 0]
    tgt_acc = tgt_all[(tgt_all["eta_eic"] > 2) & (tgt_all["eta_eic"] < 3.5) & (tgt_all["zlc"] > 0.1)]
    lead = cur[cur["is_lead"] == 1][["ev", "charge", "px", "py", "pz"]]
    axis = cur.groupby("ev")[["px", "py", "pz"]].sum().reset_index()
    axis["charge"] = 1                                     # direction only
    photon = pd.DataFrame({"ev": np.arange(n), "px": 0.0, "py": 0.0, "pz": -1.0, "charge": 1})
    print(f"\n=== {name}: {n:,} tagged events ===")
    th_axis = np.degrees(np.arctan2(np.hypot(axis["px"], axis["py"]), -axis["pz"]))
    th_lead = np.degrees(np.arctan2(np.hypot(lead["px"], lead["py"]), -lead["pz"]))
    print(f"  median polar angle from -z: tagging hadron {np.median(th_lead):.1f} deg, "
          f"charged current-hemisphere sum {np.median(th_axis):.1f} deg")
    for lab, a, b in (("tagging hadron x target acceptance (as in the plots)", lead, tgt_acc),
                      ("tagging hadron x all target-hemisphere charged hadrons", lead, tgt_all),
                      ("charged current-hemisphere axis x target acceptance", axis, tgt_acc),
                      ("charged current-hemisphere axis x all target-hemisphere", axis, tgt_all),
                      ("photon direction (-z) x all target-hemisphere", photon, tgt_all)):
        z, w = zvals(a, b[["ev", "charge", "px", "py", "pz"]])
        report(lab, n, z, w)
