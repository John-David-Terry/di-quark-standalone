#!/usr/bin/env python3
"""
Compact the original Figure 3 observables CSV (Track A, surgery v2) into a small Parquet file
with only the columns the shower study compares against. The original CSV is read, never modified.

Default input: the 9.6M-event 5×41 GeV, 4 < Q < 5 GeV production used for ``signals.pdf``.
"""

from __future__ import annotations

import os
import argparse
from pathlib import Path

import pandas as pd

DEFAULT_CSV = (
    Path.home()
    / "Data/dis_q_bins_e5_ep41_480each/Q2_16_25__Q_4_5/merged/jet_hadron_pi_pm_observables.csv"
)
DEFAULT_OUT = Path(os.environ.get("STUDY_ROOT", Path.home() / "Data/fig3_shower_study")) / "original/original_fig3_compact.parquet"

COLUMNS = [
    "event_id",
    "sample",
    "pion",
    "ok",
    "reason",
    "Q2",
    "xB",
    "zLC",
    "phi_hJ",
    "jet_pjt_GeV",
    "pion_pht_GeV",
    "pion_lab_eta",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    df = pd.read_csv(args.csv, usecols=COLUMNS, low_memory=False)
    df["ok"] = df["ok"].astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    for c in ("sample", "pion", "reason"):
        df[c] = df[c].fillna("").astype("category")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    print(f"wrote {args.out}  rows={len(df):,}  ok={int(df['ok'].sum()):,}")
    print(df[df["ok"]].groupby(["sample", "pion"], observed=True).size())


if __name__ == "__main__":
    main()
