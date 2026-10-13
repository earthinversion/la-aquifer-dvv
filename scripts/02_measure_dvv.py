"""Measure dv/v for every day and band by coda stretching.

Each day's correlation function is normalized; the reference is the mean of all
normalized days in the record. For each day, the stretch factor that best matches the
reference over the band's coda window gives dv/v (percent) and a correlation
coefficient cc. No cc threshold is applied here; step 6 keeps days with cc >= 0.6.

Usage:
    python scripts/02_measure_dvv.py                       # uses data/corr/
    python scripts/02_measure_dvv.py --corr-dir DIR1 DIR2  # other caches (first wins)
    python scripts/02_measure_dvv.py --station OLI         # -> results/dvv_oli_2004_2026.csv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C
from lib.dvv import coda_mask, lag_axis, stretch_dvv


def corr_files(dirs, since):
    """{date: file} on the 3-day grid from `since`; earlier directories take priority."""
    origin, out = pd.Timestamp(C.GRID_ORIGIN), {}
    for d in dirs:
        for f in sorted(Path(d).glob(f"{C.STATION}_*.npz")):
            day = pd.Timestamp(f.stem.split("_")[1])
            if day >= since and (day - origin).days % C.STRIDE_DAYS == 0:
                out.setdefault(day, f)
    return dict(sorted(out.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--corr-dir", nargs="+", default=[str(C.CORR_DIR)])
    ap.add_argument("--since", default=C.START)
    ap.add_argument("--station", default=C.STATION)
    ap.add_argument("--out", help="default results/dvv_<station>_2004_2026.csv")
    args = ap.parse_args()
    C.STATION = args.station
    args.out = args.out or str(C.dvv_file(args.station))
    files = corr_files(args.corr_dir, pd.Timestamp(args.since))
    print(f"{len(files)} correlation files", flush=True)
    lags, rows = lag_axis(), []
    for band, (lo, hi) in C.CODA.items():
        mask = coda_mask(lags, lo, hi)
        days, ccfs = [], []
        for day, f in files.items():
            with np.load(f) as z:
                if band not in z:
                    continue
                c = z[band]
            if np.all(np.isfinite(c)) and not np.allclose(c, 0):
                days.append(day)
                ccfs.append(c)
        X = np.vstack(ccfs)
        norms = np.linalg.norm(X, axis=1, keepdims=True)
        X /= np.where(norms == 0, 1, norms)
        ref = X.mean(axis=0)
        for day, cur in zip(days, X):
            dvv, cc = stretch_dvv(ref, cur, lags, mask)
            rows.append({"date": day, "station": C.STATION, "band": band,
                         "dvv_pct": dvv, "cc": cc})
        print(f"  {band:14s} {len(days)} days", flush=True)
    df = pd.DataFrame(rows).sort_values(["station", "band", "date"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    good = (df.cc >= C.CC_MIN).groupby(df.band).sum()
    print(f"wrote {args.out}; days with cc >= {C.CC_MIN}: {good.to_dict()}")


if __name__ == "__main__":
    main()
