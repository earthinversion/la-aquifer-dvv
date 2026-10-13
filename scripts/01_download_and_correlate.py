"""Download waveforms from SCEDC and compute daily correlation functions.

Every third day from START to END (one fixed 3-day phase anchored at GRID_ORIGIN) is
fetched through the SCEDC FDSN web service and reduced to one single-station
cross-component correlation function per frequency band. Raw waveforms are not kept;
only data/corr/<STATION>_YYYYMMDD.npz is written (station LAF unless --station).

The run can be resumed: complete days already on disk are skipped. It stops with exit
code 3 if the data center fails repeatedly, so an outage is not mistaken for missing data.

Usage:
    python scripts/01_download_and_correlate.py                       # full record
    python scripts/01_download_and_correlate.py --start 2023-01-01 --end 2023-01-10
    python scripts/01_download_and_correlate.py --station OLI                # another station
"""
import argparse
import sys
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C
from lib.correlation import ServiceDown, day_correlations, fetch_day, preprocess

warnings.filterwarnings("ignore")
_local = threading.local()


def client():
    if not hasattr(_local, "c"):
        _local.c = Client(C.FDSN_SERVER, timeout=180)
    return _local.c


def grid_days(start, end):
    """Days on the fixed 3-day grid within [start, end)."""
    origin, step = UTCDateTime(C.GRID_ORIGIN), C.STRIDE_DAYS * 86400
    t0, t1 = UTCDateTime(start), UTCDateTime(end)
    k = int(np.ceil((t0 - origin) / step))
    day, out = origin + k * step, []
    while day < t1:
        out.append(day)
        day += step
    return out


def out_file(day, station=None):
    return C.CORR_DIR / f"{station or C.STATION}_{day.strftime('%Y%m%d')}.npz"


def complete(day):
    f = out_file(day)
    if not f.exists():
        return False
    try:
        with np.load(f) as z:
            return all(b in z for b in C.BANDS)
    except Exception:
        return False


class Breaker:
    """Trip after many consecutive data-center failures."""
    LIMIT = 60

    def __init__(self):
        self.lock, self.streak, self.tripped = threading.Lock(), 0, False

    def record(self, ok):
        with self.lock:
            self.streak = 0 if ok else self.streak + 1
            if self.streak >= self.LIMIT:
                self.tripped = True


breaker = Breaker()


def run_day(day):
    """Process one day; return 1 if a file was written, else 0."""
    try:
        st = fetch_day(client(), day)
    except ServiceDown:
        breaker.record(False)
        return 0
    breaker.record(True)
    if st is None:
        return 0
    comps = preprocess(st)
    res = day_correlations(comps) if comps else None
    if res:
        np.savez_compressed(out_file(day), **res)
        return 1
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--start", default=C.START)
    ap.add_argument("--end", default=C.END)
    ap.add_argument("--workers", type=int, default=C.N_WORKERS)
    ap.add_argument("--station", default=C.STATION)
    ap.add_argument("--network", default=C.NETWORK)
    args = ap.parse_args()
    C.STATION, C.NETWORK = args.station, args.network
    C.CORR_DIR.mkdir(parents=True, exist_ok=True)
    days = grid_days(args.start, args.end)
    todo = [d for d in days if not complete(d)]
    print(f"{len(days)} days on the grid, {len(todo)} to process", flush=True)
    t0, n = time.time(), 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for i, got in enumerate(ex.map(run_day, todo), 1):
            n += got
            if i % 50 == 0:
                print(f"  {i}/{len(todo)} processed, {n} written, "
                      f"{(time.time() - t0) / 60:.1f} min", flush=True)
            if breaker.tripped:
                break
    if breaker.tripped:
        print(f"STOPPED: the data center failed {Breaker.LIMIT} times in a row; "
              f"rerun later to resume ({n} days written)", flush=True)
        sys.exit(3)
    print(f"done: {n} days written in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
