"""Monthly rainfall at Los Angeles International Airport (LAX) from NOAA NCEI.

Downloads daily precipitation (GHCN-Daily station USW00023174) from the NCEI Access
Data Service and sums it by calendar month; a month needs at least 20 reported days.

Usage:
    python scripts/03_fetch_rainfall.py
"""
import io
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

URL = ("https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries"
       "&stations={sta}&dataTypes=PRCP&startDate={t0}&endDate={t1}&units=metric&format=csv")


def main():
    url = URL.format(sta=C.RAIN_STATION, t0=C.RAIN_START, t1=C.RAIN_END)
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    daily = pd.read_csv(io.StringIO(r.text), parse_dates=["DATE"]).set_index("DATE").PRCP
    monthly = daily.astype(float).resample("MS").sum(min_count=C.RAIN_MIN_DAYS)
    out = monthly.rename("precip_mm").rename_axis("date").to_frame()
    out.to_csv(C.RAIN_FILE)
    (C.DATA / "rain_lax_monthly_MANIFEST.json").write_text(json.dumps(
        {"source": "NOAA NCEI Access Data Service, GHCN-Daily PRCP", "url": url,
         "station": C.RAIN_STATION, "accessed": date.today().isoformat(),
         "days": int(daily.notna().sum()), "months": int(out.precip_mm.notna().sum())},
        indent=1) + "\n")
    print(f"{len(daily)} days -> {out.precip_mm.notna().sum()} months; wrote {C.RAIN_FILE}")


if __name__ == "__main__":
    main()
