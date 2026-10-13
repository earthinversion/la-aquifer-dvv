"""Coordinates of the sites on the map, and the equipment history of the seismic stations.

  seismic stations   SCEDC FDSN station service: coordinates, and the dates where the
                     Z-channel metadata epoch changes (sensor or sampling rate)
  weather stations   NOAA NCEI Access Data Service (station location columns)
  nested wells       Water Replenishment District hydrograph pages (coordinates embedded)
  stream gauge       USGS Water Data API, monitoring location USGS-11101250 (Rio Hondo)

Writes data/sites/{seismic_stations,station_equipment_changes,weather_stations,
wrd_nested_wells,stream_gauge}.csv and data/sites/MANIFEST.json.

Usage:
    python scripts/04_fetch_site_metadata.py
"""
import io
import json
import re
import sys
import time
from datetime import date
from pathlib import Path

import pandas as pd
import requests
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C
from lib.stations import epochs_table, equipment_changes

UA = {"User-Agent": "la-aquifer-dvv research script"}


def seismic():
    """Coordinates and equipment changes of the stations on the map."""
    t0, t1 = UTCDateTime(C.START), UTCDateTime(C.END)
    inv = Client(C.FDSN_SERVER).get_stations(network=C.NETWORK,
                                             station=",".join(C.MAP_STATIONS),
                                             channel=C.CHANNELS, level="channel",
                                             starttime=t0, endtime=t1)
    rows = [dict(sta=s.code, lat=s.latitude, lon=s.longitude, elevation_m=s.elevation)
            for net in inv for s in net]
    coords = pd.DataFrame(rows).drop_duplicates("sta", keep="last").sort_values("sta")
    return coords, equipment_changes(epochs_table(inv, t0, t1), t0)


def weather():
    rows = []
    for sid, label in C.WEATHER_STATIONS.items():
        url = ("https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries"
               f"&stations={sid}&dataTypes=PRCP&startDate=2024-01-01&endDate=2024-01-31"
               "&includeStationLocation=1&format=csv")
        df = pd.read_csv(io.StringIO(requests.get(url, headers=UA, timeout=120).text))
        r = df.iloc[0]
        rows.append(dict(station=sid, label=label, lat=r.LATITUDE, lon=r.LONGITUDE,
                         elevation_m=r.ELEVATION))
    return pd.DataFrame(rows)


def wells():
    rows = []
    for site, wid in C.WRD_SITES.items():
        url = f"https://hydrographs.wrd.org/{wid}/2"
        html = requests.get(url, headers=UA, timeout=60).text
        lat = re.findall(r"3[34]\.\d{4,}", html)[0]
        lon = re.findall(r"-11[78]\.\d{4,}", html)[0]
        rows.append(dict(site=site, wrd_id=wid, lat=float(lat), lon=float(lon), url=url))
        time.sleep(1)
    return pd.DataFrame(rows)


def gauge():
    url = (f"https://api.waterdata.usgs.gov/ogcapi/v0/collections/monitoring-locations/"
           f"items/{C.GAUGE_ID}?f=json")
    j = requests.get(url, headers=UA, timeout=60).json()
    lon, lat = j["geometry"]["coordinates"][:2]
    return pd.DataFrame([dict(id=C.GAUGE_ID, name=j["properties"]["monitoring_location_name"],
                              lat=lat, lon=lon)]), url


def main():
    C.SITES.mkdir(parents=True, exist_ok=True)
    s, ch = seismic()
    s.to_csv(C.SITES / "seismic_stations.csv", index=False)
    ch.to_csv(C.SITES / "station_equipment_changes.csv", index=False)
    w = weather()
    w.to_csv(C.SITES / "weather_stations.csv", index=False)
    k = wells()
    k.to_csv(C.SITES / "wrd_nested_wells.csv", index=False)
    g, gurl = gauge()
    g.to_csv(C.SITES / "stream_gauge.csv", index=False)
    (C.SITES / "MANIFEST.json").write_text(json.dumps({
        "accessed": date.today().isoformat(),
        "seismic_stations": f"FDSN station service, {C.FDSN_SERVER}, network {C.NETWORK}",
        "station_equipment_changes": "Z-channel metadata epochs from the same service",
        "weather_stations": "NOAA NCEI Access Data Service, includeStationLocation=1",
        "wrd_nested_wells": "https://hydrographs.wrd.org/<id>/2 (coordinates in page)",
        "stream_gauge": gurl}, indent=1) + "\n")
    print(f"{len(s)} seismic stations ({len(ch)} equipment changes), {len(w)} weather "
          f"stations, {len(k)} wells, 1 gauge -> {C.SITES}")


if __name__ == "__main__":
    main()
