"""Station metadata from the FDSN station service: channel epochs, equipment
changes and coordinates."""
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

import config as C

FAR_FUTURE = UTCDateTime("2100-01-01")


def km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km."""
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(np.asarray(lon2) - lon1)
    return 6371.0 * np.arccos(np.clip(np.sin(p1) * np.sin(p2) +
                                      np.cos(p1) * np.cos(p2) * np.cos(dl), -1, 1))


def epochs_table(inv, t0, t1):
    """One row per channel epoch, clipped to [t0, t1]."""
    rows = []
    for net in inv:
        for sta in net:
            for ch in sta:
                start = max(ch.start_date, t0)
                end = min(ch.end_date or FAR_FUTURE, t1)
                if end <= start:
                    continue
                rows.append(dict(network=net.code, station=sta.code, lat=sta.latitude,
                                 lon=sta.longitude, location=ch.location_code,
                                 channel=ch.code, start=start.datetime, end=end.datetime,
                                 sample_rate=ch.sample_rate,
                                 sensor=(ch.sensor.description if ch.sensor and
                                         ch.sensor.description else "")))
    return pd.DataFrame(rows)


def equipment_changes(ep, t0):
    """Dates after t0 where a station's Z-channel metadata epoch changes.

    Each new epoch is listed, with what changed (sensor, sampling rate) or "new metadata
    epoch" when neither did; any of these can coincide with a step in dv/v.
    """
    out = []
    for (net, sta), g in ep.groupby(["network", "station"]):
        prev = None
        for _, r in g[g.channel.str.endswith("Z")].sort_values("start").iterrows():
            key = (r.sensor, r.sample_rate)
            if prev is not None and r.start > t0.datetime:
                what = []
                if key[0] != prev[0]:
                    what.append(f"sensor {prev[0].split(',')[0]} -> {key[0].split(',')[0]}")
                if key[1] != prev[1]:
                    what.append(f"{prev[1]:g} -> {key[1]:g} sps")
                out.append(dict(network=net, station=sta, date=r.start.date(),
                                change="; ".join(what) or "new metadata epoch"))
            prev = key
    return pd.DataFrame(out, columns=["network", "station", "date", "change"])


def station_info(station, network=None, start=None, end=None):
    """Coordinates and equipment changes of one station, from the FDSN station service."""
    t0, t1 = UTCDateTime(start or C.START), UTCDateTime(end or C.END)
    inv = Client(C.FDSN_SERVER).get_stations(network=network or C.NETWORK, station=station,
                                             channel=C.CHANNELS, level="channel",
                                             starttime=t0, endtime=t1)
    ep = epochs_table(inv, t0, t1)
    return (float(ep.lat.iloc[0]), float(ep.lon.iloc[0])), equipment_changes(ep, t0)
