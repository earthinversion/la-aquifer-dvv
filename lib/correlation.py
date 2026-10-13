"""Daily single-station cross-component correlation functions (step 1).

For one station-day: download three-component waveforms, resample, cut 30-minute
windows, whiten and one-bit normalize each window within a frequency band, correlate
the component pairs E-N, E-Z and N-Z, and stack all windows and pairs into one
correlation function per band.
"""
import time

import numpy as np
from scipy.signal import fftconvolve

import config as C


class ServiceDown(Exception):
    """The data center is failing (as opposed to having no data for a request)."""


def fetch_day(client, day, tries=4, network=None, station=None):
    """Return an obspy Stream for one station-day, or None if the archive has no data.

    "No data" (HTTP 204) is permanent and is not retried. Server errors are retried
    with exponential back-off; persistent failure raises ServiceDown.
    """
    delay = 2.0
    for k in range(tries):
        try:
            return client.get_waveforms(network or C.NETWORK, station or C.STATION, "*",
                                        C.CHANNELS, day, day + 86400)
        except Exception as e:
            name = type(e).__name__
            if "NoData" in name:
                return None
            if k == tries - 1:
                raise ServiceDown(f"{station or C.STATION} {day.date}: {name}")
            time.sleep(delay)
            delay *= 2
    return None


def whiten_onebit(x, fmin, fmax, sps):
    """Spectral whitening within [fmin, fmax], then one-bit (sign) normalization."""
    x = x - x.mean()
    n = len(x)
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / sps)
    amp = np.abs(X)
    amp[amp == 0] = 1.0
    X = X / amp
    X[~((f >= fmin) & (f <= fmax))] = 0.0
    return np.sign(np.fft.irfft(X, n))


def preprocess(st):
    """Merge, detrend, taper and resample; return {'E'|'N'|'Z': Trace} or None."""
    st.merge(fill_value=0)
    comps = {}
    for tr in st:
        c = tr.stats.channel[-1]
        if c in "ENZ" and c not in comps:
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(0.01)
            if tr.stats.sampling_rate > C.SPS:
                tr.filter("lowpass", freq=0.4 * C.SPS, zerophase=True)
                tr.resample(C.SPS, no_filter=True)
            comps[c] = tr
    return comps if len(comps) == 3 else None


def day_correlations(comps):
    """Return {band: stacked correlation, band + '__nwin': [n windows]} or None."""
    nlag = int(C.MAXLAG_S * C.SPS)
    npts_win = int(C.WIN_S * C.SPS)
    length = min(len(tr.data) for tr in comps.values())
    out = {}
    for band, (fmin, fmax) in C.BANDS.items():
        stack, nwin = np.zeros(2 * nlag + 1), 0
        for s0 in range(0, length - npts_win, npts_win):
            segs = {}
            for c, tr in comps.items():
                seg = tr.data[s0:s0 + npts_win].astype(float)
                if np.count_nonzero(seg) < C.MIN_FILL * npts_win:
                    break
                segs[c] = whiten_onebit(seg, fmin, fmax, C.SPS)
            if len(segs) < 3:
                continue
            for a, b in C.PAIRS:
                cc = fftconvolve(segs[a], segs[b][::-1], mode="full")
                mid = len(cc) // 2
                stack += cc[mid - nlag: mid + nlag + 1]
            nwin += 1
        if nwin >= C.MIN_WINDOWS:
            out[band] = stack / nwin
            out[band + "__nwin"] = np.array([nwin])
    return out or None
