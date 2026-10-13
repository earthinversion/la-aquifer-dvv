"""Velocity change by coda stretching (step 2) and the seasonal model (step 6)."""
import numpy as np
import pandas as pd

import config as C


def lag_axis():
    n = int(C.MAXLAG_S * C.SPS)
    return np.arange(-n, n + 1) / C.SPS


def coda_mask(lags, lo, hi):
    """True where lo <= |lag| <= hi (both causal and acausal coda)."""
    a = np.abs(lags)
    return (a >= lo) & (a <= hi)


def stretch_dvv(ref, cur, lags, mask):
    """Grid-search the stretch factor that best maps cur onto ref over the coda.

    Returns (dv/v in percent, correlation coefficient). A velocity increase
    compresses arrival times, so dv/v = -eps.
    """
    best_cc, best_eps = -2.0, 0.0
    r = ref[mask]
    for e in C.EPS:
        s = np.interp(lags * (1.0 + e), lags, cur)[mask]
        n = np.linalg.norm(r) * np.linalg.norm(s)
        cc = np.dot(r, s) / n if n else -2.0
        if cc > best_cc:
            best_cc, best_eps = cc, e
    return -best_eps * 100.0, best_cc


def seasonal_design(index):
    """Design matrix with a constant and annual and semiannual harmonics of the day of year."""
    tau = index.dayofyear / 365.25
    return np.column_stack([np.ones(len(index))] +
                           [f(2 * np.pi * k * tau) for k in (1, 2) for f in (np.cos, np.sin)])


def fit_seasonal(series, start, end):
    """Least-squares harmonic coefficients from the samples between start and end."""
    m = (series.index >= pd.Timestamp(start)) & (series.index <= pd.Timestamp(end)) \
        & series.notna().to_numpy()
    coef = np.linalg.lstsq(seasonal_design(series.index)[m], series.to_numpy()[m],
                           rcond=None)[0]
    return coef, int(m.sum())
