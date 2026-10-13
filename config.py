"""Settings for the LAF dv/v workflow. Every script reads its parameters from here."""
from pathlib import Path

import numpy as np

# paths
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
CORR_DIR = DATA / "corr"                 # daily correlation functions (step 1, not tracked)
SITES = DATA / "sites"                   # station, well and gauge coordinates (step 4)
RESULTS = ROOT / "results"

# seismic data (step 1)
NETWORK, STATION = "CI", "LAF"           # Southern California Seismic Network station LAF
CHANNELS = "BH?"                         # broadband, three components
FDSN_SERVER = "SCEDC"                    # Southern California Earthquake Data Center
START, END = "2003-10-01", "2026-09-16"  # water year 2004 through September 2026
GRID_ORIGIN, STRIDE_DAYS = "2012-01-01", 3   # every third day, on one fixed phase
N_WORKERS = 8                            # parallel station-day downloads

SPS = 20.0                               # resampling rate (Hz)
WIN_S = 1800                             # correlation window length (s)
MIN_FILL = 0.8                           # a window needs >= 80 % non-zero samples
MIN_WINDOWS = 12                         # >= 12 windows (6 h) per day and band
MAXLAG_S = 60.0                          # correlation lag kept on each side (s)
PAIRS = [("E", "N"), ("E", "Z"), ("N", "Z")]   # single-station cross-component pairs
BANDS = {"shallow_1_4Hz": (1.0, 4.0),    # nominal sensitivity ~50 m
         "mid_05_2Hz": (0.5, 2.0),       # ~200 m
         "deep_02_08Hz": (0.2, 0.8)}     # ~450 m

# dv/v (step 2)
CODA = {"shallow_1_4Hz": (2.0, 12.0),    # |lag| window (s) used for stretching
        "mid_05_2Hz": (2.0, 18.0),
        "deep_02_08Hz": (3.0, 18.0)}
EPS = np.linspace(-0.03, 0.03, 601)      # stretch factors, +/-3 % in 0.01 % steps


def dvv_file(station):
    """dv/v table of a station, e.g. results/dvv_laf_2004_2026.csv."""
    return RESULTS / f"dvv_{station.lower()}_2004_2026.csv"


# rainfall (step 3)
RAIN_STATION = "USW00023174"             # Los Angeles International Airport (LAX)
RAIN_START, RAIN_END = "2003-01-01", "2026-09-15"
RAIN_MIN_DAYS = 20                       # a monthly total needs >= 20 reported days
RAIN_FILE = DATA / "rain_lax_monthly.csv"

# site metadata (step 4)
MAP_STATIONS = ["BLC", "CHN", "LAF", "LLS", "OLI", "PDR"]   # seismic stations on the map
WEATHER_STATIONS = {"USW00023174": "LAX", "USW00023152": "Burbank",
                    "USW00023129": "Long Beach", "USW00003171": "Riverside"}
WRD_SITES = {"GAR1": 100020, "PM-3": 100034, "PM-4": 100038, "ING1": 100091,
             "SFS1": 100096, "LMR1": 100876, "HAW1": 100887, "WHT1": 101735,
             "GAR2": 101804}             # Water Replenishment District nested wells: page id
GAUGE_ID = "USGS-11101250"               # Rio Hondo above Whittier Narrows Dam

# map background (step 5)
MAP_REGION = [-119.02, -117.38, 33.50, 34.50]   # lon_min, lon_max, lat_min, lat_max
MAP_WIDTH_IN = 1.80556                   # GMT Mercator width, matches the plotted size
MAP_FILE = DATA / "map_background.png"
RELIEF_GRID = DATA / "la_relief.nc"
RELIEF_CPT = DATA / "relief_muted.cpt"

# figure (step 6)
CC_MIN = 0.6                             # keep days whose stretching correlation >= 0.6
PLOT_BAND = "shallow_1_4Hz"
SMOOTH_SAMPLES, SMOOTH_MIN = 20, 7       # rolling median, ~60 days at the 3-day stride
GAP_DAYS = 30                            # break the plotted line at gaps longer than this
FIT_START, FIT_END = "2003-10-01", "2011-09-30"   # seasonal model fit, water years 2004-2011
WY_MEDIAN_YEARS = (2005, 2025)           # reference for dry and wet water years
DRY_FRAC, WET_FRAC = 0.75, 1.25          # dry < 75 %, wet > 125 % of the median
T0, T1 = "2003-10-01", "2026-10-01"      # plotted time range

# Axis limits of panels A and B. Stations not listed get limits from their data.
PANEL_LIMITS = {"LAF": dict(a=(-.5, .5), a_ticks=[-.4, 0, .4],
                            b=(-.3, .38), b_ticks=[-.2, 0, .2])}


def figure_dir(station):
    """Figure and panel tables of a station, e.g. results/figure_laf/."""
    return RESULTS / f"figure_{station.lower()}"
