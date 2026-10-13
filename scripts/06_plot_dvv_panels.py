"""Draw the dv/v panels (A, B and the site map E) and save every plotted series.

A  1-4 Hz dv/v (cc >= 0.6; 20-sample rolling median, about 60 days), monthly LAX
   rainfall, dry and wet water years, and the seasonal model (annual + semiannual
   harmonics fitted over water years 2004-2011; dashed where extrapolated).
B  dv/v after removing the seasonal model (monthly mean, 3-month centered mean),
   with the fit window and the station's equipment changes.
E  Site map: the station and LAX highlighted, the other seismic and weather stations
   in gray, the nine Water Replenishment District nested wells, and the Rio Hondo gauge.

Inputs: results/dvv_<sta>_2004_2026.csv (step 2), data/rain_lax_monthly.csv (step 3),
data/sites/*.csv (step 4) and data/map_background.png (step 5). A station missing from
the step 4 tables is looked up in the FDSN station service. Axis limits come from
config.PANEL_LIMITS when the station is listed there, otherwise from the data.

Outputs go to results/figure_<sta>/: dvv_panels_<sta>.pdf, .png and the panel tables.

Usage:
    python scripts/06_plot_dvv_panels.py                  # LAF -> results/figure_laf/
    python scripts/06_plot_dvv_panels.py --station OLI    # -> results/figure_oli/
    python scripts/06_plot_dvv_panels.py --station OLI --fit-start 2005-10-01
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C
from lib.dvv import fit_seasonal, seasonal_design

INK, RAIN, MODEL = "#16324F", "#1B9E77", "#C0392B"
DRY, WET = "#EADCC0", "#D5E4F0"
STA_C, LAX_C, WELL_C = "#1F4E79", "#D95F02", "#2B8CBE"


def equipment_changes(sta):
    """Equipment-change dates of a station (step 4 table, else FDSN)."""
    f = C.SITES / "station_equipment_changes.csv"
    if f.exists():
        ch = pd.read_csv(f)
        ch = ch[(ch.network == C.NETWORK) & (ch.station == sta)]
        if len(ch):
            return pd.to_datetime(ch.date).tolist()
    from lib.stations import station_info
    return pd.to_datetime(station_info(sta)[1].date).tolist()


def coordinates(sta):
    """(lat, lon) of a station (step 4 table, else FDSN)."""
    f = C.SITES / "seismic_stations.csv"
    if f.exists():
        t = pd.read_csv(f)
        t = t[t.sta == sta]
        if len(t):
            return float(t.lat.iloc[0]), float(t.lon.iloc[0])
    from lib.stations import station_info
    return station_info(sta)[0]


def nice_limits(lo, hi, symmetric=False):
    """Axis limits padded to the next 0.1 and three ticks, from the data range."""
    if symmetric:
        lo, hi = -max(abs(lo), abs(hi)), max(abs(lo), abs(hi))
    lo, hi = np.floor((lo - .08) * 10) / 10, np.ceil((hi + .08) * 10) / 10
    if lo < 0 < hi:
        t = max(round(.8 * min(-lo, hi), 1), .1)
        return (float(lo), float(hi)), [-t, 0, t]
    return (float(lo), float(hi)), [round(lo + (hi - lo) * f, 2) for f in (.2, .5, .8)]


def prepare(sta):
    """Return a dict of every series the figure plots."""
    T0, T1 = pd.Timestamp(C.T0), pd.Timestamp(C.T1)
    d = pd.read_csv(C.dvv_file(sta), parse_dates=["date"])
    d = d[(d.station == sta) & (d.cc >= C.CC_MIN) & (d.band == C.PLOT_BAND)]
    raw = d.set_index("date").dvv_pct.sort_index()
    raw = raw[~raw.index.duplicated()]
    smooth = raw.rolling(C.SMOOTH_SAMPLES, center=True, min_periods=C.SMOOTH_MIN).median()
    gap = raw.index.to_series().diff().dt.days.rolling(3, center=True).max().fillna(0)
    smooth = smooth.where(gap.lt(C.GAP_DAYS))
    coef, n_fit = fit_seasonal(raw, C.FIT_START, C.FIT_END)
    grid = pd.date_range(T0, T1, freq="D")
    model = pd.Series(seasonal_design(grid) @ coef, index=grid)
    residual = (raw - seasonal_design(raw.index) @ coef).resample("MS").mean() \
        .rolling(3, center=True, min_periods=2).mean()
    rain = pd.read_csv(C.RAIN_FILE, parse_dates=["date"]).set_index("date").precip_mm.loc[T0:T1]
    wy = rain.index.year + (rain.index.month >= 10)            # water year ends 30 September
    total = rain.groupby(wy).sum(min_count=11)
    median = total.loc[C.WY_MEDIAN_YEARS[0]:C.WY_MEDIAN_YEARS[1]].median()
    classes = pd.Series(np.where(total < C.DRY_FRAC * median, "dry",
                                 np.where(total > C.WET_FRAC * median, "wet", "normal")),
                        index=total.index)
    if sta in C.PANEL_LIMITS:
        lim = C.PANEL_LIMITS[sta]
    else:
        a, a_ticks = nice_limits(min(smooth.min(), model.min()), max(smooth.max(), model.max()),
                                 symmetric=True)
        b, b_ticks = nice_limits(residual.min(), residual.max())
        lim = dict(a=a, a_ticks=a_ticks, b=b, b_ticks=b_ticks)
    return dict(sta=sta, raw=raw, smooth=smooth, coef=coef, n_fit=n_fit, grid=grid,
                model=model, residual=residual, rain=rain, classes=classes, median=median,
                lim=lim, changes=equipment_changes(sta), coords=coordinates(sta))


def save_tables(s, R):
    s["smooth"].rename("smoothed_dvv_pct").to_csv(R / "panelA_dvv_smoothed.csv")
    s["rain"].rename("rain_mm").to_csv(R / "panelA_rain_monthly.csv")
    s["model"].rename("model_dvv_pct").to_csv(R / "panelA_seasonal_model.csv")
    s["classes"].rename("class").rename_axis("water_year").to_csv(R / "panelA_water_year_classes.csv")
    s["residual"].rename("residual_dvv_pct").to_csv(R / "panelB_residual.csv")
    (R / "seasonal_fit.json").write_text(json.dumps(dict(
        coefficients=dict(zip(["a0", "a1", "b1", "a2", "b2"], np.round(s["coef"], 6).tolist())),
        station=f"{C.NETWORK}.{s['sta']}", band=C.PLOT_BAND,
        fit_window=[C.FIT_START, C.FIT_END], fit_samples=s["n_fit"],
        equipment_changes=[str(d.date()) for d in s["changes"]],
        water_year_median_rain_mm=round(float(s["median"]), 1)), indent=1) + "\n")


def merc(lat):
    return np.degrees(np.log(np.tan(np.pi / 4 + np.radians(lat) / 2)))


def draw(s, out):
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 8, "legend.fontsize": 8, "pdf.fonttype": 42,
                         "axes.linewidth": .6, "axes.spines.top": False,
                         "axes.spines.right": False, "xtick.major.width": .6,
                         "ytick.major.width": .6})
    W, H = 460.8, 262                                          # points; 6.4 in wide
    fig = plt.figure(figsize=(W / 72, H / 72))
    ax = lambda x, top, w, h: fig.add_axes([x / W, (H - top - h) / H, w / W, h / H])
    label = lambda x, top, t: fig.text(x / W, 1 - top / H, t, fontsize=9, weight="bold", va="top")
    T0, T1 = pd.Timestamp(C.T0), pd.Timestamp(C.T1)
    F0, F1 = pd.Timestamp(C.FIT_START), pd.Timestamp(C.FIT_END)

    def shade(a):
        for yr, c in s["classes"].items():
            if c != "normal":
                a.axvspan(pd.Timestamp(f"{yr - 1}-10-01"), pd.Timestamp(f"{yr}-09-30"),
                          color=DRY if c == "dry" else WET, lw=0, zorder=0)

    # ---- A: dv/v, rainfall, water years, seasonal model
    a = ax(69.12, 16, 352.512, 58)
    label(5.5, 16, "A")
    ar = a.twinx()
    shade(ar)
    ar.bar(s["rain"].index + pd.Timedelta(days=15), s["rain"].values, width=25, color=RAIN,
           alpha=.55, lw=0)
    ar.set_ylim(0, s["rain"].max() * 2.6)
    ar.set_yticks([0, 200])
    ar.set_ylabel("Rain (mm)", color=RAIN, labelpad=2)
    ar.tick_params(axis="y", colors=RAIN, pad=2, length=2)
    ar.spines["right"].set_visible(True)
    ar.spines["right"].set_color(RAIN)
    ar.set_zorder(a.get_zorder() - 1)
    a.patch.set_visible(False)
    a.plot(s["smooth"].index, s["smooth"], color=INK, lw=.9)
    g, m = s["grid"], s["model"]
    inside = (g >= F0) & (g <= F1)
    a.plot(g[inside], m[inside], color=MODEL, lw=.8)
    a.plot(g[g > F1], m[g > F1], color=MODEL, lw=.8, ls=(0, (2.2, 1.6)))
    a.axhline(0, color=".55", lw=.6)
    a.set_ylim(*s["lim"]["a"])
    a.set_yticks(s["lim"]["a_ticks"])
    a.set_xlim(T0, T1)
    a.set_ylabel(f"{s['sta']} dv/v (%)", color=INK, labelpad=2)
    a.tick_params(axis="x", labelbottom=False, length=0)
    a.tick_params(axis="y", pad=2, length=2)
    a.legend([Patch(facecolor=DRY), Patch(facecolor=WET), Line2D([], [], color=MODEL, lw=1)],
             ["Dry year", "Wet year", "Seasonal model"], loc="lower left",
             bbox_to_anchor=(0, 1.03), ncol=3, frameon=False, handlelength=1,
             handletextpad=.35, columnspacing=1, borderaxespad=0)

    # ---- B: residual, fit window, equipment changes
    b = ax(69.12, 92, 352.512, 44)
    label(5.5, 92, "B")
    shade(b)
    b.plot(s["residual"].index, s["residual"], color=INK, lw=1.2)
    b.add_patch(Rectangle((mdates.date2num(F0), .015), mdates.date2num(F1) - mdates.date2num(F0),
                          .97, transform=b.get_xaxis_transform(), facecolor="none",
                          edgecolor=".3", lw=.8, ls=(0, (3, 2))))
    for d in s["changes"]:
        b.axvline(d, color=".43", lw=.5, alpha=.8, ls=(0, (1.2, 1.4)))
    b.axhline(0, color=".55", lw=.6)
    b.set_ylim(*s["lim"]["b"])
    b.set_yticks(s["lim"]["b_ticks"])
    b.set_xlim(T0, T1)
    b.set_ylabel("Residual\n(%)", color=INK, labelpad=2, linespacing=1)
    b.xaxis.set_major_locator(mdates.YearLocator(4))
    b.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    b.tick_params(pad=2, length=2)
    b.legend([Line2D([], [], color=INK, lw=1.2), Patch(facecolor="none", edgecolor=".3", ls="--"),
              Line2D([], [], color=".43", ls=":")], ["Residual", "Fit window", "Equipment changes"],
             loc="lower left", bbox_to_anchor=(0, 1.03), ncol=3, frameon=False, handlelength=1,
             handletextpad=.35, columnspacing=.7, borderaxespad=0)

    # ---- E: site map
    R = C.MAP_REGION
    mapw = 130
    maph = mapw * (merc(R[3]) - merc(R[2])) / (R[1] - R[0])
    e = ax(69.12, 160, mapw, maph)
    label(5.5, 156, "E")
    e.imshow(plt.imread(C.MAP_FILE), extent=[R[0], R[1], merc(R[2]), merc(R[3])], aspect="auto")
    st = pd.read_csv(C.SITES / "seismic_stations.csv")
    wx = pd.read_csv(C.SITES / "weather_stations.csv")
    wells = pd.read_csv(C.SITES / "wrd_nested_wells.csv")
    gauge = pd.read_csv(C.SITES / "stream_gauge.csv").iloc[0]
    e.scatter(wells.lon, merc(wells.lat), s=16, facecolors="white", edgecolors=WELL_C, linewidths=.8)
    other_s, other_w = st[st.sta != s["sta"]], wx[wx.label != "LAX"]
    e.scatter(other_s.lon, merc(other_s.lat), s=20, marker="^", color=".6",
              edgecolors="white", linewidths=.4)
    e.scatter(other_w.lon, merc(other_w.lat), s=16, marker="s", facecolors="white",
              edgecolors=".6", linewidths=.8)
    (slat, slon), lax = s["coords"], wx[wx.label == "LAX"].iloc[0]
    e.scatter([slon], [merc(slat)], s=32, marker="^", color=STA_C, edgecolors="white",
              linewidths=.4)
    e.scatter([lax.lon], [merc(lax.lat)], s=24, marker="s", facecolors="white", edgecolors=LAX_C,
              linewidths=.9)
    e.scatter([gauge.lon], [merc(gauge.lat)], s=24, marker="v", color=RAIN, edgecolors="white",
              linewidths=.4)
    south = slat < lax.lat                                   # label away from LAX
    e.annotate(s["sta"], (slon, merc(slat)), xytext=(-4, -9) if south else (5, 4),
               textcoords="offset points", ha="right" if south else "left",
               va="top" if south else "bottom", color=STA_C)
    e.annotate("LAX", (lax.lon, merc(lax.lat)), xytext=(-6, 4), textcoords="offset points",
               ha="right", va="bottom", color=LAX_C)
    e.text(.50, .94, "Los Angeles", transform=e.transAxes, ha="center", va="top",
           bbox=dict(facecolor="white", alpha=.65, edgecolor="none", pad=.5))
    e.set_xlim(R[:2])
    e.set_ylim(merc(R[2]), merc(R[3]))
    e.axis("off")
    handles = [Line2D([], [], marker="^", color=STA_C, ls="none", ms=5),
               Line2D([], [], marker="s", mfc="white", mec=LAX_C, ls="none", ms=4),
               Line2D([], [], marker="o", mfc="white", mec=WELL_C, ls="none", ms=4),
               Line2D([], [], marker="v", color=RAIN, ls="none", ms=4),
               Line2D([], [], marker="^", color=".6", ls="none", ms=4)]
    fig.legend(handles, [f"{s['sta']} (seismic, dv/v)", "LAX (weather, rain)", "Nested monitoring wells",
                         "Rio Hondo stream gauge", "Other seismic and weather stations"],
               loc="upper left", bbox_to_anchor=((69.12 + mapw + 12) / W, 1 - 166 / H),
               frameon=False, handlelength=.8, handletextpad=.4, labelspacing=.5,
               borderpad=0, borderaxespad=0)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 450})):
        fig.savefig(out.with_suffix("." + ext), **kw)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--station", default=C.STATION)
    ap.add_argument("--network", default=C.NETWORK)
    ap.add_argument("--fit-start", default=C.FIT_START, help="seasonal-model fit window start")
    ap.add_argument("--fit-end", default=C.FIT_END, help="seasonal-model fit window end")
    args = ap.parse_args()
    C.NETWORK, C.FIT_START, C.FIT_END = args.network, args.fit_start, args.fit_end
    sta = args.station.upper()
    out_dir = C.figure_dir(sta)
    out_dir.mkdir(parents=True, exist_ok=True)
    s = prepare(sta)
    save_tables(s, out_dir)
    out = out_dir / f"dvv_panels_{sta.lower()}"
    draw(s, out)
    print(f"wrote {out}.pdf/.png and the panel tables in {out_dir}")
    print(f"equipment changes: {', '.join(str(d.date()) for d in s['changes']) or 'none'}")
    print(f"seasonal fit: {s['n_fit']} samples; median water-year rain {s['median']:.0f} mm")


if __name__ == "__main__":
    main()
