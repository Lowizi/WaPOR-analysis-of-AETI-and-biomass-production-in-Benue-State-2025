"""Analyze 2025 WaPOR monthly AETI and NPP for Benue State, Nigeria.

Outputs a forum-ready figure (with author caption), a monthly CSV, a text
summary, and copy/paste-ready discussion-forum content. Raster I/O uses GDAL
because rasterio is blocked on this computer.
"""

import calendar
import json
import os
import warnings

warnings.filterwarnings("ignore")

os.environ.setdefault("GDAL_DATA", r"C:\Program Files\QGIS 3.44.15\share\gdal")
os.environ.setdefault("PROJ_DATA", r"C:\Program Files\QGIS 3.44.15\share\proj")
os.environ.setdefault("PROJ_LIB", r"C:\Program Files\QGIS 3.44.15\share\proj")

from osgeo import gdal

gdal.UseExceptions()

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BOUNDARY = os.path.join(HERE, "Benue_State.geojson")
NIGERIA_ADM1 = os.path.join(HERE, "Nigeria_ADM1_geoBoundaries.geojson")
AETI_FP = os.path.join(HERE, "wapor", "Benue_WaPORv3_L1_AETI_monthly_2025.tif")
NPP_FP = os.path.join(HERE, "wapor", "Benue_WaPORv3_L1_NPP_monthly_2025.tif")
AUTHOR = "Louis Micheal Obam"
YEAR = 2025


def read_monthly_stack(fp, fallback_scale):
    ds = gdal.Open(fp)
    if ds is None:
        raise FileNotFoundError(fp)
    if ds.RasterCount != 12:
        raise ValueError(f"Expected 12 monthly bands in {fp}, got {ds.RasterCount}")

    arrays, metadata = [], []
    for i in range(1, 13):
        b = ds.GetRasterBand(i)
        raw = b.ReadAsArray().astype("float64")
        nd = b.GetNoDataValue()
        mask = ~np.isfinite(raw)
        if nd is not None:
            mask |= raw == nd

        scale = b.GetScale()
        if scale is None or scale == 1.0:
            # WaPOR L1 scale factors: AETI 0.1, NPP 0.001. Some GDAL
            # operations do not copy scale metadata, so retain a documented
            # fallback and report the scale used.
            scale = fallback_scale
        offset = b.GetOffset() or 0.0
        a = raw * scale + offset
        a[mask] = np.nan
        arrays.append(a)
        metadata.append(dict(b.GetMetadata()))

    out = np.stack(arrays)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    nx, ny = ds.RasterXSize, ds.RasterYSize
    dtype = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
    ds = None
    return out, gt, proj, metadata, float(scale), dtype, nx, ny


def annual_sum(stack):
    valid_count = np.isfinite(stack).sum(axis=0)
    summed = np.nansum(stack, axis=0)
    summed[valid_count != 12] = np.nan
    return summed


def weighted_mean(a, weights):
    m = np.isfinite(a)
    return float(np.sum(a[m] * weights[m]) / np.sum(weights[m]))


def percentiles(a):
    v = a[np.isfinite(a)]
    return tuple(float(x) for x in np.percentile(v, [5, 50, 95]))


def draw_north_arrow(ax):
    ax.annotate(
        "N",
        xy=(0.93, 0.92),
        xytext=(0.93, 0.77),
        xycoords="axes fraction",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        arrowprops=dict(facecolor="#222222", width=2.5, headwidth=8),
    )


print("Reading WaPOR rasters...")
aeti, gt_a, proj_a, md_a, scale_a, dtype_a, nx, ny = read_monthly_stack(AETI_FP, 0.1)
npp, gt_n, proj_n, md_n, scale_n, dtype_n, nx_n, ny_n = read_monthly_stack(NPP_FP, 0.001)

if aeti.shape != npp.shape or gt_a != gt_n:
    raise ValueError("AETI and NPP grids do not match")

# Pixel-centre latitude weights provide an area-weighted state mean on the
# geographic WaPOR grid. The longitude factor cancels in the weighted mean.
ys = gt_a[3] + (np.arange(ny) + 0.5) * gt_a[5]
latitude_weights = np.cos(np.deg2rad(ys))[:, None] * np.ones((1, nx))

months = [calendar.month_abbr[i] for i in range(1, 13)]
monthly_aeti = [weighted_mean(aeti[i], latitude_weights) for i in range(12)]
monthly_npp = [weighted_mean(npp[i], latitude_weights) for i in range(12)]

annual_aeti = annual_sum(aeti)
annual_npp = annual_sum(npp)
annual_tbp = annual_npp * 22.222 / 1000.0  # ton/ha/year, course formula
annual_wpb = annual_tbp / annual_aeti * 100.0  # kg/m3

mean_aeti = weighted_mean(annual_aeti, latitude_weights)
mean_npp = weighted_mean(annual_npp, latitude_weights)
mean_tbp = weighted_mean(annual_tbp, latitude_weights)
mean_wpb = weighted_mean(annual_wpb, latitude_weights)

valid = np.isfinite(annual_aeti) & np.isfinite(annual_npp)
spatial_r = float(np.corrcoef(annual_aeti[valid], annual_npp[valid])[0, 1])

aeti_p05, aeti_median, aeti_p95 = percentiles(annual_aeti)
npp_p05, npp_median, npp_p95 = percentiles(annual_npp)
tbp_p05, tbp_median, tbp_p95 = percentiles(annual_tbp)

peak_aeti_i = int(np.argmax(monthly_aeti))
peak_npp_i = int(np.argmax(monthly_npp))
low_aeti_i = int(np.argmin(monthly_aeti))

benue = gpd.read_file(BOUNDARY)
nigeria = gpd.read_file(NIGERIA_ADM1)
if benue.crs.to_epsg() != 4326:
    benue = benue.to_crs(4326)
if nigeria.crs.to_epsg() != 4326:
    nigeria = nigeria.to_crs(4326)
area_km2 = float(benue.to_crs(32632).area.sum() / 1e6)

xmin = gt_a[0]
xmax = gt_a[0] + nx * gt_a[1]
ymax = gt_a[3]
ymin = gt_a[3] + ny * gt_a[5]
extent = (xmin, xmax, ymin, ymax)

# Approximate ground resolution at the state's centre.
lat0 = float(benue.geometry.union_all().centroid.y)
pixel_x_m = abs(gt_a[1]) * 111_320.0 * np.cos(np.deg2rad(lat0))
pixel_y_m = abs(gt_a[5]) * 110_574.0
resolution_m = int(round(np.sqrt(pixel_x_m * pixel_y_m) / 10.0) * 10)

monthly = pd.DataFrame(
    {
        "month": months,
        "AETI_mm_month": monthly_aeti,
        "NPP_gC_m2_month": monthly_npp,
        "TBP_ton_ha_month": np.array(monthly_npp) * 22.222 / 1000.0,
    }
)
monthly.to_csv(os.path.join(HERE, "Benue_2025_monthly_summary.csv"), index=False)

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "bold"})
fig = plt.figure(figsize=(13.2, 8.9), facecolor="white")
gs = fig.add_gridspec(2, 3, width_ratios=[0.88, 1.15, 1.15], height_ratios=[1, 1],
                      hspace=0.26, wspace=0.22)

# A. Location within Nigeria
ax0 = fig.add_subplot(gs[:, 0])
nigeria.boundary.plot(ax=ax0, color="#8b8b8b", linewidth=0.55)
nigeria.plot(ax=ax0, color="#eeeeee", edgecolor="#8b8b8b", linewidth=0.4)
benue.plot(ax=ax0, color="#e69f00", edgecolor="#7a4f00", linewidth=1.2)
centroid = benue.geometry.union_all().centroid
ax0.annotate("Benue State", (centroid.x, centroid.y), xytext=(6, 5),
             textcoords="offset points", fontsize=9, fontweight="bold", color="#4b3100")
ax0.set_title("A  Study area", loc="left", fontsize=12)
ax0.set_xlabel("Longitude (°E)")
ax0.set_ylabel("Latitude (°N)")
ax0.grid(color="#dddddd", linewidth=0.5)
ax0.set_aspect("equal")
ax0.text(0.04, 0.035, f"Area: {area_km2:,.0f} km²", transform=ax0.transAxes,
         fontsize=9, bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#bbbbbb", alpha=.9))
draw_north_arrow(ax0)

# Robust map colour limits prevent isolated outliers from hiding patterns.
av = annual_aeti[np.isfinite(annual_aeti)]
tv = annual_tbp[np.isfinite(annual_tbp)]
a_lo, a_hi = np.percentile(av, [2, 98])
t_lo, t_hi = np.percentile(tv, [2, 98])

# B. AETI map
ax1 = fig.add_subplot(gs[0, 1])
im1 = ax1.imshow(annual_aeti, extent=extent, origin="upper", cmap="YlGnBu",
                 vmin=a_lo, vmax=a_hi, interpolation="nearest")
benue.boundary.plot(ax=ax1, color="#252525", linewidth=0.75)
ax1.set_title("B  Annual AETI (2025)", loc="left", fontsize=12)
ax1.set_xlabel("Longitude (°E)")
ax1.set_ylabel("Latitude (°N)")
cb1 = fig.colorbar(im1, ax=ax1, fraction=0.045, pad=0.018)
cb1.set_label("mm/year", labelpad=2)
cb1.set_ticks([600, 800, 1000, 1200])
ax1.tick_params(axis="y", pad=1)
ax1.text(0.03, 0.04, f"Area-weighted mean: {mean_aeti:,.0f} mm/year",
         transform=ax1.transAxes, fontsize=8.7,
         bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#bbbbbb", alpha=.9))
draw_north_arrow(ax1)

# C. Derived TBP map
ax2 = fig.add_subplot(gs[0, 2])
im2 = ax2.imshow(annual_tbp, extent=extent, origin="upper", cmap="YlOrBr",
                 vmin=t_lo, vmax=t_hi, interpolation="nearest")
benue.boundary.plot(ax=ax2, color="#252525", linewidth=0.75)
ax2.set_title("C  Total biomass production", loc="left", fontsize=12)
ax2.set_xlabel("Longitude (°E)")
ax2.tick_params(axis="y", labelleft=False)
cb2 = fig.colorbar(im2, ax=ax2, fraction=0.045, pad=0.018)
cb2.set_label("ton/ha/year", labelpad=2)
cb2.set_ticks([6, 8, 10, 12, 14, 16])
ax2.tick_params(axis="y", pad=1)
ax2.text(0.03, 0.04, f"Area-weighted mean: {mean_tbp:.1f} ton/ha/year",
         transform=ax2.transAxes, fontsize=8.7,
         bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#bbbbbb", alpha=.9))
draw_north_arrow(ax2)

# D. Monthly area-weighted time series
ax3 = fig.add_subplot(gs[1, 1:])
x = np.arange(12)
bars = ax3.bar(x, monthly_aeti, width=0.68, color="#56b4e9", edgecolor="#1c6688",
               linewidth=0.65, label="AETI")
ax3.set_xticks(x, months)
ax3.set_ylabel("AETI (mm/month)", color="#1c6688")
ax3.tick_params(axis="y", labelcolor="#1c6688")
ax3.grid(axis="y", color="#dddddd", linewidth=0.65)
ax3.set_axisbelow(True)
ax3.set_title("D  Monthly state-average pattern", loc="left", fontsize=12)

ax3b = ax3.twinx()
line = ax3b.plot(x, monthly_npp, color="#007f5f", marker="o", markersize=5,
                 linewidth=2.1, label="NPP")[0]
ax3b.fill_between(x, monthly_npp, color="#009e73", alpha=0.09)
ax3b.set_ylabel("NPP (gC/m²/month)", color="#006c50")
ax3b.tick_params(axis="y", labelcolor="#006c50")
ax3.legend([bars, line], ["AETI", "NPP"], loc="upper left", frameon=True, bbox_to_anchor=(0.01, 0.99))

stats_text = (
    f"Peak AETI: {months[peak_aeti_i]} ({monthly_aeti[peak_aeti_i]:.0f} mm)\n"
    f"Peak NPP: {months[peak_npp_i]} ({monthly_npp[peak_npp_i]:.1f} gC/m²)\n"
    f"Annual AETI–NPP spatial r: {spatial_r:.2f}"
)
ax3.text(0.20, 0.97, stats_text, transform=ax3.transAxes, ha="left", va="top",
         fontsize=8.7, bbox=dict(boxstyle="round,pad=.38", fc="white", ec="#bbbbbb", alpha=.95))

fig.suptitle("WaPOR analysis of Benue State, Nigeria — 2025", fontsize=17,
             fontweight="bold", y=0.985)
fig.text(0.5, 0.024,
         f"Created by {AUTHOR}  |  FAO WaPOR v3 Level 1 (~{resolution_m} m): monthly AETI and NPP  |  "
         "Boundary: geoBoundaries/GRID3 (CC BY 4.0)",
         ha="center", fontsize=9.1, color="#3f3f3f", fontweight="bold")
fig.text(0.5, 0.007,
         "TBP derived from NPP using TBP = NPP × 22.222 / 1000. Map colour ranges use the 2nd–98th percentiles.",
         ha="center", fontsize=7.9, color="#626262")

figure_fp = os.path.join(HERE, "Benue_WaPOR_analysis_2025_Louis_Micheal_Obam.png")
fig.savefig(figure_fp, dpi=145, bbox_inches="tight", facecolor="white", pil_kwargs={"optimize": True})
plt.close(fig)

# ------------------------------------------------------------------ reports
summary = {
    "author": AUTHOR,
    "area": "Benue State, Nigeria",
    "area_km2": area_km2,
    "period": "2025-01-01 to 2025-12-31",
    "products": ["WaPOR v3 L1-AETI-M", "WaPOR v3 L1-NPP-M"],
    "approx_resolution_m": resolution_m,
    "valid_pixels": int(valid.sum()),
    "mean_annual_AETI_mm": mean_aeti,
    "mean_annual_NPP_gC_m2": mean_npp,
    "mean_annual_TBP_ton_ha": mean_tbp,
    "mean_biomass_water_productivity_kg_m3": mean_wpb,
    "AETI_p05_median_p95_mm": [aeti_p05, aeti_median, aeti_p95],
    "NPP_p05_median_p95_gC_m2": [npp_p05, npp_median, npp_p95],
    "TBP_p05_median_p95_ton_ha": [tbp_p05, tbp_median, tbp_p95],
    "peak_AETI_month": months[peak_aeti_i],
    "peak_AETI_mm_month": monthly_aeti[peak_aeti_i],
    "lowest_AETI_month": months[low_aeti_i],
    "lowest_AETI_mm_month": monthly_aeti[low_aeti_i],
    "peak_NPP_month": months[peak_npp_i],
    "peak_NPP_gC_m2_month": monthly_npp[peak_npp_i],
    "annual_AETI_NPP_spatial_correlation_r": spatial_r,
    "raster_shape": [12, ny, nx],
    "AETI_scale_used": scale_a,
    "NPP_scale_used": scale_n,
    "AETI_dtype": dtype_a,
    "NPP_dtype": dtype_n,
}

with open(os.path.join(HERE, "Benue_2025_analysis_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

forum_subject = "WaPOR analysis of AETI and biomass production in Benue State, Nigeria (2025)"
forum_message = f"""1. Study area

I studied Benue State in north-central Nigeria. The boundary used in this analysis covers approximately {area_km2:,.0f} km². I selected 2025 because it is the latest complete calendar year available for both WaPOR layers used.

2. Data layers

• FAO WaPOR v3 Level 1 monthly Actual Evapotranspiration and Interception (L1-AETI-M), in mm/month.
• FAO WaPOR v3 Level 1 monthly Net Primary Production (L1-NPP-M), in gC/m²/month.
• Benue State administrative boundary from geoBoundaries/GRID3 (CC BY 4.0).

The WaPOR rasters have an approximate ground resolution of {resolution_m} m in Benue State.

3. Python spatial analysis

I used Python with GDAL, GeoPandas, NumPy, pandas and Matplotlib. I clipped 12 monthly AETI and NPP rasters to the Benue State boundary, removed NoData cells, applied the WaPOR scale factors and summed the monthly layers to annual totals. I calculated latitude-area-weighted state means, derived total biomass production from NPP, examined the spatial correlation between annual AETI and NPP, and plotted annual maps and a monthly time series.

4. Main findings

The area-weighted mean annual AETI was {mean_aeti:,.0f} mm/year. Mean annual NPP was {mean_npp:,.1f} gC/m²/year, equivalent to approximately {mean_tbp:.1f} ton/ha/year of total biomass production using the course conversion formula. Monthly mean AETI was highest in {months[peak_aeti_i]} ({monthly_aeti[peak_aeti_i]:.0f} mm/month) and lowest in {months[low_aeti_i]} ({monthly_aeti[low_aeti_i]:.0f} mm/month). NPP peaked in {months[peak_npp_i]} ({monthly_npp[peak_npp_i]:.1f} gC/m²/month). Annual AETI and NPP had a spatial correlation of r = {spatial_r:.2f}, showing that locations with greater water consumption generally also had greater vegetation production.

The attached image was created with my Python script and includes my name in the caption.

— {AUTHOR}
"""

forum_fp = os.path.join(HERE, "FORUM_SUBMISSION.txt")
with open(forum_fp, "w", encoding="utf-8") as f:
    f.write("SUBJECT\n" + forum_subject + "\n\nMESSAGE\n" + forum_message)

print("\nRESULTS")
print(f"Area                           : {area_km2:,.0f} km2")
print("Raster                         : 12 x %d x %d; ~%d m" % (ny, nx, resolution_m))
print("Valid annual pixels            :", int(valid.sum()))
print("Mean annual AETI               : %.2f mm/year" % mean_aeti)
print("Mean annual NPP                : %.3f gC/m2/year" % mean_npp)
print("Mean annual TBP                : %.3f ton/ha/year" % mean_tbp)
print("Mean biomass WP                : %.3f kg/m3" % mean_wpb)
print("Peak / low monthly AETI        : %s %.2f / %s %.2f mm" %
      (months[peak_aeti_i], monthly_aeti[peak_aeti_i], months[low_aeti_i], monthly_aeti[low_aeti_i]))
print("Peak monthly NPP               : %s %.3f gC/m2" %
      (months[peak_npp_i], monthly_npp[peak_npp_i]))
print("Annual AETI-NPP spatial r      : %.4f" % spatial_r)
print("AETI 5th / median / 95th       : %.1f / %.1f / %.1f mm" %
      (aeti_p05, aeti_median, aeti_p95))
print("TBP 5th / median / 95th        : %.2f / %.2f / %.2f ton/ha" %
      (tbp_p05, tbp_median, tbp_p95))
print("Figure                         :", figure_fp)
print("Figure size                    : %.1f KB" % (os.path.getsize(figure_fp) / 1024))
print("Forum text                     :", forum_fp)
print("DONE")
