"""Download 2025 monthly WaPOR AETI or NPP clipped to Benue State."""

import os
import sys

os.environ["GDAL_HTTP_UNSAFESSL"] = "YES"
os.environ.setdefault("GDAL_DATA", r"C:\Program Files\QGIS 3.44.15\share\gdal")
os.environ.setdefault("PROJ_DATA", r"C:\Program Files\QGIS 3.44.15\share\proj")
os.environ.setdefault("PROJ_LIB", r"C:\Program Files\QGIS 3.44.15\share\proj")

from osgeo import gdal

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_HTTP_UNSAFESSL", "YES")

from wapordl import wapor_map

HERE = os.path.dirname(os.path.abspath(__file__))
REGION = os.path.join(HERE, "Benue_State.geojson")
OUT = os.path.join(HERE, "wapor")
os.makedirs(OUT, exist_ok=True)

VARIABLE = sys.argv[1] if len(sys.argv) > 1 else "L1-AETI-M"
if VARIABLE not in {"L1-AETI-M", "L1-NPP-M"}:
    raise SystemExit("Use L1-AETI-M or L1-NPP-M")

short = "AETI" if "AETI" in VARIABLE else "NPP"
filename = f"Benue_WaPORv3_L1_{short}_monthly_2025"

print("Region  :", REGION)
print("Variable:", VARIABLE)
print("Period  : 2025-01-01 to 2025-12-31")
print("Output  :", OUT)

fp = wapor_map(
    REGION,
    VARIABLE,
    ["2025-01-01", "2025-12-31"],
    OUT,
    unit_conversion="none",
    overview="NONE",
    separate_unscale=False,
    filename=filename,
)

ds = gdal.Open(fp)
print("Created :", fp)
print("Bands   :", ds.RasterCount)
print("Size    :", ds.RasterXSize, "x", ds.RasterYSize)
print("Pixel   :", ds.GetGeoTransform()[1], "x", abs(ds.GetGeoTransform()[5]))
print("DONE")
