#!/usr/bin/env python3
"""lrm_3dep.py - Local Relief Model (LRM) from USGS 3DEP 1 m bare-earth DTM.

LRM = DTM - GaussianSmooth(DTM, sigma_m). Deterministic: no randomness; tile
list is sorted before mosaicking.

Sources (in order):
  1. USGS TNM Access API -> public 3DEP 1 m COGs (prd-tnm S3 / National Map),
     read by window over HTTP (/vsicurl/).
  2. OpenTopography USGS 1m API (needs env OPENTOPO_API_KEY).
  Or --input LOCAL.tif to skip the network.
"""
from __future__ import annotations
import argparse, json, os, sys, tempfile, urllib.parse, urllib.request
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from scipy.ndimage import gaussian_filter

TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
OT_API = "https://portal.opentopography.org/API/usgsdem"
NODATA_OUT = -9999.0


def local_relief(dtm: np.ndarray, nodata, sigma_px: float) -> np.ndarray:
    """Nodata-aware Gaussian LRM (normalized convolution). Returns float32."""
    a = dtm.astype("float64")
    valid = np.isfinite(a)
    if nodata is not None:
        valid &= a != nodata
    filled = np.where(valid, a, 0.0)
    w = gaussian_filter(valid.astype("float64"), sigma_px, mode="nearest")
    s = gaussian_filter(filled, sigma_px, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        smooth = s / w
    out = np.where(valid & (w > 1e-6), a - smooth, NODATA_OUT)
    return out.astype("float32")


def tnm_cog_urls(bbox):
    q = urllib.parse.urlencode({"datasets": "Digital Elevation Model (DEM) 1 meter",
                                "bbox": ",".join(map(str, bbox)),
                                "prodFormats": "GeoTIFF", "max": 50, "outputFormat": "JSON"})
    with urllib.request.urlopen(f"{TNM_API}?{q}", timeout=60) as r:
        items = json.load(r).get("items", [])
    return sorted({i["downloadURL"] for i in items if i.get("downloadURL", "").endswith(".tif")})


def read_tnm(bbox):
    urls = tnm_cog_urls(bbox)
    if not urls:
        raise RuntimeError("TNM returned no 1 m tiles for bbox")
    srcs = [rasterio.open("/vsicurl/" + u) for u in urls]
    crs = srcs[0].crs
    srcs = [s for s in srcs if s.crs == crs]  # 1 m tiles of one UTM zone
    b = transform_bounds("EPSG:4326", crs, *bbox, densify_pts=21)
    arr, tr = merge(srcs, bounds=b, nodata=srcs[0].nodata if srcs[0].nodata is not None else NODATA_OUT)
    nd = srcs[0].nodata if srcs[0].nodata is not None else NODATA_OUT
    for s in srcs: s.close()
    return arr[0], tr, crs, nd, urls


def read_opentopo(bbox):
    key = os.environ.get("OPENTOPO_API_KEY")
    if not key:
        raise RuntimeError("OPENTOPO_API_KEY not set")
    w, s, e, n = bbox
    q = urllib.parse.urlencode({"datasetName": "USGS1m", "south": s, "north": n, "west": w,
                                "east": e, "outputFormat": "GTiff", "API_Key": key})
    tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False).name
    urllib.request.urlretrieve(f"{OT_API}?{q}", tmp)
    return read_local(tmp, None) + ([OT_API],)


def read_local(path, bbox):
    with rasterio.open(path) as src:
        if bbox:
            b = transform_bounds("EPSG:4326", src.crs, *bbox, densify_pts=21)
            win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
            arr = src.read(1, window=win, boundless=True, fill_value=src.nodata or NODATA_OUT)
            tr = src.window_transform(win)
        else:
            arr, tr = src.read(1), src.transform
        return arr, tr, src.crs, src.nodata


def write_lrm(path, lrm, transform, crs):
    prof = dict(driver="GTiff", height=lrm.shape[0], width=lrm.shape[1], count=1,
                dtype="float32", crs=crs, transform=transform, nodata=NODATA_OUT,
                compress="lzw", predictor=3, tiled=True, blockxsize=256, blockysize=256)
    if lrm.shape[0] < 256 or lrm.shape[1] < 256:
        prof.update(tiled=False); prof.pop("blockxsize"); prof.pop("blockysize")
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(lrm, 1)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bbox", nargs=4, type=float, metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"))
    p.add_argument("--sigma-m", type=float, default=10.0, help="Gaussian sigma in meters (default 10)")
    p.add_argument("--source", choices=["auto", "tnm", "opentopo"], default="auto")
    p.add_argument("--input", help="local DTM GeoTIFF (skips network)")
    p.add_argument("--save-dtm", help="also write the cropped DTM here")
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args(argv)
    if not a.input and not a.bbox:
        p.error("--bbox or --input required")
    if a.input:
        arr, tr, crs, nd = read_local(a.input, a.bbox); src = a.input
    else:
        errs = []
        for name, fn in (("tnm", read_tnm), ("opentopo", read_opentopo)):
            if a.source not in ("auto", name):
                continue
            try:
                arr, tr, crs, nd, srcs = fn(a.bbox); src = f"{name}: {len(srcs)} file(s)"; break
            except Exception as ex:  # noqa: BLE001
                errs.append(f"{name}: {ex}")
        else:
            sys.exit("All sources failed:\n  " + "\n  ".join(errs))
    if not crs or not crs.is_projected:
        sys.exit(f"DTM CRS {crs} is not projected; sigma in meters needs a metric CRS")
    px = abs(tr.a)
    lrm = local_relief(arr, nd, a.sigma_m / px)
    write_lrm(a.output, lrm, tr, crs)
    if a.save_dtm:
        write_lrm(a.save_dtm, np.where(lrm == NODATA_OUT, NODATA_OUT, arr).astype("float32"), tr, crs)
    v = lrm[lrm != NODATA_OUT]
    print(f"source={src} crs={crs} px={px}m shape={lrm.shape} sigma_px={a.sigma_m/px:.2f} "
          f"valid={v.size} lrm_min={v.min() if v.size else 'na'} lrm_max={v.max() if v.size else 'na'} -> {a.output}")


if __name__ == "__main__":
    main()
