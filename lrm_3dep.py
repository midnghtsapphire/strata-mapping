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
import argparse, json, math, os, sys, tempfile, urllib.parse, urllib.request
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.merge import merge
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds
from scipy.ndimage import gaussian_filter

TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
OT_API = "https://portal.opentopography.org/API/usgsdem"
NODATA_OUT = -9999.0
TNM_PAGE_SIZE = 50
OT_TIMEOUT_S = 60
HALO_SIGMA_MULT = 4.0


def _gaussian_sigma(sigma_px):
    """Scalar, or (sigma_row, sigma_col) for scipy.ndimage (y, then x).

    A float still means the same sigma on every axis so existing callers keep working.
    """
    if isinstance(sigma_px, (str, bytes)):
        raise ValueError(
            f"sigma_px must be a positive float or (sigma_row, sigma_col), got {sigma_px!r}"
        )
    arr = np.asarray(sigma_px, dtype="float64")
    if arr.ndim == 0:
        values = (float(arr),)
        sigma_arg = values[0]
    elif arr.ndim == 1 and arr.shape == (2,):
        values = (float(arr[0]), float(arr[1]))
        sigma_arg = values
    else:
        raise ValueError(
            f"sigma_px must be a positive float or (sigma_row, sigma_col), got {sigma_px!r}"
        )
    for s in values:
        if not math.isfinite(s) or s <= 0:
            raise ValueError(f"sigma_px must be finite and > 0, got {sigma_px!r}")
    return sigma_arg


def local_relief(dtm: np.ndarray, nodata, sigma_px) -> np.ndarray:
    """Nodata-aware Gaussian LRM (normalized convolution). Returns float32.

    `sigma_px` is a float (isotropic) or `(sigma_row, sigma_col)` in pixels.
    """
    sigma_arg = _gaussian_sigma(sigma_px)
    a = dtm.astype("float64")
    valid = np.isfinite(a)
    if nodata is not None:
        valid &= a != nodata
    filled = np.where(valid, a, 0.0)
    w = gaussian_filter(valid.astype("float64"), sigma_arg, mode="nearest")
    s = gaussian_filter(filled, sigma_arg, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        smooth = s / w
    out = np.where(valid & (w > 1e-6), a - smooth, NODATA_OUT)
    return out.astype("float32")


def effective_nodata(src_nodata) -> float:
    """Single nodata used for fill AND masking (handles None and 0 correctly)."""
    return float(src_nodata) if src_nodata is not None else NODATA_OUT


def meters_per_crs_unit(crs) -> float:
    """Factor: 1 CRS linear unit * factor = meters. Hard-fail geographic / unknown."""
    if crs is None or not crs.is_projected:
        raise ValueError(
            f"DTM CRS {crs} is not projected; sigma in meters needs a projected CRS "
            "with linear units (not geographic lon/lat)"
        )
    try:
        _name, to_m = crs.linear_units_factor
    except Exception as ex:  # noqa: BLE001
        raise ValueError(f"DTM CRS {crs} has no linear units factor: {ex}") from ex
    if to_m is None or not math.isfinite(to_m) or to_m <= 0:
        raise ValueError(f"DTM CRS {crs} linear units factor is invalid: {to_m!r}")
    return float(to_m)


def _reject_rotated_or_sheared(transform) -> None:
    """|a| and |e| are not the ground pixel size once the axes are rotated or sheared."""
    b = float(transform.b)
    d = float(transform.d)
    scale = max(abs(float(transform.a)), abs(float(transform.e)), 1.0)
    tol = 1e-9 * scale
    if abs(b) > tol or abs(d) > tol:
        raise ValueError(
            "rotated or sheared affine is not supported for Gaussian sigma "
            f"(b={b}, d={d}; axis-aligned tolerance {tol:g}). "
            "Refusing a silent wrong radius; a rotated-pixel kernel is not implemented."
        )


def pixel_size_meters(crs, transform) -> tuple[float, float]:
    """Per-axis pixel size in meters from the affine, not mean(|a|, |e|).

    Returns `(px_x, px_y)` = `(|a|, |e|)` times the CRS linear-unit factor
    (column, then row). North-up rasters have `e < 0`; the sign is dropped.
    Raises ValueError when `b` or `d` is non-zero (rotation or shear), or when
    either axis is not a finite positive length.
    """
    _reject_rotated_or_sheared(transform)
    to_m = meters_per_crs_unit(crs)
    px_x = abs(float(transform.a)) * to_m
    px_y = abs(float(transform.e)) * to_m
    if px_x <= 0 or px_y <= 0 or not math.isfinite(px_x) or not math.isfinite(px_y):
        raise ValueError(f"invalid pixel size from transform {transform}")
    return (px_x, px_y)


def gaussian_sigma_px(sigma_m: float, px_x: float, px_y: float) -> tuple[float, float]:
    """Gaussian sigma in pixels as `(sigma_row, sigma_col)` for scipy.ndimage.

    Row is northing (`sigma_m / px_y`); column is easting (`sigma_m / px_x`).
    """
    if not math.isfinite(sigma_m) or sigma_m <= 0:
        raise ValueError(f"sigma_m must be finite and > 0, got {sigma_m!r}")
    if min(px_x, px_y) <= 0 or not math.isfinite(px_x) or not math.isfinite(px_y):
        raise ValueError(f"pixel sizes must be finite and > 0, got px_x={px_x!r} px_y={px_y!r}")
    return (sigma_m / px_y, sigma_m / px_x)


def halo_pixels(sigma_m: float, px_x: float, px_y: float) -> tuple[int, int]:
    """Kernel halo in `(rows, cols)`, about 4 sigma on each axis independently."""
    sigma_row, sigma_col = gaussian_sigma_px(sigma_m, px_x, px_y)
    return (
        int(math.ceil(HALO_SIGMA_MULT * sigma_row)),
        int(math.ceil(HALO_SIGMA_MULT * sigma_col)),
    )


def tnm_fetch_page(bbox, offset: int, max_n: int = TNM_PAGE_SIZE) -> dict:
    q = urllib.parse.urlencode({
        "datasets": "Digital Elevation Model (DEM) 1 meter",
        "bbox": ",".join(map(str, bbox)),
        "prodFormats": "GeoTIFF",
        "max": max_n,
        "offset": offset,
        "outputFormat": "JSON",
    })
    with urllib.request.urlopen(f"{TNM_API}?{q}", timeout=60) as r:
        return json.load(r)


def tnm_cog_urls(bbox, page_fetcher=tnm_fetch_page) -> list[str]:
    """Paginate TNM products until all items fetched. Hard-fail on silent truncation."""
    items: list = []
    offset = 0
    total = None
    while True:
        payload = page_fetcher(bbox, offset, TNM_PAGE_SIZE)
        page = payload.get("items") or []
        if total is None:
            total = payload.get("total")
            if total is None:
                total = payload.get("totalCount")
        items.extend(page)
        if total is not None:
            try:
                total_i = int(total)
            except (TypeError, ValueError):
                total_i = None
            if total_i is not None:
                if len(items) >= total_i:
                    break
                if not page:
                    raise RuntimeError(
                        f"TNM pagination stalled: got {len(items)} items but total={total_i}"
                    )
                offset = len(items)
                continue
        if len(page) < TNM_PAGE_SIZE:
            break
        offset = len(items)
        if offset > 5000:
            raise RuntimeError(
                f"TNM returned >= {offset} items without a total; refusing unbounded pagination"
            )
    if total is not None:
        try:
            total_i = int(total)
        except (TypeError, ValueError):
            total_i = None
        if total_i is not None and len(items) < total_i:
            raise RuntimeError(
                f"TNM results truncated: fetched {len(items)} of total={total_i}"
            )
    return sorted({
        i["downloadURL"] for i in items
        if str(i.get("downloadURL", "")).endswith(".tif")
    })


def open_aligned_tiles(urls: list[str]):
    """Open COGs; reproject CRS mismatches into the first tile's CRS via WarpedVRT.

    Never silently drops tiles. Returns (aligned_datasets, target_crs, close_handles).
    """
    if not urls:
        raise RuntimeError("TNM returned no 1 m tiles for bbox")
    raw = [rasterio.open("/vsicurl/" + u) for u in urls]
    target_crs = raw[0].crs
    if target_crs is None:
        for s in raw:
            s.close()
        raise RuntimeError("first TNM tile has no CRS")
    aligned = []
    wrappers = []
    try:
        for s in raw:
            if s.crs == target_crs:
                aligned.append(s)
            else:
                # GDAL's default destination nodata is 0 when the source tag is
                # missing. That writes a 0 seam around the warped tile. Use the
                # same effective nodata as read_local, and nearest resampling so
                # nodata does not bleed into valid elevations.
                vrt = WarpedVRT(
                    s,
                    crs=target_crs,
                    src_nodata=s.nodata,
                    nodata=effective_nodata(s.nodata),
                    resampling=Resampling.nearest,
                )
                wrappers.append(vrt)
                aligned.append(vrt)
        return aligned, target_crs, raw + wrappers
    except Exception:
        for s in raw:
            s.close()
        for v in wrappers:
            v.close()
        raise


def merge_aligned(srcs, bounds=None):
    """Merge datasets that already share a CRS.

    Output nodata is `effective_nodata` of the first dataset. Pixels outside
    every tile (a bounds window larger than the mosaic, the boundless case)
    stay that nodata value. Each source is masked with its own nodata, so a
    nodata sample on one tile does not punch a hole where a later tile has data.
    Resolution and alignment follow the first dataset. Nearest resampling keeps
    nodata from blending into neighboring elevations.
    """
    if not srcs:
        raise RuntimeError("no tiles to merge")
    nd = effective_nodata(srcs[0].nodata)
    kwargs = {"nodata": nd, "resampling": Resampling.nearest, "method": "first"}
    if bounds is not None:
        kwargs["bounds"] = bounds
    arr, tr = merge(srcs, **kwargs)
    return arr[0], tr, nd


def read_tnm(bbox, sigma_m: float | None = None):
    urls = tnm_cog_urls(bbox)
    srcs, crs, handles = open_aligned_tiles(urls)
    try:
        b = transform_bounds("EPSG:4326", crs, *bbox, densify_pts=21)
        if sigma_m is not None and sigma_m > 0:
            to_m = meters_per_crs_unit(crs)
            pad = HALO_SIGMA_MULT * sigma_m / to_m
            b = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
        arr, tr, nd = merge_aligned(srcs, bounds=b)
        # Crop halo after merge so LRM edges use real terrain
        core = transform_bounds("EPSG:4326", crs, *bbox, densify_pts=21)
        return arr, tr, crs, nd, urls, core
    finally:
        for s in handles:
            s.close()


def _download_to(path: str, url: str, timeout: float = OT_TIMEOUT_S) -> None:
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(path, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def read_opentopo(bbox, sigma_m: float | None = None):
    key = os.environ.get("OPENTOPO_API_KEY")
    if not key:
        raise RuntimeError("OPENTOPO_API_KEY not set")
    w, s, e, n = bbox
    req_bbox = (w, s, e, n)
    if sigma_m is not None and sigma_m > 0:
        dlat = (HALO_SIGMA_MULT * sigma_m) / 111320.0
        dlon = dlat / max(math.cos(math.radians((s + n) / 2.0)), 0.1)
        req_bbox = (w - dlon, s - dlat, e + dlon, n + dlat)
    q = urllib.parse.urlencode({
        "datasetName": "USGS1m",
        "south": req_bbox[1], "north": req_bbox[3],
        "west": req_bbox[0], "east": req_bbox[2],
        "outputFormat": "GTiff", "API_Key": key,
    })
    tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False).name
    try:
        _download_to(tmp, f"{OT_API}?{q}", timeout=OT_TIMEOUT_S)
        arr, tr, crs, nd = read_local(tmp, None)
        core = transform_bounds("EPSG:4326", crs, *bbox, densify_pts=21)
        return arr, tr, crs, nd, [OT_API], core
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def _halo_rc(halo_px) -> tuple[int, int]:
    """Accept an int (both axes) or `(rows, cols)` so square callers stay valid."""
    if isinstance(halo_px, (int, np.integer)):
        h = int(halo_px)
        if h < 0:
            raise ValueError(f"halo_px must be >= 0, got {halo_px!r}")
        return h, h
    if isinstance(halo_px, (tuple, list)) and len(halo_px) == 2:
        hy, hx = int(halo_px[0]), int(halo_px[1])
        if hy < 0 or hx < 0:
            raise ValueError(f"halo_px must be >= 0, got {halo_px!r}")
        return hy, hx
    raise ValueError(f"halo_px must be an int or (rows, cols), got {halo_px!r}")


def read_local(path, bbox, halo_px: int | tuple[int, int] = 0):
    """Read local GeoTIFF. Boundless fill value == returned nodata (never mismatched).

    `halo_px` is an int, or `(rows, cols)` when the pixel axes differ.
    """
    with rasterio.open(path) as src:
        nd = effective_nodata(src.nodata)
        if bbox:
            b = transform_bounds("EPSG:4326", src.crs, *bbox, densify_pts=21)
            win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
            hy, hx = _halo_rc(halo_px)
            if hy or hx:
                win = Window(
                    win.col_off - hx,
                    win.row_off - hy,
                    win.width + 2 * hx,
                    win.height + 2 * hy,
                )
            arr = src.read(1, window=win, boundless=True, fill_value=nd)
            tr = src.window_transform(win)
        else:
            arr, tr = src.read(1), src.transform
        return arr, tr, src.crs, nd


def crop_to_bounds(arr, transform, bounds):
    """Crop array to CRS bounds; return (cropped, new_transform)."""
    win = from_bounds(*bounds, transform=transform).round_offsets().round_lengths()
    # Clip window to array
    row_off = max(0, int(win.row_off))
    col_off = max(0, int(win.col_off))
    row_end = min(arr.shape[0], int(win.row_off + win.height))
    col_end = min(arr.shape[1], int(win.col_off + win.width))
    if row_end <= row_off or col_end <= col_off:
        raise ValueError(f"crop window empty for bounds={bounds}")
    cropped = arr[row_off:row_end, col_off:col_end]
    tr = rasterio.windows.transform(
        Window(col_off, row_off, cropped.shape[1], cropped.shape[0]),
        transform,
    )
    return cropped, tr


def crop_center(arr, transform, halo_px):
    """Drop a halo that is an int or `(rows, cols)`."""
    hy, hx = _halo_rc(halo_px)
    if hy == 0 and hx == 0:
        return arr, transform
    if arr.shape[0] <= 2 * hy or arr.shape[1] <= 2 * hx:
        raise ValueError(f"array {arr.shape} too small to crop halo_px={halo_px!r}")
    cropped = arr[hy:arr.shape[0] - hy, hx:arr.shape[1] - hx]
    tr = rasterio.windows.transform(
        Window(hx, hy, cropped.shape[1], cropped.shape[0]),
        transform,
    )
    return cropped, tr


def write_lrm(path, lrm, transform, crs):
    prof = dict(
        driver="GTiff", height=lrm.shape[0], width=lrm.shape[1], count=1,
        dtype="float32", crs=crs, transform=transform, nodata=NODATA_OUT,
        compress="lzw", predictor=3, tiled=True, blockxsize=256, blockysize=256,
    )
    if lrm.shape[0] < 256 or lrm.shape[1] < 256:
        prof.update(tiled=False)
        prof.pop("blockxsize")
        prof.pop("blockysize")
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
    if not math.isfinite(a.sigma_m) or a.sigma_m <= 0:
        p.error(f"--sigma-m must be finite and > 0, got {a.sigma_m!r}")

    core_bounds = None
    halo_rc = (0, 0)
    if a.input:
        with rasterio.open(a.input) as src:
            peek_crs, peek_tr = src.crs, src.transform
        try:
            px_x_peek, px_y_peek = pixel_size_meters(peek_crs, peek_tr)
            halo_rc = halo_pixels(a.sigma_m, px_x_peek, px_y_peek) if a.bbox else (0, 0)
        except ValueError as ex:
            sys.exit(str(ex))
        arr, tr, crs, nd = read_local(a.input, a.bbox, halo_px=halo_rc)
        src_label = a.input
    else:
        errs = []
        arr = tr = crs = nd = None
        src_label = None
        for name, fn in (("tnm", read_tnm), ("opentopo", read_opentopo)):
            if a.source not in ("auto", name):
                continue
            try:
                arr, tr, crs, nd, srcs, core_bounds = fn(a.bbox, sigma_m=a.sigma_m)
                src_label = f"{name}: {len(srcs)} file(s)"
                break
            except Exception as ex:  # noqa: BLE001
                errs.append(f"{name}: {ex}")
        else:
            sys.exit("All sources failed:\n  " + "\n  ".join(errs))

    # Keep padded copies for halo crop of both LRM and optional DTM
    arr_pad, tr_pad = arr, tr

    try:
        px_x, px_y = pixel_size_meters(crs, tr_pad)
        sigma_px = gaussian_sigma_px(a.sigma_m, px_x, px_y)
    except ValueError as ex:
        sys.exit(str(ex))
    lrm_pad = local_relief(arr_pad, nd, sigma_px)

    if a.input and a.bbox and halo_rc != (0, 0):
        lrm, tr = crop_center(lrm_pad, tr_pad, halo_rc)
        arr_core, _ = crop_center(arr_pad, tr_pad, halo_rc)
    elif core_bounds is not None:
        lrm, tr = crop_to_bounds(lrm_pad, tr_pad, core_bounds)
        arr_core, _ = crop_to_bounds(arr_pad, tr_pad, core_bounds)
    else:
        lrm, tr, arr_core = lrm_pad, tr_pad, arr_pad

    write_lrm(a.output, lrm, tr, crs)
    if a.save_dtm:
        dtm_out = np.where(
            (arr_core == nd) | ~np.isfinite(arr_core), NODATA_OUT, arr_core
        ).astype("float32")
        write_lrm(a.save_dtm, dtm_out, tr, crs)

    v = lrm[lrm != NODATA_OUT]
    sy, sx = sigma_px
    print(
        f"source={src_label} crs={crs} px_x={px_x}m px_y={px_y}m shape={lrm.shape} "
        f"sigma_px=({sy:.2f},{sx:.2f}) "
        f"valid={v.size} lrm_min={v.min() if v.size else 'na'} lrm_max={v.max() if v.size else 'na'} "
        f"-> {a.output}"
    )


if __name__ == "__main__":
    main()
