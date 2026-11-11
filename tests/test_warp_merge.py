"""SYNTHETIC mixed-CRS warp/merge tests. No network."""
from __future__ import annotations
import math
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import reproject, transform, transform_bounds
import lrm_3dep as m


def _open_local(monkeypatch, paths):
    real_open = rasterio.open

    def local_open(path, *args, **kwargs):
        if isinstance(path, str) and path.startswith("/vsicurl/"):
            path = path[len("/vsicurl/"):]
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(m.rasterio, "open", local_open)
    aligned, crs, handles = m.open_aligned_tiles([str(p) for p in paths])
    return aligned, crs, handles


def _write(path, epsg, west, south, east, north, fill, nodata=-9999.0, px=20.0, paint=None):
    minx, miny, maxx, maxy = transform_bounds(
        "EPSG:4326", f"EPSG:{epsg}", west, south, east, north, densify_pts=21
    )
    minx = math.floor(minx / px) * px
    miny = math.floor(miny / px) * px
    maxx = math.ceil(maxx / px) * px
    maxy = math.ceil(maxy / px) * px
    width = int(round((maxx - minx) / px))
    height = int(round((maxy - miny) / px))
    tr = from_origin(minx, maxy, px, px)
    arr = np.full((height, width), fill, dtype="float32")
    if paint is not None:
        paint(arr, tr, f"EPSG:{epsg}")
    profile = dict(
        driver="GTiff", height=height, width=width, count=1, dtype="float32",
        crs=f"EPSG:{epsg}", transform=tr,
    )
    if nodata is not None:
        profile["nodata"] = nodata
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr, 1)
    return tr


def _chip(arr, tr, crs, lon, lat, n=3):
    xs, ys = transform("EPSG:4326", crs, [lon], [lat])
    row, col = rasterio.transform.rowcol(tr, xs[0], ys[0])
    row, col = int(row), int(col)
    h = n // 2
    chip = arr[row - h:row + h + 1, col - h:col + h + 1]
    assert chip.shape == (n, n), (chip.shape, row, col, arr.shape)
    return chip


def test_warped_vrt_nearest_matches_independent_reproject(tmp_path, monkeypatch):
    """SYNTHETIC step edge in another UTM zone. Nearest warp keeps the two values."""
    target = tmp_path / "z12.tif"
    other = tmp_path / "z13.tif"
    _write(target, 26912, -108.03, 40.00, -108.005, 40.015, 100.0)
    def paint_step(arr, _tr, _crs):
        arr[:, : arr.shape[1] // 2] = 10.0
        arr[:, arr.shape[1] // 2 :] = 90.0
        arr[:2, :2] = -9999.0
    _write(other, 26913, -108.02, 40.00, -107.995, 40.015, 90.0, paint=paint_step)

    aligned, crs, handles = _open_local(monkeypatch, [target, other])
    try:
        assert crs.to_epsg() == 26912
        vrt = aligned[1]
        assert isinstance(vrt, WarpedVRT)
        assert vrt.resampling == Resampling.nearest
        assert vrt.nodata == -9999.0
        got = vrt.read(1)
        with rasterio.open(other) as src:
            dest = np.empty(got.shape, dtype="float32")
            reproject(
                source=rasterio.band(src, 1),
                destination=dest,
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=src.nodata,
                dst_transform=vrt.transform,
                dst_crs=vrt.crs,
                dst_nodata=vrt.nodata,
                resampling=Resampling.nearest,
            )
        np.testing.assert_array_equal(got, dest)
        vals = set(np.unique(got).tolist())
        assert vals <= {-9999.0, 10.0, 90.0}
        assert 10.0 in vals and 90.0 in vals
        assert 0.0 not in vals
    finally:
        for h in handles:
            h.close()


def test_mixed_crs_nodata_seam_and_boundless_fill(tmp_path, monkeypatch):
    """SYNTHETIC: hole on the west tile is filled from the east tile; outside stays nodata."""
    west = tmp_path / "west.tif"
    east = tmp_path / "east.tif"

    def hole(arr, tr, crs):
        rows, cols = np.indices(arr.shape)
        xs, ys = rasterio.transform.xy(tr, rows, cols, offset="center")
        lons, lats = transform(crs, "EPSG:4326", np.ravel(xs), np.ravel(ys))
        lons = np.asarray(lons).reshape(arr.shape)
        lats = np.asarray(lats).reshape(arr.shape)
        arr[(lons > -108.012) & (lons < -108.004) & (lats > 40.005) & (lats < 40.015)] = -9999.0

    _write(west, 26912, -108.04, 40.000, -108.000, 40.020, 111.0, paint=hole)
    _write(east, 26913, -108.016, 40.000, -107.980, 40.020, 222.0)
    aligned, crs, handles = _open_local(monkeypatch, [west, east])
    try:
        arr, tr, nd = m.merge_aligned(aligned)
        assert nd == -9999.0
        assert abs(tr.a - aligned[0].res[0]) < 1e-6
        assert abs(abs(tr.e) - aligned[0].res[1]) < 1e-6
        assert abs(tr.b) < 1e-9 and abs(tr.d) < 1e-9
        assert np.all(_chip(arr, tr, crs, -108.030, 40.010) == 111.0)
        assert np.all(_chip(arr, tr, crs, -108.008, 40.010) == 222.0)
        assert np.all(_chip(arr, tr, crs, -107.990, 40.010) == 222.0)
        assert int(np.sum(arr == 0)) == 0

        height, width = arr.shape
        left, top = tr.c, tr.f
        right = left + width * tr.a
        bottom = top + height * tr.e
        pad = 400.0
        wide, wtr, wnd = m.merge_aligned(
            aligned, bounds=(left - pad, bottom - pad, right + pad, top + pad)
        )
        assert wnd == nd
        assert wide.shape[0] > height and wide.shape[1] > width
        assert wide[0, 0] == nd and wide[-1, -1] == nd
        assert int(np.sum(wide == 0)) == 0
        # Interior samples survive the expanded window (alignment, not a shifted grid of zeros).
        assert np.all(_chip(wide, wtr, crs, -108.030, 40.010) == 111.0)
        assert np.all(_chip(wide, wtr, crs, -108.008, 40.010) == 222.0)
    finally:
        for h in handles:
            h.close()


def test_missing_nodata_tag_does_not_zero_fill_across_crs(tmp_path, monkeypatch):
    """SYNTHETIC: source with no nodata tag must not become a 0 halo after WarpedVRT."""
    a = tmp_path / "tagged.tif"
    b = tmp_path / "untagged.tif"
    _write(a, 26912, -108.025, 40.000, -108.005, 40.012, 100.0, nodata=-9999.0)
    _write(b, 26913, -108.012, 40.000, -107.992, 40.012, 50.0, nodata=None)
    aligned, _crs, handles = _open_local(monkeypatch, [a, b])
    try:
        assert aligned[1].nodata == m.NODATA_OUT
        got = aligned[1].read(1)
        assert set(np.unique(got).tolist()) <= {m.NODATA_OUT, 50.0}
        assert 0.0 not in set(np.unique(got).tolist())
        arr, _tr, nd = m.merge_aligned(aligned)
        assert nd == -9999.0
        assert set(np.unique(arr).tolist()) <= {-9999.0, 50.0, 100.0}
        assert int(np.sum(arr == 0)) == 0
    finally:
        for h in handles:
            h.close()


def test_nodata_zero_stays_masked_across_crs(tmp_path, monkeypatch):
    """SYNTHETIC: nodata 0 on the warped tile must not be rewritten as a valid elevation."""
    a = tmp_path / "meters.tif"
    b = tmp_path / "zero_nd.tif"
    _write(a, 26912, -108.025, 40.000, -108.005, 40.012, 100.0, nodata=-9999.0)

    def corner_zero(arr, _tr, _crs):
        arr[0, 0] = 0.0

    _write(b, 26913, -108.012, 40.000, -107.992, 40.012, 80.0, nodata=0.0, paint=corner_zero)
    aligned, _crs, handles = _open_local(monkeypatch, [a, b])
    try:
        assert aligned[1].nodata == 0.0
        warped = aligned[1].read(1)
        assert 0.0 in set(np.unique(warped).tolist())
        assert m.NODATA_OUT not in set(np.unique(warped).tolist())
        arr, _tr, nd = m.merge_aligned(aligned)
        assert nd == -9999.0
        # Output nodata follows the first tile. The warped tile's 0s are masked,
        # so they do not show up as a valid elevation of 0.
        assert int(np.sum(arr == 0)) == 0
        assert 80.0 in set(np.unique(arr).tolist())
        assert 100.0 in set(np.unique(arr).tolist())
    finally:
        for h in handles:
            h.close()


def test_mixed_crs_seam_nodata_tag_and_no_zero_fill(tmp_path, monkeypatch):
    """SYNTHETIC nodata and 0-fill guard, not a per-row continuity check.

    The east tile has no nodata tag and overlaps the west tile by about one pixel.
    The old WarpedVRT (no explicit nodata) leaves nodata unset and initialises
    uncovered dest pixels to 0, so those zeros land on the seam columns. This
    test proves the warp's nodata is NODATA_OUT and that the seam band has no 0.
    A few rows where 111 meets 222 only show the tiles touch. They do not claim
    every row is continuous.
    """
    px = 20.0
    lat0, lat1 = 40.00, 40.01
    # One 20 m pixel of longitude at this latitude.
    overlap_deg = px / (111320.0 * math.cos(math.radians(40.0)))
    west_path = tmp_path / "west.tif"
    east_path = tmp_path / "east.tif"
    _write(west_path, 26912, -108.01, lat0, -108.0, lat1, 111.0, nodata=-9999.0, px=px)
    _write(east_path, 26913, -108.0 - overlap_deg, lat0, -107.99, lat1, 222.0, nodata=None, px=px)
    with rasterio.open(east_path) as east_src:
        assert east_src.nodata is None

    aligned, crs, handles = _open_local(monkeypatch, [west_path, east_path])
    try:
        assert crs.to_epsg() == 26912
        assert isinstance(aligned[1], WarpedVRT)
        assert aligned[1].nodata == m.NODATA_OUT
        arr, _tr, nd = m.merge_aligned(aligned)
        assert nd == m.NODATA_OUT
        c111 = np.flatnonzero(np.any(arr == 111.0, axis=0))
        c222 = np.flatnonzero(np.any(arr == 222.0, axis=0))
        assert c111.size and c222.size
        # About one pixel of overlap or a shared edge: the values meet.
        assert abs(int(c222.min()) - int(c111.max())) <= 2
        seam = arr[:, int(min(c111.max(), c222.min())):int(c222.max()) + 1]
        assert not np.any(seam == 0)
        met = 0
        for row in arr:
            a = np.flatnonzero(row == 111.0)
            b = np.flatnonzero(row == 222.0)
            if a.size and b.size and abs(int(b.min()) - int(a.max())) == 1:
                met += 1
        assert met > 10
    finally:
        for h in handles:
            h.close()


def test_nearest_resample_slope_error_is_bounded(tmp_path, monkeypatch):
    """SYNTHETIC plane sampled in two UTM zones.

    z = 1500 + 0.05*(x-x0) + 0.02*(y-y0) metres in EPSG:26912.
    Nearest can land anywhere in the half-pixel box, so one stage errs by at
    most |0.05|*10 + |0.02|*10 = 0.7 m on a 20 m pixel. The value is resampled
    twice (WarpedVRT nearest, then merge nearest): about 1.4 m. The cap is
    that worst case plus 0.1 m, i.e. 1.5 m. The real guard against bilinear
    is `assert aligned[1].resampling == Resampling.nearest`. The 0.1 m mean
    lower bound is only a weak secondary check.
    """
    x0, y0 = 756000.0, 4433200.0
    px = 20.0
    height, width = 8, 8
    seam = x0 + width * px

    def z_of(x, y):
        return 1500.0 + 0.05 * (x - x0) + 0.02 * (y - y0)

    tr_w = from_origin(x0, y0, px, px)
    west = np.empty((height, width), np.float32)
    for r in range(height):
        for c in range(width):
            west[r, c] = z_of(x0 + (c + 0.5) * px, y0 - (r + 0.5) * px)

    # East grid lives in the neighbouring zone. Each pixel is the analytic
    # value at its own centre, so the only resample is WarpedVRT nearest.
    corners_x = [seam, seam + width * px, seam, seam + width * px]
    corners_y = [y0, y0, y0 - height * px, y0 - height * px]
    lons, lats = transform("EPSG:26912", "EPSG:4326", corners_x, corners_y)
    ex, ey = transform("EPSG:4326", "EPSG:26913", lons, lats)
    minx = math.floor((min(ex) - px) / px) * px
    maxx = math.ceil((max(ex) + px) / px) * px
    miny = math.floor((min(ey) - px) / px) * px
    maxy = math.ceil((max(ey) + px) / px) * px
    width_e = int(round((maxx - minx) / px))
    height_e = int(round((maxy - miny) / px))
    tr_e = from_origin(minx, maxy, px, px)
    east = np.empty((height_e, width_e), np.float32)
    for r in range(height_e):
        for c in range(width_e):
            x13 = minx + (c + 0.5) * px
            y13 = maxy - (r + 0.5) * px
            lon, lat = transform("EPSG:26913", "EPSG:4326", [x13], [y13])
            x12, y12 = transform("EPSG:4326", "EPSG:26912", lon, lat)
            east[r, c] = np.float32(z_of(x12[0], y12[0]))

    west_path = tmp_path / "west.tif"
    east_path = tmp_path / "east.tif"
    for path, arr, tr, epsg in (
        (west_path, west, tr_w, "EPSG:26912"),
        (east_path, east, tr_e, "EPSG:26913"),
    ):
        with rasterio.open(
            path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
            dtype="float32", crs=epsg, transform=tr, nodata=-9999.0,
        ) as dst:
            dst.write(arr, 1)

    aligned, crs, handles = _open_local(monkeypatch, [west_path, east_path])
    try:
        assert isinstance(aligned[1], WarpedVRT)
        assert aligned[1].resampling == Resampling.nearest
        arr, tr, nd = m.merge_aligned(aligned)
        errs = []
        east_errs = []
        for r in range(arr.shape[0]):
            for c in range(arr.shape[1]):
                if arr[r, c] == nd:
                    continue
                x = tr.c + (c + 0.5) * tr.a
                y = tr.f + (r + 0.5) * tr.e
                err = abs(float(arr[r, c]) - z_of(x, y))
                errs.append(err)
                if x >= seam:
                    east_errs.append(err)
        errs = np.asarray(errs)
        east_errs = np.asarray(east_errs)
        assert errs.size > 50 and east_errs.size > 20
        # Derivation of the published cap, not the same expression as limit_m:
        # one stage is |0.05|*(20/2) + |0.02|*(20/2) = 0.7 m, two nearest
        # stages are 1.4 m, plus 0.1 m is the literal 1.5 m.
        per_stage_m = abs(0.05) * (px / 2.0) + abs(0.02) * (px / 2.0)
        limit_m = 2.0 * per_stage_m + 0.1
        assert per_stage_m == 0.7
        assert limit_m == 1.5
        assert float(east_errs.max()) < limit_m
        assert float(errs.max()) < limit_m
        assert float(east_errs.mean()) > 0.1
        assert float(east_errs.mean()) < 0.8
    finally:
        for h in handles:
            h.close()


def test_abutting_tiles_align_and_nodata_edge_is_filled(tmp_path):
    """SYNTHETIC same-grid seam: no gap between abutting tiles; overlap nodata is filled."""
    left_path = tmp_path / "left.tif"
    right_path = tmp_path / "right.tif"
    tr_l = from_origin(500000, 4400010, 1.0, 1.0)
    tr_r = from_origin(500010, 4400010, 1.0, 1.0)
    left = np.full((10, 10), 10.0, dtype="float32")
    right = np.full((10, 10), 20.0, dtype="float32")
    for path, arr, tr in ((left_path, left, tr_l), (right_path, right, tr_r)):
        with rasterio.open(
            path, "w", driver="GTiff", height=10, width=10, count=1, dtype="float32",
            crs="EPSG:26913", transform=tr, nodata=-9999.0,
        ) as dst:
            dst.write(arr, 1)
    srcs = [rasterio.open(left_path), rasterio.open(right_path)]
    try:
        arr, tr, nd = m.merge_aligned(srcs)
        assert nd == -9999.0
        assert arr.shape == (10, 20)
        assert tr.c == 500000 and tr.a == 1.0 and tr.e == -1.0
        assert np.all(arr[:, :10] == 10.0)
        assert np.all(arr[:, 10:] == 20.0)
        assert int(np.sum(arr == nd)) == 0
    finally:
        for s in srcs:
            s.close()

    # One-pixel overlap: left's east column is nodata and right covers it.
    tr_r2 = from_origin(500009, 4400010, 1.0, 1.0)
    left[:, -1] = -9999.0
    overlap = tmp_path / "overlap.tif"
    with rasterio.open(
        left_path, "w", driver="GTiff", height=10, width=10, count=1, dtype="float32",
        crs="EPSG:26913", transform=tr_l, nodata=-9999.0,
    ) as dst:
        dst.write(left, 1)
    with rasterio.open(
        overlap, "w", driver="GTiff", height=10, width=10, count=1, dtype="float32",
        crs="EPSG:26913", transform=tr_r2, nodata=-9999.0,
    ) as dst:
        dst.write(right, 1)
    srcs = [rasterio.open(left_path), rasterio.open(overlap)]
    try:
        arr, _tr, nd = m.merge_aligned(srcs)
        assert arr.shape == (10, 19)
        assert np.all(arr[:, :9] == 10.0)
        assert np.all(arr[:, 9:] == 20.0)
        assert int(np.sum(arr == nd)) == 0
    finally:
        for s in srcs:
            s.close()
