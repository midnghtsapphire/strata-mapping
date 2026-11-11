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
