"""Tests use SYNTHETIC TEST FIXTURES (generated rasters), NOT measured elevation data."""
from __future__ import annotations
import math
import numpy as np
import pytest
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
from rasterio.warp import transform_bounds
import lrm_3dep as m


def make_fixture(path, nodata=True, crs="EPSG:26913", px=1.0, nodata_tag=-9999.0, write_nodata_tag=True):
    """SYNTHETIC: tilted plane + 2 m-deep 5x5 pit."""
    y, x = np.mgrid[0:100, 0:100].astype("float32")
    dtm = 1500 + 0.1 * x + 0.05 * y
    dtm[48:53, 48:53] -= 2.0
    if nodata:
        dtm[0:5, 0:5] = -9999.0 if nodata_tag is None else nodata_tag
    tr = from_origin(500000, 4400000, px, px)
    profile = dict(
        driver="GTiff", height=100, width=100, count=1, dtype="float32",
        crs=crs, transform=tr,
    )
    if write_nodata_tag and nodata_tag is not None:
        profile["nodata"] = nodata_tag
    with rasterio.open(path, "w", **profile) as d:
        d.write(dtm, 1)
    return dtm, tr


def test_plane_gives_zero_lrm():
    y, x = np.mgrid[0:60, 0:60].astype("float64")
    lrm = m.local_relief(3 * x + 2 * y, None, 5.0)
    assert np.abs(lrm[21:39, 21:39]).max() < 1e-4


def test_pit_negative_and_formula(tmp_path):
    from scipy.ndimage import gaussian_filter
    dtm, _ = make_fixture(tmp_path / "s.tif", nodata=False)
    lrm = m.local_relief(dtm, None, 4.0)
    ref = dtm - gaussian_filter(dtm.astype("float64"), 4.0, mode="nearest")
    np.testing.assert_allclose(lrm, ref, atol=1e-3)
    assert lrm[50, 50] < -1.0 and abs(lrm[20, 20]) < 1e-3


def test_cli_georef_nodata_compression_determinism(tmp_path):
    make_fixture(tmp_path / "s.tif")
    for o in ("a.tif", "b.tif"):
        m.main(["--input", str(tmp_path / "s.tif"), "--sigma-m", "4", "-o", str(tmp_path / o)])
    with rasterio.open(tmp_path / "a.tif") as a, rasterio.open(tmp_path / "s.tif") as s, rasterio.open(tmp_path / "b.tif") as b:
        assert a.crs == s.crs and a.transform == s.transform
        assert a.dtypes[0] == "float32" and a.nodata == -9999.0
        assert a.profile["compress"].lower() == "lzw"
        A = a.read(1)
        assert (A[0:5, 0:5] == -9999.0).all() and np.isfinite(A).all()
        assert np.array_equal(A, b.read(1))
        assert A[50, 50] < -1.0
        # Valid pit center must not be the nodata sentinel (guards all-nodata regression)
        assert A[50, 50] != -9999.0


def test_sigma_meters_scales_with_pixel_size(tmp_path):
    """Equivalent --sigma-m on 1 m vs 2 m pixels must produce matching sigma_px effect."""
    make_fixture(tmp_path / "p1.tif", nodata=False, px=1.0)
    make_fixture(tmp_path / "p2.tif", nodata=False, px=2.0)
    m.main(["--input", str(tmp_path / "p1.tif"), "--sigma-m", "4", "-o", str(tmp_path / "o1.tif")])
    m.main(["--input", str(tmp_path / "p2.tif"), "--sigma-m", "8", "-o", str(tmp_path / "o2.tif")])
    with rasterio.open(tmp_path / "o1.tif") as a, rasterio.open(tmp_path / "o2.tif") as b:
        # Both use sigma_px = 4; pit center should both be strongly negative
        assert a.read(1)[50, 50] < -1.0
        assert b.read(1)[50, 50] < -1.0
        # Distinct from a half-sigma run on 1 m pixels
    m.main(["--input", str(tmp_path / "p1.tif"), "--sigma-m", "2", "-o", str(tmp_path / "o1b.tif")])
    with rasterio.open(tmp_path / "o1.tif") as a, rasterio.open(tmp_path / "o1b.tif") as c:
        assert a.read(1)[50, 50] != c.read(1)[50, 50]


# ---- P1: read_local fill/nodata match ----
def test_read_local_fill_matches_nodata_when_tag_missing(tmp_path):
    """SYNTHETIC: no nodata tag; bbox past edge; fill must equal returned nodata."""
    path = tmp_path / "nonodata.tif"
    y, x = np.mgrid[0:40, 0:40].astype("float32")
    dtm = (1000 + x + y).astype("float32")
    tr = from_origin(500000, 4400000, 1.0, 1.0)
    with rasterio.open(
        path, "w", driver="GTiff", height=40, width=40, count=1, dtype="float32",
        crs="EPSG:26913", transform=tr,  # no nodata tag
    ) as d:
        d.write(dtm, 1)
    # Bbox in lon/lat covering area that extends beyond raster (west of origin)
    # Origin 500000,4400000 EPSG:26913 ≈ lon/lat near -105. something
    # Use a bbox that transform_bounds into a window past the left edge
    west, south, east, north = transform_bounds("EPSG:26913", "EPSG:4326", 499980, 4399980, 500020, 4400020, densify_pts=21)
    arr, out_tr, crs, nd = m.read_local(str(path), (west, south, east, north))
    assert nd == m.NODATA_OUT
    # Padded (out-of-bounds) pixels must equal nd so local_relief masks them
    assert (arr == nd).any()
    lrm = m.local_relief(arr, nd, 2.0)
    # Edge fill must not appear as extreme valid elevation spikes in LRM valid set
    valid = lrm[lrm != m.NODATA_OUT]
    assert valid.size > 0
    assert np.nanmax(np.abs(valid)) < 500  # plane residual is small; spike would be ~1000s


def test_read_local_nodata_zero_fill_consistent(tmp_path):
    """SYNTHETIC: nodata=0 must fill with 0, not -9999."""
    path = tmp_path / "zero_nd.tif"
    dtm = np.full((30, 30), 100.0, dtype="float32")
    assert np.all(dtm != 0.0)  # source cells all nonzero; zeros may come only from OOB fill
    tr = from_origin(500000, 4400000, 1.0, 1.0)
    with rasterio.open(
        path, "w", driver="GTiff", height=30, width=30, count=1, dtype="float32",
        crs="EPSG:26913", transform=tr, nodata=0.0,
    ) as d:
        d.write(dtm, 1)
    west, south, east, north = transform_bounds(
        "EPSG:26913", "EPSG:4326", 499990, 4399990, 500015, 4400015, densify_pts=21
    )
    arr, out_tr, _, nd = m.read_local(str(path), (west, south, east, north))
    assert nd == 0.0
    assert (arr == 0.0).any()
    rr, cc = np.meshgrid(np.arange(arr.shape[0]), np.arange(arr.shape[1]), indexing="ij")
    xs, ys = rasterio.transform.xy(out_tr, rr.ravel(), cc.ravel(), offset="center")
    src_rows, src_cols = rasterio.transform.rowcol(tr, xs, ys)
    src_rows = np.asarray(src_rows).reshape(arr.shape)
    src_cols = np.asarray(src_cols).reshape(arr.shape)
    in_bounds = (src_rows >= 0) & (src_rows < 30) & (src_cols >= 0) & (src_cols < 30)
    assert np.all(arr[in_bounds] == 100.0)
    assert np.all(arr[~in_bounds] == 0.0)
    assert (~in_bounds).any()


# ---- P1: TNM pagination ----
def test_tnm_paginates_until_total():
    """SYNTHETIC mock: two full pages + remainder; no silent truncation."""
    calls = []

    def fetcher(bbox, offset, max_n):
        calls.append(offset)
        total = 120
        start = offset
        end = min(offset + max_n, total)
        items = [
            {"downloadURL": f"https://example.invalid/tile_{i}.tif"}
            for i in range(start, end)
        ]
        return {"total": total, "items": items}

    urls = m.tnm_cog_urls((-105.3, 40.0, -105.2, 40.1), page_fetcher=fetcher)
    assert len(urls) == 120
    assert calls == [0, 50, 100]


def test_tnm_hard_fails_when_truncated():
    def fetcher(bbox, offset, max_n):
        # Claims total=80 but only ever returns first page
        if offset == 0:
            return {
                "total": 80,
                "items": [{"downloadURL": f"https://example.invalid/t_{i}.tif"} for i in range(50)],
            }
        return {"total": 80, "items": []}

    with pytest.raises(RuntimeError, match="truncat|stalled"):
        m.tnm_cog_urls((-105.3, 40.0, -105.2, 40.1), page_fetcher=fetcher)


# ---- P1: mixed CRS — WarpedVRT alignment, no silent drop ----
def test_open_aligned_tiles_reprojects_mixed_crs(tmp_path, monkeypatch):
    """SYNTHETIC: two local GeoTIFFs in different projected CRS; both kept via WarpedVRT."""
    p1 = tmp_path / "utm13.tif"
    p2 = tmp_path / "utm12.tif"
    for path, epsg, x0 in ((p1, 26913, 500000), (p2, 26912, 400000)):
        dtm = np.full((20, 20), 1500.0, dtype="float32")
        tr = from_origin(x0, 4400000, 1.0, 1.0)
        with rasterio.open(
            path, "w", driver="GTiff", height=20, width=20, count=1, dtype="float32",
            crs=f"EPSG:{epsg}", transform=tr, nodata=-9999.0,
        ) as d:
            d.write(dtm, 1)

    # Bypass /vsicurl by opening local paths: patch rasterio.open used inside open_aligned_tiles
    urls = [str(p1), str(p2)]
    real_open = rasterio.open

    def local_open(path, *args, **kwargs):
        # strip /vsicurl/ prefix if present
        if isinstance(path, str) and path.startswith("/vsicurl/"):
            path = path[len("/vsicurl/"):]
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(m.rasterio, "open", local_open)
    aligned, crs, handles = m.open_aligned_tiles(urls)
    try:
        assert len(aligned) == 2
        assert crs.to_epsg() == 26913
        assert all(s.crs == crs for s in aligned)
    finally:
        for h in handles:
            h.close()


# ---- P1: projected meters / reject geographic / convert feet ----
def test_geographic_crs_hard_fails(tmp_path):
    path = tmp_path / "geo.tif"
    dtm = np.zeros((10, 10), dtype="float32")
    tr = from_origin(-105.3, 40.1, 0.0001, 0.0001)
    with rasterio.open(
        path, "w", driver="GTiff", height=10, width=10, count=1, dtype="float32",
        crs="EPSG:4326", transform=tr,
    ) as d:
        d.write(dtm, 1)
    with pytest.raises(SystemExit):
        m.main(["--input", str(path), "--sigma-m", "10", "-o", str(tmp_path / "out.tif")])


def test_feet_crs_converts_sigma(tmp_path):
    """SYNTHETIC EPSG:2263 (US survey foot): sigma_m must convert via linear_units_factor."""
    path = tmp_path / "ft.tif"
    # 1 US survey foot ≈ 0.3048006 m; use 1-ft pixels
    dtm = np.zeros((40, 40), dtype="float32")
    dtm[20, 20] = 5.0
    tr = from_origin(1000000, 200000, 1.0, 1.0)
    with rasterio.open(
        path, "w", driver="GTiff", height=40, width=40, count=1, dtype="float32",
        crs="EPSG:2263", transform=tr, nodata=-9999.0,
    ) as d:
        d.write(dtm, 1)
    to_m = m.meters_per_crs_unit(CRS.from_epsg(2263))
    assert abs(to_m - 0.30480060960121924) < 1e-9
    px = m.pixel_size_meters(CRS.from_epsg(2263), tr)
    assert abs(px - to_m) < 1e-9
    m.main(["--input", str(path), "--sigma-m", "3", "-o", str(tmp_path / "out.tif")])
    # sigma_px = 3 / px ≈ 9.84; CLI path is local_relief(arr, nd, sigma_px)
    expected_sigma_px = 3.0 / px
    with rasterio.open(tmp_path / "out.tif") as out:
        assert out.crs.to_epsg() == 2263
        # Fixture has no nodata pixels; nodata=None matches the (dtm, nodata, sigma_px) signature
        np.testing.assert_allclose(
            out.read(1), m.local_relief(dtm, None, expected_sigma_px), atol=1e-5
        )


def test_sigma_m_must_be_positive(tmp_path):
    make_fixture(tmp_path / "s.tif", nodata=False)
    with pytest.raises(SystemExit):
        m.main(["--input", str(tmp_path / "s.tif"), "--sigma-m", "0", "-o", str(tmp_path / "o.tif")])
    with pytest.raises(SystemExit):
        m.main(["--input", str(tmp_path / "s.tif"), "--sigma-m", "-1", "-o", str(tmp_path / "o.tif")])


def test_effective_nodata_helpers():
    assert m.effective_nodata(None) == m.NODATA_OUT
    assert m.effective_nodata(0) == 0.0
    assert m.effective_nodata(-9999) == -9999.0
