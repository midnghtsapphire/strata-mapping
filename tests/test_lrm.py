"""Tests use a SYNTHETIC TEST FIXTURE (generated plane + pit), NOT measured elevation data."""
import numpy as np, rasterio, pytest
from rasterio.transform import from_origin
import lrm_3dep as m

def make_fixture(path, nodata=True):
    # SYNTHETIC: tilted plane plus a 2 m-deep 5x5 pit, 1 m pixels, EPSG:26913
    y, x = np.mgrid[0:100, 0:100].astype("float32")
    dtm = 1500 + 0.1 * x + 0.05 * y
    dtm[48:53, 48:53] -= 2.0
    if nodata:
        dtm[0:5, 0:5] = -9999.0
    tr = from_origin(500000, 4400000, 1.0, 1.0)
    with rasterio.open(path, "w", driver="GTiff", height=100, width=100, count=1,
                       dtype="float32", crs="EPSG:26913", transform=tr, nodata=-9999.0) as d:
        d.write(dtm, 1)
    return dtm, tr

def test_plane_gives_zero_lrm():
    y, x = np.mgrid[0:60, 0:60].astype("float64")
    lrm = m.local_relief(3 * x + 2 * y, None, 5.0)
    assert np.abs(lrm[21:39, 21:39]).max() < 1e-4  # interior beyond 4-sigma kernel radius: Gaussian preserves linear surfaces

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

def test_sigma_meters_scales_with_pixel():
    dtm = np.zeros((40, 40)); dtm[20, 20] = 1
    assert m.local_relief(dtm, None, 2.0)[20, 20] != m.local_relief(dtm, None, 4.0)[20, 20]
