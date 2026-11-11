"""Opt-in live network smoke tests for TNM and OpenTopography.

Skipped unless STRATA_LIVE_SMOKE=1. OpenTopography reads OPENTOPO_API_KEY from
the environment and skips when it is missing. This file must not contain a key.
"""
from __future__ import annotations
import os
import pytest
import rasterio
from rasterio.windows import Window
import lrm_3dep as m

# README example: a few hundred metres of Boulder 3DEP, small enough for a smoke read.
BBOX = (-105.2800, 40.0100, -105.2750, 40.0140)

_LIVE = pytest.mark.skipif(
    os.environ.get("STRATA_LIVE_SMOKE") != "1",
    reason="live smoke disabled; set STRATA_LIVE_SMOKE=1 to contact TNM and OpenTopography",
)


@_LIVE
@pytest.mark.live_smoke
def test_tnm_live_smoke_lists_and_reads_one_pixel():
    urls = m.tnm_cog_urls(BBOX)
    assert urls, "TNM returned no 1 m GeoTIFF for the smoke bbox"
    assert all(u.startswith("https://") and u.endswith(".tif") for u in urls)
    with rasterio.open("/vsicurl/" + urls[0]) as src:
        assert src.crs is not None
        chip = src.read(1, window=Window(0, 0, 1, 1))
    assert chip.shape == (1, 1)


@_LIVE
@pytest.mark.live_smoke
def test_opentopo_live_smoke_reads_key_from_env_only():
    if not os.environ.get("OPENTOPO_API_KEY"):
        pytest.skip(
            "OPENTOPO_API_KEY is not set; OpenTopography smoke reads it from the environment only"
        )
    arr, _tr, crs, _nd, srcs, _core = m.read_opentopo(BBOX, sigma_m=None)
    assert arr.ndim == 2 and arr.size > 0
    assert crs is not None and crs.is_projected
    assert srcs == [m.OT_API]
    assert all("API_Key" not in str(s) and "OPENTOPO_API_KEY" not in str(s) for s in srcs)
