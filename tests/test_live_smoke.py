"""Opt-in live network smoke tests for TNM and OpenTopography.

Skipped unless STRATA_LIVE_SMOKE=1. OpenTopography reads OPENTOPO_API_KEY from
the environment and skips when it is missing. This file must not contain a key.

The TNM smoke calls read_tnm, so it lists products and merges them. The smoke
bbox sits in one UTM zone, so that live call does not open a WarpedVRT.
Mixed-CRS WarpedVRT resampling is covered by the synthetic tests, not here.
"""
from __future__ import annotations
import os
import numpy as np
import pytest
import lrm_3dep as m
from rasterio.warp import transform_bounds

# README example: a few hundred metres of Boulder 3DEP, small enough for a smoke read.
BBOX = (-105.2800, 40.0100, -105.2750, 40.0140)

_LIVE = pytest.mark.skipif(
    os.environ.get("STRATA_LIVE_SMOKE") != "1",
    reason="live smoke disabled; set STRATA_LIVE_SMOKE=1 to contact TNM and OpenTopography",
)


def _finite_elevations(arr, nodata):
    valid = np.isfinite(arr)
    if nodata is not None:
        valid &= arr != nodata
    return arr[valid]


def _intersects_bbox(arr, transform, crs):
    height, width = arr.shape
    left = transform.c
    top = transform.f
    right = left + width * transform.a
    bottom = top + height * transform.e
    west, south, east, north = transform_bounds(
        crs, "EPSG:4326",
        min(left, right), min(bottom, top), max(left, right), max(bottom, top),
        densify_pts=21,
    )
    return east > BBOX[0] and west < BBOX[2] and north > BBOX[1] and south < BBOX[3]


@_LIVE
@pytest.mark.live_smoke
def test_tnm_live_smoke_reads_merged_mosaic():
    """read_tnm merges the listed tiles. This bbox does not force a live WarpedVRT."""
    arr, transform, crs, nd, urls, _core = m.read_tnm(BBOX, sigma_m=None)
    assert urls, "TNM returned no 1 m GeoTIFF for the smoke bbox"
    assert all(u.startswith("https://") and u.endswith(".tif") for u in urls)
    assert arr.ndim == 2 and arr.shape[0] >= 2 and arr.shape[1] >= 2
    assert crs is not None and crs.is_projected
    vals = _finite_elevations(arr, nd)
    assert vals.size > 10
    assert 500.0 < float(vals.min()) and float(vals.max()) < 5000.0


@_LIVE
@pytest.mark.live_smoke
def test_opentopo_live_smoke_reads_key_from_env_only():
    key = os.environ.get("OPENTOPO_API_KEY")
    if not key:
        pytest.skip(
            "OPENTOPO_API_KEY is not set; OpenTopography smoke reads it from the environment only"
        )
    arr, transform, crs, nd, srcs, _core = m.read_opentopo(BBOX, sigma_m=None)
    assert arr.ndim == 2 and arr.shape[0] >= 2 and arr.shape[1] >= 2
    assert crs is not None and crs.is_projected
    vals = _finite_elevations(arr, nd)
    assert vals.size > 10
    assert 500.0 < float(vals.min()) and float(vals.max()) < 5000.0
    assert _intersects_bbox(arr, transform, crs)
