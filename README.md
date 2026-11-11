# Strata Mapping

Permanent repository for the Strata / Global Mapping sequential build.

Daily automation advances one step at a time and commits here so the work is not vapor.

## Sequential Plan
1. Streamed topography foundation
2. Mapping driver + gap fill
3. First residual site
4. Noticing threshold
5. Global layer stays current

## Status
Step 1 marked complete in prior runs. Next work lands in this repo.

**2026-11-05 (MechSensorForge):** Step 1 evidence landed as `lrm_3dep.py` — Local Relief Model from USGS 3DEP 1 m DTM (port of local R&D work for revvel-standards#18207). Steps 2–5 still open; see `PLAN.md`.

## Local Relief Model (`lrm_3dep.py`)

Computes LRM = DTM − Gaussian(DTM, σ) from USGS 3DEP 1 m bare-earth DTM for a lon/lat (EPSG:4326) bounding box.

### Setup and tests
```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python3 -m pytest -v tests            # synthetic-fixture tests; live smoke skipped
# optional network (OpenTopography also needs OPENTOPO_API_KEY in the environment; never commit it):
# STRATA_LIVE_SMOKE=1 python3 -m pytest -v tests -m live_smoke
```

### Manual run
```bash
python3 lrm_3dep.py --bbox -105.2800 40.0100 -105.2750 40.0140 --sigma-m 10 -o lrm.tif
# optional fallback: export OPENTOPO_API_KEY=...; add --source opentopo
# offline full file:
python3 lrm_3dep.py --input my_dtm.tif --sigma-m 10 -o lrm.tif
# offline crop (four lon/lat floats required):
python3 lrm_3dep.py --input my_dtm.tif --bbox -105.2800 40.0100 -105.2750 40.0140 --sigma-m 10 -o lrm.tif
```

### Behaviour
- Source 1 (TNM): USGS TNM Access API lists 1 m DEM COGs; **paginated** until complete (hard-fail if truncated). Windowed over HTTP via GDAL `/vsicurl/` (no full-tile download). Mixed-CRS tiles are reprojected into the first tile's CRS via `WarpedVRT` (never silently dropped). Resampling is nearest in the warp and again when the tiles are merged. On a 20 m pixel the half-pixel box of the synthetic slope (0.05 per metre east, 0.02 per metre north) is 0.7 m, so two stages are about 1.4 m; the test allows 1.5 m. Deterministic mosaics for a fixed TNM response. The opt-in TNM smoke calls `read_tnm` (list and merge). Its bbox is one UTM zone, so that live call does not open a `WarpedVRT`; mixed-CRS warp is covered by the synthetic tests.
- Source 2 (OpenTopography): downloads the requested bbox GeoTIFF (not streamed), with request timeout; temp file always unlinked. Needs `OPENTOPO_API_KEY`.
- σ in meters is applied per axis from the absolute pixel sizes `|a|` and `|e|` (after converting feet). The axes are not averaged. A rotated or sheared affine (`b` or `d` nonzero) is an error: a rotated-pixel kernel is not implemented. CRS must be projected with known linear units. Geographic CRS hard-fails.
- Edge halo (~4σ) is read before smoothing, then cropped back to the requested bbox.
- Nodata-aware smoothing; boundless fill value always matches the nodata mask; output float32, nodata −9999, LZW.

Sources listed under `Library/`.
