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
python3 -m pytest -v tests            # synthetic-fixture tests
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
- Source 1 (TNM): USGS TNM Access API lists 1 m DEM COGs; **paginated** until complete (hard-fail if truncated). Windowed over HTTP via GDAL `/vsicurl/` (no full-tile download). Mixed-CRS tiles are reprojected into the first tile's CRS via `WarpedVRT` (never silently dropped). Deterministic mosaics for a fixed TNM response.
- Source 2 (OpenTopography): downloads the requested bbox GeoTIFF (not streamed), with request timeout; temp file always unlinked. Needs `OPENTOPO_API_KEY`.
- σ in meters; CRS must be projected with known linear units (meters or feet — feet are converted). Geographic CRS hard-fails.
- Edge halo (~4σ) is read before smoothing, then cropped back to the requested bbox.
- Nodata-aware smoothing; boundless fill value always matches the nodata mask; output float32, nodata −9999, LZW.

Sources listed under `Library/`.
