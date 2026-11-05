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
python3 -m pytest -v tests            # synthetic-fixture tests (4)
```

### Manual run
```bash
python3 lrm_3dep.py --bbox -105.2800 40.0100 -105.2750 40.0140 --sigma-m 10 -o lrm.tif
# optional fallback: export OPENTOPO_API_KEY=...; add --source opentopo
# offline: python3 lrm_3dep.py --input my_dtm.tif [--bbox ...] -o lrm.tif
```

### Behaviour
- Source 1: USGS TNM Access API lists 1 m DEM COGs; windowed over HTTP via GDAL `/vsicurl/` (no full-tile download). Deterministic mosaics for a fixed TNM response.
- Source 2: OpenTopography `API/usgsdem?datasetName=USGS1m` with `OPENTOPO_API_KEY`.
- σ in meters, converted to pixels from DTM pixel size (projected metric CRS required).
- Nodata-aware smoothing; output float32, nodata −9999, LZW, source CRS/transform preserved.

Sources listed under `Library/`.
