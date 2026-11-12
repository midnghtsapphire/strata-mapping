# Library — strata-mapping sources

Source notes the LRM work relies on (no large binaries in-repo).

| Author / org | Year | Title | Notes / filename |
|---|---|---|---|
| USGS National Map | — | TNM Access API (`tnmaccess.nationalmap.gov`) | Lists public 3DEP 1 m COGs; primary stream source for `lrm_3dep.py` |
| USGS 3DEP | — | 1-meter bare-earth DEM (Cloud-Optimized GeoTIFF) | Read by window via GDAL `/vsicurl/` |
| OpenTopography | — | USGS 1m DEM API (`portal.opentopography.org/API/usgsdem`) | Fallback when `OPENTOPO_API_KEY` is set |
| midnghtsapphire/revvel-standards | 2026 | Work request #18207 (LRM GeoTIFF script) | Governing WR for this port |

Do not store measured elevation GeoTIFFs here unless Audrey asks; tests use a synthetic fixture only.
