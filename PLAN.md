# Sequential Build Order

Finish each item before opening the next.

1. Streamed topography foundation  
   Done when a player can move across one real-world patch from public data only.  
   **Evidence (2026-11-05):** `lrm_3dep.py` CLI streams USGS 3DEP 1 m bare-earth DTM from TNM COGs (HTTP windowed `/vsicurl/`) for a lon/lat bbox and writes a nodata-aware Gaussian Local Relief Model GeoTIFF. OpenTopography is a fallback that **downloads** the requested bbox GeoTIFF (not streamed). Synthetic pytest suite in `tests/test_lrm.py`. Refs midnghtsapphire/revvel-standards#18207.

2. Mapping driver + gap fill  
   Done when a gap becomes navigable after the driver pass and the player can ignore the driver.

3. First residual site  
   Done when a single site (clandestine-grave or equivalent) can be spotted, marked, and resolved.

4. Noticing threshold on the residual  
   Done when a player who was not looking still registers the mismatch.

5. Global layer stays current  
   Done when the driver or tile stream advances while the player is only doing the residual.
