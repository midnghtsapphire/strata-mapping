# Changelog

## [Unreleased] - 2026-11-11
### Fixed
- Rotation/shear rejection uses `not (abs(b) <= tol)` (and the same for `d`), so a NaN component fails closed instead of passing the comparison.
### Added
- Tolerance-edge and NaN affine tests; anisotropic halo checked through the CLI `--input` + `--bbox` path with expectations computed in the test.
- Mixed-CRS abutment of two non-overlapping tiles (no nodata or zero at the seam) and a bound on nearest-resample error for a plane that crosses UTM zones.
- Live smoke checks real shape, finite elevations, projected CRS, and bbox intersection. The TNM smoke calls `read_tnm` (merge). That bbox is one UTM zone, so it does not open a live `WarpedVRT`.

## [Unreleased] - 2026-11-11
### Changed
- Renamed the mixed-CRS seam test to `test_mixed_crs_seam_nodata_tag_and_no_zero_fill`. It is a nodata and 0-fill guard, not a per-row continuity check.
- Slope docstring and README: `resampling == nearest` is the guard against bilinear. The 0.1 m mean lower bound is only a weak secondary check. Thresholds unchanged.
### Fixed
- Mixed-CRS seam test writes the east tile with no nodata tag and about one pixel of overlap, so the old `WarpedVRT` 0-fill fails it (`nodata` stays unset and the seam band contains zeros).
### Added
- Nearest-slope check requires east-of-seam mean error above 0.1 m, which bilinear warp misses. The 1.5 m cap is two half-pixel stages (0.7 m each, about 1.4 m) plus 0.1 m.
### Changed
- Gaussian sigma uses per-axis pixel size from the affine (`|a|` and `|e|` independently, in metres) instead of `mean(|a|, |e|)`. Rotation or shear (`b` or `d` nonzero) raises `ValueError`.
- Mixed-CRS `WarpedVRT` sets source and destination nodata explicitly and resamples with nearest, so a missing nodata tag cannot fill a seam with 0.
### Added
- Synthetic tests for mixed-CRS warp/merge resampling, nodata seams at tile edges, merge alignment, and boundless/nodata consistency.
- Opt-in live TNM and OpenTopography smoke tests (`live_smoke` marker and `STRATA_LIVE_SMOKE=1`). OpenTopography reads `OPENTOPO_API_KEY` from the environment only and skips when it is unset.

## [Unreleased] - 2026-11-06
### Fixed
- P1 `read_local`: boundless `fill_value` and returned nodata now use the same `effective_nodata()` (handles `None` and `0`).
- P1 TNM: paginate `offset`/`max` until `total` items fetched; hard-fail on truncation/stall (no silent max=50 cut).
- P1 Mixed CRS: reproject tiles into the first tile's CRS via `WarpedVRT` (never silently drop).
- P1 Units: require projected CRS; convert feet via `linear_units_factor`; hard-fail geographic CRS before `sigma_m / px`.
- Cubic follow-up on merged PR #2: append-only ledger correction restoring `3712474` file/symbol/branch identifiers and recording `.gitignore` (added); zero-nodata test keeps source cells nonzero so only out-of-bounds fill is 0; feet-CRS test compares CLI output to `local_relief(dtm, None, expected_sigma_px)`.
### Added
- P2: OpenTopography temp unlink + request timeout; ~4σ edge halo then crop; `--sigma-m` must be finite and > 0; synthetic P1 tests; README/PLAN nits (offline example; TNM vs OT streaming wording).

## [Unreleased] - 2026-11-05
### Added
- `lrm_3dep.py`, `tests/test_lrm.py`, `requirements.txt`: MechSensorForge R&D port of USGS 3DEP Local Relief Model CLI (synthetic pytest suite, 4 tests). Refs midnghtsapphire/revvel-standards#18207.
- `LEDGER.md`, `Library/README.md`: write-as-you-go ledger and source notes (USGS TNM, OpenTopography).
- `PLAN.md` Step 1 evidence note for LRM CLI; steps 2–5 unchanged.
