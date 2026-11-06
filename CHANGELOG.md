# Changelog

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
