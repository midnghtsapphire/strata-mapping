# LEDGER — strata-mapping

## 2026-11-11 07:10 America/Denver — Sense-Bedrock — anisotropic Gaussian sigma, warp/merge seams, opt-in live smoke (#3)

- **Actor:** Sense-Bedrock
- **Requested by:** issue #3 (follow-up to #2; TriageChief routing 2026-11-11). OpenHands comment 6724206545 ignored as out of scope.
- **Role:** executed
- **Files:** `lrm_3dep.py` (changed); `tests/test_lrm.py` (changed); `tests/test_warp_merge.py` (added); `tests/test_live_smoke.py` (added); `tests/conftest.py` (added); `README.md` (changed); `CHANGELOG.md` (changed); `LEDGER.md` (changed)
- **What was done:** Gaussian sigma is per-axis from affine `|a|` and `|e|` (metres), not `mean(|a|, |e|)`. `local_relief` still accepts a scalar sigma; equal axes use `(sigma, sigma)`, which matches that scalar kernel. Rotation or shear (`b` or `d` nonzero) raises `ValueError` (no rotated kernel, no silent wrong radius). Mixed-CRS `WarpedVRT` sets nodata explicitly and uses nearest resampling so a missing nodata tag cannot fill seams with 0. Synthetic tests cover WarpedVRT resampling, tile-edge nodata seams, merge alignment, and boundless/nodata fill across CRS boundaries. Live TNM and OpenTopography smoke tests are marked `live_smoke` and skipped unless `STRATA_LIVE_SMOKE=1`; OpenTopography reads `OPENTOPO_API_KEY` from the environment only and skips if it is unset. Pytest from repo root: 23 passed, 2 skipped. No merge; no push to `daily/*` or `main`; no workflow files; nothing deleted.
- **Decisions:** Reject rotation/shear instead of inventing a rotated kernel — separable `|a|`/`|e|` sigma is the wrong radius when `b` or `d` is nonzero. Keep the scalar `local_relief` path so existing square-pixel callers keep working.
- **Sources:** Issue #3; scipy.ndimage `gaussian_filter` sigma axis order (row, col); rasterio `WarpedVRT` destination nodata defaults to 0 when the source tag is missing.
- **Commits / PRs:** branch `sensecrew/bedrock/issue-3-aniso-sigma` (this commit) → draft PR into `daily/2026-11-05`. Closes #3.
- **Open / handoff:** Draft PR left open for TriageChief. Live smoke was not executed in this run (`STRATA_LIVE_SMOKE` unset; no `OPENTOPO_API_KEY`). Next owner: TriageChief / Audrey
- **Credit note:** Sense-Bedrock on issue #3

## 2026-11-06 10:30 America/Denver — MechSensorForge — correction to 2026-11-06 08:20 review-fix identifiers and 2026-11-05 Files inventory

Prior 2026-11-06 08:20 review-fix ledger entry lost backtick identifiers (shell backtick expansion). Restoring the concrete identifiers here per the append-only / never-rewrite rule; the original entry below is left intact. The 2026-11-06 correction note to the 2026-11-05 entry also lost its filename.

- **Actor:** MechSensorForge
- **Requested by:** Audrey Evans (Cubic review on PR #2; Sense-Strata follow-up)
- **Role:** executed
- **Restored Files for `3712474`:** `CHANGELOG.md` (changed); `LEDGER.md` (changed); `PLAN.md` (changed); `README.md` (changed); `lrm_3dep.py` (changed); `tests/test_lrm.py` (changed)
- **Restored symbols:** `effective_nodata`; `WarpedVRT`; TNM pagination (`offset`/`max`); `pixel_size_meters` (interim mean of `|a|` and `|e|` in metres)
- **Restored branch / commit:** `mech/strata-lrm-2026-11-05` @ `37124744890ca78a2774833b25d3ac520e024fd5` (short `3712474`); PR #2 (later squash-merged into `daily/2026-11-05` as `ee6e980bea854a89a6ace310a3d84cafec2a3894` by Proofy; MechSensorForge did not self-merge)
- **2026-11-05 inventory omission:** `.gitignore` (added in the port commit `0aa94eb96ebf7db9fe26f1305a85f248d2748c55`). Recorded here; original 2026-11-05 Files list left intact.
- **This follow-up (same commit as this entry):** `LEDGER.md` (changed); `CHANGELOG.md` (changed); `tests/test_lrm.py` (changed) — Cubic zero-nodata fixture (source cells all nonzero; only out-of-bounds fill is 0) and feet-CRS CLI vs `local_relief(dtm, None, expected_sigma_px)` check. Branch `mech/strata-cubic-fixes-2026-11-06` targeting `daily/2026-11-05`.
- **What was done:** Append-only identifier restoration and two synthetic test tightenings. No merge; no main edits; no Actions; nothing deleted.
- **Open / handoff:** Follow-up WR for anisotropic/rotated-pixel sigma, warp/merge nodata-seam tests, and opt-in live TNM/OpenTopo smoke: #3. Next owner: Sense-Strata / Audrey
- **Credit note:** MechSensorForge cubic-follow-up on own authored work


## 2026-11-06 08:20 America/Denver — MechSensorForge — Sense-Strata / gitar / cubic P1 review fixes

- **Actor:** MechSensorForge
- **Requested by:** Audrey Evans (Sense-Strata review feedback on PR #2)
- **Role:** executed
- **Files:**  (changed);  (changed);  (changed);  (changed);  (changed);  (changed)
- **What was done:** Addressed P1 findings on draft PR #2 (): (1) boundless fill/nodata match via ; (2) TNM pagination with hard-fail on truncation; (3) mixed-CRS  alignment (no silent drop); (4) projected meter/feet unit check before sigma scaling. Cheap P2: OT temp unlink + timeout, edge halo, sigma>0 validation, tighter synthetic tests, README/PLAN nits. Full pytest: 13 passed. No merge; no main/daily edits; no Actions; nothing deleted.
- **Decisions:** Reproject mixed CRS with WarpedVRT rather than hard-fail-only — because TNM 1 m mosaics often straddle UTM zones and dropping tiles silently was the bug being fixed. Paginate TNM instead of hard-failing at exactly 50 when more exist.
- **Sources:** PR #2 review threads (gitar-bot boundless-fill fix shape; cubic P1/P2); USGS TNM Access pagination (//); rasterio  / .
- **Commits / PRs:** same branch  / PR #2
- **Open / handoff:** Sense-Strata CLEAR path; evening daily→main still Audrey/Proofy. Next owner: Sense-Strata / Audrey
- **Credit note:** MechSensorForge review-fix pass on own authored PR

## 2026-11-06 — correction note to 2026-11-05 entry

Prior entry Files list omitted  (added in the same port commit). Correcting here per LEDGER never-rewrite rule; original entry left intact below.

## 2026-11-05 08:25 America/Denver — MechSensorForge — Port lrm_3dep (standards#18207)

- **Actor:** MechSensorForge
- **Requested by:** Audrey Evans (silent R&D; draft PR into daily/2026-11-05)
- **Role:** executed
- **Files:** `lrm_3dep.py` (added); `tests/test_lrm.py` (added); `requirements.txt` (added); `README.md` (changed); `PLAN.md` (changed); `CHANGELOG.md` (added); `LEDGER.md` (added); `Library/README.md` (added)
- **What was done:** Created `daily/2026-11-05` from `main`. Ported Local Relief Model CLI and synthetic tests from local `/workspace/revvel-mechsensor` (excluded `real_lrm_boulder.tif`). Extended README with run/test steps; kept prior Status text. Marked PLAN Step 1 evidence = LRM CLI; left steps 2–5 open. No Scitor issue edits; no GitHub Actions; nothing deleted.
- **Decisions:** Keep port in repo root to match local R&D layout and pytest import path — because strata-mapping is still a thin sequential-build home and Step 1 is the streamed topography foundation.
- **Sources:** USGS National Map TNM Access API; OpenTopography USGS 1m DEM API — see `Library/README.md`. Work request midnghtsapphire/revvel-standards#18207.
- **Commits / PRs:** task branch `mech/strata-lrm-2026-11-05` → draft PR into `daily/2026-11-05`
- **Open / handoff:** Steps 2–5 (mapping driver, residual site, noticing threshold, global layer). Next owner: MechSensorForge / Audrey
- **Credit note:** MechSensorForge R&D port of prior local implementation tied to #18207
