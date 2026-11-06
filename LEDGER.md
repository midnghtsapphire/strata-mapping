# LEDGER — strata-mapping

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
