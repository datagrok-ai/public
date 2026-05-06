# Edit — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Reload the tree in Browser | 3s | PASS | PASSED | Clicked `[name="icon-sync"]` in the Postgres connections view; tree refreshed without error |
| 2 | Right-click test_postgres connection | 2s | PASS | PASSED | `contextmenu` dispatched on `[data-link="/db/null.test_postgres"]`; full context menu appeared (Browse, New Query…, Delete…, Edit…, Rename…, Clone…, Clear cache, Browse queries, Test connection, Share…, Copy → ID/Grok name/Markup/URL) |
| 3 | Select Edit from context menu | 2s | PASS | PASSED | Clicked `Edit…` in `.d4-menu-item`; "Edit Connection" dialog opened with Name=test_postgres, Server=db.datagrok.ai, Port=54322, Db=northwind, Login=datagrok |
| 4 | Change name to `new_test_postgres` and click OK | 6s | PASS | PASSED | Set `[name="input-Name"]` via native value setter + change event; OK clicked; verified via JS API that `friendlyName` became `new_test_postgres` (slug `name` stays `test_postgres`) |
| 5 | Change login/password with arbitrary data and save | 6s | PASS | PASSED | Right-clicked → Edit; set `[name="input-Login"]`=wronguser, `[name="input-Password"]`=wrongpassword; OK clicked; saved without error |
| 6 | Test the connection — should return error | 33s | PASS | PASSED | UI balloon disappeared too fast to capture, so verified via `await conn.test()`: returned `org.postgresql.util.PSQLException: FATAL: password authentication failed for user "wronguser"` — matches expected format |
| 7 | Set right login/password and test — should be fine | n/a | SKIP | SKIPPED | Skipped — real DB password not available; mechanism verified in step 6 |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 0s |
| grok-browser execution (scenario steps) | 1m 0s |
| Execute via grok-browser (total) | 2m 0s |
| Spec file generation | 4m 0s |
| Spec script execution | 1m 6s |
| **Total scenario run (with model)** | 7m 30s |

## Summary

6 of 7 steps fully passed; step 7 skipped (real DB password unavailable). The connection rename, credential modification, and Test connection error all behave correctly on dev.datagrok.ai. Total scenario run (with model): ~7m 30s. Two pre-existing connections (`Agolovko:NewTestPostgres` and `null:test_postgres`) shared the friendlyName `new_test_postgres` after rename — UI must be disambiguated by `data-link` / `name` attribute rather than label text.

## Retrospective

### What worked well
- `[data-link="/db/<namespace>.<slug>"]` and `[name="span-<slug>"]` are unique per connection, even when friendlyName collides — the right way to scope right-clicks
- `await conn.test()` returns the same `PSQLException` string the UI balloon displays, so verification works headlessly
- Setting `[name="input-Name"]` etc. via native setter + `change` event is the only way to make Datagrok's Dart-side form pick up the new value

### What did not work
- Step 7 cannot be tested without DevOps-supplied password — same blocker as the prior public.datagrok.ai run
- `expect(page.locator('text=new_test_postgres')).toBeVisible()` is unreliable: the existing `Agolovko:NewTestPostgres` connection also displays as `new_test_postgres`, so the locator passes regardless of whether our rename succeeded — must verify via `grok.dapi.connections.find(id).friendlyName` instead
- The Edit Connection dialog renders with `inputCount: 0` for ~1s after open. Reading inputs immediately after `Edit…` click sees an empty form and silently no-ops — must `waitForFunction(() => dialog.querySelectorAll('input').length >= 5)` before reading/writing
- The Test connection toast is shown only briefly after a failure; polling balloons every 500ms missed it consistently. Use the JS API for verification

### Suggestions for the platform
- Persist the most recent Test-connection result on the entity so the Details panel shows pass/fail without rerunning the test
- Add a unique `data-id` (connection UUID) on connection cards so automation never has to disambiguate by namespace-prefixed slug
- Validate friendlyName uniqueness on save — the platform happily allows two connections in different namespaces to share a display name, which is confusing in the UI

### Suggestions for the scenario
- Step 4 should clarify that what is changed is the connection's display name (`friendlyName`), not the slug `name` — current wording suggests a full rename
- Step 7 needs a concrete credential source ("use the password from `~/.grok/config.yaml dev`" or "ask DevOps for the `datagrok@db.datagrok.ai:54322` password") so it isn't routinely skipped
- Add a precondition: "If a previous run left `friendlyName=new_test_postgres`, restore it first" — otherwise step 4's "rename" is a no-op
