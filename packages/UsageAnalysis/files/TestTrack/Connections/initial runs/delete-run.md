# Delete — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Find connection `new_test_postgres` in Browse → Connections | 4s | PASS | PASSED | Pre-existing `Agolovko:NewTestPostgres` (friendlyName `new_test_postgres`) located in Postgres connections via `[data-link="/db/Agolovko.NewTestPostgres"]` |
| 2a | Right-click and select Delete from context menu | 2s | PASS | PASSED | `contextmenu` dispatched on the data-link element; menu showed Browse, New Query…, New Visual Query…, Delete…, Edit…, Rename…, Clone…, Clear cache, Browse queries, Test connection, Share…, Copy → ID/Grok name/Markup/URL — clicked Delete… |
| 2b | In confirmation dialog, click YES | 2s | PASS | PASSED | Dialog showed `DELETE` / `CANCEL` (not `YES`) — clicked `[name="button-DELETE"]` |
| 2c | Check connection disappeared | 2s | PASS | PASSED | After click + sync refresh, `[data-link="/db/Agolovko.NewTestPostgres"]` gone from DOM; `grok.dapi.connections.find(id)` returns null |
| 3 | Find connection `test_postgres_2` in Browse → Databases | 1s | PASS | PASSED | Connection had to be created via grok s before the run — adding scenario hadn't been re-run on dev. Located via `[data-link="/db/null.test_postgres_2"]` |
| 4a | Right-click and select Delete from context menu | 1s | PASS | PASSED | Same right-click flow as 2a; Delete… clicked |
| 4b | In confirmation dialog, click YES | 1s | PASS | PASSED | Clicked DELETE in confirmation dialog |
| 4c | Check connection is no longer present | 2s | PASS | PASSED | After refresh, `[data-link="/db/null.test_postgres_2"]` gone from DOM and JS API returns null |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 30s |
| grok-browser execution (scenario steps) | 25s |
| Execute via grok-browser (total) | 1m 0s |
| Spec file generation | 2m 0s |
| Spec script execution | 27s |
| **Total scenario run (with model)** | 4m 30s |

## Summary

All 8 sub-steps passed. Both connections (`new_test_postgres` and `test_postgres_2`) were deleted successfully through right-click → Delete… → DELETE-button confirmation. The confirmation dialog uses `DELETE` / `CANCEL` buttons — the scenario's "click YES" wording is wrong. `test_postgres_2` did not exist on dev before the run (Adding-spec hadn't been executed); created it via `grok s connections save --json` as part of preconditions, then verified deletion. Total scenario run (with model): ~4m 30s.

## Retrospective

### What worked well
- `[data-link="/db/<namespace>.<slug>"]` again proved unique and resilient — much more reliable than label-text matching
- `[name="button-DELETE"]` is the canonical confirm button selector inside the d4-dialog
- Verifying deletion both in DOM (data-link gone) AND via `grok.dapi.connections.find(id)` returning null catches both UI staleness and server-side failures
- The spec is now self-contained: `ensureConnection()` creates the targets via `DG.DataConnection.create + name/friendlyName override + save()` if they're missing, so it can run independently of the Adding/Edit specs

### What did not work
- Initial spec used `DG.DataConnection.fromInfo` — that method does not exist on the JS API. The working constructor is `DG.DataConnection.create(slug, {dataSource: 'Postgres', ...params})`, but it auto-PascalCases the slug to a friendlyName so the slug must be reassigned via `conn.name = slug` before saving
- `grok s connections save --json` capitalizes the auto-derived friendlyName (e.g. typed `test_postgres_2` → friendlyName `Test_postgres_2`) — must reset friendlyName explicitly after the save

### Suggestions for the platform
- Confirmation dialog's primary button text — keep `DELETE` (more explicit than `YES`) but update Test Track scenario wording to match
- Add a brief success toast after deletion (`Connection deleted`) — currently the only feedback is the row disappearing
- `DG.DataConnection.create()` should NOT auto-capitalize the `name` slug — only the `friendlyName`. As-is, plugin authors who want a lowercase slug must manually reassign `name` after construction

### Suggestions for the scenario
- Steps 2b and 4b say "click YES" — the actual button label is `DELETE`. Update the scenario
- Add a precondition note: "Run Adding scenario first; if `test_postgres_2` is missing, create it before starting"
- Clarify that "Browse → Platform → Connections" and "Browse → Databases" both list the same connections — currently the scenario implies they are different
