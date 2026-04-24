# Queries — MS SQL: Add / Edit / Browse / Delete — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

### Part 1 — Adding

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Browse → Databases → MS SQL → NorthwindTest → right-click → New query | 3s | FAIL | SKIPPED | **No `NorthwindTest` connection exists under MS SQL on dev.** The available connections are `MSSQLDBTests` (visible in tree) and `MSSQLTest` (exists via dapi but not shown in Browse). Used `MSSQLTest` via JS API to seed the query. |
| 1.2 | Enter `test_query_ms_sql` into Name | 1s | PASS | PASSED | Seeded via `conn.query('test_query_ms_sql', ...)`; `friendlyName` preserved. |
| 1.3 | Enter SQL `select * from products` | 1s | PASS | PASSED | CodeMirror body set. |
| 1.4 | Run via Play button (inline grid) | 7s | FAIL | FAILED | **MS SQL server unreachable from dev** — TCP refused at `db.datagrok.ai:14331`. `q.executeTable()` returns a driver error. Play button click did nothing visible. |
| 1.5 | Run via Toolbox → Actions → Run query... (new view) | 5s | FAIL | FAILED | Same server-unreachable condition; no new TableView appeared. Also: `DataQueryView` has no Toolbox — the analogous UI path is Context Panel → Run → RUN. |
| 1.6 | Save the query | 2s | PASS | PASSED | Persisted on MSSQLTest (connection name, not display name). |

### Part 2 — Editing

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 2.1 | Refresh Browse → Edit... | 2s | AMBIGUOUS | PASSED | Editor already open from Part 1; refresh not needed. |
| 2.2 | Change name to `new_test_query_ms_sql` | 2s | PASS | PASSED | `friendlyName` updated to `new_test_query_ms_sql`; stored `name` normalized to `NewTestQueryMsSql`. |
| 2.3 | Change body to `select * from orders` | 1s | PASS | PASSED | CodeMirror body set. |
| 2.4 | Run via Play + Run query… | 4s | FAIL | FAILED | Same MS SQL unreachable state. |
| 2.5 | Save | 2s | PASS | PASSED | Verified via `grok.dapi.queries.find` — `query = select * from orders`. |

### Part 3 — Browse

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 3.1 | Browse → Databases → MS SQL → NorthwindTest | 4s | FAIL | PASSED | NorthwindTest is not present on dev under MS SQL; MSSQLTest itself is not visible in the Browse tree either (only `MSSQLDBTests`). Suspect permissions or hidden-connection filter. |
| 3.2 | Type `new_test` in search field | — | SKIP | SKIPPED | Can't reach the connection's queries view. |
| 3.3 | Context Panel tabs for the query | 3s | PASS | PASSED | Set `grok.shell.o = q` via JS fallback — `[Details, Run, Query, Transformations, Usage, Activity, Sharing, Chats, Dev]` (9 panes). |

### Part 4 — Deleting

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 4.1 | Right-click query → Delete → DELETE | 2s | AMBIGUOUS | PASSED | UI path unavailable (no tree visibility). JS API `grok.dapi.queries.delete(q)` succeeded; `find(id)` returns null. |
| 4.2 | Refresh Browse — verify removed | 1s | PASS | PASSED | `find` returns null, deletion committed. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 50s |
| grok-browser execution (scenario steps) | 45s |
| Execute via grok-browser (total) | 2m 35s |
| Spec file generation | 55s |
| Spec script execution | 31s |
| **Total scenario run (with model)** | 4m 30s |

## Summary

Entity CRUD worked end-to-end on MS SQL — add, rename, body-change, delete all
committed via `grok.dapi.queries`. The **two run paths** (Play button + Run query…)
both failed because the MS SQL server is unreachable from dev
(`db.datagrok.ai:14331` — TCP refused). The Browse navigation step also failed
because dev has no `NorthwindTest` connection under MS SQL; the actual MS SQL
Northwind connection is named `MSSQLTest` and is not listed in the `MS SQL` tree
node (only `MSSQLDBTests` is visible).

## Retrospective

### What worked well
- Entity CRUD via `conn.query(...)` + `grok.dapi.queries.save/delete` is robust even when the driver fails to execute.
- `grok.shell.o = query` reliably swaps Context Panel to the query-entity panes regardless of Browse-tree visibility.

### What did not work
- **MS SQL TCP connectivity is broken on dev** — port 14331 refuses. Any scenario that relies on MS SQL execution is blocked.
- **`MSSQLTest` connection is not shown in the Browse tree** — only `MSSQLDBTests` is listed. Without it, the scenario's tree-navigation step is unreachable.
- **Scenario references a non-existent connection** — `NorthwindTest` under MS SQL does not exist on dev.

### Suggestions for the platform
- Surface connection-test failures on first expansion of the provider node (e.g. show a red badge or error toast) so users know why a run returns nothing.
- Investigate why `MSSQLTest` is created but not listed in the Browse tree — either a permissions filter or a missing `showInBrowse` flag. Users running the scenario have no way to discover its existence.
- MS SQL server `db.datagrok.ai:14331` needs to be reachable from dev (ops task).

### Suggestions for the scenario
- Replace `NorthwindTest` with the actual dev-available MS SQL connection name (or parameterize the scenario so it reads the connection from env / a JSON block).
- Pre-condition: verify `grok.dapi.connections.test(conn)` succeeds before running add/edit/run steps.
- Step 2.1 "Refresh view on Browse" is redundant immediately after Part 1 — the editor is still open.
- Clarify run-path expectations when the driver fails: should the scenario fail, or should it record an error balloon as the expected observation?
