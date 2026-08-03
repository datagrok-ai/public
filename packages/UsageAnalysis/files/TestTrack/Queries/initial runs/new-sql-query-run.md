# Queries — New SQL Query from Products Table — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → Databases → Postgres → NorthwindTest → Schemas → public — list of tables opens | 7s | PASS | PASSED | `Databases` was already expanded after `/browse`; walked `Postgres → NorthwindTest → Schemas → public` and verified `products` is among the listed tables (also: `categories, customers, employees, orders, …`). |
| 2 | Right-click `products` → `New SQL Query…` → editor opens | 4s | PASS | PASSED | `contextmenu` on `.d4-tree-view-node` exposes `Get All / Get Top 100 / New SQL Query… / New Visual Query…`. Clicking `New SQL Query…` opens a `DataQueryView` with CodeMirror seeded `select * from public.products`. |
| 3a | Run via Menu Ribbon Play button — inline result | 8s | PASS | PASSED | `[name="icon-play"]` ribbon-item click flaky on first attempt; the inline `[name="viewer-Grid"]` (1145×423) appeared after retrying both the icon and the wrapping `.d4-ribbon-item` for ~3s. Result = full products table in a docked grid below the editor. |
| 3b | Run via Toolbox → Actions → Run query… — new view | 5s | PASS | PASSED | The `DataQueryView` Toolbox does have an `Actions` section with a `Run query...` link (`label.d4-link-action`). Click opens a new `TableView` with 77 rows × 10 cols (no parameter dialog — query is parameter-free). View count went 3 → 4. |
| 4 | Save the query | 5s | PASS | PASSED | Renamed editor `Name` from `products` to `tt_new_sql_query_<ts>` via `execCommand('insertText')` so Dart change listener fires; `[name="button-Save"]` persisted; `grok.dapi.queries.filter('name = "<ts>"')` confirmed it was saved on the server; cleanup deleted the entity. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 5s |
| grok-browser execution (scenario steps) | 30s |
| Execute via grok-browser (total) | 1m 35s |
| Spec file generation | 30s |
| Spec script execution | 28s |
| **Total scenario run (with model)** | 2m 33s |

## Summary

End-to-end New-SQL-Query-from-table scenario passes via the UI-first path.
Right-click → context menu → editor → Play (inline) → Toolbox `Run query…` (new view) → Save — all four steps reproduced cleanly in both the MCP run and the standalone Playwright spec (28s, all softSteps PASSED). The previous run-doc claim that `DataQueryView` has no Toolbox no longer holds — the `Actions` section with the `Run query...` link is present on dev.

## Retrospective

### What worked well
- CodeMirror seeding (`select * from public.products`) is reliable and inspectable via `document.querySelector('.CodeMirror').CodeMirror.getValue()`.
- Toolbox `Actions → Run query...` is a `label.d4-link-action` that `.click()` activates without needing pointer coordinates.
- `grok.dapi.queries.filter('name = "..."').first()` confirms saved entities round-trip to the server.
- Tree-node right-click via synthetic `contextmenu` on `.d4-tree-view-node` (not the label) — works first try.

### What did not work
- The ribbon Play button (`[name="icon-play"]`) needs multiple click attempts on the first run — the click handler is wired up shortly after editor mount, and a single `.click()` on either the icon or the wrapping `.d4-ribbon-item` immediately after the editor opens is silently dropped. Spec retries the click for up to 60s before failing.

### Suggestions for the platform
- Wire the Play-button click handler synchronously during editor construction so the first click reliably triggers the run (avoids retry-loop heuristics in tests).
- Add a `name=` attribute to the Toolbox `Run query...` link (`label.d4-link-action`) so specs can target it without a text match.

### Suggestions for the scenario
- Step numbering is `1, 3, 1, 8` — renumber monotonically `1..4` for clarity.
- Spell `Schemas` and `public` exactly as displayed in the tree (currently fine).
- Add an explicit "Save with a unique name" sub-step under step 4 — overwriting `products` collides with any pre-existing query of that name on shared servers.
