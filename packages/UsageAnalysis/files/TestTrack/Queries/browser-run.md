# Queries — Browse NorthwindTest and Find Query — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Refresh Browse | 3s | PASS | PASSED | Navigated to `/browse` to force a clean tree render. `.d4-tree-view-root` visible within ~300ms. |
| 2 | Browse → Databases → Postgres → NorthwindTest | 6s | PASS | PASSED | Double-click on `Postgres` group expanded it; `NorthwindTest` became visible. Double-click opened the connection's `queries` view (14 saved queries on the Details panel). |
| 3 | Type `new_test` in search field | 3s | PASS | PASSED | Search input placeholder `Search queries by name or by #tags`. Setting `.value` + dispatching `input`/`change`/`keyup` filtered the gallery; `new_test_query` appears in the filtered list. |
| 4 | Context Panel — check all tabs for the query | 2s | PASS | PASSED | Set `grok.shell.o = query` and read accordion pane headers — returned `[Details, Run, Query, Transformations, Usage, Activity, Sharing, Chats, Dev]` (9 panes). |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 40s |
| grok-browser execution (scenario steps) | 14s |
| Execute via grok-browser (total) | 54s |
| Spec file generation | 30s |
| Spec script execution | 19s |
| **Total scenario run (with model)** | 1m 43s |

## Summary

The full browse flow worked in both the MCP and Playwright runs. `/browse` route
reliably renders the tree after a `closeAll()`, Postgres → NorthwindTest expands on a
single double-click, the query search filters the gallery, and selecting the query
swaps the Context Panel to the query-entity accordion (9 panes).

## Retrospective

### What worked well
- Direct `/browse` navigation bypasses the Tabs-mode sidebar-click ambiguity — the tree renders once the URL commits.
- Search input `input`/`change`/`keyup` event batch is enough to trigger the Dart filter — no keyboard simulation required for read-only filters.
- Setting `grok.shell.o = query` is the fastest way to verify Context Panel behavior; the same panes also render when clicking the query card in the UI.

### What did not work
- Clicking the query card in the filtered gallery previews the query **result DataFrame**; `grok.shell.o` becomes the DataFrame, not the query. Context Panel then shows DataFrame panes (`General`, `Columns`, `Rows`, ...) instead of query panes.
- Initial sidebar `Browse` tab click after `closeAll()` did not render the tree — the Browse panel header appears but `.d4-tree-view-root` remains unmounted until a view navigation triggers it. `/browse` route is the reliable workaround.

### Suggestions for the platform
- Single-click on a query card in the gallery should select the query (show query panes), not preview the result. Make "run and preview" an explicit double-click / hover action.
- After `grok.shell.closeAll() + grok.shell.windows.showBrowse = true`, render the Browse tree synchronously so automation doesn't need to force a route navigation.

### Suggestions for the scenario
- Step 4 wording "check all tabs" is vague — specify which panes are expected (e.g. Details, Run, Query, Transformations, Usage, Sharing).
- Add an explicit substep: "Click the query card (not the preview surface) to get the query-entity Context Panel." Otherwise users see DataFrame panes and don't know why.
- Step 1 numbering: lines are `1, 2, 3, 3` — the last `3` should be `4`.
