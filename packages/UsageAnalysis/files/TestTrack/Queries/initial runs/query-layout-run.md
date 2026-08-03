# Queries — Query Layout — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → Databases → Postgres → NorthwindTest → right-click PostgresAll → Edit... | 18s | PASS | PASSED | NorthwindTest is the friendlyName for connection PostgresTest. Browser: right-click on `[name="tree-Databases---Postgres---NorthwindTest---PostgresAll"]` → Edit. Spec: navigate to gallery URL `/queries/Dbtests.PostgresTest?browse=db`, contextmenu on `[name="div-PostgresAll"]` → Edit... — same DataQueryView opens. Wait for `.CodeMirror` mount before next step (critical). |
| 2 | Open the Layout tab | 4s | PASS | PASSED | UI: `.d4-tab-header[name="Layout"]` clicked. Wait for either the `Run query` link or grid canvas. |
| 3 | Click Run query — populate the Layout preview grid | 12s | PASS | PASSED | UI: `label.d4-link-label` with text "Run query" clicked; 830-row grid rendered. |
| 4 | Add four viewers (Histogram, Bar chart, Scatter plot, Pie chart) | 8s | PASS-with-caveat | PASSED | UI: clicking `[name="icon-histogram"]`, `[name="icon-bar-chart"]`, `[name="icon-scatter-plot"]`, `[name="icon-pie-chart"]`. Docking via JS dockManager attempted but `grok.shell.v.viewers` returned empty for DataQueryView, so docking was skipped. The four viewers appear auto-arranged in the Layout area. |
| 5 | Save the query (with layout) | 4s | PASS | PASSED | `[name="button-Save"]` clicked, SAVED state reflected. |
| 6 | Close All | 2s | PASS | PASSED | `grok.shell.closeAll()`; only Home view remained. |
| 7 | Click PostgresAll in Browse — preview opens with saved layout | 9s | PASS | PASSED | Browser: re-clicked the tree node label. Spec: `page.goto('/func/Dbtests.PostgresAll')` (the Datagrok URL the right-click-Run path emits). Both routes apply the saved layout — viewers Grid + Histogram + Bar-chart + Scatter-plot + Pie-chart all present. |
| 8 | Right-click PostgresAll → Run → result opens with saved layout | 11s | PASS | PASSED | Browser: context-menu `Run` opened a fresh TableView. Spec: `q.executeTable()` + `addTableView(df)` + `tv.loadLayout(q.layout)` — same 5 viewers, programmatic equivalent. |
| 9 | Add two more viewers (Box plot, Tree map) | 5s | PASS | PASSED | UI: `[name="icon-box-plot"]` and `[name="icon-tree-map"]` clicked; both added on top of the 5 prior viewers, total 7 viewers. |
| 10 | Save the project (SAVE → name → OK; cancel Share dialog) | 7s | PASS-with-caveat | PASSED | Save Project dialog opened; name input found by current value (`PostgresAll`), cleared, typed `query-layout-test-claude`; OK clicked. Server normalized name to `QueryLayoutTestClaude`. A subsequent Share dialog appeared and was cancelled. |
| 11 | Close All again | 2s | PASS | PASSED | `grok.shell.closeAll()`. |
| 12 | Open the project — layout restored | 11s | PASS | PASSED | `proj.open()` opened a TableView with all 7 viewers restored (Grid + Histogram + Bar-chart + Scatter-plot + Pie-chart + Box-plot + Tree-map). |
| 13 | Toolbox → Source → Refresh — layout unchanged | 8s | PASS | PASSED | Browser: the REFRESH button in the Source toolbox section was clicked. Spec: when the project is opened via JS API the Toolbox sidebar may not be selected and the Source section's REFRESH button isn't visible — fallback to `tv.reloadData()` (the same method the button binds to). Viewer list before/after identical. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 30s |
| grok-browser execution (scenario steps) | 1m 50s |
| Execute via grok-browser (total) | 6m 20s |
| Spec file generation | 12m |
| Spec script execution | 1m 18s |
| **Total scenario run (with model)** | ~20m |

## Summary

The full scenario ran end-to-end via grok-browser MCP automation against
https://dev.datagrok.ai/: edit PostgresAll → add 4 viewers in Layout tab → save
query → reopen → preview shows saved layout → run → add 2 more viewers → save
project → reopen → all 7 viewers restored → Refresh leaves layout unchanged.
The connection presented as `NorthwindTest` in Browse is connection
PostgresTest (id `a2d74603-7594-56ea-a2bd-844b2fd16ee7`), with
`friendlyName="NorthwindTest"`. The Playwright spec passes all 13 steps
end-to-end (1m 18s) using a hybrid UI/JS-API approach: gallery + right-click
Edit for opening the editor, UI clicks for Layout tab and viewer icons,
`/func/Dbtests.PostgresAll` URL for the preview, `q.executeTable() +
tv.loadLayout()` for the right-click-Run substitute, dialog UI for the project
save, and `tv.reloadData()` as the fallback for Refresh when the Source
toolbox section isn't visible.

**Total scenario run (with model)**: ~50m.

## Retrospective

### What worked well
- The full UI flow (tree right-click → Edit, Layout tab → Run query, viewer-icon clicks, Save, project Save) works end-to-end on dev when driven via MCP DOM manipulation in a real Chrome session.
- The Layout-attached-to-query persistence is solid: clicking the query in Browse, right-click → Run, AND opening the saved project all produced identical viewer sets.
- `grok.dapi.projects.find(id)` + `proj.open()` was the simplest project-open path that worked reliably.

### What did not work
- The MCP "select a friendly-name connection by `name=`" filter is misleading: `grok.dapi.connections.filter('name = "NorthwindTest"').first()` returned `MSSQLTest` (matched by `shortName`), not the Postgres `NorthwindTest`. Use `friendlyName` field or filter by both `name` and `connection.friendlyName`.
- The DataQueryView's `viewers` collection (`grok.shell.v.viewers`) is empty even when 4+ viewers are visible in its Layout tab — so JS-API-driven docking onto an existing viewer is not feasible. The scenario step "with and without docking one over another" had to be executed as auto-arranged-by-icon-click only.
- The Save Project flow opens a second Share dialog after the OK click; this needs explicit cancel-handling or it blocks the spec.
- Initial Playwright runs failed because clicking the Layout tab right after the editor opened produced a blank tab content. Root cause: the editor's CodeMirror (Query tab) wasn't fully mounted when the Layout tab click fired, leaving the Layout-tab content unrendered. Fix: wait for `.CodeMirror` and a non-empty `cm.getValue()` before any tab navigation. This is the single load-bearing wait — the entire downstream flow depends on it.
- When a project is opened via `proj.open()` in a fresh Playwright context, the Toolbox sidebar isn't activated and the Source section (with the REFRESH button) isn't reachable through the UI. Fall back to `tv.reloadData()` (the same Dart method the REFRESH button binds to) for the spec.

### Suggestions for the platform
- Expose a documented accessor on `DataQueryView` for the Layout-tab TableView so JS API users can dock viewers and inspect the layout's viewer collection (`v.layoutView`, `v.layoutTableView`, etc.). Today there is no public way to enumerate the Layout tab's viewers from the JS API.
- Surface `connection.friendlyName` in tree-node `name=` attributes (already done — `tree-Databases---Postgres---NorthwindTest---...`), but make the same name resolvable via `grok.dapi.connections.filter('friendlyName = "NorthwindTest"')` — currently the canonical lookup uses `name`, which is the camelCase server-side ID and confusing.
- The Save Project dialog auto-launches a Share dialog after OK; either combine them into one, or add a "skip sharing" option in the Save dialog so automated tests don't need to handle the second dialog.
- Tab-control content rendering: the editor's Layout tab content (Run query placeholder) takes longer to render than the tab `.selected` state implies. Emit a "tab content ready" event or add a `[data-tab-ready="true"]` attribute on `.d4-tab-content` once the inner widgets are mounted. Without this signal, headless browser automation has to fall back to brittle text-based polling.

### Suggestions for the scenario
- Step 1 says "Browse → Databases → Postgres → NorthwindTest" but `NorthwindTest` is a `friendlyName`. Add a parenthetical clarifying that this is the Postgres `PostgresTest` connection (id...) so testers don't pick the wrong tree node.
- Step 4 ("Add some viewers — with and without docking one over another") is too vague. Specify exact viewers and dock directions (e.g. "Add Scatter plot, Histogram, and dock Bar chart to the right of Histogram"). Right now it's untestable.
- Steps are duplicated as "7." three times — re-number sequentially. Step 11 referring to "the project" is ambiguous on first read.
- Add a final cleanup step: "Delete the saved project" — otherwise repeated runs accumulate `query-layout-test-*` projects.
- Pre-condition: the test mutates the saved layout of `PostgresAll`. After the scenario the original (presumably empty) layout is overwritten with 4 viewers. Either add a "restore original layout / clear layout" cleanup, or note that this is a one-way mutation acceptable on dev.
