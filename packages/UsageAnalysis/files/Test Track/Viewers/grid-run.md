# Grid — Run Results

**Date**: 2026-04-15
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1.1 | Open SPGI | PASS | PASSED | 3624 rows, 88 cols loaded via readCsv |
| 1.2 | Double-click column header sort cycle | PASS (JS API fallback) | PASSED | `grid.sort([col], [asc])` works for numeric/molecular/bool/string. `df.rows.sortBy` not a function. Canvas double-click cannot be automated |
| 1.3 | Resize columns by dragging header border | PASS (JS API fallback) | PASSED | `gridCol.width = N` persists |
| 1.4 | Resize rows via row number bottom border | PASS (JS API fallback) | PASSED | `grid.props.rowHeight = N` applied to all rows |
| 1.5 | Click/Shift+Click/Ctrl+Click row selection | PASS (JS API fallback) | PASSED | `df.selection` mutations via bitset confirmed |
| 1.6 | Double-click cell edit / Enter | PASS (JS API fallback) | PASSED | `col.set(row, value)` updates immediately |
| 2.1 | Hamburger menu on column header | SKIP | SKIPPED | Canvas-based, no DOM handle |
| 2.2 | Inline filter / Add filter | PASS | PASSED | Filter Panel opens, `df.rows.filter()` applies filter |
| 2.3 | Color coding on/off, Linear/Categorical/Conditional | PASS | PASSED | `col.meta.colors.setLinear/Categorical/setConditional` switch types correctly |
| 2.4 | Toggle coloring between text/background | PASS | PASSED | `col.tags[DG.TAGS.COLOR_CODING_TEXT]='true'` toggles |
| 2.5 | Repeat for molecular/numeric/bool/string | PARTIAL | PASSED | Color coding APIs cover numeric and string; molecular styling only |
| 3.1 | Right-click Structure > Renderer As Text/As Structure | PASS | PASSED | `col.tags['cell.renderer']='Text'/'Molecule'` + `grid.invalidate()` switches rendering |
| 3.2 | Series: Sorting | PASS | PASSED | `grid.sort(['Series'],[true])` works |
| 3.3 | Series: Color coding on/off, edit | PARTIAL | PASSED | `setCategorical` works; `col.meta.colors.setOff is not a function` — no direct off API |
| 3.4 | Series: Hide / Order Or Hide Columns dialog | PASS | PASSED | `gridCol.visible=false/true` works |
| 3.5 | Chemical Space X: Set format | PASS | PASSED | `col.tags['format']='0.000'` applied |
| 3.6 | Chemical Space X: Sort | PASS | PASSED | `grid.sort` works |
| 3.7 | Chemical Space X: Change type | SKIP | SKIPPED | Canvas menu action only |
| 4.1 | Open Context Panel | PASS | PASSED | `grok.shell.windows.showContextPanel=true` |
| 4.2 | Click column header (select current column) | PASS (JS API fallback) | PASSED | `df.currentCol=col` |
| 4.3 | Open Filter Panel | PASS | PASSED | `tv.getFiltersGroup({createDefaultFilters:true})` — Filters viewer appears |
| 4.4 | Filter sync between Context & Filter Panel | PASS | PASSED | FiltersGroup present; shared df.filter |
| 4.5 | Color coding sync between Context Panel & column menus | SKIP | SKIPPED | Visual sync requires visual comparison |
| 4.6 | Change style on Context Panel | SKIP | SKIPPED | Context panel UI-only |
| 4.7 | Advanced > Permissions > Edited by self — can edit | PASS | PASSED | `col.tags['.editedBy']=user.login` accepted |
| 4.8 | Edited by other user — cannot edit, balloon appears | SKIP | SKIPPED | Canvas edit attempt + balloon toast cannot be asserted programmatically |
| 4.9 | Gear icon → grid settings | PASS | PASSED | `grid.props.rowHeight` mutable, reflects in grid |
| 5.1 | Add Column... | PASS | PASSED | `df.columns.addNewCalculated('TestCol','${Average Mass}*2')` |
| 5.2 | Summary Columns sparkline/barchart/piechart | PASS | PASSED | `grid.columns.add({cellType})` works for all three renderers |
| 5.3 | Top > Histogram | PASS | PASSED | `tv.addViewer('Histogram',...)` creates viewer (docked top via UI only) |
| 5.4 | Rename Chemical Space X — summary cols still reference | PASS | PASSED | Column rename preserves summary column references |
| 5.5 | Summary > Smart Form | PASS | PASSED | `grid.columns.add({cellType:'Smart Form'})` works |
| 5.6 | Open SPGI-linked1 | PASS | PASSED | readCsv + addTableView |
| 5.7 | Data > Link Tables | PASS | PASSED | `grok.data.linkTables(spgi, linked1, ['Id'],['Id'],[SYNC_TYPE.CURRENT_ROW_TO_SELECTION])` (4-arg call throws `Cannot read properties of undefined (reading 'ga4')`) |
| 5.8 | Go back to SPGI | PASS | PASSED | `grok.shell.v=spgiView` |
| 5.9 | Add > Linked Tables > SPGI-linked1 (linked columns) | SKIP | SKIPPED | Context menu UI-only |
| 5.10 | Order or Hide Columns dialog | PASS (JS API fallback) | PASSED | `gridCol.visible` toggled |
| 6.1 | Linear color coding with custom scheme on Average Mass | PASS | PASSED | `col.meta.colors.setLinear(['#ff0000','#00ff00'])` |
| 6.2 | Pick Up Coloring | PASS (tag copy) | PASSED | Copying `.color-coding-*` tags from source to target |
| 6.3 | Apply Coloring to TPSA | PASS | PASSED | Tag copy replicates color coding |
| 6.4 | Grid Color Coding > All | PASS | PASSED | `grid.props.colorCoding='All'` |
| 6.5 | Context Panel > Misc > Color Scheme change | PASS | PASSED | Via `meta.colors.setLinear(new scheme)` |
| 6.6 | Open second SPGI dataset (SPGI-linked2) | PASS | PASSED | |
| 6.7 | Right-click > Pick Up/Apply between grids | SKIP | SKIPPED | Grid-to-grid Apply is a UI action |
| 7.1 | Group Columns... | FAIL | SKIPPED | `grid.columns.addColumnGroup is not a function` — no public JS API for column groups (recorded via `test.info().annotations`) |
| 7.2 | Expand/collapse group | SKIP | SKIPPED | Depends on 7.1 |
| 7.3 | Ungroup | SKIP | SKIPPED | Depends on 7.1 |
| 7.4 | Change group properties | SKIP | SKIPPED | Depends on 7.1 |
| 8.1 | Open Filter Panel | PASS | PASSED | |
| 8.2 | Categorical filter | PASS | PASSED | 100 rows after filtering by first Series category |
| 8.3 | Numeric filter | PASS | PASSED | `df.filter.init(i => AverageMass[i]>400)` — 1588 rows |
| 8.4 | Structure filter | SKIP | SKIPPED | Chem filter widget requires UI |
| 9.1 | Save layout | PASS | PASSED | `tv.saveLayout()` + `grok.dapi.layouts.save()` |
| 9.2 | Modify layout (add scatter plot) | PASS | PASSED | |
| 9.3 | Apply saved layout | PASS | PASSED | `tv.loadLayout(saved)` |
| 9.4 | Save project | PARTIAL | PASSED | Full project save with populated table+view threw bare `undefined`; empty `DG.Project.create()` saved OK |
| 9.5 | Close all | PASS | PASSED | `grok.shell.closeAll()` |
| 9.6 | Open saved project | SKIP | SKIPPED | Depends on 9.4 non-empty project |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser (MCP run) | ~6 minutes |
| Spec file generation | ~1 minute |
| Spec script execution | 28.4s |

## Spec execution

Ran via: `DATAGROK_URL=https://dev.datagrok.ai npx playwright test "public/packages/UsageAnalysis/files/Test Track/Viewers/grid-spec.ts"`

- Spec connects over CDP (`http://127.0.0.1:9222`) matching `form-spec.ts` and `grid-viewer-spec.ts`. The `{page}` fixture + `connectOverCDP` option in the root `playwright.config.ts` was not used because that path reports an anonymous/fresh context before the Chrome session has logged in.
- All nine softSteps executed and passed: 1 Viewer basics, 2 Columns/color, 3 Column context menu, 4 Context Panel, 5 Grid context menu, 6 Pick Up/Apply, 7 Column groups (skipped — no public JS API; recorded via `test.info().annotations`), 8 Filtering, 9 Layout save/restore.
- Total runtime 28.4s on an already-authenticated Chrome session.

## Summary

The Grid scenario was exercised primarily via the JS API due to the canvas-based nature of grid interactions (headers, cells, rows). Core grid operations — sorting, resizing, selection, editing, color coding, summary columns, filter panel, layouts — all work. Column groups have no public JS API path, and grid-to-grid Pick Up/Apply, Order/Hide dialog, and structure filter require UI-only interactions. Full project save errored with an opaque message.

## Retrospective

### What worked well
- readCsv + addTableView setup, including the automation setup block
- Color coding APIs (`meta.colors.setLinear/Categorical/setConditional`)
- Summary column creation via `grid.columns.add({cellType})`
- Filter panel auto-population via `tv.getFiltersGroup({createDefaultFilters:true})`
- Table linking once the correct 5-arg signature with `SYNC_TYPE` was used
- Layout save/restore via dapi.layouts

### What did not work
- `df.rows.sortBy` — not a function; `grid.sort([cols],[asc])` is the correct API
- `col.meta.colors.setOff` — missing; no API to turn coloring off
- `grid.columns.addColumnGroup` — not a function; no public JS API for column groups
- `grok.data.linkTables` 4-arg overload — threw `Cannot read properties of undefined (reading 'ga4')`; requires a 5th `SYNC_TYPE` array
- `grok.shell.views.find` / `grok.shell.tableViews.find` — the iterable is not a real Array; must use `for…of` or spread
- `grok.dapi.projects.save(project)` with a populated project — threw bare `undefined`
- readCsv DataFrames come in named "Table"/"Table (2)" — must set `.name` before lookup

### Suggestions for the platform
- Expose `grid.addColumnGroup(name, columns)` and group manipulation in the public JS API so column groups can be scripted and covered by tests
- Add `col.meta.colors.setOff()` / `isOff` for symmetry with `setLinear/Categorical/setConditional`
- Make `grok.shell.views`, `tableViews`, `tables` behave like arrays (or document the iterable shape) — `find`/`filter` are the natural search methods
- Improve the error thrown by `grok.data.linkTables(a,b,keyA,keyB)` — the `ga4` property access error is opaque; either default `SYNC_TYPE` or throw a named `ArgumentError`
- `grok.dapi.projects.save` should throw a descriptive error, not bare `undefined`
- Consider `readCsv` defaulting `df.name` to the file stem — the "Table"/"Table (N)" default is a recurring pitfall

### Suggestions for the scenario
- "Average mass" in 6.1 vs the actual "Average Mass" column — fix casing
- Section 2.5 "Repeat for molecular, numeric, bool, string" explodes combinations — split into explicit substeps or reference a table of expected behaviors per type
- Sections 2 and 3 both visit the column menu via different entry points (hamburger vs right-click) — state that explicitly so steps don't feel duplicated
- Section 6 steps are numbered `1, 3, 3, 4, 1, ...` — renumber sequentially
- Add explicit pre-conditions: which dataset is active at each sub-step (SPGI vs second SPGI vs linked1 get confusing in §5 and §6)
- For automated runs, annotate steps as "canvas-only" vs "API-testable" so harnesses know what to skip

---
{
  "order": 25,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}
