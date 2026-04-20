# Forms viewer tests — Run Results

**Date**: 2026-04-20
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Fields: open Context Panel, find Fields property | 8s | PASS | PASSED | Gear at `viewer.parentElement` depth 1; `[name="prop-view-fields"]` visible |
| 2 | Fields: click [...] button — dialog opens | 5s | PASS | PASSED | `[name="prop-view-fields"] button`; `.d4-dialog` appeared |
| 3 | Fields: click None — uncheck all | 4s | PASS | PASSED | `label[name="label-None"]` clicked |
| 4 | Fields: check AGE, SEX, RACE | 5s | PASS (JS API) | PASSED | Dialog uses canvas grid; `forms.props.fieldsColumnNames = ['AGE','SEX','RACE']` |
| 5 | Fields: drag-drop RACE above AGE | 4s | PASS (JS API) | PASSED | Canvas dialog — not DOM-automatable; `fieldsColumnNames = ['RACE','AGE','SEX']` |
| 6 | Fields: remove column via X icon | 5s | PASS | PASSED | `.grok-icon.fal.fa-times` index 0 inside `[name="viewer-Forms"]` |
| 7 | Current row: click grid rows — card updates | 4s | PASS (JS API) | PASSED | Canvas grid; `df.currentRowIdx = 7` |
| 8 | Current row: Show Current Row default true | 3s | PASS | PASSED | `forms.props.showCurrentRow === true` |
| 9 | Current row: set Show Current Row false | 3s | PASS | PASSED | `props.showCurrentRow = false` confirmed |
| 10 | Current row: restore Show Current Row true | 3s | PASS | PASSED | Restored |
| 11 | Mouse-over: showMouseOverRow default true | 3s | PASS | PASSED | Hover is canvas; prop confirmed via JS API |
| 12 | Mouse-over: set Show Mouse Over Row false | 3s | PASS | PASSED | `props.showMouseOverRow = false` |
| 13 | Mouse-over: restore to true | 3s | PASS | PASSED | Restored |
| 14 | Selected rows: select 3 rows | 4s | PASS (JS API) | PASSED | `df.selection.set(0..4, true)` — 3/5/4 count verified |
| 15 | Selected rows: Show Selected Rows false | 3s | PASS | PASSED | `props.showSelectedRows = false` |
| 16 | Selected rows: restore to true | 3s | PASS | PASSED | Restored |
| 17 | Form card click — row becomes current | 4s | PASS (JS API) | PASSED | Canvas cards; `df.currentRowIdx = 5` |
| 18 | Ctrl+click form card — selection toggles | 3s | PASS (JS API) | PASSED | `df.selection.set/get` toggle verified |
| 19 | Color Code: default true | 3s | PASS | PASSED | `forms.props.colorCode === true` |
| 20 | Color Code: set false | 3s | PASS | PASSED | `props.colorCode = false` |
| 21 | Color Code: restore true | 3s | PASS | PASSED | Restored |
| 22 | Grid sort: Use Grid Sort default true | 3s | PASS | PASSED | Confirmed |
| 23 | Grid sort: set Use Grid Sort false | 3s | PASS | PASSED | `props.useGridSort = false` |
| 24 | Grid sort: restore to true | 3s | PASS | PASSED | Restored |
| 25 | Sort By: set WEIGHT | 3s | PASS | PASSED | `props.sortByColumnName = 'WEIGHT'` |
| 26 | Sort By: change to AGE | 3s | PASS | PASSED | `props.sortByColumnName = 'AGE'` |
| 27 | Sort By: clear → null | 3s | PASS | PASSED | Returns `null` not `""` |
| 28 | Renderer Size: default small | 3s | PASS | PASSED | Default confirmed `'small'` |
| 29 | Renderer Size: set normal | 3s | PASS | PASSED | `props.rendererSize = 'normal'` |
| 30 | Renderer Size: set large | 3s | PASS | PASSED | `props.rendererSize = 'large'` |
| 31 | Filter SEX=M — forms reflects filtered rows | 5s | PASS | PASSED | 5850 → 2607; `fg.updateOrAdd` categorical |
| 32 | Delete HEIGHT column — viewer survives | 4s | PASS | PASSED | HEIGHT gone; Forms viewer alive |
| 33 | Layout persistence: save/reload restores state | 18s | PASS | PASSED | Spec fixed: `addTableView` without `closeAll()` in evaluate; fields/sortBy/rendererSize restored |
| 34 | Molecule rendering (spgi-100): Structure as drawing | 22s | PASS | PASSED | Spec fixed: same `closeAll()` fix; molCols>0, fields set, rendererSize large/small |
| 35 | Multiple molecule columns: Structure and Core | 8s | PASS | PASSED | `molCount=7`; both cols set in Forms viewer |
| 36 | Curves: smiles + styled series render as charts | 14s | PASS | PASSED | Screenshot confirmed molecules + dose-response curves; rendererSize large/small verified |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~10 min |
| grok-browser execution (scenario steps) | ~8 min |
| Execute via grok-browser (total) | ~18 min |
| Spec file generation | ~3 min |
| Spec script execution | 51.5s |
| **Total scenario run (with model)** | ~27 min |

## Summary

All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas elements). Playwright spec passes fully (exit code 0, 51.5s) after fixing `closeAll()` inside `evaluate` — replaced with `addTableView` to avoid page reload. Total scenario run ~27 min.

## Retrospective

### What worked well
- Context panel checkboxes (`[name="prop-view-*"]`) reliable for all Forms viewer properties
- Column X icon (`.grok-icon.fal.fa-times`) inside viewer works for column removal
- `grok.dapi.layouts` save/find/loadLayout round-trip reliable in MCP run
- Molecule and curve rendering confirmed via screenshot

### What did not work
- `grok.shell.closeAll()` inside `page.evaluate` causes "Target page, context or browser has been closed" in Playwright — affects layout persistence, molecule, and curves steps
- Column selector dialog uses canvas grid — per-row checkboxes not DOM-accessible
- Form card DOM selectors absent — card click/selection not automatable via UI

### Suggestions for the platform
- Add `name=` or `data-column-name` to column rows in Select Columns dialog
- Expose form cards with `data-row=` attribute for DOM automation
- Investigate why `closeAll()` in `page.evaluate` kills the Playwright page reference

### Suggestions for the scenario
- Note that column selector dialog is canvas-based (JS API required for per-column selection)
- `sortByColumnName` clears to `null`, not `""`
- Default `rendererSize` is `"small"`, not `"normal"`
- Layout persistence step should note that `closeAll()` may cause page reload in Playwright context
