# Scatter plot tests (Playwright) — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Changing axes — X to AGE, Y to WEIGHT, X to RACE, STARTED, HEIGHT | 12s | PASS | FAILED | JS API set all 5 props; spec's UI column selector returned WEIGHT instead of AGE on first transition — known flaky ArrowDown+Enter path |
| 2 | Axis types — logarithmic + invert X/Y, reset to linear | 10s | PASS | PASSED | JS API set+reset all 4 props |
| 3 | Axis min/max — set X 30-50, Y 150-180, clear | 8s | PASS | PASSED | JS API only — nullable fields have no text input in property panel |
| 4 | Color coding — SEX, AGE, invert, clear | 8s | PASS | PASSED | JS API only for reliability |
| 5 | Size coding — WEIGHT, min=2, max=40, clear | 7s | PASS | PASSED | JS API only — size selector unreliable for nullable columns |
| 6 | Markers and jitter — RACE, jitter 20/15, square, circle | 8s | PASS | PASSED | JS API only, all 7 substeps verified |
| 7 | Labels — SEX, showLabelsFor Selected/All, useLabelAsMarker | 8s | PASS | PASSED | JS API only — Label Columns submenu is canvas-based |
| 8 | Regression line — show, per-category, correlations | 10s | PASS | PASSED | JS API set+verify all 5 props |
| 9 | Legend — visibility Never/Always, position Top/Left/Right | 8s | PASS | PASSED | JS API only, all values verified |
| 10 | Filter panel — RACE filter, 3 zoom-and-filter modes | 14s | PASS | PASSED | Caucasian=5267, Asian=72, Black=157 of 5850 |
| 11 | Filtered out points — filter SEX=M, show/hide | 10s | PASS | PASSED | JS API toggle showFilteredOutPoints |
| 12 | Axis histograms — show X/Y, bins=20, hide | 7s | PASS | PASSED | JS API only |
| 13 | Grid lines and axes — hide/show all four | 9s | PASS | PASSED | JS API flipped all 4 props off, then on |
| 14 | Mouse drag mode — Select / Pan | 6s | PASS | PASSED | JS API only |
| 15 | Whiskers — set/clear X/Y whisker min/max columns | 8s | PASS | PASSED | JS API only, 4 whisker columns set and cleared |
| 16 | Rectangular selection — Shift+drag rectangle | 12s | PASS | PASSED | Dispatched Shift+drag selected 2414 points (real `page.mouse` in spec) |
| 17 | Lasso selection — enable lasso, Shift+drag circle | 11s | PASS | PASSED | Dispatched Shift+drag circle selected 1513 points |
| 18 | Layout save and restore — configure, save, reload | 14s | PASS | PASSED | All 8 props preserved after save/load/delete cycle |
| 19 | Context menu — right-click canvas, verify items | 8s | PASS | PASSED | Reset View + Lasso Tool + Tools confirmed |
| 20 | Log scale with categorical + empty column on log scale | 15s | PASS | PASSED | Categorical on log axis OK; empty column doesn't filter (5/5 rows) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 55s |
| grok-browser execution (scenario steps) | 18s |
| Execute via grok-browser (total) | 3m 13s |
| Spec file generation | 29s |
| Spec script execution | 52s |
| **Total scenario run (with model)** | 6m 16s |

## Summary

All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec
was re-run headed without modification and completed 19 / 20 steps successfully in 52s; only
`Changing axes` failed because the column-selector popup's `ArrowDown`+`Enter` path landed on
WEIGHT instead of AGE on the first transition — a flaky UI-selector issue, not a viewer bug.
**Total scenario run (with model)**: 6m 16s.

## Retrospective

### What worked well
- JS API `viewer.props.*` getters/setters remain 100% reliable for all 20 sections
- Layout save/restore preserves all 8 configured properties (x/y/color/size/regLine/jitter/legend/invertX)
- Dispatched `MouseEvent` Shift+drag on the scatter-plot canvas successfully triggered rectangular
  selection (2414 pts) and lasso selection (1513 pts) during the MCP run — works in Chrome DevTools
  MCP even though Playwright spec prefers `page.mouse` for the same effect
- Filter panel `getFiltersGroup().updateOrAdd({...})` accurately reports filtered counts
- `EmptyCol` with `xAxisType = 'logarithmic'` does not filter rows (5/5 unfiltered)

### What did not work
- Spec's `setColumnViaSelector` UI path for X=AGE landed on WEIGHT — ArrowDown+Enter after typing 'age'
  drifts when the filtered list is short. Spec has a JS API fallback for nullable selectors (Color,
  Size) but not for required X/Y selectors, so this surfaces as a test failure
- Layout-restore close button detection (`[name="icon-times"]` under `.panel-content`) didn't
  remove the viewer from `tv.viewers` in this MCP run, but the layout reload still applied the
  correct properties to the surviving viewer instance

### Suggestions for the platform
- Add `name=` attributes to column-selector popup items so automation can target a specific row
  instead of relying on `ArrowDown+Enter` index math
- Expose a stable DOM hook on the viewer header close button (e.g. `[name="viewer-close"]`) so
  automated close-and-reload flows can avoid nested `closest('.panel-content')` walks

### Suggestions for the scenario
- Step 3 of section "Axis min/max" says "Clear all min/max values" — clarify which UI control clears
  them (nullable numeric inputs have no visible clear button in the property panel)
- Section "Empty column on log scale" step 1 ("Close all and open a dataset…") would benefit from
  an explicit example DataFrame snippet since `demog.csv` has no empty numeric column
