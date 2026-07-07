# Line chart tests (Playwright) — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Chart types: Area, Stacked Area, Stacked Bar, Line | PASS | 4s | PASSED | Context menu `div-Chart-Type---*` selectors work |
| 2 | Whiskers: Avg\|Min,Max; Med\|Q1,Q3; +/-StDev; +/-StError; None | PASS | 3s | PASSED | JS API `whiskersType` prop |
| 3 | Markers: visibility, shape, size, opacity, columns | PASS | 5s | PASSED | Context menu for visibility, JS API for shape/size props |
| 4 | Axis config: log, invert, grid lines, label orient | PASS | 4s | PASSED | JS API props, all 12 sub-steps |
| 5 | Interpolation: Spline, tension 1.0, None | PASS | 2s | PASSED | JS API `interpolation`, `splineTension` |
| 6 | Aggregation types: avg, min, max, med, sum, stdev | PASS | 3s | PASSED | JS API `aggrType` |
| 7 | Left panel histogram: Histogram, None | PASS | 2s | PASSED | JS API `leftPanel` |
| 8 | Controls visibility: 6 controls off then on | PASS | 3s | PASSED | JS API boolean props |
| 9 | Y global scale: multi axis + global scale toggle | PASS | 3s | PASSED | JS API `yGlobalScale` |
| 10 | Split by column: SEX, RACE, both, clear | PASS | 3s | PASSED | `splitColumnName` and `splitColumnNames` |
| 11 | Multi-axis mode: 3 Y cols, toggle on/off | PASS | 2s | PASSED | JS API `multiAxis` |
| 12 | Title and description: show, set, position, vis | PASS | 3s | PASSED | JS API string props |
| 13 | Custom axis min/max: set and clear | PASS | 3s | PASSED | `xMin`, `xMax`, `yMin`, `yMax` |
| 14 | Date/time X axis: STARTED, Year/Month/DoW/None | PASS | 3s | PASSED | JS API `xMap` |
| 15 | Line styling: width, transparency, coloring type | PASS | 2s | PASSED | JS API line props |
| 16 | Axis tickmarks: MinMax, Auto for X and Y | PASS | 2s | PASSED | JS API `xAxisTickmarksMode` |
| 17 | Overview chart: Line, Area, Stacked Bar, None | PASS | 4s | PASSED | Context menu `div-Overview---*` |
| 18 | Legend: visibility, position (L/T/B/R) | PASS | 3s | PASSED | JS API `legendVisibility`, `legendPosition` |
| 19 | Axes follow filter: default true, toggle, filter | PASS | 4s | PASSED | `axesFollowFilter` + histogram filter |
| 20 | Context menu items present | PASS | 2s | PASSED | All 8 expected items found |
| 21 | Layout save and restore | PASS | 8s | PASSED | `grok.dapi.layouts` save/find/load/delete |
| 22 | Table switching and row source (SPGI) | PASS | 10s | PASSED | Table switch, RowSource Selected/Filtered |
| 23 | Filter expression and collaborative filtering | PASS | 4s | PASSED | In-viewer filter + panel filter collaborative |
| 24 | Split and Y-columns sync with Context Panel | PASS | 5s | PASSED | Split add/remove, Y cols add/remove, multi-axis |
| 25 | Selection checkboxes | PASS | 6s | PASSED | Shift+drag selection, layout save/restore |
| 26 | Data panel checkboxes | PASS | 5s | PASSED | Layout save/restore preserves settings |
| 27 | GROK-17835 regression (SPGI) | PASS | 6s | PASSED | Multi axis + split + hover = no errors |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 50s |
| grok-browser execution (scenario steps) | 7s |
| Execute via grok-browser (total) | 57s |
| Spec file generation | 8s |
| Spec script execution | 1m 43s |
| **Total scenario run (with model)** | 2m 48s |

All rows are full-phase wall-clock (incl. model thinking and retries). The two `scenario steps`
rows sum to `Execute via grok-browser (total)`. Spec was not rewritten (existing spec is correct
per the user's instruction); the 8s reflects verification/sanity-pass only. Spec execution used
the existing `playwright.config.files-and-sharing.ts` (testMatch: `/-spec\.ts$/`).

## Summary

All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout save/restore, table switching, and GROK-17835 regression all work correctly. Formula lines and Legend filtering sections moved to line-chart-tests-ui.md (require complex dialog/canvas interaction).

## Retrospective

### What worked well
- JS API `lc.props.*` for property manipulation is reliable and fast
- Context menu `name=` attribute selectors (`div-Chart-Type---*`, `div-Overview---*`, `div-Markers---*`) work consistently
- Layout save/restore via `grok.dapi.layouts` is stable with full property round-trip
- Table switching and row source filtering work correctly
- In-viewer `filter` property for expression filtering works

### What did not work
- `showSelectedOnly` property not directly accessible on Line chart viewer

### Suggestions for the platform
- Expose legend click filtering via JS API (e.g. `lc.legendFilter(['category1'])`)

### Suggestions for the scenario
- Selection checkboxes section references "Context Panel > Selection" but doesn't specify exact checkbox names
