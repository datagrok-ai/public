# Scatter plot tests (Playwright) — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Changing axes — X to AGE, Y to WEIGHT, X to RACE, STARTED, HEIGHT | PASS | 5s | PASSED | UI column selector: mousedown `.d4-column-selector-column`, type, ArrowDown, Enter |
| 2 | Axis types — logarithmic + invert X/Y via context menu | PASS | 6s | PASSED | Real right-click context menu for log/invert. JS API fallback for dropdown reset |
| 3 | Axis min/max — set X 30-50, Y 150-180, clear | PASS | 2s | PASSED | JS API only — nullable fields have no text input in property panel |
| 4 | Color coding — SEX, AGE, invert, clear | PASS | 4s | PASSED | UI selector with auto-fallback to JS API when popup search fails |
| 5 | Size coding — WEIGHT, min=2, max=40, clear | PASS | 2s | PASSED | JS API only — size selector unreliable for nullable columns |
| 6 | Markers and jitter — RACE, jitter 20/15, square, circle | PASS | 2s | PASSED | JS API only, all 7 substeps verified |
| 7 | Labels — SEX, showLabelsFor Selected/All, useLabelAsMarker | PASS | 2s | PASSED | JS API only — Label Columns submenu is canvas-based |
| 8 | Regression line — show via ctx menu, per-category, correlations | PASS | 4s | PASSED | Real right-click for Show Regression Line toggle |
| 9 | Legend — visibility Never/Always, position Top/Left/Right | PASS | 2s | PASSED | JS API only, all values verified |
| 10 | Filter panel — RACE filter, 3 zoom-and-filter modes | PASS | 4s | PASSED | Caucasian=5267, Asian=72, Black=157 of 5850 |
| 11 | Filtered out points — filter SEX=M, show/hide | PASS | 3s | PASSED | JS API toggle showFilteredOutPoints |
| 12 | Axis histograms — show X/Y, bins=20, hide | PASS | 2s | PASSED | JS API only |
| 13 | Grid lines and axes — hide/show all four | PASS | 8s | PASSED | 8 real right-click context menu clicks to toggle off then on |
| 14 | Mouse drag mode — Select / Pan | PASS | 1s | PASSED | JS API only |
| 15 | Whiskers — set/clear X/Y whisker min/max columns | PASS | 2s | PASSED | JS API only, 4 whisker columns set and cleared |
| 16 | Rectangular selection — Shift+drag rectangle | PASS | 4s | PASSED | Real `page.mouse` Shift+drag, >0 points selected |
| 17 | Lasso selection — enable lasso, Shift+drag circle | PASS | 4s | PASSED | Real `page.mouse` circular drag, >0 points selected |
| 18 | Layout save and restore — configure, save, close, restore | PASS | 8s | PASSED | All 8 properties preserved after save/load/delete cycle |
| 19 | Context menu — right-click canvas, verify items | PASS | 2s | PASSED | Real right-click, Reset View + Lasso Tool + Tools confirmed |
| 20 | Log scale with categorical + empty column on log scale | PASS | 5s | PASSED | Categorical on log axis OK; empty column doesn't filter (5/5) |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~8 min |
| Spec file generation | ~5s |
| Spec script execution | 1.1 min |

## Summary

All 20 test sections passed on dev.datagrok.ai. Playwright spec executed in 1.1 min (headed mode).
Context menu, rectangular selection, and lasso selection use real Playwright mouse events — visually
verified in browser. Color selector has auto-fallback to JS API when popup search fails.

## Retrospective

### What worked well
- Column selector popup: `mousedown` on `.d4-column-selector-column` + type to search + `ArrowDown` + `Enter`
- Context menu via `page.mouse.click(x, y, {button: 'right'})` — reliable, visually confirmed
- Rectangular and lasso selection via real `page.mouse` Shift+drag — points selected correctly
- All JS API `viewer.props.*` getters/setters reliable
- Layout save/restore preserved all 8 configured properties
- Setup wait `grok.shell.settings != null` handles slow Dart interop loading in headed mode

### What did not work
- Color/Size selector popup search intermittently fails — auto-fallback to JS API added
- Property panel dropdowns (X Axis Type): custom `.d4-combo-popup` doesn't apply via click — JS API fallback needed
- Dispatched `MouseEvent` via `evaluate` doesn't trigger scatter plot selection — must use real `page.mouse`

### Suggestions for the platform
- Add `name=` attributes to context menu items for reliable automation
- Expose property panel dropdown controls as standard `<select>` elements

### Suggestions for the scenario
- Clarify "click the Color selector" means clicking `.d4-column-selector-column`, not the `+` icon
- Consider splitting 20 sections into smaller independent test files for better isolation
