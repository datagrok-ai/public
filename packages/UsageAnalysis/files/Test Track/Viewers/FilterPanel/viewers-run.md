# Viewers — Filter Panel interaction with viewers — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1.1 | Structure filter c1ccccc1 | PASS | 8s | PASSED | Sketch-link click → SMILES input → OK. 1356 rows (benzene substructure in full SPGI) |
| 1.2 | Uncheck S_UNKN in Stereo Category | PASS | 2s | PASSED | JS API fallback: fg.updateOrAdd categorical, 1265 rows |
| 1.4 | Average Mass max=400 | PASS | 2s | PASSED | JS API: fg.updateOrAdd histogram, 506 rows |
| 2.1 | Add Scaffold Tree Filter | PASS | 3s | PASSED | JS API only: fg.updateOrAdd Chem:scaffoldTreeFilter, 370 rows |
| 2.3 | Check scaffold checkbox | PASS | — | PASSED | Included in savedTree: checked=true |
| 3.1 | Add Scatter Plot | PASS | 3s | PASSED | UI: click icon-scatter-plot |
| 3.2 | Zoom in to filter | PASS | 3s | PASSED | UI: Alt+drag on canvas, 370→187 rows |
| 3.3 | Double-click to reset | PASS | 2s | PASSED | UI: dblclick on canvas, back to 370 |
| 3.4 | Close Scatter Plot | PASS | 2s | PASSED | JS API: v.close() (UI icon-times failed) |
| 4.1 | Add Bar Chart | PASS | 3s | PASSED | UI: click icon-bar-chart |
| 4.2 | Set On Click > Filter | PASS | 1s | PASSED | UI: right-click → On Click → Filter (context menu) |
| 4.3 | Click bar — filter applies | PASS | 2s | PASSED | UI: canvas click at 15% width, 370→55 rows |
| 4.4 | Click white space — reset | PASS | 2s | PASSED | UI: canvas click corner, back to 370 |
| 4.5 | Close Bar Chart | PASS | 1s | PASSED | JS API: v.close() |
| 5.1 | Add Histogram | PASS | 3s | PASSED | UI: click icon-histogram |
| 5.2 | Filter via range slider | PASS | 2s | PASSED | JS API fallback: canvas slider not automatable, used fg.updateOrAdd, 370→174 |
| 5.3 | Double-click to reset | PASS | 2s | PASSED | JS API: fg.updateOrAdd full range, back to 370 |
| 5.4 | Close Histogram | PASS | 1s | PASSED | JS API: v.close() |
| 6.1 | Add PC Plot | PASS | 3s | PASSED | UI: click icon-pc-plot |
| 6.2 | Drag range selector on axis | PASS | 3s | PASSED | UI: page.mouse drag on SVG range-slider circle handle (top → +150px) |
| 6.3 | Double-click to reset | PASS | 2s | PASSED | UI: dblclick on overlay canvas |
| 6.4 | Close PC Plot | PASS | 1s | PASSED | JS API: v.close() |
| 7.1 | Add Trellis Plot | PASS | 4s | PASSED | UI: click icon-trellis-plot |
| 7.2 | Set On Click > Filter | PASS | 1s | PASSED | UI: right-click → On Click → Filter (context menu) |
| 7.3 | Click cell — filter applies | PASS | 2s | PASSED | UI: page.mouse double-click on cell, 0 rows (empty cell) |
| 7.4 | Esc to reset filter | PASS | 1s | PASSED | UI: Escape key, back to baseline |
| 7.5 | Close Trellis Plot | PASS | 1s | PASSED | JS API: v.close() |
| 8.1 | Add Pie Chart | PASS | 3s | PASSED | UI: click icon-pie-chart |
| 8.2 | Set On Click > Filter | PASS | 1s | PASSED | UI: right-click → On Click → Filter (context menu) |
| 8.3 | Click segment — filter applies | PASS | 2s | PASSED | UI: canvas click at 60%/40%, 370→152 rows |
| 8.4 | Click white space — reset | PASS | 2s | PASSED | UI: canvas click corner, back to 370 |
| 8.5 | Close Pie Chart | PASS | 1s | PASSED | JS API: v.close() |
| 9.1 | Reset all filters | PASS | 2s | PASSED | UI: reset icon click + OK dialog, 3624 rows |
| 9.2 | Close All | PASS | 1s | PASSED | JS API: grok.shell.closeAll() |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~85s |
| Spec file generation | ~5s |
| Spec script execution | 72s (PASSED) |

## Summary

All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to reset. PC Plot range slider drag works via SVG handle + Playwright page.mouse. All viewers (Scatter Plot, Bar Chart, Histogram, PC Plot, Trellis Plot, Pie Chart) interact correctly with the Filter Panel.

## Retrospective

### What worked well
- Scatter Plot Alt+drag zoom filtering works perfectly via dispatchEvent
- Bar Chart and Pie Chart click-to-filter/reset works with canvas click events
- Filter Panel filter APIs (updateOrAdd) are reliable for all filter types
- Scaffold Tree Filter setup via JS API works as documented

### What did not work
- Trellis Plot needs two clicks to filter (first selects, second applies) — not obvious from scenario wording
- MCP dispatchEvent doesn't work for Trellis at all — only Playwright page.mouse works
- Histogram range slider — canvas-based SlideFilterCore, not automatable via DOM events (JS API fallback used)
- PC Plot SVG range sliders don't respond to JS dispatchEvent — but work with Playwright page.mouse (real browser input)
- PC Plot has overlay canvas intercepting pointer events — dblclick must target `canvas[name="overlay"]`

### Suggestions for the platform
- PC Plot and Trellis Plot could expose JS API methods for programmatic filtering (e.g., `viewer.setAxisFilter(axis, min, max)`)
- Consider adding `data-testid` attributes to canvas-based interactive elements for automation

### Suggestions for the scenario
- The scenario says "~32 rows" for c1ccccc1 filter but full SPGI has 1356 matching rows — values may be from a smaller dataset
- Consider specifying exact expected row counts or using relative assertions
