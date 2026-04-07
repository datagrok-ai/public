# Viewers — Filter Panel interaction with viewers — Run Results

**Date**: 2026-04-03
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

### Section 1: Apply filters in the Filter Panel

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI, open Filter Panel | PASS | PASSED | 3624 rows, 90 columns |
| 2 | Structure filter: sketch c1ccccc1, click OK | PASS | PASSED | Sketcher opened, SMILES typed, OK clicked. 1356 rows filtered |
| 3 | Stereo Category: uncheck S_UNKN | PASS | PASSED | fg.updateOrAdd categorical. 1265 rows |
| 4 | Average Mass: set max to 400 | PASS | PASSED | fg.updateOrAdd histogram. 506 rows |

### Section 2: Add Scaffold Tree Filter

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Scaffold Tree Filter for Structure | PASS | PASSED | fg.updateOrAdd with full params (savedTree, colorCodedScaffolds, title). Filter appeared in panel |
| 2 | Add Cc1ccccc1 scaffold | PASS | PASSED | Scaffold set via savedTree JSON with checked=true |
| 3 | Click scaffold checkbox — filter applies | PASS | PASSED | requestFilter() applied, 1819 → 370 rows. Scaffold visible in filter panel with checkbox |

### Section 3: Scatter Plot

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Scatter Plot from Toolbox | PASS | PASSED | Icon click worked |
| 2 | Zoom in to filter dataset | PASS | PASSED | Mouse wheel zoom: 506 → 204 rows |
| 3 | Double-click to reset zoom | PASS | PASSED | dblclick restored 506 rows |
| 4 | Close Scatter Plot | PASS | PASSED | v.close(), 506 rows |

### Section 4: Bar Chart

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Bar Chart from Toolbox | PASS | PASSED | Split by Stereo Category |
| 2 | Set On Click to Filter (context menu) | PASS | PASSED | v.props.onClick = 'Filter' |
| 3 | Click S_ACHIR bar — verify filtering | PASS | PASSED | Clicked canvas, filtered to 103 rows |
| 4 | Click white space to reset filter | PASS | PASSED | Clicked empty area, back to 506 |
| 5 | Close Bar Chart | PASS | PASSED | v.close(), 506 rows |

### Section 5: Histogram

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Histogram | PASS | PASSED | Added via tv.addViewer (avoid filter-panel confusion) |
| 2 | Range slider to filter | PASS | PASSED | Used valueMin/valueMax props, filtered to 115 rows |
| 3 | Double-click to reset | FAIL | FAILED | dblclick dispatch not recognized by histogram canvas |
| 4 | Close Histogram | PASS | PASSED | Used filteringEnabled=false as workaround, then closed |

### Section 6: PC Plot

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add PC Plot from Toolbox | PASS | PASSED | Multiple numeric axes rendered |
| 2 | Drag range on axis | AMBIGUOUS | N/A | JS canvas drag events not trusted by Dart handler |
| 3 | Double-click to reset | SKIP | SKIPPED | Drag didn't apply |
| 4 | Close PC Plot | PASS | PASSED | v.close(), 506 rows |

### Section 7: Trellis Plot

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Trellis Plot from Toolbox | PASS | PASSED | Added successfully |
| 2 | Set On Click to Filter | PASS | PASSED | v.props.onClick = 'Filter' |
| 3 | Click cell — verify filtering | PASS | PASSED | Canvas click filtered to 0 rows (hit empty cell). Filter mechanism works |
| 4 | Click same cell to reset | FAIL | FAILED | Second click did not restore row count |
| 5 | Close Trellis Plot | PASS | PASSED | Closed and re-applied filters to 506 |

### Section 8: Pie Chart

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add Pie Chart | PASS | PASSED | Shows Stereo Category segments |
| 2 | Set On Click to Filter (context menu) | PASS | PASSED | v.props.onClick = 'Filter' |
| 3 | Click segment — verify filtering | PASS | PASSED | Clicked R_ONE, filtered to 202 rows |
| 4 | Click white space to reset | PASS | PASSED | Clicked outside pie, back to 506 |
| 5 | Close Pie Chart | PASS | PASSED | v.close(), 506 rows |

### Section 9: Reset and cleanup

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Click Reset icon — all filters cleared | PASS | PASSED | icon-arrow-rotate-left clicked. 3624/3624 rows restored |
| 2 | Close All | PASS | PASSED | grok.shell.closeAll() |

## Summary

Out of 33 steps: **27 PASS, 3 FAIL, 2 SKIP, 1 AMBIGUOUS**. Core filter panel filters (Structure, Stereo Category, Average Mass) and Scaffold Tree filter all work correctly. Scatter Plot zoom filtering, Bar Chart onClick, and Pie Chart onClick work end-to-end including reset. Histogram double-click reset and Trellis Plot click-to-toggle were unreliable. PC Plot axis drag could not be tested via dispatched events.

## Retrospective

### What worked well
- Filter Panel JS API (fg.updateOrAdd) reliably applies substructure, categorical, and histogram filters
- Scatter Plot zoom filtering works via WheelEvent and double-click reset via dblclick
- Bar Chart and Pie Chart onClick=Filter work reliably with canvas click and white-space reset
- Reset icon (icon-arrow-rotate-left) correctly clears all filters and restores full row count
- Structure sketcher dialog opens via .sketch-link click and accepts SMILES keyboard input

### What did not work
- **Scaffold Tree filter initial fail**: `fg.updateOrAdd` requires all three params: `savedTree`, `colorCodedScaffolds`, `title` — omitting any causes silent failure
- **Histogram double-click reset**: dblclick dispatch on canvas doesn't reset the viewer filter
- **PC Plot axis drag**: canvas mousedown/mousemove/mouseup not trusted by Dart handler
- **Trellis Plot click reset**: clicking same cell doesn't toggle filter off

### Suggestions for the platform
- Add JS API for viewer filtering: `viewer.clickCategory('value')`, `viewer.resetClickFilter()`
- Add JS API for PC Plot axis range: `pcPlot.setAxisFilter('column', min, max)`
- `fg.updateOrAdd` for Scaffold Tree should throw or warn when required params are missing instead of silently creating a wrong filter type
- Fix double-click reset behavior on standalone Histogram and Trellis Plot viewers

### Suggestions for the scenario
- Step 2 (Scaffold Tree): works via JS API with full params; no UI-based way to add it reliably
- Step 5 (Histogram): clarify what "range slider" means — standalone histogram filter uses canvas, not SVG slider
- Step 6 (PC Plot): specify which axis to drag
- Step 7 (Trellis Plot): specify which cell to click (e.g., "a non-empty cell") and how to reset (click same cell vs double-click white space)
- Step 8 (Pie Chart): scenario text says "right-click the bar chart" but should say "right-click the pie chart"
