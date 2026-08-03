# 3D Scatter Plot — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 0 | Setup: closeAll, open demog.csv, add 3D Scatter Plot via UI | 8s | PASS | PASSED | `[name="icon-3d-scatter-plot"]` click added the viewer; `grok.shell.tv.viewers` → `["Grid","3d scatter plot"]`; demog.csv loaded 5850 rows, 11 columns. ~1.8s of MCP latency. |
| 1 | Axis column assignment: AGE/HEIGHT/WEIGHT → WEIGHT/AGE/HEIGHT → back | 3s | PASS | PASSED | JS API `v.props.{x,y,z}ColumnName = ...`. All three permutations confirmed via `getProps`. |
| 2 | Axis types: X/Y/Z → logarithmic, then back to linear | 2s | PASS | PASSED | JS API `v.props.{x,y,z}AxisType`. Each axis flipped individually; final state all `linear`. |
| 3 | Color coding — categorical (SEX, RACE, clear) | 2s | PASS | PASSED | `colorColumnName` set then cleared (empty string). |
| 4 | Color coding — numerical (AGE, clear) | 1s | PASS | PASSED | Same pattern; final empty. |
| 5 | Size coding (WEIGHT, AGE, clear) | 2s | PASS | PASSED | `sizeColumnName` set then cleared. |
| 6 | Labels (SEX, clear) | 1s | PASS | PASSED | `labelColumnName` set then cleared. |
| 7 | Marker type (sphere/box/cylinder/tetrahedron/dodecahedron/octahedron) | 2s | PASS | PASSED | All six values round-tripped via `markerType`. |
| 8 | Marker opacity 20/100/69 + Random Rotation on/off | 2s | PASS | PASSED | `markerOpacity` and `markerRandomRotation` toggled cleanly. |
| 9 | Filtered out points: filter AGE 20–40, toggle, clear filter | 5s | PASS | PASSED | `fg.updateOrAdd({type:'histogram', column:'AGE', min:20, max:40})` → 2071 of 5850 rows passed filter; `showFilteredOutPoints` true→false; cleared by widening to `col.stats.min/max` → 5850 rows. |
| 10 | Axes visibility and grid lines | 2s | PASS | PASSED | `showAxes`, `showVerticalGridLines`, `showHorizontalGridLines` all toggled and restored. |
| 11 | Background and colors (JS API) | 1s | PASS | PASSED | `backColor=0xFF000000`, `axisLineColor=0xFFFFFFFF`, restored to white + original (`0xFF777878`). Scenario explicitly JS API only. |
| 12 | Dynamic camera movement (toggle on, toggle off) | 1s | PASS | PASSED | `dynamicCameraMovement` true → false. |
| 13 | Zoom and navigation (wheel up x5, wheel down x5, right-click Reset View) | 4s | PASS | PASSED | Dispatched `WheelEvent` + `MouseEvent('contextmenu')` on the canvas; context menu showed "Reset View" as the first item; clicked, menu closed. |
| 14 | Mouse-over row group highlight (add Bar Chart, toggle, close Bar Chart) | 5s | PASS | PASSED | `[name="icon-bar-chart"]` added Bar Chart; default `showMouseOverRowGroup=true`; toggled false then true; closed via `panel.closest('.panel-base').querySelector('.panel-titlebar-button-close')`. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 25s |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 40s |
| Spec file generation | 10s |
| Spec script execution | 24s |
| **Total scenario run (with model)** | 1m 30s |

All rows are full-phase wall-clock (incl. model thinking and retries), not just tool latency. The two `scenario steps` rows sum to `Execute via grok-browser (total)`. Reported Playwright runtime was 20.7s; the 24s figure includes `npx`/Node startup overhead. Spec generation is short because the existing spec already encoded the patterns the MCP run reproduced — no rewrite was needed beyond a sanity pass.

## Summary

All 14 steps (setup + 13 scenario sections) passed in both the MCP browser run against `https://dev.datagrok.ai` and the standalone Playwright replay. The MCP run executed in three batched JS API calls (sections 1–8, sections 9–12, section 13, section 14) for a total of ~15s of script latency; the spec replay completed in 20.7s. The spec was authored yesterday (2026-04-22) and required no edits — every selector, JS API call, and timing it encodes still matches the live behaviour on dev. Filter-AGE-20–40 still yields 2071 of 5850 rows; default `showMouseOverRowGroup` is still `true`; bar-chart close-button is still `.panel-titlebar-button-close`. **Total scenario run (with model)**: ~1m 30s.

## Retrospective

### What worked well
- JS API property patches (`v.props[k] = ...`) are the fastest and most reliable way to drive the property panel — sections 1–12 ran in ~10s combined.
- `fg.updateOrAdd({type: 'histogram', column, min, max})` is a clean, idempotent way to create and clear a numeric filter without having to interact with the filter panel UI.
- Canvas wheel-zoom and right-click-context-menu via dispatched `WheelEvent`/`MouseEvent('contextmenu')` worked reliably; Reset View was found by text on `.d4-menu-popup .d4-menu-item-label`.
- Bar Chart added via `[name="icon-bar-chart"]` Toolbox click and closed via `.panel-titlebar-button-close` (scoped to `.panel-base` ancestor).
- The previously-debugged spec passed on the first replay attempt — confirming the prior run's fixes (split setup phases, JS API property patches, `.first()` on generic locators) are still load-bearing.

### What did not work
- Nothing failed this run. Same caveats as the prior run still apply (documented below).

### Suggestions for the platform
- (Carried over from prior run) Make `filterType` string case consistent across `updateOrAdd`/`ApplyState` inputs and `filter.filterType` outputs (`'histogram'` vs `'Histogram'`).
- (Carried over) Add `[name="icon-times"]` to viewer title bars to match the documented convention (currently only `.panel-titlebar-button-close` exists on `.panel-base`).
- (Carried over) Provide `fg.clearFilter(columnName)` (or `filter.reset()`) so tests can clear a histogram filter without inspecting `col.stats`.
- (Carried over) Normalize viewer `.type` strings: `3d scatter plot` (lowercase) vs `Scatter plot` / `Bar chart` (title case) is inconsistent and forces test code to compare against the wrong form half the time.

### Suggestions for the scenario
- Step 1 "Set X to AGE, Y to HEIGHT, Z to WEIGHT" is the default for demog.csv — instruct testers to first set X/Y/Z to other columns so the assignment is a visible state change.
- "Background and colors" section is explicitly JS-only — note in the scenario that the colour-swatch UI is not automatable, so manual testers know what to click instead.
- "Move mouse over a bar in the Bar Chart" is hard to verify — the canvas-rendered bars have no DOM targets. Replace with an observable check like "confirm `df.mouseOverRowGroup.trueCount > 0`" or "tooltip appears".
- Add a final cleanup step (remove 3D Scatter Plot, close filter panel) so the next scenario starts from a clean view.
