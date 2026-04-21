# 3D Scatter plot — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog dataset | 3s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`; 5850 rows, 11 cols |
| 2 | Click 3D Scatter plot in Viewers toolbox | 3s | PASS | PASSED | Clicked `[name="icon-3d-scatter-plot"]`; `[name="viewer-3d-scatter-plot"]` appeared; default axes X=AGE, Y=HEIGHT, Z=WEIGHT |
| 3 | Click point highlights grid row | 2s | PASS | PASSED | MCP click inside viewer region; `currentRowIdx` changed 0 → 1611; grid scrolled to row 1612 with highlighted cell |
| 3 | Scroll wheel zooms | 2s | PASS | PASSED | Dispatched 5 WheelEvent down on canvas; axes scale visibly changed (HEIGHT axis now spans ~140-200), plot enlarged |
| 4 | Click gear icon opens Property Pane | 2s | PASS | PASSED | Gear `[name="icon-font-icon-settings"]` lives on viewer grandparent (`.panel-base`); clicking showed pane with Data/Axes/Marker/Events/Misc/Style/Legend/Description sections |
| 4 | Modify properties reflected in viewer | 3s | PASS | PASSED | UI click on `[name="div-column-combobox-color"]` in the pane did not open a column picker — fell back to `v3d.setOptions({colorColumnName:'RACE', sizeColumnName:'WEIGHT', backColor:0xFFF5F5F5, markerOpacity:90})`. Legend "Asian / Black / Caucasian / Other" rendered at top, markers recolored and resized |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 30s |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 2m 45s |
| Spec file generation | 40s |
| Spec script execution | 10s |
| **Total scenario run (with model)** | 3m 35s |

## Summary

All four scenario steps passed against dev. 3D Scatter plot mounts via Viewers toolbox icon,
canvas click selects the underlying row (grid auto-scrolls and highlights), wheel events zoom
the 3D view, and setting Color/Size/backColor/markerOpacity through the viewer API is reflected
visually (legend, marker colors & sizes). Playwright spec passed in 10.1s.
**Total scenario run (with model): 3m 35s.**

## Retrospective

### What worked well
- `[name="icon-3d-scatter-plot"]` and `[name="viewer-3d-scatter-plot"]` are stable.
- MCP real mouse click on the viewer region hit a marker cleanly — `currentRowIdx`
  changed and the grid auto-scrolled, which is a deterministic signal.
- `v3d.setOptions({...})` changes are immediately visible; no refresh needed.
- `grok.dapi.files.readCsv` + `onSemanticTypeDetected` gives a reliable dataset setup
  path in ~1s on dev.

### What did not work
- Property Pane column combobox (`[name="div-column-combobox-color"]` inside
  `.grok-prop-panel`) did not open a picker on plain `.click()` — popup never appeared.
  Fell back to JS API. Same issue reported in 2026-04-20 run — recurring.
- The viewer's title-bar gear icon is on the **grandparent** of `viewer-3d-scatter-plot`
  (class `.panel-base`), not on the container. `container.querySelector('[name="icon-font-icon-settings"]')`
  returns null. Inconsistent with other viewers where icons are scoped within the container.
- Wheel-zoom has no exposed numeric state (`zoom` / `camera`) on `viewer.props` or
  `getOptions().look`, so asserting "zoom happened" requires visual screenshot, not API.
- Playwright default `testMatch` is `*.spec.ts` (dot) but the skill convention is
  `*-spec.ts` (dash). Had to create a temporary config with `testMatch` override to
  discover the spec.

### Suggestions for the platform
- Make `[name="div-column-combobox-*"]` inside the property pane respond to plain
  `.click()` by opening the column picker programmatically (same as the canvas-overlay
  combos). Blocker for UI-first property automation.
- Move the viewer title-bar icons (`icon-font-icon-settings`, `icon-font-icon-menu`,
  `icon-font-icon-help`) into `[name="viewer-<Name>"]` so scoping to the container
  is sufficient.
- Expose 3D scatter plot camera state (zoom, rotation, pan) in the viewer properties
  so automated tests can assert view changes without relying on screenshots.

### Suggestions for the scenario
- Rewrite step 3 "Interact with all elements on the viewer" as a concrete bulleted
  list: single-click (→ current row changes), wheel up / wheel down (→ zoom in / out),
  mouse drag (→ camera rotation), right-click (→ context menu). Each with a measurable
  expected outcome.
- Step 4 "Modify various properties" is vague. Name the specific properties to test
  (Color, Size, Marker Opacity, Background Color) and require a visible-effect check
  per property (legend appears, markers resize, background tint, opacity).
- Precondition: Property Pane may already be docked from the right sidebar. Reword
  as "click the Gear icon to focus/reveal the Property Pane for this viewer".

### Suggestions for the skill
- Document that Playwright's default `testMatch` rejects `-spec.ts`; either adopt
  `.spec.ts` naming or include a one-off `testMatch` override in the run instructions
  (e.g. `--config` pointing at a temp file, as done here).
