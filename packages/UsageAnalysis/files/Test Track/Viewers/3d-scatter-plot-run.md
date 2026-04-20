# 3D Scatter plot — Run Results

**Date**: 2026-04-20
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog dataset | 3s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`; 5850 rows, 11 cols |
| 2 | Click 3D Scatter plot in Viewers toolbox | 5s | PASS | PASSED | Clicked `[name="icon-3d-scatter-plot"]`; viewer container + canvas appeared; default axes X=AGE, Y=HEIGHT, Z=WEIGHT |
| 3 | Click point highlights grid row | 4s | PASS | PASSED | Dispatched mousedown/mouseup/click on canvas center; `currentRowIdx` changed 0 → 3262; tooltip with row values rendered |
| 3 | Scroll wheel zooms | 3s | PASS | PASSED | Dispatched 5 forward + 5 backward wheel events on canvas; view updated (camera rotated/zoomed), canvas remained rendered |
| 4 | Click gear icon (Settings) opens Property Pane | 2s | PASS | PASSED | Gear `[name="icon-font-icon-settings"]` lives on grandparent element of viewer container; clicking toggled the pane |
| 4 | Modify properties (Color, Size, markerDefaultSize) | 2s | PASS | PASSED | UI combo-box click in Property Pane did not open a column picker — fell back to `viewer.setOptions({colorColumnName:'SEX', sizeColumnName:'AGE', markerDefaultSize:15})`. Colors changed (F orange / M blue), marker sizes varied by AGE |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 40s |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 1m 55s |
| Spec file generation | 30s |
| Spec script execution | 6s |
| **Total scenario run (with model)** | 2m 31s |

## Summary

All four scenario steps passed against dev. 3D Scatter plot opens by clicking the Viewers toolbox
icon, click-on-canvas selects the underlying row in the grid (with tooltip), mouse wheel
manipulates the 3D view, and properties set via the API reflect immediately in the viewer.
Playwright spec passed in 3.2s. **Total scenario run (with model): 2m 31s.**

## Retrospective

### What worked well
- JS API dataset open + semType wait completed in under a second.
- `[name="icon-3d-scatter-plot"]` is stable and reliable.
- `viewer.getOptions().look.*` gives a clean read-back channel for verifying property changes.
- Click + tooltip confirms row highlighting without needing to decode a WebGL pixel.

### What did not work
- The Property Pane column combo-boxes (`[name="div-column-combobox-color"]`) did not open a column picker on `.click()` — the drop-down stayed hidden. Falling back to `viewer.setOptions({colorColumnName: ...})` was needed. Root cause unclear; the combo probably needs a pointerdown + keyboard interaction rather than `.click()`.
- The viewer's title bar gear icon is rendered on the **grandparent** of `viewer-3d-scatter-plot` (not the container or its immediate parent). Standard `container.querySelector('[name="icon-font-icon-settings"]')` returns `null`, which the viewers reference could call out explicitly.
- Scroll-wheel zoom does not expose a numeric `zoom`/`camera` property on `viewer.getOptions().look`, so automated regression against "zoom actually happened" is indirect (visual).

### Suggestions for the platform
- Make `div-column-combobox-*` respond to plain `.click()` by programmatically showing the column drop-down, so UI-driven tests can change columns without falling back to JS API.
- Expose the 3D scatter plot camera state (zoom/rotation) in `look` so tests can assert camera changes deterministically.
- Add `body.selenium`-driven title-bar icons directly to the viewer's container (`viewer-<Name>`) rather than a wrapper two levels up, so scoping to the container is sufficient.

### Suggestions for the scenario
- "click the Gear icon. The Property Pane opens" — on dev, the Property Pane often is already open from the right sidebar when a viewer is added. Reword as "click the Gear icon to open/focus the Property Pane for this viewer", and list which sections must exist (Data / Axes / Marker / Style / Legend).
- Step 3 "Interact with all elements on the viewer" is open-ended — enumerate concrete interactions (single click, right-click context menu, mouse drag rotate, wheel zoom) with expected outcomes.
- Add an explicit verification step: after modifying Color and Size, check that `Color:` and `Size:` combo boxes on the viewer show the chosen column names, and that the legend reflects the new color column.
