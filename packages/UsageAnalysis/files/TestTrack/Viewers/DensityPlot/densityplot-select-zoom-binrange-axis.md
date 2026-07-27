---
feature: densityplot
realizes_atlas:
  - densityplot.cp.select-zoom-binrange-axis
  - densityplot.bin-to-range-viewport
realizes:
  - viewers.density-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - densityplot-select-zoom-binrange-axis-spec.ts
related_bugs: []
expected_results:
  - anchor: "Bin click selection, tooltip, Esc round-trip"
    expectation: "df.selection.trueCount rises above zero after a click at the
      center of the plot canvas (recon reference: deterministically 18 rows at
      bins=50 hexagon on a 430x748 viewer — asserted as > 0, the exact count is
      a reference only); hovering a populated bin shows the DOM tooltip with
      both axis bin ranges and the bin's row count; Esc with focus inside the
      viewer returns df.selection.trueCount to zero (round-trip)"
  - anchor: "Wheel zoom and context-menu Reset View"
    expectation: "A mouse-wheel zoom produces a settle-gated canvas repaint
      (reference non-white fraction 0.27 -> 0.48); Reset View from the canvas
      context menu restores the render DIRECTIONALLY toward the pre-zoom
      baseline — not pixel-exact, since the bins are recomputed on restore"
  - anchor: "Bin To Range toggled while zoomed"
    expectation: "With the plot zoomed in, toggling Bin To Range produces a
      settle-gated canvas repaint in each direction — the bins recompute against
      the full data range instead of the zoomed viewport and back"
  - anchor: "Explicit axis bounds with Bin To Range"
    expectation: "With explicit X/Y Min/Max bounds set, toggling Bin To Range
      repaints or holds an error-free floor; clearing all four bounds restores
      the render directionally toward the full-range baseline"
  - anchor: "Axis invert and logarithmic round-trip"
    expectation: "Inverting the Y axis and switching it to logarithmic on HEIGHT (an
      all-positive column) and reverting both raises no page or console errors;
      setting a logarithmic Y axis on a column whose minimum is at or below zero
      also raises no page or console errors — the viewer renders an on-canvas
      warning instead of failing (the warning is canvas-only, its text is not
      asserted)"
---

# Density Plot — Bin Selection, Zoom, Bin To Range, Axis Configuration

## Purpose

Verifies the Density Plot's interactive surface on the demog dataset:
selecting rows by clicking a bin (with the tooltip identifying the bin and the
Esc round-trip clearing the selection), zooming with the mouse wheel and
restoring the view from the context menu, the coupling between the Bin To
Range setting and the zoomed viewport, explicit axis bounds, and axis
inversion and logarithmic scale including the non-positive-column edge.
Selection is read back from the table state, so the click gesture is verified
by its product effect; rendered outcomes are checked as settle-gated pixel
deltas or an honest error-free floor.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Density Plot viewer via the Toolbox icon.
4. Set the X column to AGE and the Y column to HEIGHT explicitly (the auto-pick
   already matches on demog, but the auto-selection order is not contractual).
5. Confirm the defaults: Bins 50, hexagon bin shape.

## Scenarios

### Scenario 1: Bin click selection, tooltip, and the Esc round-trip

Steps:
1. Read the selected-row count — record the baseline value (expected: 0).
2. Click the center of the plot canvas (fraction 0.5, 0.5 of the canvas) —
   with X=AGE / Y=HEIGHT, bins 50, hexagon, the center lands in a populated
   bin.
3. Verify the selected-row count rises above zero.
4. Hover the same populated bin and wait for the tooltip.
5. Verify the tooltip lists the AGE bin range, the HEIGHT bin range, and the
   bin's row count.
6. Press Esc while focus is still inside the viewer.
7. Verify the selected-row count returns to zero (round-trip).

Expected:
- The center-of-canvas click selects the clicked bin's rows: the selected-row count is above zero
- The hover tooltip shows both bin ranges and the row count
- Esc clears the selection back to zero

### Scenario 2: Mouse-wheel zoom and Reset View from the context menu

Steps:
1. Measure the canvas ink at the full data range (baseline).
2. Zoom in one step with the mouse wheel over the canvas center.
3. Verify a render change against the baseline (the zoomed viewport repaints
   the bins).
4. Right-click the plot canvas and click **Reset View** in the context menu
   (a text menu item).
5. Verify the render restores DIRECTIONALLY toward the baseline — the restore
   is not pixel-exact because the bins are recomputed; a double-click does
   NOT reset the view on this viewer and must not be used as the restore
   path.

Expected:
- Wheel zoom produces a settle-gated canvas repaint
- Reset View restores the render directionally toward the pre-zoom baseline (not pixel-exact)

### Scenario 3: Bin To Range toggled while zoomed

The Bin To Range setting couples binning to the viewport: with it off, bins
are recomputed against the current zoom viewport; with it on, bins span the
full data range and ignore the zoom. The toggle is only distinguishable while
zoomed.

Steps:
1. Zoom in one step with the mouse wheel and measure the canvas ink in the
   zoomed state.
2. Open the viewer's settings and, in the **Misc** section, turn
   **Bin To Range** on.
3. Verify a render change — the bins now span the full data range instead of
   the zoomed viewport.
4. Turn **Bin To Range** back off.
5. Verify a render change back (round-trip).

Expected:
- Toggling Bin To Range while zoomed repaints the canvas in each direction

### Scenario 4: Explicit axis bounds with Bin To Range

Steps:
1. Right-click the plot canvas and click **Reset View** to return to the full
   data range; measure the canvas ink (full-range baseline).
2. In the property panel's **X** section set **X Min** to 30 and **X Max** to
   60; in the **Y** section set **Y Min** to 150 and **Y Max** to 190 — the
   plot re-renders inside the explicit bounds.
3. Turn **Bin To Range** on, then off again.
4. Verify a render change from the toggle, or that the sequence completes
   with no page or console errors if the live-calibrated delta is too modest.
5. Clear all four bounds (X Min, X Max, Y Min, Y Max).
6. Verify the render restores directionally toward the full-range baseline.

Expected:
- With explicit bounds set, the Bin To Range toggle repaints or holds an error-free floor
- Clearing the bounds restores the render toward the full-range baseline

### Scenario 5: Axis invert and logarithmic scale, including the non-positive edge

Steps:
1. Note the current page-error and console-error counts.
2. In the property panel's **Y** section, turn **Invert Y Axis** on — the Y
   direction flips.
3. Set **Y Axis Type** to logarithmic — HEIGHT is an all-positive column, so
   the log scale applies cleanly.
4. Revert both: **Y Axis Type** back to linear, **Invert Y Axis** off
   (round-trip).
5. Verify no new page or console errors appeared across the round-trip.
6. Add a computed column whose minimum is at or below zero (for example
   HEIGHT − 200) and set it as the Y column.
7. Set **Y Axis Type** to logarithmic.
8. Verify the viewer stays alive with no new page or console errors — the
   viewer renders an on-canvas warning for the non-positive values instead of
   failing; the warning is drawn on the canvas only, so its text is not
   asserted.
9. Restore the Y column to HEIGHT with a linear axis type and remove the
   computed column.

Expected:
- The invert + logarithmic round-trip on HEIGHT raises no page or console errors
- A logarithmic Y axis on a column with minimum at or below zero raises no page or console errors (the on-canvas warning has no DOM signal and is not asserted)

## Automation notes

Setup: the viewer handle is
`const dp = grok.shell.tv.viewers.find(v => v.type === 'Density plot');` the
viewer is added via `[name="icon-density-plot"]`; X/Y are set explicitly via
`dp.props.xColumnName` / `yColumnName` (the UI selector path is owned by the
p0 scenario). The selection count is read as
`grok.shell.tv.dataFrame.selection.trueCount`, never hard-coded.

Scenario 1: both a real CDP click and the synthetic
`mousedown`+`mouseup`+`click` chain on `canvas[name="canvas"]` trigger the
bin selection (a rare canvas gesture where synthetic events work). The recon
reference is 18 rows selected, twice identically, at the canvas center on a
430x748 viewer — the assert is `trueCount > 0` because the exact count
depends on viewer size and bin geometry. Esc must be pressed while focus is
still inside the viewer (right after the click); Esc with focus outside can
trigger the platform's beforeunload dialog. The tooltip is the document-level
`.d4-tooltip`; its text carries both bin ranges and a `N rows` line
(regex `(\d+) rows`).

Scenarios 2-4 — pixels: ink is measured on the single
`canvas[name="canvas"]` as a settle-gated non-white pixel fraction (settle on
`d4-viewer-rendered`). Zoom is a synthetic `WheelEvent` (`deltaY: -300`) at
the canvas center (recon reference 0.27 -> 0.48; re-calibrate live). Reset
View: dispatch a synthetic `contextmenu` on the canvas, then find the menu
item by exact trimmed text `Reset View` among `[role="menuitem"]` elements
(they carry no `name=` attribute). The restore is NOT pixel-exact (recon:
0.27 pre-zoom vs 0.39 after reset with identical props) — assert a
directional delta toward the baseline, not equality. `props.viewport` stays
`null` through zoom, so it is not a readable zoom signal. Bin To Range is the
`[name="prop-bin-to-range"]` checkbox row (caption "Bin To Range", category
Misc).

Spec floor calibration (re-calibrate if the viewer size changes): wheel zoom —
settledPx basePx ~189310 -> zoomPx ~318536 (delta ~129226), spec floor 50000
(~2.5x margin); Bin To Range toggled while zoomed — summed per-color delta
~527302 in each direction, spec floor 50000 (wide margin).

Scenario 4 — bounds: the X Min and Y Min rows BOTH carry `name="prop-min"`
(and the Max rows `name="prop-max"`), so the rows are scoped via the unique
view cells `[name="prop-view-x-min"]`, `[name="prop-view-x-max"]`,
`[name="prop-view-y-min"]`, `[name="prop-view-y-max"]` and
`.closest('.property-grid-item')`; the editor is
`input.property-grid-item-editor-textbox` (click the view cell, focus the
editor, type the value, Enter).

Scenario 5: demog's numeric columns are all positive, so the non-positive
edge uses a helper column added via the JS API (e.g.
`df.columns.addNewCalculated('LOG_EDGE', '${HEIGHT} - 200')`), assigned via
`dp.props.yColumnName`, and removed in cleanup. The axis rows are
`[name="prop-invert-y-axis"]` (checkbox) and `[name="prop-y-axis-type"]`
(select editor, linear/logarithmic). The on-canvas log warning has no DOM
signal, so the edge check is an honest no-error floor. All console-error
guards filter infra noise (the shared dev server logs unrelated
`WebSocket ... 503` reconnect errors).

---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv"]
}
