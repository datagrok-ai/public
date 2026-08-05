# 3D Scatter Plot (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add 3D scatter plot from the Viewers toolbox | PASS | On-viewer selectors read X: AGE, Y: HEIGHT, Z: WEIGHT |
| 2 | Reassign X and Z with the on-viewer selectors | PASS | Trusted click + typing; selector text and scene both change |
| 3 | Color by SEX shows a categorical legend | PASS | Legend rendered with exactly F and M |
| 4 | Color by AGE switches the legend to a gradient | PASS | Legend no longer categorical, scene repaints |
| 5 | Marker type redraws the markers | PASS | box / sphere / cylinder each paint a different scene |
| 6 | Marker opacity redraws the markers | PASS | 100 → 25 → 100 |
| 7 | Show Axes hides and restores the axes | PASS | Checkbox round-trip, scene repaints |
| 8 | X axis type switches to logarithmic | PASS | logarithmic → linear, scene repaints |
| 9 | Drag rotates the scene and Reset View restores it | PASS | Real drag; **Reset View** from the context menu |
| 10 | Mouse wheel zooms the scene | PASS | Real `mouse.wheel`, in and out |
| 11 | Click makes a row current, Shift+click selects it | PASS | currentRowIdx changes, selection count grows |
| 12 | Show Filtered Out Points repaints the filtered-away rows | PASS | `SEX = F`, ghost markers toggle |
| 13 | Hovering a bar chart bin highlights the matching 3D points | PASS | Real hover over the Bar Chart repaints the 3D scene |
| 14 | Legend position moves the legend | PASS | Auto → Left → Auto |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 1m 06s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The 3D scene is WebGL, so `getImageData` is unavailable — "did it repaint" is
  proved by a screenshot hash of the viewer region, and the meaning of each change
  by the value the UI shows (selector text, property-grid value, legend items).
* Properties are driven through the Context Panel property grid, which is the same
  `.property-grid` widget the Map viewer uses. The gear icon must be clicked first:
  the panel is not open by default in Tabs mode.
* Synthetic `dispatchEvent(new WheelEvent(...))` / `new MouseEvent('mousemove')`,
  used by the previous version of this spec, does not drive the viewer — d4 tracks
  pointer input through its own handlers and ignores untrusted events. Every
  gesture here uses real Playwright input.

## Deliberately not asserted

* "With **Show Mouse Over Row Group** off the plot does not change" is asserted on
  the setting, not on pixels. An unchanged screenshot would not prove the group
  highlight is gone — the single mouse-over row and the Bar Chart's own tooltip
  repaint that region too.

## Not automated

* Ctrl+click deselection, colour-scheme pickers, Style colours, Dynamic Camera
  Movement, Marker Random Rotation, grid-line toggles — kept as a manual list at
  the end of `3d-scatter-plot.md` and in `3d-scatter-plot-ui.md`.
