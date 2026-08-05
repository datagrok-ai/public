# Map viewer (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Map viewer from the Viewers toolbox | PASS | `[name="icon-Map"]`; Latitude/Longitude auto-detected and shown in the property grid; Markers GL visible, Heatmap hidden |
| 2 | Color by Magnitude | PASS | Column selector driven with real click + typing; property grid shows *Magnitude*, map repaints |
| 3 | Size by Depth | PASS | Same path; property grid shows *Depth*, map repaints |
| 4 | Marker Min Size redraws the markers | PASS | Property-grid slider editor: 2 → 12 → 2, map repaints |
| 5 | Layers panel toggles Heatmap and Markers GL | PASS | Canvas-drawn checkboxes hit by geometry; layer visibility verified through the OpenLayers stack; UI stays responsive after toggling Markers GL |
| 6 | Render Type cycles markers / heatmap / both | PASS | `<select>` in the property grid; each mode paints a different map |
| 7 | Zoom in and out with the map buttons | PASS | `.ol-zoom-in` / `.ol-zoom-out`; zoom 1.97 → 2.97 → 1.97 |
| 8 | Ctrl+drag selects points, Escape clears | PASS | ~166 of 2426 rows selected by a real drag; Escape clears to 0 |
| 9 | Filtering the table narrows the map | PASS | `MagType = Mw`; map repaints, filtered count below total |
| 10 | Show Tooltip reveals a tooltip over a point | PASS | Hover targets projected from feature geometry; tooltip lists the row's columns |
| 11 | Close the viewer after moving the pointer over it | PASS | Title-bar **Close**; viewer removed, no page error |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 25s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* Every action goes through the UI. The JS API is read-only here: OpenLayers layer visibility, zoom level, selected-row count, and the feature coordinates used to aim the tooltip hover.
* Verification is user-visible: the value the property grid displays, the layer stack, the zoom level, the selection count, and a pixel signature of the viewer region (screenshot hash) that proves the map actually repainted. The marker layer is WebGL, so `getImageData` is not available on the OL canvas.
* Layer checkboxes live on the layers mini-grid canvas — no DOM node — so they are hit by geometry (24px rows, checkbox column at the panel's left edge).
* Acting on the table (filtering) replaces the Context Panel content; steps that need the property grid afterwards reopen it through the viewer gear.

## Not automated

* Mouse-drag panning, **Selected Color**, and drag and drop of KML/KMZ/GeoJSON/TopoJSON layers — kept as a manual checklist at the end of `map.md`.
