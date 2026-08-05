# Tree Map (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Tree map from the Viewers toolbox | PASS | Split selector pre-filled; trailing placeholder empty |
| 2 | Split by RACE draws one leaf per race | PASS | Leaf tooltip name and count match the data |
| 3 | Adding SEX as a second level nests the leaves | PASS | New placeholder appears; leaves shrink below the largest race |
| 4 | Colour by AGE and switch the aggregation | PASS | avg → max, canvas recoloured both times |
| 5 | Size by WEIGHT rescales the rectangles | PASS | Column set, then `max` aggregation |
| 6 | Clicking a leaf selects exactly that group | PASS | Selection count equals the tooltip's row count |
| 7 | Filtering the table reshapes the map | PASS | Leaf counts follow the filter |
| 8 | Show Column Selection Panel hides the on-viewer selectors | PASS | Selectors hidden, then restored |
| 9 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 21s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The leaf tooltip (`"RA\n2550 rows"`) is the strongest read this viewer offers: it
  names the group and counts its rows, so the split, the filter and the click
  selection can all be asserted against real numbers instead of property
  round-trips.
* Split levels are plain `<select class="d4-column-selector-tree-map">` elements on
  the viewer, driven with `selectOption`. The trailing one is always the empty
  placeholder — picking a column in it appends a level and a new placeholder.
* Two things that cost a debugging cycle and are now encoded in the spec:
  `sum` is already the default size aggregation (selecting it changes nothing, so
  the step uses `max`), and hiding the column selection panel leaves the selectors
  in the DOM — only their visibility changes.

## Not automated

* Default Color, nested-border and caption-centring visuals, outer margins — kept
  as a manual list in `tree-map-viewer.md` and `tree-map-viewer-ui.md`.
