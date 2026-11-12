# Scatter plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `scatterplot/CLAUDE.md` "Automation
surface", the TestTrack `scatterplot-*.md` and `scatter-plot-ui.md`. The "(no value)" legend click is
in `scatter-plot-legend-no-value.feature`; the rejected line formula (GROK-17282) in
PowerPack's `formula-lines-dialog.feature`.

## Still open

- Structures on the X axis (`scatter-plot-ui.md`, items 7-8): needs an `x axis labels drawn as` reading (`text` /
  `structure`), claimed after the pick and after a saved layout is applied.
- public#3091 — a layout applied after its Y column was removed draws `y: Column "HEIGHT" does not exist.`, but the
  status carries no `error` (`scatterplot_status.dart` should pass `error: p.error`).
- Markers the renderer drew (public#2108, public#1411; Size on HEIGHT hides the 128 empty rows while `rows shown` reads
  1000): needs a `markers drawn` reading counted by `_ScatterPlotRenderer`.
- A wheel zoom's end: the zoom animation is not in `isRenderPending`, so a claim right after a wheel zoom reads a window
  still narrowing.
- core#1339 — a viewer's picture (General > Save as PNG, Control+C) logs three 404s for `/font/MaterialIcons-Regular.*`;
  once fixed, Control+C copying the picture (public#3052) can be claimed in `viewer-chrome.feature`.

## Won't do

- Moving Average Window typed past its bounds: core#1200's clamp tests (`../MISSING.md`).
- A plot that filters by its zoom, closed: `filter-panel/compose-with-viewers.feature` claims the rows come back.
- A layout naming a missing Color column, a removed Color or Size column, a renamed X or Color column, no rows, a one-row
  table: the rules every viewer shares (`../MISSING.md`).
- Scale (1M-row hover and drag): timing belongs to the Dart benchmarks, not a feature.
