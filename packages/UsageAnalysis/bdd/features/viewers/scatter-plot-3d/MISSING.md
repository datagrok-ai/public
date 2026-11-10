# 3D scatter plot — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against
the viewer's "Automation surface" (`core/client/d4/lib/src/viewers/scatterplot3d/CLAUDE.md`) and the
library's rules (`public/libraries/bdd/CLAUDE.md`). `scatter-plot-3d-look.feature` and
`scatter-plot-3d-interaction.feature` hold what could be said; `scatter-plot-3d.feature` is the TestTrack
translation and is not touched here.

## Still open

- The current-row and mouse-over boxes (GROK-17780, core#129): needs `current row drawn` / `mouse-over row drawn`
  readings from `JsScatterPlot3dCore`'s two boxes.
- `current row` is 0-based where every other viewer counts from 1: `currentRow + 1` in `scatterplot_3d_core.dart`.
- The log-scale error is a balloon only: needs `'error': validateState() ?? ''` in the status values.
- Labels drawn: needs a `labels shown` reading (`labelsList`, and `getPendingLabels` for labels still loading).
- Axis categories (public#3594): needs `x categories` / `y categories` / `z categories` readings.
- core#929 — Show Vertical/Horizontal Grid Lines, Grid Line Color, Axis Line Color, Axis Text Color and Filtered Out Rows
  Color never reach the renderer (`fillLook` in `utils/js_interop/scatterplot_3d.dart`).
- Linear Color Scheme: needs a `color scale` hit area to open the scale's context menu.

## Won't do

- The first frame taller than the plot area, the hovered row kept after the pointer leaves, a tall narrow viewer not
  fitting the cube: the lead, 2026-10-19 (not a bug, not reproduced, fine for now).
- Dynamic Camera Movement: an endless rotation timer, nothing settles.
- Marker size as a number: the `point` area is a fixed 6×6 box; Size and Marker Size Scaling are claimed as frame changes.
- Removed, renamed and re-added columns, no rows, rows removed and added, a table swap, WebGL contexts released
  (core#1027): xamgle `d4 | Viewers | Column removal | 3D scatter plot…`, `d4 | Viewers | Rows change | 3D scatter
  plot…`, `d4 | Viewers | Table swap | 3D scatter plot…`, `d4 | Viewers | Data Manipulation | 3D scatter plot`; the
  rename is ApiTests `AI: Viewers: Shared rules`.
