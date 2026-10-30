# 3D scatter plot — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against
the viewer's "Automation surface" (`core/client/d4/lib/src/viewers/scatterplot3d/CLAUDE.md`) and the
library's rules (`public/libraries/bdd/CLAUDE.md`). `scatter-plot-3d-look.feature` and
`scatter-plot-3d-interaction.feature` hold what could be said; `scatter-plot-3d.feature` is the TestTrack
translation and is not touched here.

## Signals the status does not carry

- **The current-row and mouse-over highlight boxes (GROK-17780, core#129).** The fix makes the box appear only
  round a point the plot draws (row source, filter, Show Filtered Out Points), but the status says nothing about
  the boxes: `currentRowBox.visible` / `mouseOverBox.visible` and the row each surrounds live in
  `JsScatterPlot3dCore` only. A feature can claim the current row of the table, not that the plot drew its box.
  Needs readings such as `current row drawn` (true/false) and `mouse-over row drawn`, read from the two boxes
  in `getStatus()`. Then: make a row current, filter it out → `current row drawn` false; Show Filtered Out
  Points on → true.
- **`current row` is 0-based.** `getWidgetStatus` reports `dataFrame.currentRow` as is, while every other
  viewer's `current row` reading counts from 1 with 0 for none (library CLAUDE.md, "Readings"). On a table
  whose row 0 is current (what `user opens … dataset` leaves) it reads 0, i.e. "no current row". The features
  claim the table's current row (`"USUBJID" of the current row should be …`) instead. Fix: `currentRow + 1`
  in `scatterplot_3d_core.dart` (`scatter-plot-3d.feature` only claims it "differs from before", so nothing
  breaks).
- **The log-scale error is not in the status.** `validateState()` returns "'turbo' contains negative or zero
  values and can't be visualized in the log scale", `updateAll` shows it as a balloon and keeps drawing the
  previous frame, and `getWidgetStatus` reports no `error` (the histogram, PC plot, pie chart and others put
  `validateState()` into `values['error']`). So `3d scatter plot viewer should report the error …` cannot be
  written, and nothing on the plot itself says why it did not change; the feature claims the balloon once.
  Needs `'error': validateState() ?? ''` in the values (and, for a person, the message on the canvas).
- **Labels drawn.** No reading counts the label meshes (`labelsList`, visible ones), so "a label column draws
  labels" is a signature change and "clearing takes them away" is the signature coming back to the
  remembered one. A `labels shown` reading would claim it directly, and would also separate "the labels are
  still loading" (`getPendingLabels`) from "there are none".
- **Axis categories (public#3594).** A string column on an axis is claimed only as a frame change with every
  row still shown; nothing reports how the axis lays the column out — how many categories it has, their
  order, or where each one's markers land — so "one column of markers per category" cannot be said. Needs an
  `x categories` / `y categories` / `z categories` reading (the category names in axis order, empty for a
  numeric column) in `getStatus()`. Then: X = RACE → `x categories` lists RACE's four categories; back to
  AGE → "".

## Look properties the renderer never receives

Measured on the stand (headless runner, the frame's `scene signature` before, after and back) and by a
screenshot diff of the viewer in a GPU browser: these properties are in the property grid but change nothing,
because `fillLook` (`utils/js_interop/scatterplot_3d.dart`) passes only backColor, markerType,
markerRandomRotation, dynamicCameraMovement, showAxes, the column names and axisFont, and the axis helper's
colors are hard-coded in `MyAxesHelper`:

- Show Vertical Grid Lines, Show Horizontal Grid Lines, Grid Line Color — the scene has no grid at all;
- Axis Line Color, Axis Text Color (the TestTrack "Background and colors" case sets Axis Line Color to white);
- Filtered Out Rows Color — the ghosts Show Filtered Out Points brings back are painted in Missing Value
  Color (`getColor`: `!filter[i] ? look.missingValueColor`), which does change the frame.

Filed as core#929 (2026-10-19). Not written as `@known-failure`: whether the fix is to draw them or to drop
the properties is open there, and a scenario would then claim the wrong thing.

## Steps the features did without

- **Linear Color Scheme.** The property is a color list; `user sets {string} property` cannot type one, and the
  scheme picker in the property grid has no step. The color-scale context menu (`addColorCodingMenu`) is
  reachable only by a right-click on the scale, which has no hit area (`_colorScaleBounds` is not reported).
  Needs a `color scale` hit area; then `user picks "…" from the context menu of the "color scale" area`.

## Answered

- **The first frame is taller than the plot area** (lead, 2026-10-19: leave as is, not a bug). Until the
  first layout change the WebGL canvas keeps the size it was created with, 36 px taller than the area
  between the selectors (996 against 960 on a 1000 px tall viewer); toggling any Show … Selector fixes it.
  No scenario claims either size.
- **The pointer leaving the plot keeps the hovered row** (lead, 2026-10-19: not reproduced by hand). In
  the harness `user moves the pointer away from …` left the table's mouse-over row on the hovered point;
  a person's mouse leaving the viewer clears it. The tooltip scenario claims the hover only.
- **A tall, narrow viewer does not fit the cube** (lead, 2026-10-19: fine for now). With X = RACE the
  Asian and Other columns and the start of the HEIGHT title fall outside the camera; not claimed.

## Won't do

- **Save as PNG (public#859).** `viewer-chrome.feature:13` skips it on purpose for every viewer; left to the lead.
- **Dynamic Camera Movement.** An endless rotation timer: nothing settles, a signature claim has no fixed frame.
- **Marker size as a number.** The `point` area is a fixed 6×6 box, so Size and Marker Size Scaling are
  claimed as frame changes, not as marker sizes.

## Edge cases (`scatter-plot-3d-edge-cases.feature`)

Written for `scatter-plot-3d-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

1. **Markers drawn.** `rows shown` is the filter count, not what the scene holds, so two core#1027
   defects have no claim: a one-row table draws nothing (the marker scale `1/(10*ln 1)` is
   Infinity) and rows removed while the plot is open stay drawn (200 markers after 50 of 200 rows
   are removed). Needs a `markers drawn` reading from the instanced mesh (`JsScatterPlot3dCore`);
   then a one-row table → 1, and remove 50 of 200 rows → 150, both tagged `@realizes:core#1027` (closed).

#### 2. Decided

1. **No column takes over** after the Size, Color or Label column is removed: the encoding is
   dropped (the removal rule of `viewer-rules-edge-cases.feature`). The two removal scenarios stay in
   `scatter-plot-3d-edge-cases.feature` (core#1027), since they claim the X, Y and Z kept and the plot
   still answering a property change, which the outline does not.

#### 3. Questions for owners

None. A null label is empty, as on the scatter plot (2026-10-30, not reproduced).

#### 3. API or Dart tests to write

1. Dart regTest: WebGL contexts released after 17 plots are opened and closed (core#1027; a loop no
   user claim measures).
2. Dart regTest: options or a layout naming a Label or Size column the table lacks (core#1027, same
   throw path as the removal).

### Edge cases: Won't do

Nothing.
