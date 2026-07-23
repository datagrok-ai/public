---
feature: linechart
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes:
  - viewers.line-chart
realized_as:
  - line-chart-spec.ts
related_bugs: []
expected_results:
  - anchor: "Axes follow filter"
    expectation: >-
      On a non-aggregated chart (unique-x scratch column rowIdx), with Axes
      Follow Filter enabled (default true, confirmed by read), a rowIdx
      2000-3500 filter synchronizes the X axis with the filter's range: the
      axis bounds read back through screenToWorld land in bands around the
      filter range (min in (1600, 2200), max in (3300, 3900)) and the axis
      span shrinks below 70% of the unfiltered span. With Axes Follow Filter
      disabled, the same filter drops rows (filter.trueCount < 5850) but the
      X-axis min/max stays within 5% of its pre-filter values, and the fixed
      axis span stays more than twice the synced span. Flipping the flag back
      to true with the filter active re-ranges the axis to the filtered span
      again, then the full range is restored. Canvas repaint deltas are
      logged for observability but are not load-bearing. Teardown removes the
      rowIdx scratch column (the filters group drops its histogram filter
      together with the column), the full row set is confirmed restored
      (filter.trueCount back to 5850), and X/Y return to AGE/HEIGHT.
  - anchor: "Table switching and row source"
    expectation: >-
      After binding the line chart to spgi-100 its bound dataFrame.rowCount reads
      100, and after binding back to demog it reads 5850 — the row source
      really moves between tables. Under Row Source = Selected, selecting rows
      changes the canvas (per-color diff deltaPx > 1000 against the pre-selection
      snapshot), and a Filtered row source with an AGE 20-40 filter drops
      df.filter.trueCount below 5850.
  - anchor: "Filter expression and collaborative filtering"
    expectation: >-
      Setting the in-viewer Filter expression drops the chart's own
      lc.filter.trueCount below the full spgi-100 row count (points are filtered on
      the chart) while leaving it above zero; a Filter Panel numeric filter then
      narrows the shared df.filter.trueCount below the full count on top of the
      expression filter.
  - anchor: "Selection checkboxes"
    expectation: >-
      A Shift-drag rubber-band selects rows (selection.trueCount > 0) and the
      chart paints the orange selection overlay (SELECTION_HUE_RANGE pixel count
      > 0). Unchecking Show Selected Rows removes the overlay (hue count drops
      below half); rechecking restores it. Toggling the remaining Selection
      checkboxes raises no errors, and a non-default Selection combo survives a
      layout save/close/restore round-trip unchanged.
  - anchor: "Lasso selection"
    expectation: >-
      Checking Tools > Lasso Tool in the context menu flips the viewer's
      lassoTool property to true (menu-to-prop). A freeform Shift+drag over the
      chart selects rows (df.selection.trueCount > 0), and pressing Escape
      clears the selection (trueCount back to 0). Unchecking Lasso Tool through
      the same menu item flips the property back to false, and the whole step
      raises no page or console errors.
  - anchor: "Data panel checkboxes"
    expectation: >-
      Setting the real viewer checkboxes (Pack Categories from Data, Multi Axis
      from Controls) to a non-default combo (packCategories=false,
      multiAxis=true) raises no errors, and the combo survives a layout
      save/close/restore round-trip unchanged.
  - anchor: "Chart types"
    expectation: >-
      Cycling Chart Type through Area Chart, Stacked Area Chart, Stacked Bar
      Chart, and back to Line Chart via the context menu lands each selected
      entry in the viewer's chart-type state after the menu click
      (menu-to-prop), with no errors.
  - anchor: "Axis configuration"
    expectation: >-
      Driving the full axis-configuration cycle (logarithmic/linear X and Y,
      Invert X Axis, grid-line toggles, label orientation) and reverting to
      defaults raises no console or page errors and leaves the chart alive
      (no-error floor).
  - anchor: "Interpolation"
    expectation: >-
      Switching Interpolation to Spline with tension 1.0 and back to None
      raises no console or page errors and leaves the chart alive (no-error
      floor).
  - anchor: "Left panel histogram"
    expectation: >-
      Enabling and disabling the left-bar Histogram raises no console or page
      errors and leaves the chart alive (no-error floor).
  - anchor: "Controls visibility"
    expectation: >-
      Unchecking all six Controls entries and rechecking them raises no console
      or page errors and leaves the chart alive (no-error floor).
  - anchor: "Y global scale"
    expectation: >-
      Toggling Y Global Scale on and off under Multi Axis with two Y columns
      raises no console or page errors and leaves the chart alive (no-error
      floor).
  - anchor: "Title and description"
    expectation: >-
      Setting the title and description, moving the description position, and
      hiding the description raises no console or page errors and leaves the
      chart alive (no-error floor).
  - anchor: "Date/time X axis"
    expectation: >-
      Setting X to STARTED and cycling X Map through Year, Month, Day of week,
      and None raises no console or page errors and leaves the chart alive
      (no-error floor).
  - anchor: "Line styling"
    expectation: >-
      Driving Line Width, Line Transparency, and Line Coloring Type through
      non-default values and back raises no console or page errors and leaves
      the chart alive (no-error floor).
  - anchor: "Axis tickmarks modes"
    expectation: >-
      Switching the X and Y Axis Tickmarks Mode to MinMax and back to Auto
      raises no console or page errors and leaves the chart alive (no-error
      floor).
  - anchor: "Overview chart"
    expectation: >-
      Selecting Line, Area Chart, Stacked Bar Chart, and None from the Overview
      menu updates the viewer's overview-type property after each selection
      (menu-to-prop round-trip).
  - anchor: "Legend"
    expectation: >-
      With Split = SEX the legend lists exactly the two categories F and M;
      Legend Visibility Never removes the legend and Auto restores it with the
      same two categories; position changes keep it rendered.
  - anchor: "Context menu — chart area"
    expectation: >-
      Right-clicking the chart area opens a context menu containing Reset View,
      Tools, Data, Markers, Chart Type, Overview, Selection, and Controls.
  - anchor: "Layout save and restore"
    expectation: >-
      After save → close → apply layout, the restored viewer reads back
      X = STARTED, Y = AGE + HEIGHT, Split = SEX, Multi Axis on, Line Width 3,
      and Interpolation Spline.
---

# Line chart tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a line chart by clicking the Line Chart icon in the Toolbox > Viewers section

## Chart types

1. Right-click the viewer and go to Chart Type > Area Chart
2. Right-click the viewer and go to Chart Type > Stacked Area Chart
3. Right-click the viewer and go to Chart Type > Stacked Bar Chart
4. Right-click the viewer and go to Chart Type > Line Chart

## Axis configuration

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to WEIGHT
2. On the viewer, right-click the X Axis and set X Axis Type to logarithmic
3. From the context menu, check Invert X Axis
4. On the viewer, right-click the Y Axis and set Y Axis Type to logarithmic
5. Right-click the X Axis and uncheck Show Vertical Grid Lines
6. Right-click the Y Axis and uncheck Show Horizontal Grid Lines
7. Go to the Context Panel > X, set X Axis Label Orientation to Vert
8. Right-click the X Axis and set X Axis Type to linear
9. Right-click the Y Axis and set Y Axis Type to linear
10. Right-click the X Axis and uncheck Invert X Axis
11. Right-click the X Axis and recheck Show Vertical Grid Lines; right-click the Y Axis and recheck Show Horizontal Grid Lines
12. Go to the Context Panel > X, set X Axis Label Orientation to Auto

## Interpolation

1. Go to the Context Panel > Lines, set Interpolation to Spline
2. Set Spline Tension to 1.0
3. Set Interpolation to None

## Left panel histogram

1. Right-click the viewer and go to Misc > Left Bar Chart > Histogram
2. Right-click the viewer and go to Misc > Left Bar Chart > None

## Controls visibility

1. Right-click the viewer and go to Controls — uncheck Show X Selector
2. From the Controls submenu, uncheck Show Y Selectors
3. Uncheck Show Aggr Selectors
4. Uncheck Show Split Selector
5. Uncheck Show X Axis
6. Uncheck Show Y Axis
7. Recheck all six controls

## Y global scale

1. Go to the Context Panel > Y, set Y to AGE and HEIGHT
2. Right-click the viewer and go to Controls — check Multi Axis
3. Go to the Context Panel > Y, set Y Global Scale to true
4. Set Y Global Scale to false
5. Right-click the viewer and go to Controls — uncheck Multi Axis

## Title and description

1. Go to the Context Panel > General, set Show Title to true
2. Set Title to "My Line Chart"
3. Set Description to "Test description"
4. Set Description Position to Top
5. Set Description Position to Bottom
6. Set Description Visibility Mode to Never

## Date/time X axis

1. Go to the Context Panel > X, set X to STARTED
2. Go to the Context Panel > X, set X Map to Year
3. Set X Map to Month
4. Set X Map to Day of week
5. Set X Map to None

## Line styling

1. Go to the Context Panel > Lines, set Line Width to 3
2. Set Line Transparency to 0.5
3. Set Line Coloring Type to Custom
4. Set Line Width to 1
5. Set Line Transparency to 0
6. Set Line Coloring Type to Auto

## Axis tickmarks modes

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to HEIGHT
2. Go to the Context Panel > X, set X Axis Tickmarks Mode to MinMax
3. Set X Axis Tickmarks Mode to Auto
4. Go to the Context Panel > Y, set Y Axis Tickmarks Mode to MinMax
5. Set Y Axis Tickmarks Mode to Auto

## Overview chart

1. Right-click the viewer and go to Overview — select Line
2. Right-click the viewer and go to Overview — select Area Chart
3. Right-click the viewer and go to Overview — select Stacked Bar Chart
4. Right-click the viewer and go to Overview — select None

## Legend

1. Right-click the viewer and go to Data > Split Columns and set Split to SEX
2. Go to the Context Panel > Legend, set Legend Visibility to Always
3. Set Legend Position to Left
4. Set Legend Position to Top
5. Set Legend Position to Bottom
6. Set Legend Position to Right
7. Set Legend Visibility to Never -- the legend disappears
8. Set Legend Visibility to Auto -- the legend reappears listing the same categories (F, M)
9. Clear the split column

## Axes follow filter

1. Add a unique-valued integer column rowIdx (0..5849) to demog; set X to
   rowIdx and Y to HEIGHT — a unique X keeps the chart non-aggregated, the
   regime where the flag governs the axis range
2. Go to the Context Panel > Data, verify Axes Follow Filter is enabled by default
3. With Axes Follow Filter enabled, open the filter panel by clicking the Filters
   section in the Toolbox and narrow the rowIdx range to 2000-3500 — the X axis
   synchronizes with the filter's range values, so its min/max pulls in to
   roughly [2000, 3500]
4. Restore the full range, then set Axes Follow Filter to false
5. Apply the same rowIdx 2000-3500 filter — points are filtered out on the
   chart but the X-axis min/max stays where it was
6. Set Axes Follow Filter back to true — the axis re-ranges to the still-active
   filtered span — then restore the full range
7. Remove the rowIdx scratch column — its histogram filter disappears from the
   filter panel together with the column — and put X back on AGE and Y back on
   HEIGHT

## Context menu

1. Right-click the center of the line chart viewer — a context menu should appear
2. Verify the menu contains items: Reset View, Tools, Data, Markers, Chart Type, Overview, Selection, Controls
3. Close the menu by pressing Escape

## Layout save and restore

1. Go to the Context Panel > X, set X to STARTED; go to the Context Panel > Y, set Y to AGE and HEIGHT
2. Right-click the viewer and go to Data > Split Columns and set Split to SEX
3. Right-click the viewer and go to Controls — check Multi Axis
4. Go to the Context Panel > Lines, set Line Width to 3, set Interpolation to Spline
5. Save the layout (top menu **View > Layout > Save to Gallery**)
6. Close the line chart viewer by clicking the X icon on the viewer title bar
7. Apply the saved layout (**View > Layout > Open Gallery**, click the saved layout)
8. Verify all properties match the configuration from steps 1-4
9. Delete the saved layout

## Table switching and row source

Setup: Close all, open demog, then also open spgi-100

1. Add a line chart on demog by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > Data, switch the table to spgi-100 — the line chart should update to show spgi-100 columns
3. Switch back to demog — the line chart should restore demog columns
4. Set Row Source to Selected
5. Select some rows in the grid (Ctrl+click a row to toggle its selection, Shift+drag to select a range) — the line chart should show only selected rows
6. Set Row Source to Filtered
7. Open the filter panel and apply a filter — the line chart should reflect the filter

## Filter expression and collaborative filtering

1. Close all and open spgi-100
2. Add a line chart
3. Go to the Context Panel > Data, set Filter to `${CAST Idea ID} <634834` — verify points are filtered on the chart
4. Open the filter panel
5. Apply a categorical filter — verify the chart reflects both the expression filter and the panel filter together
6. Apply a numeric filter — verify collaborative filtering between the Filter Panel and the in-viewer Filter value
7. Remove all filters

## Selection checkboxes

1. Close all and open demog
2. On the viewer, set X to AGE, Y to HEIGHT
3. Hold Shift and drag a rectangle to select points — verify selection in the grid
4. Go to the Context Panel > Selection
5. Toggle each Selection checkbox — **Show Selected Rows**, **Show Current Row Line**, **Show Mouse-over Category**, **Show Mouse-over Row Line** — verify the line chart interaction updates accordingly
6. Save the layout, close the line chart, apply the layout — verify selection settings are preserved

## Lasso selection

1. On the viewer, set X to AGE, Y to HEIGHT
2. Right-click the viewer and go to Tools — check Lasso Tool
3. Hold Shift and drag a freeform path over the chart
4. Verify some rows are selected
5. Press Escape — the selection is cleared
6. Right-click the viewer and go to Tools — uncheck Lasso Tool

## Data panel checkboxes

1. Go to the Context Panel > Data
2. Set the checkboxes to a non-default combination: uncheck **Pack Categories** (Context Panel > Data) and check **Multi Axis** (right-click > Controls)
3. Save the layout, close the line chart, apply the layout — verify the checkbox states are preserved

## Automation notes

Axes follow filter: the X-axis bounds are read as product state, not from the
canvas. `LineChartViewer.screenToWorld` (published JS API backed by the
viewer's X range slider) is sampled at two canvas points and the linear
screen-to-world map is extrapolated to the canvas edges, giving the axis
min/max (slightly widened by the axis margins, which the assert bands absorb).
With the flag enabled the bounds must land in bands around the filter range
and the span must shrink below 70% of the unfiltered span; with the flag
disabled the rows drop (filter.trueCount) while both bounds stay within 5% of
their pre-filter values, and the fixed span stays more than twice the synced
span. Canvas repaint deltas are still logged (console) for observability, but
they are not load-bearing — in a headless canvas the per-color histogram delta
of a re-range is indistinguishable from the point-drop repaint. The contrast
deliberately runs on a unique X: on an aggregated chart (duplicated x values,
e.g. x=AGE) the aggregate frame is rebuilt from the filtered rows and the
axis range is taken from it, so the axes follow the filter regardless of the
flag.

Table switching and row source: the table switch is driven by assigning `lc.props.table` (the
Context Panel > Data > Table dropdown backs the same property); the observable is
the viewer's bound `dataFrame.rowCount` moving between the two tables, which is
the "shows the other table's columns" outcome. "Shows only selected rows" is
measured as a canvas change: under Row Source = Selected an empty selection still
paints the axes/lines, so a blank-vs-drawn threshold is meaningless — instead the
per-color canvas histogram is snapshotted (after a settle-precheck diff) and the
selection is applied, and the resulting `deltaPx` shows the selected subset
re-rendered.

Filter expression and collaborative filtering: "points are filtered on the chart" is measured as the viewer's
own `lc.filter.trueCount` dropping below the full row
count (the expression filter is viewer-local — it narrows the plotted points
without touching the shared `df.filter`). The Filter Panel numeric filter is
applied through `getFiltersGroup().updateOrAdd` (the range strip is a
canvas-drawn embedded Histogram with no scriptable DOM handle), and
collaboration is confirmed by the shared `df.filter.trueCount` narrowing on top
of the expression filter.

Selection checkboxes: the Selection section on a line chart has no
"Show Selected Only" or "Invert Selection" checkbox — those were illustrative
names in the source scenario. The real Selection checkboxes are Show Selected
Rows, Show Current Row Line, Show Mouse-over Category, and Show Mouse-over Row
Line. Show Selected Rows has a measurable canvas effect: with rows selected and
no split column (so the only orange is the selection overlay), unchecking it
drops the SELECTION_HUE_RANGE pixel count and rechecking restores it. The three
mouse-over checkboxes need a live hovering pointer that headless automation
cannot drive, so they are exercised through their stored value plus a no-error
guard; persistence of a non-default combo is verified by the layout round-trip.

Lasso selection: on a line chart the Lasso Tool checkbox governs the shape of
annotation-region drawing (the Shift+drag row selection itself is X-range
based in the source), so the assertable signals are the menu-to-prop
round-trip of the checkbox, rows getting selected by the freeform Shift+drag
(selection.trueCount > 0), and Escape clearing the selection back to zero
(the app-level Reset Selection command). The freeform gesture's visual (lasso
outline while drawing annotation regions) remains a human-side check.

Data panel checkboxes: a line chart's Data section has no "Show Null
Values" or "Show Missing Values" checkbox — those were illustrative names in the
source scenario (they are grid properties). The real viewer booleans set
here are Pack Categories (Data) and Multi Axis (Controls). Multi Axis behavior
(enable/disable with canvas signals) is owned by multi-axis-and-split.md; this
step's own signal is the non-default combo surviving a layout
save/close/restore round-trip with a no-error guard.
---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv,System:AppData/Chem/tests/spgi-100.csv"]
}
