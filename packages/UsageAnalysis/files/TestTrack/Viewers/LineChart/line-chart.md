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
related_bugs:
  - id: GROK-17835
    status: fixed
expected_results:
  - anchor: "Table switching and row source"
    expectation: >-
      After binding the line chart to SPGI its bound dataFrame.rowCount reads
      3624, and after binding back to demog it reads 5850 — the row source
      really moves between tables. Under Row Source = Selected, selecting rows
      changes the canvas (per-color diff deltaPx > 1000 against the pre-selection
      snapshot), and a Filtered row source with an AGE 20-40 filter drops
      df.filter.trueCount below 5850.
  - anchor: "Filter expression and collaborative filtering"
    expectation: >-
      Setting the in-viewer Filter expression drops the chart's own
      lc.filter.trueCount below the full SPGI row count (points are filtered on
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
  - anchor: "Data panel checkboxes"
    expectation: >-
      Toggling real viewer checkboxes (Pack Categories from Data, Multi Axis
      from Controls) off their defaults reads the new value back and raises no
      errors, and the toggled combo (packCategories=false, multiAxis=true)
      survives a layout save/close/restore round-trip unchanged.
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

1. Go to the Context Panel > X, set X to AGE
2. Go to the Context Panel > Data, verify Axes Follow Filter is enabled by default
3. Set Axes Follow Filter to false
4. Open the filter panel by clicking the Filters section in the Toolbox and narrow the AGE range
5. Go to the Context Panel > Data, set Axes Follow Filter to true

## Context menu

1. Right-click the center of the line chart viewer — a context menu should appear
2. Verify the menu contains items: Reset View, Tools, Data, Markers, Chart Type, Overview, Selection, Controls
3. Close the menu by pressing Escape

## Layout save and restore

1. Go to the Context Panel > X, set X to STARTED; go to the Context Panel > Y, set Y to AGE and HEIGHT
2. Right-click the viewer and go to Data > Split Columns and set Split to SEX
3. Right-click the viewer and go to Controls — check Multi Axis
4. Go to the Context Panel > Lines, set Line Width to 3, set Interpolation to Spline
5. Save the layout
6. Close the line chart viewer by clicking the X icon on the viewer title bar
7. Apply the saved layout
8. Verify all properties match the configuration from steps 1-4
9. Delete the saved layout

## Table switching and row source

Setup: Close all, open demog, then also open SPGI

1. Add a line chart on demog by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > Data, switch the table to SPGI — the line chart should update to show SPGI columns
3. Switch back to demog — the line chart should restore demog columns
4. Set Row Source to Selected
5. Select some rows in the grid — the line chart should show only selected rows
6. Set Row Source to Filtered
7. Open the filter panel and apply a filter — the line chart should reflect the filter

Actuation note (recon 2026-07-21): the table switch is driven by assigning
`lc.props.table` (the Context Panel > Data > Table dropdown backs the same
property); the observable is the viewer's bound `dataFrame.rowCount` moving
between the two tables (3624 ↔ 5850), which is the "shows the other table's
columns" outcome. "Shows only selected rows" is measured as a canvas change:
under Row Source = Selected an empty selection still paints the axes/lines, so a
blank-vs-drawn threshold is meaningless — instead the per-color canvas histogram
is snapshotted (after a settle-precheck diff) and the selection is applied, and
the resulting `deltaPx` proves the selected subset re-rendered.

## GROK-17835 regression (uses SPGI dataset)

Setup: Close all, open SPGI (not demog)

1. Add a line chart by clicking the Line Chart icon in the Toolbox > Viewers section
2. Right-click the viewer and go to Controls — check Multi Axis
3. Right-click the viewer and go to Data > Split Columns and set Split to Series and Scaffold Names
4. Go to the Context Panel > X, set X to Chemist 521
5. Hover over the viewer — no errors should occur

## Filter expression and collaborative filtering

1. Close all and open SPGI
2. Add a line chart
3. Go to the Context Panel > Data, set Filter to `${CAST Idea ID} <636500` — verify points are filtered on the chart
4. Open the filter panel
5. Apply a categorical filter — verify the chart reflects both the expression filter and the panel filter together
6. Apply a numeric filter — verify collaborative filtering between the Filter Panel and the in-viewer Filter value
7. Remove all filters

Actuation note (recon 2026-07-21): "points are filtered on the chart" is
measured as the viewer's own `lc.filter.trueCount` dropping below the full row
count (the expression filter is viewer-local — it narrows the plotted points
without touching the shared `df.filter`). The Filter Panel numeric filter is
applied through `getFiltersGroup().updateOrAdd` (the range strip is a
canvas-drawn embedded Histogram with no scriptable DOM handle), and
collaboration is confirmed by the shared `df.filter.trueCount` narrowing on top
of the expression filter.

## Selection checkboxes

1. Close all and open demog
2. On the viewer, set X to AGE, Y to HEIGHT
3. Hold Shift and drag a rectangle to select points — verify selection in the grid
4. Go to the Context Panel > Selection
5. Toggle the available checkboxes — verify the line chart interaction updates accordingly
6. Save the layout, close the line chart, apply the layout — verify selection settings are preserved

Actuation note (recon 2026-07-21): the Selection section on a line chart has no
"Show Selected Only" or "Invert Selection" checkbox — those were illustrative
names in the source scenario. The real Selection checkboxes are Show Selected
Rows, Show Current Row Line, Show Mouse-over Category, and Show Mouse-over Row
Line. Show Selected Rows has a measurable canvas effect: with rows selected and
no split column (so the only orange is the selection overlay), unchecking it
drops the SELECTION_HUE_RANGE pixel count and rechecking restores it. The three
mouse-over checkboxes need a live hovering pointer that headless automation
cannot drive, so they are exercised through their stored value plus a no-error
guard; persistence of a non-default combo is verified by the layout round-trip.

## Data panel checkboxes

1. Go to the Context Panel > Data
2. Toggle the available checkboxes one by one — verify the line chart updates after each toggle
3. Save the layout, close the line chart, apply the layout — verify the checkbox states are preserved

Actuation note (recon 2026-07-21): a line chart's Data section has no "Show Null
Values" or "Show Missing Values" checkbox — those were illustrative names in the
source scenario (they are grid properties). The real viewer booleans toggled
here are Pack Categories (Data) and Multi Axis (Controls); each is set off its
default with the new value read back and a no-error guard, and the toggled combo
is confirmed to survive a layout save/close/restore round-trip.
---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
