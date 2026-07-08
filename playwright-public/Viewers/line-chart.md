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

## Whiskers (error bars)

1. Go to the Context Panel > X, set X to RACE; go to the Context Panel > Y, set Y to AGE
2. Go to the Context Panel > Data, set Whiskers Type to Avg | Min, Max
3. Set Whiskers Type to Med | Q1, Q3
4. Set Whiskers Type to Avg | +/-StDev
5. Set Whiskers Type to Avg | +/-StError
6. Set Whiskers Type to None

## Markers

1. Right-click the viewer and go to Markers — set visibility to Always
2. From the Markers submenu, select the square marker shape
3. Select the triangle up marker shape
4. Drag the Size slider to 10
5. Drag the Opacity slider to 50
6. From the Markers submenu, go to Marker Column and set it to SEX
7. Go to the Context Panel > Marker, set Size Column to WEIGHT
8. Right-click the viewer, go to Markers > Marker Column and select the empty first row to clear it
9. Right-click the viewer, go to Markers and set visibility to Auto

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

## Aggregation types

1. Go to the Context Panel > X, set X to RACE; go to the Context Panel > Y, set Y to AGE
2. Right-click the viewer and go to Data > Aggregation > avg
3. Right-click the viewer and go to Data > Aggregation > min
4. Right-click the viewer and go to Data > Aggregation > max
5. Right-click the viewer and go to Data > Aggregation > med
6. Right-click the viewer and go to Data > Aggregation > sum
7. Right-click the viewer and go to Data > Aggregation > stdev
8. Right-click the viewer and go to Data > Aggregation > avg

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

## Split by column

1. Right-click the viewer and go to Data > Split Columns and set Split to SEX
2. Right-click the viewer and go to Data > Split Columns and set Split to RACE
3. Set Split to SEX and RACE
4. Remove all split columns

## Multi-axis mode

1. Go to the Context Panel > Y, set Y to AGE, HEIGHT, and WEIGHT
2. Right-click the viewer and go to Controls — check Multi Axis
3. Right-click the viewer and go to Controls — uncheck Multi Axis

## Title and description

1. Go to the Context Panel > General, set Show Title to true
2. Set Title to "My Line Chart"
3. Set Description to "Test description"
4. Set Description Position to Top
5. Set Description Position to Bottom
6. Set Description Visibility Mode to Never

## Custom axis min/max

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to HEIGHT
2. Go to the Context Panel > X, set Min to 30, Max to 60
3. Go to the Context Panel > Y, set Min to 150, Max to 190
4. Clear all custom min/max values

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
7. Set Legend Visibility to Never
8. Set Legend Visibility to Auto

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

## Table switching and row source (uses SPGI dataset)

Setup: Close all, open demog, then also open SPGI

1. Add a line chart on demog by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > Data, switch the table to SPGI — the line chart should update to show SPGI columns
3. Switch back to demog — the line chart should restore demog columns
4. Set Row Source to Selected
5. Select some rows in the grid — the line chart should show only selected rows
6. Set Row Source to Filtered
7. Open the filter panel and apply a filter — the line chart should reflect the filter

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

## Split and Y-columns sync with Context Panel

1. Close all and open demog
2. Add a line chart
3. On the viewer, use the in-plot Split selector to add SEX as a split column — verify the Context Panel > Data > Split shows 1 split column
4. Use the in-plot Split selector to add RACE — verify the Context Panel shows 2 split columns
5. Remove RACE via the in-plot selector — verify the Context Panel shows 1 split column
6. Go to the Context Panel > Y and set Y to AGE, HEIGHT, and WEIGHT
7. On the viewer, click the 'x' button on one of the Y lines to remove it — verify the Context Panel > Y immediately shows 2 Y columns
8. Remove another Y line via the 'x' button — verify Context Panel shows 1 Y column
9. On the viewer, right-click the canvas and go to Controls — check Multi Axis
10. Verify that legend categories' names start with the corresponding Y column name
11. Use the in-plot Y column selector to change one of the Y columns — verify the legend updates accordingly

## Selection checkboxes

1. Close all and open demog
1. On the viewer, set X to AGE, Y to HEIGHT
2. Hold Shift and drag a rectangle to select points — verify selection in the grid
3. Go to the Context Panel > Selection
4. Toggle the available checkboxes (e.g., Show Selected Only, Invert Selection) — verify the line chart interaction updates accordingly
5. Save the layout, close the line chart, apply the layout — verify selection settings are preserved

## Data panel checkboxes

1. Go to the Context Panel > Data
2. Toggle the available checkboxes (e.g., Show Null Values, Show Missing Values) one by one — verify the line chart updates after each toggle
3. Save the layout, close the line chart, apply the layout — verify the checkbox states are preserved
---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
