# Scatter plot tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add scatter plot

## Changing axes

1. Set X to AGE
2. Set Y to WEIGHT
3. Set X to RACE — axes should switch to categorical
4. Set X to STARTED — axes should switch to date
5. Set X back to HEIGHT

## Axis types and inversion

1. Set X to AGE, Y to WEIGHT
2. Open properties, set X Axis Type to logarithmic
3. Set Y Axis Type to logarithmic
4. Check Invert X Axis — X axis should flip
5. Check Invert Y Axis — Y axis should flip
6. Uncheck both inversions
7. Set both axis types back to linear

## Axis min/max

1. Set X to AGE, Y to HEIGHT
2. Open properties, set X Min to 30 and X Max to 50
3. Plot should show only the 30–50 range on X
4. Set Y Min to 150 and Y Max to 180
5. Clear all min/max values — plot should reset to full range

## Zooming

1. Zoom in using Alt+drag over a group of points
2. Zoom in using mouse wheel
3. Zoom out using mouse wheel
4. Press H to reset view
5. Zoom in using Alt+drag again
6. Double-click empty area to reset view

## Panning

1. Zoom in using Alt+drag
2. Press Left arrow — plot should pan left
3. Press Right arrow — plot should pan right
4. Press Up arrow — plot should pan up
5. Press Down arrow — plot should pan down
6. Press H to reset view

## Range sliders

1. Zoom in using Alt+drag
2. Drag horizontal range slider all the way to the left
3. Drag horizontal range slider all the way to the right
4. Double-click horizontal range slider to reset
5. Drag vertical range slider to the top
6. Double-click vertical range slider to reset

## Rectangular selection

1. Shift+drag to select a group of points — selection should appear in the grid
2. Shift+drag another area — previous selection should be replaced
3. Ctrl+click a single point — it should toggle selection
4. Ctrl+Shift+drag over selected points — they should be deselected
5. Click on empty background — selection should clear

## Lasso selection

1. Press L to enable lasso tool
2. Shift+drag to draw a freeform selection — enclosed points should be selected
3. Check selection in the grid
4. Press L to disable lasso tool
5. Shift+drag — selection should be rectangular again

## Color coding

1. Set Color to SEX — points should be colored by category
2. Set Color to AGE — points should use a linear color scheme
3. Open properties, check Invert Color Scheme — colors should reverse
4. Uncheck Invert Color Scheme
5. Clear Color — all points should return to default color

## Size coding

1. Set Size to WEIGHT — point sizes should vary
2. Open properties, change Marker Min Size to 2
3. Change Marker Max Size to 40
4. Clear Size — points should return to uniform size

## Markers and jitter

1. Set Marker column to RACE — marker shapes should vary by category
2. Open properties, set Jitter Size to 20 — overlapping points should spread
3. Set Jitter Size Y to 15
4. Clear Marker column
5. Change Marker Type to square
6. Change Marker Type back to circle
7. Set Jitter Size and Jitter Size Y back to 0

## Labels

1. Open properties, set Label Columns to SEX
2. Labels should appear next to markers
3. Set Show Labels For to Selected
4. Select a few points — only selected points should show labels
5. Set Show Labels For to All
6. Check Use Label as Marker — labels should replace markers
7. Uncheck Use Label as Marker
8. Clear Label Columns

## Regression line

1. Press R — regression line should appear with equation
2. Set Color to RACE
3. Open properties, check Regression Per Category — separate lines per category
4. Check Show Spearman Correlation — ρ value should appear
5. Check Show Pearson Correlation — r value should appear
6. Uncheck Show Regression Line Equation — equation should hide
7. Press R to hide regression line

## Tooltips

1. Hover over a point — tooltip should show
2. Right-click, open Tooltip settings
3. Set Show Tooltip to "show custom tooltip"
4. Add AGE and RACE to Row Tooltip
5. Hover over a point — tooltip should show only AGE and RACE
6. Set Show Tooltip to "do not show"
7. Hover over a point — no tooltip should appear
8. Set Show Tooltip back to "inherit from table"

## Legend

1. Set Color to RACE — legend should appear
2. Open properties, set Legend Visibility to Never — legend should hide
3. Set Legend Visibility to Always — legend should reappear
4. Change Legend Position to top
5. Change Legend Position to left
6. Set Legend Position back to right

## Filter panel interaction

1. Click the Filters icon to open the filter panel
2. Filter RACE to show only one category
3. Scatter plot should update to show filtered points only
4. Open properties, set Zoom And Filter to "zoom by filter"
5. Change the filter — plot should auto-zoom to filtered points
6. Set Zoom And Filter to "no action"
7. Change the filter — plot should not zoom
8. Set Zoom And Filter back to "filter by zoom"

## Filtered out points

1. Click the Filters icon and filter some rows
2. Open properties, check Show Filtered Out Points — filtered out points should appear in gray
3. Uncheck Show Filtered Out Points

## Context menu

1. Right-click on the plot area — context menu should appear
2. Check that Reset View, Lasso Tool, and Tools submenu are present
3. Close the menu
4. Right-click on the X axis — axis-specific menu should appear
5. Right-click on the Y axis — axis-specific menu should appear
6. Close the menu

## Drop lines

1. Right-click the plot, go to Tools, check Show Drop Lines
2. Move the mouse over points — crosshair lines should follow the cursor showing coordinates
3. Right-click the plot, go to Tools, uncheck Show Drop Lines

## Formula lines

1. Right-click the plot, select Tools > Formula Lines — editor should open
2. Add a new horizontal line
3. Close the editor — line should appear on the plot
4. Right-click, select Tools > Formula Lines
5. Delete the line
6. Close the editor

## Annotation regions

1. Set X and Y to numerical columns (AGE, HEIGHT)
2. Right-click, select Tools > Draw Annotation Region
3. Draw a rectangular region on the plot
4. The region should appear with a label
5. Right-click on the region — region context menu should appear

## Axis histograms

1. Open properties, check Show X Histogram — histogram should appear above X axis
2. Check Show Y Histogram — histogram should appear next to Y axis
3. Change Histogram Bins to 20
4. Uncheck both histograms

## Grid lines and axes visibility

1. Open properties, uncheck Show Vertical Grid Lines — vertical lines should disappear
2. Uncheck Show Horizontal Grid Lines — horizontal lines should disappear
3. Uncheck Show X Axis — X axis should hide
4. Uncheck Show Y Axis — Y axis should hide
5. Recheck all four options

## Layout save and restore

1. Set X to AGE, Y to WEIGHT, Color to RACE, Size to HEIGHT
2. Enable regression line, set Jitter Size to 10
3. Save layout
4. Close all
5. Open demog
6. Apply the saved layout — all settings should restore

## Whiskers (error bars)

1. Open a dataset with error bar columns (or manually set whisker columns)
2. Open properties, set X Whisker Min and X Whisker Max to numerical columns
3. Horizontal error bars should appear on each point
4. Set Y Whisker Min and Y Whisker Max to numerical columns
5. Vertical error bars should appear

## Mouse drag mode

1. Open properties, set Mouse Drag to Select
2. Drag on the plot — it should select points instead of panning
3. Set Mouse Drag back to Pan
4. Drag on the plot — it should pan

## Keyboard shortcuts

1. Press Ctrl+A — all points should be selected
2. Press Ctrl+Shift+A — all points should be deselected
3. Press + (plus) — should zoom in
4. Press - (minus) — should zoom out
5. Press H — should reset view
6. Press R — regression line should toggle
7. Press L — lasso tool should toggle

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
