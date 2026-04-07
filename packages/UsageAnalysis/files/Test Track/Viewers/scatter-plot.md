# Scatter plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a scatter plot by clicking the Scatter Plot icon in the Toolbox > Viewers section

## Changing axes

1. On the viewer, click the X column selector and set X to AGE
2. On the viewer, click the Y column selector and set Y to WEIGHT
3. On the viewer, click the X column selector and set X to RACE
4. On the viewer, click the X column selector and set X to STARTED
5. On the viewer, click the X column selector and set X back to HEIGHT

## Axis types and inversion

1. On the viewer, click the X column selector and set X to AGE
2. On the viewer, click the Y column selector and set Y to WEIGHT
3. On the viewer, right-click the X Axis and select from the context menu X Axis Type > logarithmic
4. From the context menu check the Invert X Axis checkbox
5. On the viewer, right-click the Y Axis and select from the context menu Y Axis Type > logarithmic
6. From the context menu check the Invert Y Axis checkbox
7. Go to the Context Panel > X, uncheck Invert X Axis checkbox and set X Axis Type to linear
8. Go to the Context Panel > Y, uncheck Invert Y Axis checkbox and set Y Axis Type to linear

## Axis min/max

1. On the viewer, click the X column selector and set X to AGE, click the Y column selector and set Y to HEIGHT
2. Go to the Context Panel > X, set Min to 30, Max to 50
3. Go to the Context Panel > Y, set Min to 150, Max to 180
4. Clear all min/max values

## Color coding

1. Hover over the viewer — the Color selector appears — click it and set Color to SEX
2. Click the Color selector and set Color to AGE
3. Go to the Context Panel > Color, check Invert Color Scheme
4. Uncheck Invert Color Scheme
5. Click the Color selector on the viewer and select the empty first row to clear Color

## Size coding

1. Hover over the viewer — the Size selector appears — click it and set Size to WEIGHT
2. Go to the Context Panel > Marker, set Marker Min Size to 2, set Marker Max Size to 40
3. Click the Size selector on the viewer and select the empty first row to clear Size

## Markers and jitter

1. Go to the Context Panel > Marker and set Markers column to RACE
2. Set Jitter Size to 20
3. Set Jitter Size Y to 15
4. Select the empty first row in the Markers column picker to clear it
5. Set Marker Type to square
6. Set Marker Type back to circle
7. Set Jitter Size and Jitter Size Y back to 0

## Labels

1. On the viewer, right-click the canvas and go to Labels — check SEX in the column list
2. From the Labels submenu, set Show Labels For to Selected
3. Set Show Labels For to All
4. Check Use Label as Marker
5. Uncheck Use Label as Marker
6. Uncheck SEX in the Labels column list to clear labels

## Regression line

1. On the viewer, right-click the canvas and go to Tools > Show Regression Line
2. Hover over the viewer — click the Color selector and set Color to RACE
3. Go to the Context Panel > Lines, check Regression Per Category
4. Check Show Spearman Correlation
5. Check Show Pearson Correlation
6. Uncheck Show Regression Line Equation
7. Press the R key to toggle regression line off

## Legend

1. Hover over the viewer — click the Color selector and set Color to RACE
2. On the viewer, right-click the legend area and set Legend Visibility to Never
3. Set Legend Visibility to Always
4. Set Legend Position to Top
5. Set Legend Position to Left
6. Set Legend Position back to Right

## Filter panel interaction

1. Open the filter panel by clicking the Filters section in the Toolbox
2. On the Filter Panel, in the RACE filter click 'Caucasian'
3. Go to the Context Panel > Data, set Zoom and Filter to 'zoom by filter'
4. On the Filter Panel, in the RACE filter click 'Asian'
5. Go to the Context Panel > Data, set Zoom and Filter to 'no action'
6. On the Filter Panel, in the RACE filter click 'Black'
7. Go to the Context Panel > Data, set Zoom and Filter to 'filter by zoom'

## Filtered out points

1. Open the filter panel and on the Filter Panel, in the SEX filter click 'M' to filter some rows
2. On the viewer, right-click the canvas and go to Filter > Show Filtered Out Points
3. Right-click again and go to Filter > uncheck Show Filtered Out Points

## Axis histograms

1. Go to the Context Panel > Axes, check Show X Histogram
2. Check Show Y Histogram
3. Set Histogram Bins to 20
4. Uncheck both Show X Histogram and Show Y Histogram

## Grid lines and axes visibility

1. On the viewer, right-click the X Axis and uncheck Show Vertical Grid Lines
2. From the context menu uncheck Show X Axis
3. On the viewer, right-click the Y Axis and uncheck Show Horizontal Grid Lines
4. From the context menu uncheck Show Y Axis
5. Right-click the X Axis and recheck Show Vertical Grid Lines and Show X Axis
6. Right-click the Y Axis and recheck Show Horizontal Grid Lines and Show Y Axis

## Mouse drag mode

1. Go to the Context Panel > Misc, set Mouse Drag to Select
2. Set Mouse Drag to Pan

## Whiskers (error bars)

1. Go to the Context Panel > X, set X Whisker Min to AGE, X Whisker Max to WEIGHT
2. Go to the Context Panel > Y, set Y Whisker Min to HEIGHT, Y Whisker Max to WEIGHT
3. Clear all whisker columns by selecting the empty first row in each picker

## Rectangular selection

1. On the viewer, click the X column selector and set X to AGE, click the Y column selector and set Y to HEIGHT
2. Hold Shift and drag a rectangle over the center area of the plot
3. Verify some points are selected
4. Press Escape to deselect
5. Hold Shift and drag a rectangle over a different area
6. Verify some points are selected
7. Press Escape to deselect

## Lasso selection

1. On the viewer, click the X column selector and set X to AGE, click the Y column selector and set Y to HEIGHT
2. On the viewer, right-click the canvas and check Lasso Tool
3. Hold Shift and drag a freeform lasso area over the plot
4. Verify some points are selected
5. Press Escape to deselect
6. Right-click the canvas and uncheck Lasso Tool

## Layout save and restore

1. On the viewer, set X to AGE, Y to WEIGHT via column selectors
2. Set Color to RACE and Size to HEIGHT via the viewer selectors
3. Right-click the canvas > Tools > Show Regression Line; go to Context Panel > Marker, set Jitter Size to 10
4. Right-click the legend, set Legend Visibility to Always; go to Context Panel > X, check Invert X Axis
5. Save the layout
6. Close the scatter plot viewer by clicking the X icon on the viewer title bar
7. Apply the saved layout
8. Verify all properties match the configuration from steps 1-4
9. Delete the saved layout

## Context menu

1. Right-click the center of the scatter plot canvas — a context menu should appear
2. Verify the menu contains items: Reset View, Lasso Tool, Tools
3. Close the menu by pressing Escape

## Log scale with negative and zero values

1. On the viewer, click the X column selector and set X to AGE, click the Y column selector and set Y to HEIGHT
2. Right-click the X Axis and set X Axis Type to logarithmic
3. Right-click the Y Axis and set Y Axis Type to logarithmic
4. On the viewer, click the X column selector and set X to RACE (categorical)
5. Verify the scatter plot handles the categorical column on a logarithmic axis gracefully

## Empty column on log scale

1. Close all and open a dataset that contains an empty numerical column (or create a DataFrame with one)
2. Add a scatter plot by clicking the Scatter Plot icon in the Toolbox
3. On the viewer, click the X column selector and set it to the empty numerical column
4. Right-click the X Axis and set X Axis Type to logarithmic
5. Verify the table remains unfiltered

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
