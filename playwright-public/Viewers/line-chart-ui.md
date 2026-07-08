Manual checklist. Not included in Playwright automation.

# Line chart — manual tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a line chart


## Per-chart operations in multi-axis

1. Go to the Context Panel > Y and set Y to AGE, HEIGHT, and WEIGHT
2. On the viewer, right-click the canvas and go to Controls — check Multi Axis
3. Right-click on the AGE chart area, go to AGE > Chart Type > Area — only AGE chart becomes area
4. Right-click on the HEIGHT chart area, go to HEIGHT > Hide other charts — only HEIGHT remains
5. Verify Y columns reduced to HEIGHT only

## Custom tooltip (uses SPGI dataset)

Setup: Close all, open SPGI (not demog)

1. Add a line chart by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > X, set X to Chemist 521
3. Right-click the viewer and go to Tooltip > Edit — a dialog opens
4. Add categorical, numerical, and date columns — the available aggregations should match the column types
5. Close the dialog
6. Hover over the line chart — custom tooltip should appear with the configured columns

## Rectangular selection

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to HEIGHT
2. Hold Shift and drag a rectangle over the center area of the chart
3. Verify some rows are selected
4. Press Escape to deselect
5. Hold Shift and drag a rectangle over a different area
6. Verify some rows are selected
7. Press Escape to deselect

## Lasso selection

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to HEIGHT
2. Right-click the viewer and go to Tools — check Lasso Tool
3. Hold Shift and drag a freeform lasso area over the chart
4. Verify some rows are selected
5. Press Escape to deselect
6. Right-click the viewer and go to Tools — uncheck Lasso Tool

## Formula lines (uses SPGI dataset)

Setup: Close all, open SPGI (not demog)

1. Add a line chart
2. Go to the Context Panel > Data, set Split to Stereo Category
3. Go to the Context Panel > X and set X to Chemical Space X; go to the Context Panel > Y and set Y to Average Mass
4. Right-click the Y Axis and set Y Axis Type to logarithmic
5. On the viewer, right-click the canvas and select Tools > Formula Lines
6. In the Formula Lines editor, click ADD NEW and add a line with formula `${Average Mass} = 0.75* ${Chemical Space X}* ${Chemical Space X} - 4 * ${Chemical Space X} +300`
7. Change the line Color and Style — verify the line updates on the viewer
8. Right-click the Y Axis and set Y Axis Type to linear — verify the formula line adjusts
9. Switch Y Axis Type back to logarithmic — verify the formula line adjusts
10. Go to the grid, right-click the header of Chemical Space X and rename the column — verify the formula line caption updates
11. Delete the created line

## Legend filtering and color consistency (uses SPGI dataset)

Setup: Close all, open SPGI (not demog)

1. Add a line chart
2. Go to the Context Panel > Data, set Split to Primary Series Name, Series, and Core
3. Verify that legend colors are distinct and match the line colors on the chart
4. Open the filter panel, in the Primary Series Name filter toggle some checkboxes off and on — verify legend colors remain consistent with the corresponding line colors
5. On the legend, Ctrl+click a category — verify the line chart filters to show only that category
6. Ctrl+click another category — verify both categories are now shown
7. Ctrl+click the first category again to deselect it — verify only the second category remains

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
