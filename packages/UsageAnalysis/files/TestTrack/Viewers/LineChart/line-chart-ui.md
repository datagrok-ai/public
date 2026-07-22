Human-only visual/gesture checks. The automatable steps (rectangular selection, formula
lines, legend filtering and color consistency) have moved into the spec(s); what remains
here needs a human eye — freeform gestures, per-pane operations, and hover tooltip visuals.

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

## Custom tooltip (uses spgi-100 dataset)

Setup: Close all, open spgi-100 (not demog)

1. Add a line chart by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > X, set X to Chemist 521
3. Right-click the viewer and go to Tooltip > Edit — a dialog opens
4. Add categorical, numerical, and date columns — the available aggregations should match the column types
5. Close the dialog
6. Hover over the line chart — custom tooltip should appear with the configured columns

## Lasso selection (automatable — candidate)

1. Go to the Context Panel > X, set X to AGE; go to the Context Panel > Y, set Y to HEIGHT
2. Right-click the viewer and go to Tools — check Lasso Tool
3. Hold Shift and drag a freeform lasso area over the chart
4. Verify some rows are selected
5. Press Escape to deselect
6. Right-click the viewer and go to Tools — uncheck Lasso Tool

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv,System:AppData/Chem/tests/spgi-100.csv"]
}
