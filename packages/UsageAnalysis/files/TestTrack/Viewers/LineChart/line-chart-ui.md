Automated coverage lives in the paired spec(s) of this section; this file lists
only checks that need a human eye.

# Line chart — manual tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a line chart


## Per-chart operations in multi-axis

1. Go to the Context Panel > Y and set Y to AGE, HEIGHT, and WEIGHT
2. On the viewer, right-click the canvas and go to Controls — check Multi Axis
3. Right-click on the AGE chart area, go to AGE > Chart Type > Area — only AGE chart becomes area

Note: the per-chart 'Hide other charts' action (Y columns reduced to the one chart)
is automated in multi-axis-and-split-spec.ts; only the per-chart chart-type VISUAL
(one chart becomes area, the others stay lines) needs a human eye.

## Custom tooltip (uses spgi-100 dataset)

Setup: Close all, open spgi-100 (not demog)

1. Add a line chart by clicking the Line Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > X, set X to Chemist 521
3. Right-click the viewer and go to Tooltip > Edit — a dialog opens
4. Add categorical, numerical, and date columns — the available aggregations should match the column types
5. Close the dialog
6. Hover over the line chart — custom tooltip should appear with the configured columns

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv,System:AppData/Chem/tests/spgi-100.csv"]
}
