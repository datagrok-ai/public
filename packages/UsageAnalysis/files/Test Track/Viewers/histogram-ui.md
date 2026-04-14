# Histogram tests (manual checklist)

Manual checklist. Not included in Playwright automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open the appropriate dataset (demog or SPGI as noted)
3. Add Histogram

## Current and mouse-over row indicators

1. Set Value to AGE
2. Verify **Selection > Show Current Row** is enabled -- a dot should appear on the X axis at the current row's AGE value
3. Click a row in the grid -- the dot on the histogram X axis should move to the new row's AGE
4. Hover over a row in the grid -- a second indicator dot appears for the mouse-over row
5. Disable **Show Current Row** -- the current row dot disappears
6. Disable **Show Mouse Over Row** -- the mouse-over dot disappears

## Mouse-over row group

1. Add a second viewer (e.g., scatter plot with AGE vs HEIGHT)
2. Ensure **Selection > Show Mouse Over Row Group** is enabled on the histogram
3. Hover over a cluster of points in the scatter plot -- the histogram should highlight the distribution of the hovered group
4. Move the mouse away -- highlight disappears

## Grid selection reflected in histogram

1. Go to the grid and select the first 50 rows
2. Verify the selection is reflected on the histogram (selected bins highlighted)

## Context menu: Y axis area

1. Enable **Show Y Axis**
2. Right-click on the Y axis area -- menu should show "Show Y Axis", "Controls Font", "Axis Font"
3. Verify toggling items works

---
{
  "order": 105,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
