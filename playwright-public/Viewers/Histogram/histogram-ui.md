---
feature: histogram
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Automated coverage lives in the paired spec(s) of this section; this file
  lists only checks that need a human eye.
---

# Histogram tests (manual checklist)

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

## Pick Up / Apply

1. Customize the histogram: set Value to HEIGHT and change the number of bins
2. Add a second Histogram with default settings
3. Right-click the customized histogram and select Pick Up / Apply > Pick Up
4. Right-click the second histogram and select Pick Up / Apply > Apply — the second viewer now matches the first (value column, bins)

---
{
  "order": 105,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/chem/SPGI.csv"]
}
