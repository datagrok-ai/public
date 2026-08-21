---
feature: correlationplot
target_layer: manual-only
coverage_type: smoke
---

# Correlation plot tests (manual)

Manual checklist — only items not covered by the automated tests.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Correlation plot

## Color coding of cells

1. Verify cells with positive correlation are colored differently from negative correlation
2. Verify cells with stronger correlation (closer to 1 or -1) have more saturated color
3. Verify cells with near-zero correlation appear close to white

## Selection highlighting

1. Select rows in the grid where SEX = M
2. Verify the correlation plot visually reflects the selection (histograms update)
3. Clear selection -- plot should restore

## Adjust column width

1. Click and drag the edge of any column header to resize it
2. Column width adjusts smoothly without overlapping or cutting off content


---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
