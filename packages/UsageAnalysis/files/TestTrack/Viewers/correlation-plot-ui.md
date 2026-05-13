# Correlation plot tests (manual)

Manual checklist. Not included in Playwright automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Correlation plot

## Double-click and cell interaction

1. Double-click an non-diagonal cell (e.g., AGE vs HEIGHT) -- a new scatter plot viewer opens
2. Close the scatter plot
3. On the Context Panel enable **Ignore Double Click**
4. Double-click the same cell -- no scatter plot should open
5. On the Context Panel disable **Ignore Double Click**
6. Click an non-diagonal cell -- context panel updates to show correlation details for the selected pair

## Default state

1. Verify the plot shows a matrix with all numerical columns (AGE, HEIGHT, WEIGHT) on both axes
2. Verify diagonal cells show histograms for each column
3. Verify off-diagonal cells show color-coded correlation coefficients
4. Verify Pearson R values are displayed in each off-diagonal cell
5. Verify column headers are shown vertically at the top



## Color coding of cells

1. Verify cells with positive correlation are colored differently from negative correlation
2. Verify cells with stronger correlation (closer to 1 or -1) have more saturated color
3. Verify cells with near-zero correlation appear close to white

## Diagonal histogram cells

1. Verify each diagonal cell shows a histogram of that column's distribution
2. Verify the histogram reflects the current filter state

## Tooltip with scatter plot

1. Hover over an off-diagonal cell (e.g., AGE vs WEIGHT)
2. Verify a tooltip appears with an embedded scatter plot for those two columns
3. Verify the tooltip shows the correlation type and R value
4. Move to a different off-diagonal cell -- tooltip should update with the new pair

## Tooltip on row header

1. Hover over a column name in the Y axis (row header)
2. Verify a tooltip appears showing column statistics

## Selection highlighting

1. Select rows in the grid where SEX = M
2. Verify the correlation plot visually reflects the selection (histograms update)
3. Clear selection -- plot should restore

## Adjust column width

1. Click and drag the edge of any column header to resize it
2. Column width adjusts smoothly without overlapping or cutting off content

## Test Popup menu - Column hide bug

1. Right-click on the viewer > popup menu opens
2. Open Grid > Order or Hide Columns dialog
3. Select and unselect columns -- row names should not disappear (regression for #3492)

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
