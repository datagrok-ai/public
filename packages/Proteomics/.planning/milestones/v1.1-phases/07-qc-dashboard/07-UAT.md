---
status: complete
phase: 07-qc-dashboard
source: [07-01-SUMMARY.md, 07-02-SUMMARY.md, 07-03-SUMMARY.md]
started: 2026-03-06T22:10:00Z
updated: 2026-03-07T00:30:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Open QC Dashboard from Menu (re-test)
expected: Navigate to Proteomics | Visualize | QC Dashboard... in the menu. The dashboard opens without errors, showing multiple docked QC viewer panels within the current TableView. No circular JSON error.
result: pass
note: Fixed by changing colorColumnName to stackColumnName for categorical Group column

### 2. Missing Values Heatmap (re-test)
expected: A grid viewer is docked showing a missingness matrix -- samples as columns, proteins as rows, with cells showing present/missing status. Created via DG.Viewer.grid() on the missingness DataFrame.
result: pass
note: Grid present with correct structure. All cells show "1" (no missing values in test dataset). User suggests boolean columns for better UX.

### 3. Missing Values Bar Chart (re-test)
expected: A bar chart is docked showing per-sample missing value percentages, colored by group. Created via DG.Viewer.barChart() bound to the bar data DataFrame.
result: pass
note: Fixed by using stackColumnName instead of colorColumnName for categorical Group column

### 4. Intensity Box Plots (re-test)
expected: A box plot viewer is docked showing per-sample intensity distributions, allowing comparison of sample medians and spreads. Created via DG.Viewer.boxPlot() bound to the long-format DataFrame.
result: pass

### 5. MA Plot Viewer
expected: MA plot (scatter plot) showing log2 fold-change (M) on Y-axis vs average intensity (A) on X-axis, with one point per protein.
result: pass

### 6. MA Trend Line
expected: A smaller scatter plot docked showing a moving-average trend line of M vs A.
result: pass
note: MA Trend plot present but docked to the left instead of below the MA plot (cosmetic layout issue)

### 7. Sample Correlation Viewer
expected: A correlation plot viewer is docked showing pairwise sample-to-sample correlations based on intensity values.
result: pass

### 8. CV Plot
expected: A scatter plot showing coefficient of variation (CV) per protein within replicates.
result: pass

## Summary

total: 8
passed: 8
issues: 0
pending: 0
skipped: 0

## Gaps

- truth: "QC Dashboard opens without errors showing all viewers"
  status: resolved
  reason: "Fixed: colorColumnName is for numerical columns; changed to stackColumnName for categorical Group column"
  severity: major
  test: 1

- truth: "Missing values bar chart shows per-sample missing percentages colored by group"
  status: resolved
  reason: "Fixed: stackColumnName:'Group' correctly colors bars by category without triggering numeric aggregation"
  severity: major
  test: 3

## User Feedback

- Each visualization should have a visible title (user request)
