# NX (formula lines and legend chain) — what the features cannot say yet

Written during the viewer gap review of 23 Oct 2026. What could be said is in
`nx-regression-line.feature`.

## Needs a reading

- **The regression line is fitted through the selected rows** (`formula-lines-and-legend-ui.md`,
  Row Source = FilteredSelected, then Tools > Show regression line). The feature claims the line is
  drawn (`regression lines` 1, a `regression line` area) on a plot that shows only the 499 selected
  rows that have both values, but not which rows the fit used. The fit takes
  `indexesOfNotEmpty([xCol, yCol], filter: plot.combinedFilter)` (`scatterplot/features/regression_line.dart:51`),
  so it should be the same 499 — but nothing reports the count, and the line's band moves for any
  fit. Needs a `rows fitted` reading (the length of `valueIndexes` for the last fit, summed over the
  categories) beside `regression lines` in `scatterplot/features/scatterplot_status.dart`:

  ```gherkin
  Scenario: The regression line is fitted through the selected rows the plot shows
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName | HEIGHT |
      | yColumnName | WEIGHT |
    And user selects rows where "SEX" is "F"
    When user sets "Row Source" property of scatter plot viewer to "FilteredSelected"
    And user picks "Tools > Show Regression Line" from the context menu of scatter plot viewer
    And user closes the context menu
    Then the "regression lines" reading of scatter plot viewer should be 1
    And the "rows fitted" reading of scatter plot viewer should be 499
    And the "rows fitted" and "rows shown" readings of scatter plot viewer should be the same
    And no errors should have been logged
  ```
