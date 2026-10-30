@viewers @realizes:viewers.scatter-plot
Feature: NX — Tools > Show regression line on a scatter plot that shows only the selected rows
  The step of `formula-lines-and-legend-ui.md` that `nx-chain.feature` could only claim as a
  property: Row Source = FilteredSelected, then Tools > Show regression line from the context
  menu. In the chain nothing is selected at that point, so the plot draws no rows and there is
  nothing to fit. Here the same two steps run on demog-1000 with rows selected, so the plot shows
  rows when the line is asked for — the `regression lines` reading counts the fits the last frame
  drew and `regression line` is the band they sweep, both absent before the command. Which rows the
  fit was computed over is not reported, so "through the selected rows" is not claimed; see
  MISSING.md.
  Of the 553 F rows, 499 have both a HEIGHT and a WEIGHT and so are drawn (`rows shown`).
  Clearing the selection then leaves the plot with no markers and takes the fit away with them
  (core#1089, fixed 2026-10-23).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName | HEIGHT |
      | yColumnName | WEIGHT |
    And user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    When user sets "Row Source" property of scatter plot viewer to "FilteredSelected"
    Then the "rows shown" reading of scatter plot viewer should be 499

  Scenario: The regression line is drawn on a scatter plot that shows only the selected rows
    Then the "regression lines" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "regression line" area
    When user picks "Tools > Show Regression Line" from the context menu of scatter plot viewer
    And user closes the context menu
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "regression line" area
    And the "rows shown" reading of scatter plot viewer should be 499
    And no errors should have been logged

  Scenario: Select > None takes the regression line away with the rows
    When user picks "Tools > Show Regression Line" from the context menu of scatter plot viewer
    And user closes the context menu
    Then the "regression lines" reading of scatter plot viewer should be 1
    When user picks "Select > None" from the top menu
    Then no rows should be selected
    And the "rows shown" reading of scatter plot viewer should be 0
    And the "regression lines" reading of scatter plot viewer should be 0
    And no errors should have been logged
