@journey
Feature: Fit statistics of a curve cell become columns
  Translated from TestTrack Curves/curves-statistics.md, scenarios 1 and 2. The Fit pane of a curve
  cell lists the statistics of each series; the "+" next to one (its tooltip says "Calculate <name> for
  the whole column") adds that statistic for the whole column, right after the curve column. For a
  cell with several series, Series "aggregated" and an Aggregation type summarise them first.

  Not translated: toggling an outlier by clicking a point (scenario 4) — the points have no names a
  test can click; the statistics drawn on the plot (scenario 5), a picture. The series without a name
  (scenario 3) is in curve-json-columns.feature.

  Each scenario opens the table afresh and reads the panel after one click: a click on a second curve
  cell often leaves the plain table cell in the Context Panel instead of the curve panes.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And user opens curves dataset
    And the context panel is open

  Scenario: The IC50 of a one-curve column is added right after it, filled with numbers
    When user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 1 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    Given "Fit" pane in context panel is expanded
    Then "IC50" input in context panel should be visible
    And "R²" input in context panel should be visible
    When user clicks on first "Calculate ic50 for the whole column" icon in context panel
    Then the table should have the columns "smiles, multiple prefit, multiple styled series, styled proprtional with IC50, styled proprtional with IC50 1 ic50, log min/max candlesticks, custom fit function"
    And "styled proprtional with IC50 1 ic50" column should have type "double"
    And "styled proprtional with IC50 1 ic50" column should have at least 10 distinct values
    And no errors should have been logged

  Scenario: The average IC50 of a two-curve column is added right after it
    Given user closes all views
    And user opens curves dataset
    And the context panel is open
    When user scrolls the grid to the "multiple styled series" column
    And user clicks on the "cell 2 of multiple styled series" area of grid
    Then the current column should be "multiple styled series"
    Given "Fit" pane in context panel is expanded
    Then "Series" input in context panel should offer the choice "all"
    And "Series" input in context panel should offer the choice "aggregated"
    When user selects "aggregated" in "Series" input in context panel
    And user selects "avg" in "Aggregation type" input in context panel
    And user clicks on first "Calculate ic50 avg for the whole column" icon in context panel
    Then the table should have the columns "smiles, multiple prefit, multiple styled series, multiple styled series avg ic50, styled proprtional with IC50, log min/max candlesticks, custom fit function"
    And "multiple styled series avg ic50" column should have type "double"
    And "multiple styled series avg ic50" column should have at least 10 distinct values
    And no errors should have been logged
