@viewers @edge-cases @realizes:viewers.line-chart
Feature: Line chart edge cases — a log Y axis below zero, an empty row set, a removed column, SPC under a filter
  The line chart on data and settings its other features do not reach. `Chemical Space X` starts at
  −4.25 and 41 of spgi-100's values are not positive, so a logarithmic Y axis over it drops them as
  the logarithmic X axis does (`line-chart-axes-and-filter.feature`). A table filtered down to
  nothing takes a Shift-drag that selects nothing, a table whose rows are all deleted draws no
  point, and the Y column removed from the table while the chart shows it raises no error (core#1100
  lists it among the cases that render without errors). SPC computes its process average over the
  rows that pass the filter (`spc average` is "the filtered process mean"): with the first 50
  CAST Idea IDs kept (634783 to 634834) their Chemical Space X averages 2.497, against 3.483 over
  all 100 rows.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And the "markers drawn" reading of line chart viewer should be 100
    And line chart viewer should report no error

  Scenario: A logarithmic Y axis over a column that starts below zero drops the non-positive points
    When user sets "yAxisType" property of line chart viewer to "logarithmic"
    Then the "markers drawn" reading of line chart viewer should be 59
    And line chart viewer should report no error
    And line chart viewer should be painted
    And no errors should have been logged

  Scenario: With every row filtered out a Shift-drag over the plot selects nothing and logs no error
    When user filters rows where "CAST Idea ID" is between 1 and 2
    Then 0 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 0
    And the "markers drawn" reading of line chart viewer should be 0
    When user drags a selection box over the "plot" area of line chart viewer
    Then no rows should be selected
    And the "rows selected" reading of line chart viewer should be 0
    And no errors should have been logged

  Scenario: Once every row of the table is deleted the chart draws no point and logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "rows shown" reading of line chart viewer should be 0
    And the "markers drawn" reading of line chart viewer should be 0
    And no errors should have been logged

  Scenario: The Y column removed from the table while the chart shows it raises no error
    When user removes "Chemical Space X" column
    Then the table should not have a column "Chemical Space X"
    When user hovers over the "view" area of line chart viewer
    And user moves the pointer away from line chart viewer
    Then no errors should have been logged

  @known-failure:64 @realizes:core#1100
  Scenario: SPC computes its process average over the rows that pass the filter (core#1100)
    The average stays at the all-rows value (3.48) while the filter keeps 50 rows that average 2.50:
    SPC samples every row (core#1100, wrong behaviour 5). Remove the tag once the average follows
    the filter.
    When user sets "showStatisticalProcessControl" property of line chart viewer to "true"
    Then the "spc average" reading of line chart viewer should be between 3.4 and 3.6
    When user filters rows where "CAST Idea ID" is between 634783 and 634834
    Then 50 rows should pass the filter
    And the "spc average" reading of line chart viewer should be between 2.4 and 2.6

  Scenario: X Min and X Max outside every X value are kept and the data stays in view
    When user sets properties of line chart viewer:
      | X Min | 700000 |
      | X Max | 800000 |
    Then "X Min" property of line chart viewer should be "700000"
    And "X Max" property of line chart viewer should be "800000"
    And the "markers drawn" reading of line chart viewer should be 100
    And line chart viewer should report no error
    And no errors should have been logged
