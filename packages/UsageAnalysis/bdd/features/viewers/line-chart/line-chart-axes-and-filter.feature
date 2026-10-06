@journey @viewers @realizes:viewers.line-chart
Feature: Line chart axes, date bucketing and following the filter
  What the X axis spans, how a date column is bucketed into it, and what happens to both when the
  table is filtered down — including to nothing at all.
  The old spec derived the X window with two `screenToWorld` probes and a linear extrapolation off
  the widest canvas it could find, then polled the span; the chart publishes `x axis min`,
  `x axis max` and `x axis span` off the range slider it actually drew. The date-bucketing half is
  the reason the reading matters: the old scenario set `xMap` to "Year quarter", asserted the
  property read "Year quarter" back, and passed — but "Year quarter" is not one of the choices
  (`ValueFunction.dateCategorizations` is all lower case), so the chart it was measuring was in the
  "No X column selected" state for the whole GROK-18375 scenario. `x categories` and `aggregated`
  say whether the buckets were built.
  Fixture: spgi-100. `CAST Idea ID` runs 634783…634885, 100 distinct; `Competition assay Date` has
  20 blanks and ~50 distinct dates over 2017-03…2019-06, which bucket to 3 years, 12 months and 11
  year-quarters; `Chemical Space X` starts at −4.25, so a logarithmic X axis over it is the
  non-positive path.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And "axesFollowFilter" property of line chart viewer should be "true"
    And the "x axis min" reading of line chart viewer should be between 634780 and 634782
    And the "x axis max" reading of line chart viewer should be between 634886 and 634888
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And the "markers drawn" reading of line chart viewer should be 100
    And line chart viewer should report no error

  Scenario: Axes Follow Filter pulls the X axis onto the filtered rows
    When user adds a range filter on "CAST Idea ID" from 634800 to 634850
    Then 49 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 49
    And the "markers drawn" reading of line chart viewer should be 49
    And the "x axis span" reading of line chart viewer should be 52
    And the "x axis min" reading of line chart viewer should be between 634798 and 634800
    And the "x axis max" reading of line chart viewer should be between 634850 and 634852
    And line chart viewer should have repainted
    When user hovers over "CAST Idea ID" filter card
    And user clicks on close of "CAST Idea ID" filter card
    Then 100 rows should pass the filter
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And no errors should have been logged

  Scenario: With Axes Follow Filter off the same filter leaves the axis where it was
    When user sets "axesFollowFilter" property of line chart viewer to "false"
    And user adds a range filter on "CAST Idea ID" from 634800 to 634850
    Then 49 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 49
    And the "markers drawn" reading of line chart viewer should be 49
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    When user sets "axesFollowFilter" property of line chart viewer to "true"
    Then the "x axis span" reading of line chart viewer should be 52
    When user hovers over "CAST Idea ID" filter card
    And user clicks on close of "CAST Idea ID" filter card
    Then 100 rows should pass the filter
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And no errors should have been logged

  Scenario: A logarithmic X axis over a column that starts below zero drops the non-positive points
    When user sets "xColumnName" property of line chart viewer to "Chemical Space X"
    Then the "markers drawn" reading of line chart viewer should be 100
    And the "x axis min" reading of line chart viewer should be between -4.7 and -4.6
    When user sets "xAxisType" property of line chart viewer to "logarithmic"
    Then the "markers drawn" reading of line chart viewer should be 59
    And the "x axis min" reading of line chart viewer should be between 0 and 1
    And line chart viewer should report no error
    And line chart viewer should be painted
    When user sets properties of line chart viewer:
      | xAxisType   | linear       |
      | xColumnName | CAST Idea ID |
    Then the "markers drawn" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: Hovering an empty logarithmic chart raises nothing (github-2574)
    When user sets properties of line chart viewer:
      | xColumnName | Chemical Space X |
      | xAxisType   | logarithmic      |
    And user filters rows where "TPSA" is between 0 and 1
    Then 0 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 0
    And the "lines" reading of line chart viewer should be 0
    And the "markers drawn" reading of line chart viewer should be 0
    When user hovers over the "plot" area of line chart viewer
    Then no errors should have been logged
    And no error or warning balloon should have been shown
    When user resets the filter
    Then 100 rows should pass the filter
    And the "markers drawn" reading of line chart viewer should be 59
    When user sets properties of line chart viewer:
      | xAxisType   | linear       |
      | xColumnName | CAST Idea ID |
    Then no errors should have been logged

  Scenario: A date column buckets into as many X positions as the mapping asks for
    When user sets "xColumnName" property of line chart viewer to "Competition assay Date"
    Then the "aggregated" reading of line chart viewer should be "true"
    And the "x categories" reading of line chart viewer should be 50
    And the "x column" reading of line chart viewer should be "Competition assay Date"
    When user sets "xMap" property of line chart viewer to "year"
    Then the "x column" reading of line chart viewer should be "Competition assay Date year"
    And the "x categories" reading of line chart viewer should be 3
    And the "markers drawn" reading of line chart viewer should be 3
    And line chart viewer should have repainted
    When user sets "xMap" property of line chart viewer to "month"
    Then the "x categories" reading of line chart viewer should be 12
    When user sets "xMap" property of line chart viewer to "year quarter"
    Then the "x categories" reading of line chart viewer should be 11
    And the "aggregation" reading of line chart viewer should be "avg"
    And line chart viewer should report no error
    When user sets properties of line chart viewer:
      | xMap        |              |
      | xColumnName | CAST Idea ID |
    Then the "aggregated" reading of line chart viewer should be "false"
    And no errors should have been logged

  Scenario: Filtering while the X axis is year-quarter buckets keeps the chart drawing (GROK-18375)
    When user sets properties of line chart viewer:
      | xColumnName | Competition assay Date |
      | xMap        | year quarter           |
    Then the "x categories" reading of line chart viewer should be 11
    When user filters rows where "Series" is one of "Aminopiperidines, Pyrrolidines"
    Then 25 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 25
    And the "x categories" reading of line chart viewer should be lower than before
    And the "markers drawn" reading of line chart viewer should be at least 1
    And line chart viewer should report no error
    And line chart viewer should have repainted
    When user resets the filter
    Then 100 rows should pass the filter
    And the "x categories" reading of line chart viewer should be 11
    When user sets properties of line chart viewer:
      | xMap        |              |
      | xColumnName | CAST Idea ID |
    Then no errors should have been logged

  Scenario: Narrowing the X column to its middle half keeps 29 rows (GROK-20185)
    When user sets "xColumnName" property of line chart viewer to "Chemical Space X"
    Then the "rows shown" reading of line chart viewer should be 100
    When user adds a range filter on "Chemical Space X" from 0.5731 to 10.2274
    Then 29 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 29
    And the "markers drawn" reading of line chart viewer should be 29
    And the "x axis span" reading of line chart viewer should be lower than before
    And line chart viewer should report no error
    When user hovers over "Chemical Space X" filter card
    And user clicks on close of "Chemical Space X" filter card
    Then 100 rows should pass the filter
    And the "markers drawn" reading of line chart viewer should be 100
    When user sets "xColumnName" property of line chart viewer to "CAST Idea ID"
    Then no errors should have been logged

  Scenario: Inverting the axis and dropping the grid lines redraw without moving the range
    When user sets "invertXAxis" property of line chart viewer to "true"
    Then line chart viewer should have repainted
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    And the "markers drawn" reading of line chart viewer should be 100
    When user sets "showVerticalGridLines" property of line chart viewer to "false"
    Then the "plot" area of line chart viewer should have less ink than before
    When user sets "showHorizontalGridLines" property of line chart viewer to "false"
    Then the "plot" area of line chart viewer should have less ink than before
    When user sets properties of line chart viewer:
      | invertXAxis             | false |
      | showVerticalGridLines   | true  |
      | showHorizontalGridLines | true  |
    Then the "plot" area of line chart viewer should have more ink than before
    And no errors should have been logged

  Scenario: MinMax tickmarks leave the axis with fewer labels than Auto
    When user sets "xAxisTickmarksMode" property of line chart viewer to "MinMax"
    Then the "x axis" area of line chart viewer should have less ink than before
    And the "x axis span" reading of line chart viewer should be between 106 and 107
    When user sets "yAxisTickmarksMode" property of line chart viewer to "MinMax"
    Then the "y axis" area of line chart viewer should have less ink than before
    When user sets properties of line chart viewer:
      | xAxisTickmarksMode | Auto |
      | yAxisTickmarksMode | Auto |
    Then the "x axis" area of line chart viewer should have more ink than before
    And no errors should have been logged
