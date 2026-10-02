@journey @viewers @realizes:viewers.line-chart
Feature: Line chart chrome, chart types and the empty chart
  What the chart puts around the series — the two axes, the four on-viewer selectors, the overview
  strip, the left histogram, the title and the description — and what it says when it has nothing
  to draw.
  Every scenario here replaces one whose whole claim was `chartAlive()`: "a <canvas> element still
  exists under the viewer root and the error count did not go up". A canvas exists whether the
  chart drew a thing or gave up, so that check passed for a chart in an error state — the old
  date-axis scenario proves it (see line-chart-axes-and-filter.feature). The chart now reports its
  own `error`, and every piece of chrome is a hit area that is absent when it is not drawn: the
  selector readings are the auto-layout's computed result, not the look's request, so
  `x selector shown` says what the chart decided rather than what it was asked.
  Fixture: spgi-100, 100 rows, one point per row on `CAST Idea ID` (100 distinct values, so no
  aggregation) — asserted in the Background so a wrong table fails there.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then 100 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 100
    And the "x column" reading of line chart viewer should be "CAST Idea ID"
    And the "y columns" reading of line chart viewer should be "Chemical Space X"
    And the "charts" reading of line chart viewer should be 1
    And the "lines" reading of line chart viewer should be 1
    And the "markers drawn" reading of line chart viewer should be 100
    And line chart viewer should report no error
    And line chart viewer should be painted

  Scenario: With no Y column the chart says so and draws no chart box
    Then line chart viewer should have a "chart 1" area
    When user sets "yColumnNames" property of line chart viewer to ""
    Then line chart viewer should report the error "No Y columns selected"
    And line chart viewer should not have a "chart 1" area
    And the "charts" reading of line chart viewer should be 0
    And the "lines" reading of line chart viewer should be 0
    And the "markers drawn" reading of line chart viewer should be 0
    And the "y axes" reading of line chart viewer should be 0
    When user sets "yColumnNames" property of line chart viewer to "Chemical Space X"
    Then line chart viewer should report no error
    And line chart viewer should have a "chart 1" area
    And the "markers drawn" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: Chart Type redraws the same five series four ways
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the "lines" reading of line chart viewer should be 5
    When user picks "Chart Type > Area Chart" from the context menu of line chart viewer
    Then "chartTypes" property of line chart viewer should contain "Area Chart"
    And line chart viewer should have repainted
    And the "lines" reading of line chart viewer should be 5
    When user picks "Chart Type > Stacked Area Chart" from the context menu of line chart viewer
    Then "chartTypes" property of line chart viewer should contain "Stacked Area Chart"
    And line chart viewer should have repainted
    When user picks "Chart Type > Stacked Bar Chart" from the context menu of line chart viewer
    Then "chartTypes" property of line chart viewer should contain "Stacked Bar Chart"
    And line chart viewer should have repainted
    And the "lines" reading of line chart viewer should be 5
    When user picks "Chart Type > Line Chart" from the context menu of line chart viewer
    Then "chartTypes" property of line chart viewer should contain "Line Chart"
    And line chart viewer should have repainted
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "lines" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: Interpolation and line width repaint the same hundred points
    When user sets "interpolation" property of line chart viewer to "Spline"
    Then line chart viewer should have repainted
    And the "markers drawn" reading of line chart viewer should be 100
    When user sets "splineTension" property of line chart viewer to "1"
    Then line chart viewer should have repainted
    When user sets "lineWidth" property of line chart viewer to "5"
    Then line chart viewer should have more ink than before
    When user sets properties of line chart viewer:
      | lineWidth     | 1    |
      | interpolation | None |
    Then line chart viewer should have repainted
    And the "markers drawn" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: The left histogram takes a box of its own out of the plot
    Then line chart viewer should not have a "left panel" area
    And the "left panel" reading of line chart viewer should be "None"
    When user sets "leftPanel" property of line chart viewer to "Histogram"
    Then line chart viewer should have a "left panel" area
    And the "left panel" reading of line chart viewer should be "Histogram"
    And the "left panel" area of line chart viewer should be painted
    And the "plot" area of line chart viewer should be narrower than before
    When user sets "leftPanel" property of line chart viewer to "None"
    Then line chart viewer should not have a "left panel" area
    And the "left panel" reading of line chart viewer should be "None"
    And no errors should have been logged

  Scenario: The overview strip appears under the plot and Overview None takes it away
    Then line chart viewer should not have an "overview" area
    When user picks "Overview > Line Chart" from the context menu of line chart viewer
    Then "overviewType" property of line chart viewer should be "Line Chart"
    And line chart viewer should have an "overview" area
    And the "overview" area of line chart viewer should be painted
    And the "plot" area of line chart viewer should be shorter than before
    When user picks "Overview > Area Chart" from the context menu of line chart viewer
    Then "overviewType" property of line chart viewer should be "Area Chart"
    And the "overview" area of line chart viewer should have repainted
    When user picks "Overview > None" from the context menu of line chart viewer
    Then "overviewType" property of line chart viewer should be "None"
    And line chart viewer should not have an "overview" area
    And no errors should have been logged

  Scenario: Show X Axis and Show Y Axis take the axis boxes away
    Then line chart viewer should have an "x axis" area
    And line chart viewer should have a "y axis" area
    And the "y axes" reading of line chart viewer should be 1
    When user sets "showXAxis" property of line chart viewer to "false"
    Then line chart viewer should not have an "x axis" area
    And line chart viewer should have a "y axis" area
    When user sets "showYAxis" property of line chart viewer to "false"
    Then line chart viewer should not have a "y axis" area
    And the "y axes" reading of line chart viewer should be 0
    When user sets properties of line chart viewer:
      | showXAxis | true |
      | showYAxis | true |
    Then line chart viewer should have an "x axis" area
    And the "y axes" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: The four selector flags are read as the auto-layout resolved them
    Then the "x selector shown" reading of line chart viewer should be "true"
    And the "y selectors shown" reading of line chart viewer should be "true"
    And the "split selector shown" reading of line chart viewer should be "true"
    And the "aggr selector shown" reading of line chart viewer should be "false"
    When user sets "showXSelector" property of line chart viewer to "false"
    Then the "x selector shown" reading of line chart viewer should be "false"
    And the "y selectors shown" reading of line chart viewer should be "true"
    When user sets "showYSelectors" property of line chart viewer to "false"
    Then the "y selectors shown" reading of line chart viewer should be "false"
    When user sets "showSplitSelector" property of line chart viewer to "false"
    Then the "split selector shown" reading of line chart viewer should be "false"
    When user sets properties of line chart viewer:
      | showXSelector     | true |
      | showYSelectors    | true |
      | showSplitSelector | true |
    Then the "x selector shown" reading of line chart viewer should be "true"
    And the "split selector shown" reading of line chart viewer should be "true"
    And no errors should have been logged

  Scenario: Asking for the aggregation selector is not enough to get it
    Then "showAggrTypeSelector" property of line chart viewer should be "true"
    And the "aggr selector shown" reading of line chart viewer should be "false"
    And the "aggregated" reading of line chart viewer should be "false"
    When user sets "showAggrTypeSelector" property of line chart viewer to "false"
    Then the "aggr selector shown" reading of line chart viewer should be "false"
    When user sets "showAggrTypeSelector" property of line chart viewer to "true"
    Then no errors should have been logged

  Scenario: The description sits above the plot and Description Position moves it below
    When user sets "description" property of line chart viewer to "Chemical space over the idea ids"
    Then the description of line chart viewer should be above its content
    When user sets "descriptionPosition" property of line chart viewer to "Bottom"
    Then the description of line chart viewer should be below its content
    When user sets "descriptionVisibilityMode" property of line chart viewer to "Never"
    Then the description of line chart viewer should be hidden
    When user sets properties of line chart viewer:
      | descriptionVisibilityMode | Always |
      | descriptionPosition       | Top    |
      | description               |        |
    Then no errors should have been logged

  Scenario: The chart-area context menu offers the groups the chart is configured through
    When user opens the context menu of line chart viewer
    Then the open menu should list "Reset View"
    And the open menu should list "Tools"
    And the open menu should list "Data"
    And the open menu should list "Markers"
    And the open menu should list "Chart Type"
    And the open menu should list "Overview"
    And the open menu should list "Selection"
    And the open menu should list "Controls"
    When user closes the context menu
    Then no errors should have been logged
