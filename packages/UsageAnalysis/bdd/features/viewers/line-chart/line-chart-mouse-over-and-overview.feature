@journey @viewers @realizes:viewers.line-chart
Feature: Line chart mouse-over row line and the stacked bar overview
  The overview strip drawn as a stacked bar chart (Overview > Stacked Bar Chart, an overview
  type line-chart.feature does not walk) takes its box out of the plot and goes with Overview None;
  a hover over a grid row draws the vertical mouse-over row line, and Show Mouse-over Row Line
  switches it off. Same fixture as line-chart.feature: spgi-100, one point per row on CAST Idea ID.
  Not translated: Show Mouse-over Category — the line chart never sets the flag its overlay waits
  for (`LineChartCore._catDraw`), so the property has nothing to switch; see the round's report.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then the "rows shown" reading of line chart viewer should be 100
    And line chart viewer should report no error

  Scenario: Overview Stacked Bar Chart draws the overview strip under the plot
    Then line chart viewer should not have an "overview" area
    When user picks "Overview > Stacked Bar Chart" from the context menu of line chart viewer
    Then line chart viewer should have an "overview" area
    And the "overview" area of line chart viewer should be painted
    And the "plot" area of line chart viewer should be shorter than before
    When user picks "Overview > None" from the context menu of line chart viewer
    Then line chart viewer should not have an "overview" area
    And no errors should have been logged

  Scenario: A grid hover draws the mouse-over row line, and Show Mouse-over Row Line switches it off
    When user takes a snapshot of line chart viewer
    And user hovers over the "cell 3 of Id" area of grid
    Then the mouse-over row of the table should be 3
    And line chart viewer should have more ink than before
    When user moves the pointer away from grid
    And user sets "showMouseOverRowLine" property of line chart viewer to "false"
    And user takes a snapshot of line chart viewer
    And user hovers over the "cell 3 of Id" area of grid
    Then the mouse-over row of the table should be 3
    And line chart viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "showMouseOverRowLine" property of line chart viewer to "true"
    Then no errors should have been logged
