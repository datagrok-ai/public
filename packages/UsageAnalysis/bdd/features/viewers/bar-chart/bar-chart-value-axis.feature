@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart value axis range and scale
  Min and Max constrain the value axis: bars wholly below the min and bars beyond the max are
  clipped and marked by the clipped-bar indicators (GROK-19346), the value-axis scroll bar stays
  on the constrained range, a logarithmic axis re-scales the positive counts with the clipping
  intact, and clearing Min and Max restores the full range. One journey on spgi-100 with a bar
  chart counting CAST Idea ID by Primary Series Name.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a bar chart viewer with:
      | Split                       | Primary Series Name |
      | Value                       | CAST Idea ID        |
      | Value Aggr Type             | count               |
      | Show Clipped Bar Indicators | true                |
    Then the "bars" reading of bar chart viewer should be 5
    And the "clipped bars" reading of bar chart viewer should be 0

  Scenario: Min above the shortest bars clips them and the indicators show it
    When user sets "Min" property of bar chart viewer to "10"
    Then "Min" property of bar chart viewer should be "10"
    And the "clipped bars" reading of bar chart viewer should be 3
    When user sets "Show Clipped Bar Indicators" property of bar chart viewer to "false"
    Then bar chart viewer should have repainted by at least 150 pixels
    When user sets "Show Clipped Bar Indicators" property of bar chart viewer to "true"
    Then bar chart viewer should have repainted by at least 150 pixels
    And no errors should have been logged

  Scenario: Max below the tallest bar clips it too
    When user sets "Max" property of bar chart viewer to "60"
    Then "Max" property of bar chart viewer should be "60"
    And the "clipped bars" reading of bar chart viewer should be 4
    When user sets "Show Clipped Bar Indicators" property of bar chart viewer to "false"
    Then bar chart viewer should have repainted by at least 150 pixels
    When user sets "Show Clipped Bar Indicators" property of bar chart viewer to "true"
    Then bar chart viewer should have repainted by at least 150 pixels
    And no errors should have been logged

  Scenario: The value-axis scroll bar shows on the constrained range
    When user hovers over bar chart viewer
    Then x-slider range slider in bar chart viewer should be visible
    And no errors should have been logged
    When user moves the pointer away from bar chart viewer

  Scenario: A logarithmic axis keeps the clipping and raises no error
    Then "Axis Type" property of bar chart viewer should be "linear"
    When user sets "Axis Type" property of bar chart viewer to "logarithmic"
    Then "Axis Type" property of bar chart viewer should be "logarithmic"
    And bar chart viewer should have repainted
    And "Min" property of bar chart viewer should be "10"
    And "Show Clipped Bar Indicators" property of bar chart viewer should be "true"
    And the "clipped bars" reading of bar chart viewer should be 4
    And no errors should have been logged
    When user sets "Axis Type" property of bar chart viewer to "linear"
    Then bar chart viewer should have repainted
    And no errors should have been logged

  Scenario: Clearing Min and Max restores the full range
    When user sets properties of bar chart viewer:
      | Min | |
      | Max | |
    Then properties of bar chart viewer should be:
      | Min       |        |
      | Max       |        |
      | Axis Type | linear |
    And the "clipped bars" reading of bar chart viewer should be 0
    And bar chart viewer should be painted
    When user moves the pointer away from bar chart viewer
    Then x-slider range slider in bar chart viewer should be hidden
    And no errors should have been logged
