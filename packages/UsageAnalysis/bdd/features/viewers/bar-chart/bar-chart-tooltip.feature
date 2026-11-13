@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart bar tooltip
  The tooltip of a hovered bar follows what the chart is configured to aggregate (bar-chart-ui
  "Tooltip configuration"): under count it names the category, its count and its rows, and after the
  value is switched to the average of WEIGHT it reports that average (79.49) instead of the count. One
  journey on demog-1000 with a bar chart counting rows by RACE (Caucasian 896 rows).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split           | RACE  |
      | Value Aggr Type | count |

  Scenario: The bar tooltip reports the category, its count and its rows
    When user hovers over the "bar Caucasian" area of bar chart viewer
    Then tooltip should contain text "Caucasian"
    And tooltip should contain text "count: 896"
    And tooltip should contain text "896 rows"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Switched to the average of WEIGHT, the bar tooltip reports the average instead of the count
    When user sets properties of bar chart viewer:
      | Value           | WEIGHT |
      | Value Aggr Type | avg    |
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then tooltip should contain text "avg: 79.49"
    And tooltip should not contain text "count: 896"
    When user moves the pointer away from grid
    And user sets properties of bar chart viewer:
      | Value Aggr Type | count |
    Then no errors should have been logged
