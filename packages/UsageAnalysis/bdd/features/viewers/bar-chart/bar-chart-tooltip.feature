@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart bar tooltip
  The tooltip of a hovered bar follows what the chart is configured to aggregate (bar-chart-ui
  "Tooltip configuration"): under count it names the category, its count and its rows, and after the
  value is switched to the average of WEIGHT it reports that average (79.49) instead of the count. A
  category whose name holds markup characters (`a<b`, `x &amp; y`, on a table of its own) is shown in
  the tooltip as typed, not as HTML, which would show `a` and `x & y` (core#974). Closing a bar chart gives back the rows its click
  filtered out: a filter that outlives its viewer is a defect (core#233). One journey on demog-1000
  with a bar chart counting rows by RACE (Caucasian 896 rows, Asian 15).

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

  Scenario: Closing the bar chart gives back the rows its Filter click filtered out
    When user sets "onClick" property of bar chart viewer to "Filter"
    And user clicks on the "bar Asian" area of bar chart viewer
    Then 15 rows should pass the filter
    When user clicks on close icon of bar chart viewer
    Then the open tableview should have 0 bar chart viewers
    And all rows should pass the filter
    And no errors should have been logged

  @realizes:core#974
  Scenario: A category with markup characters is shown in the bar tooltip as typed (core#974)
    Given user opens a table "markup-categories" with:
      | Category  | Value |
      | a<b       | 1     |
      | a<b       | 2     |
      | x &amp; y | 3     |
    And user adds a bar chart viewer with:
      | Split | Category |
    Then the "bars" reading of bar chart viewer should be 2
    When user hovers over the "bar a<b" area of bar chart viewer
    Then tooltip should contain text "a<b"
    And tooltip should contain text "2 rows"
    When user hovers over the "bar x &amp; y" area of bar chart viewer
    Then tooltip should contain text "x &amp; y"
    And no errors should have been logged
