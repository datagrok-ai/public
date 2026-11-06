@viewers @realizes:viewers.bar-chart
Feature: Bar chart known failures
  Temporary: the known failures of the bar chart while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into its topical feature as an ordinary @edge-cases scenario: core#1174 into
  bar-chart-selection-overlays.feature.

  Each scenario is a test of its own on demog-1000 with a bar chart of average AGE by RACE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split           | RACE |
      | Value           | AGE  |
      | Value Aggr Type | avg  |
    Then the "bars" reading of bar chart viewer should be 4

  @edge-cases @known-failure:25 @realizes:core#1174
  Scenario: Removing the column of the category a click filtered by gives the rows back (core#1174)
    When user sets "onClick" property of bar chart viewer to "Filter"
    And user clicks on the "bar Asian" area of bar chart viewer
    Then 15 rows should pass the filter
    When user picks "Remove" from the context menu of the "header RACE" area of grid
    Then all rows should pass the filter
