@viewers @realizes:viewers.bar-chart
Feature: Bar chart regressions from closed issues
  Closed bar chart issues, each claimed where the bug showed. From the bar chart audit (core#974): a
  vertical bar chart with every row filtered out draws nothing and logs nothing (the datetime Stack
  filter kept across a property change is claimed in `bar-chart-edge-cases.feature`). A datetime Split with Split Map = month comes back from a layout
  with its twelve bars (public#2562). Each scenario adds its own bar chart to demog-1000 (Caucasian
  896 rows, 426 of them started in 1990).
  Not translated: public#463 (a column picked in the on-chart value selector under count) — the
  library's selector step opens the selector's aggregation half, see MISSING-regressions.md; the audit's
  negative Bar Corner Radius (a value set from JS is an API test) and its performance half.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # core#974 (a)
  Scenario: A vertical bar chart with every row filtered out logs nothing (core#974)
    When user adds a bar chart viewer with:
      | Split       | RACE     |
      | Orientation | vertical |
    Then the "bars" reading of bar chart viewer should be 4
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And bar chart viewer should show 0 rows
    And the "bars" reading of bar chart viewer should be 0
    When user hovers over the "view" area of bar chart viewer
    Then no errors should have been logged
    When user resets the filter
    Then the "bars" reading of bar chart viewer should be 4
    And bar chart viewer should show 1000 rows
    And no errors should have been logged

  # public#2562
  Scenario: A datetime Split with Split Map month comes back from a layout with its twelve bars (public#2562)
    When user adds a bar chart viewer with:
      | Split     | STARTED |
      | Split Map | month   |
    Then the "bars" reading of bar chart viewer should be 12
    And bar chart viewer should have a "bar January" area
    When user saves the layout of the current table view
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    When user loads the saved layout
    Then "Split Map" property of bar chart viewer should be "month"
    And the "bars" reading of bar chart viewer should be 12
    And bar chart viewer should have a "bar January" area
    And bar chart viewer should report no error
    And no errors should have been logged
