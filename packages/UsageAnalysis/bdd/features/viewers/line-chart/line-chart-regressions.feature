@viewers @realizes:viewers.line-chart
Feature: Line chart legend regressions
  Closed line chart issues about what the split legend lists: with two Y columns that each have
  values in one split category only, Multi Axis and Y Global Scale on, the legend lists only the
  two Y / split pairs that have data (core#825, public#3463 — the calculated H_F holds HEIGHT for
  the women, W_M holds WEIGHT for the men); and under Row Source = Selected the legend lists only the
  categories of the selected rows and follows the selection as it grows (public#1651). Each
  scenario adds its own chart to demog-1000.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # core#825, public#3463
  Scenario: The split legend lists no Y / split pair that has no data (core#825)
    When user adds a calculated column "H_F" with formula "if(${SEX} == \"F\", ${HEIGHT}, null)"
    And user adds a calculated column "W_M" with formula "if(${SEX} == \"M\", ${WEIGHT}, null)"
    And user adds a line chart viewer with:
      | xColumnName      | AGE      |
      | yColumnNames     | H_F, W_M |
      | splitColumnNames | SEX      |
      | multiAxis        | true     |
      | yGlobalScale     | true     |
    Then the legend of line chart viewer should list 2 items
    And "H_F / F" legend item in legend of line chart viewer should be visible
    And "W_M / M" legend item in legend of line chart viewer should be visible
    And "H_F / M" legend item in legend of line chart viewer should be absent
    And "W_M / F" legend item in legend of line chart viewer should be absent
    And no errors should have been logged

  # public#1651
  Scenario: Under Row Source = Selected the split legend lists the categories of the selected rows (public#1651)
    When user adds a line chart viewer with:
      | xColumnName      | AGE      |
      | yColumnNames     | WEIGHT   |
      | splitColumnNames | RACE     |
      | Row Source       | Selected |
    And user selects rows where "RACE" is "Asian"
    Then the legend of line chart viewer should list 1 item
    And "Asian" legend item in legend of line chart viewer should be visible
    When user selects rows where "RACE" is one of "Asian, Black"
    Then the legend of line chart viewer should list 2 items
    And "Black" legend item in legend of line chart viewer should be visible
    And "Caucasian" legend item in legend of line chart viewer should be absent
    And no errors should have been logged
