@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart selected and filtered rows overlays
  The Selected Rows and Filtered Rows overlays exist only for cumulative aggregations: under count,
  sum and value count a selection paints its share of each bar in the selection color and a filter
  outlines the filtered share; under min or avg neither overlay is drawn, while the selection and
  the filter themselves stay. The shares are the chart's own `selected <category>` and
  `filtered <category>` hit areas. The filtered share needs Row Source All (a chart of the
  filtered rows has nothing to outline). One journey on demog-1000 with a bar chart of AGE by RACE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split              | RACE  |
      | Value              | AGE   |
      | Value Aggr Type    | count |
      | Row Source         | All   |
      | Show Selected Rows | true  |
      | Show Filtered Rows | true  |
    Then the table should have 1000 rows
    And bar chart viewer should show no selection highlight
    And bar chart viewer should not have a "selected Asian" area

  Scenario: Under count a selection paints the selected-rows overlay on its bar
    When user selects rows where "RACE" is "Asian"
    Then only rows where "RACE" is "Asian" should be selected
    And all rows should pass the filter
    And bar chart viewer should have a "selected Asian" area
    And bar chart viewer should not have a "selected Caucasian" area
    And the "selected Asian" area of bar chart viewer should contain the color "#FF8C00"
    And bar chart viewer should show a selection highlight
    And no errors should have been logged

  Scenario: A filter shrinks the filtered share and leaves the selection alone
    Then bar chart viewer should have a "filtered Caucasian" area
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And only rows where "RACE" is "Asian" should be selected
    And bar chart viewer should have repainted by at least 100 pixels
    And the "filtered Caucasian" area of bar chart viewer should have less ink than before
    And the "filtered Caucasian" area of bar chart viewer should contain the color "#0000A0"
    And no errors should have been logged

  Scenario: Min drops both overlays and keeps the state
    When user sets "Value Aggr Type" property of bar chart viewer to "min"
    Then bar chart viewer should show no selection highlight
    And bar chart viewer should not have a "selected Asian" area
    And bar chart viewer should not have a "filtered Caucasian" area
    And only rows where "RACE" is "Asian" should be selected
    And 553 rows should pass the filter
    And no errors should have been logged

  Scenario: Sum restores the overlays for a fresh selection and filter
    When user clears the row selection
    And user resets the filter
    Then no rows should be selected
    And all rows should pass the filter
    When user sets "Value Aggr Type" property of bar chart viewer to "sum"
    Then bar chart viewer should show no selection highlight
    And bar chart viewer should not have a "selected Asian" area
    When user selects rows where "RACE" is "Asian"
    And user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And only rows where "RACE" is "Asian" should be selected
    And bar chart viewer should have a "selected Asian" area
    And the "selected Asian" area of bar chart viewer should contain the color "#FF8C00"
    And the "filtered Caucasian" area of bar chart viewer should have less ink than before
    And no errors should have been logged

  Scenario: Value count is cumulative too
    When user sets "Value Aggr Type" property of bar chart viewer to "values"
    And user selects rows where "RACE" is one of "Asian, Caucasian"
    Then only rows where "RACE" is one of "Asian, Caucasian" should be selected
    And bar chart viewer should have a "selected Asian" area
    And bar chart viewer should have a "selected Caucasian" area
    And bar chart viewer should show a selection highlight
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And bar chart viewer should have repainted by at least 100 pixels
    And only rows where "RACE" is one of "Asian, Caucasian" should be selected
    And no errors should have been logged
    When user filters rows where "SEX" is "F"

  Scenario: Average is not
    When user sets "Value Aggr Type" property of bar chart viewer to "avg"
    Then bar chart viewer should show no selection highlight
    And bar chart viewer should not have a "selected Caucasian" area
    And bar chart viewer should not have a "filtered Caucasian" area
    And only rows where "RACE" is one of "Asian, Caucasian" should be selected
    And 553 rows should pass the filter
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Show Selected Rows | false |
      | Show Filtered Rows | false |
      | Value Aggr Type    | count |
    And user clears the row selection
    And user resets the filter
    Then all rows should pass the filter
    And no rows should be selected
    And no errors should have been logged
