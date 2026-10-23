@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart's own filter and the Auto legend
  The Filter formula of the Data section narrows what the pie draws and leaves the table alone:
  under `${AGE} > 40` the pie counts 635 of demog-1000's rows (Caucasian 562, Other 48, Black 19,
  Asian 6) while all 1000 still pass the table's filter. Under Legend Visibility Auto the legend is
  shown while the wedges carry only their percentages and is put away once every wedge is labelled
  with its category, since it would only repeat them (GROK-20793). One journey on demog-1000 with
  a pie chart of RACE; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should show 1000 rows

  Scenario: The pie's own Filter narrows the wedges and leaves the table filter alone
    Then the "angle value of Caucasian" reading of pie chart viewer should be 896
    When user sets "Filter" property of pie chart viewer to "${AGE} > 40"
    Then pie chart viewer should show 635 rows
    And the "slices" reading of pie chart viewer should be 4
    And the "angle value of Caucasian" reading of pie chart viewer should be 562
    And the "angle value of Other" reading of pie chart viewer should be 48
    And the "angle value of Black" reading of pie chart viewer should be 19
    And the "angle value of Asian" reading of pie chart viewer should be 6
    And all rows should pass the filter
    When user sets "Filter" property of pie chart viewer to ""
    Then pie chart viewer should show 1000 rows
    And the "angle value of Caucasian" reading of pie chart viewer should be 896
    And no errors should have been logged

  Scenario: Legend Visibility Auto shows the legend while the wedges carry only percentages and hides it once they are labelled (GROK-20793)
    When user sets properties of pie chart viewer:
      | Legend Visibility | Auto  |
      | Show Label        | false |
    Then legend of pie chart viewer should be visible
    And the legend of pie chart viewer should list 4 items
    When user sets "Show Label" property of pie chart viewer to "true"
    Then legend of pie chart viewer should be hidden
    When user sets "Show Label" property of pie chart viewer to "false"
    Then legend of pie chart viewer should be visible
    And no errors should have been logged
