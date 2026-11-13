@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart category order
  Under Pie Sort Type custom, the Category Order property opens "Order categories", whose every
  change the slices show at once; Default gives the column's own order. Pie Sort Order reverses the
  column's order, not an order of the viewer, which carries its own direction: the Dart test
  `d4 | Axes | Custom Sort Type`. The order is read from the slices' start angles, clockwise from the
  start: a slice starts where the ones before it end.
  demog-1000, the pie chart of RACE: rows Asian 15, Black 27, Caucasian 896, Other 62, so the slices
  span 5.4, 9.72, 322.56 and 22.32 degrees. One journey; the scenario puts the sort back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4

  Scenario: Category Order orders the slices by row count, and Default gives the column's order
    Given the context panel is open
    When user sets "Pie Sort Type" property of pie chart viewer to "custom"
    And user clicks on settings icon of pie chart viewer
    And "Data" category in context panel is expanded
    And user clicks on button in "Category Order" property in context panel
    Then "Order categories: RACE" dialog should be visible
    And the "start angle of Asian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "start angle of Black" reading of pie chart viewer should be between 5.39 and 5.41
    And the "start angle of Caucasian" reading of pie chart viewer should be between 15.11 and 15.13
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "start angle of Caucasian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "start angle of Other" reading of pie chart viewer should be between 322.55 and 322.57
    And the "start angle of Black" reading of pie chart viewer should be between 344.87 and 344.89
    And the "start angle of Asian" reading of pie chart viewer should be between 354.59 and 354.61
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "start angle of Asian" reading of pie chart viewer should be between -0.01 and 0.01
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user sets "Pie Sort Type" property of pie chart viewer to "by value"
    Then no errors should have been logged
