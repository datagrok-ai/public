@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart category order
  Order > Custom Order... switches Pie Sort Type to custom and opens "Order categories", whose every
  change the slices show at once; Default gives the column's own order. Pie Sort Order reverses the
  column's order, not an order of the viewer, which carries its own direction. The order is read from
  the slices' start angles, clockwise from the start: a slice starts where the ones before it end.
  demog-1000, the pie chart of RACE: rows Asian 15, Black 27, Caucasian 896, Other 62, so the slices
  span 5.4, 9.72, 322.56 and 22.32 degrees. One journey; every scenario puts the sort back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4

  Scenario: Custom Order... orders the slices by row count in both directions, and Default gives the column's order
    When user picks "Order > Custom Order..." from the context menu of the "slice Caucasian" area of pie chart viewer
    Then "Order categories: RACE" dialog should be visible
    And "Pie Sort Type" property of pie chart viewer should be "custom"
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

  Scenario: Pie Sort Order reverses the column's order, but not an order of the viewer
    When user sets properties of pie chart viewer:
      | Pie Sort Type  | custom |
      | Pie Sort Order | desc   |
    Then the "start angle of Other" reading of pie chart viewer should be between -0.01 and 0.01
    And the "start angle of Caucasian" reading of pie chart viewer should be between 22.31 and 22.33
    When user sets "categoryOrder" property of pie chart viewer to '{"column":"RACE","aggr":"values","order":"desc"}'
    Then the "start angle of Caucasian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "start angle of Asian" reading of pie chart viewer should be between 354.59 and 354.61
    When user sets "categoryOrder" property of pie chart viewer to ""
    And user sets properties of pie chart viewer:
      | Pie Sort Type  | by value |
      | Pie Sort Order | asc      |
    Then the "slices" reading of pie chart viewer should be 4
    And no errors should have been logged
