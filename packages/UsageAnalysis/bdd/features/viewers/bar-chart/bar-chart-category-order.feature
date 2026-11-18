@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart category order
  Order > Custom Order... switches Bar Sort Type to custom and Bar Sort Order to asc (core#1516), and
  opens "Order categories", whose every change the bars show at once: row count in either direction,
  Default (the column's own order).
  Bar Sort Order reverses the column's order, not an order of the viewer, which carries its own
  direction (under every sort order, and the same for the pie chart, the Dart test
  `d4 | Axes | Custom Sort Type`). A numerical split offers no Custom Order....
  demog-1000, the bar chart split by RACE: rows Asian 15, Black 27, Caucasian 896, Other 62. The
  bars lie one under another, the first at the top, so Default shows the column's order from Asian
  down to Other. One journey; the scenario puts the sort back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
    Then the "bars" reading of bar chart viewer should be 4
    And the bars of bar chart viewer should lie one under another

  Scenario: Custom Order... orders the bars by row count in both directions, Default gives the column's order, and a numerical split offers none
    When user picks "Order > Custom Order..." from the context menu of the "bar Caucasian" area of bar chart viewer
    Then "Order categories: RACE" dialog should be visible
    And "Bar Sort Type" property of bar chart viewer should be "custom"
    And the "bar Asian" area of bar chart viewer should lie above the "bar Black" area
    And the "bar Caucasian" area of bar chart viewer should lie above the "bar Other" area
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "bar Caucasian" area of bar chart viewer should lie above the "bar Other" area
    And the "bar Other" area of bar chart viewer should lie above the "bar Black" area
    And the "bar Black" area of bar chart viewer should lie above the "bar Asian" area
    When user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    Then the "bar Asian" area of bar chart viewer should lie above the "bar Black" area
    And the "bar Other" area of bar chart viewer should lie above the "bar Caucasian" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "bar Asian" area of bar chart viewer should lie above the "bar Black" area
    And the "bar Black" area of bar chart viewer should lie above the "bar Caucasian" area
    And the "bar Caucasian" area of bar chart viewer should lie above the "bar Other" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user sets "Bar Sort Type" property of bar chart viewer to "by category"
    And user sets "Split" property of bar chart viewer to "AGE"
    And user right-clicks on the "view" area of bar chart viewer
    Then the open menu should list "Order > by value"
    And the open menu should not list "Order > Custom Order..."
    When user closes the context menu
    And user sets "Split" property of bar chart viewer to "RACE"
    Then the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged
