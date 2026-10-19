@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot category order
  Order... in the menu of a categorical X or Y axis opens "Order categories: <column>", and the
  markers move to their category's new place as the dialog changes: Row count in both directions on
  X, Row count on Y, and an order kept while Pack and Zoom by Filter packs the categories the filter
  leaves. A numerical axis offers no Order....
  The place of a category is claimed through the markers of rows the plot reports (the first 20
  drawn): on demog-1000 row 1 is Caucasian (896 rows), row 2 Other (62) and row 10 Asian (15); Black
  (27) has no row among them. The column's own order puts them Asian, Caucasian, Other; by row count
  they go Asian, Other, Caucasian. With Asian filtered out and the categories packed, the X axis
  holds the three that are left. The dialog's modes are walked in box-plot-category-order.feature;
  the ranking rules are the Dart test `d4 | Axes | Category Order`. One journey; every scenario puts
  the order back to Default.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | RACE |
      | Y | AGE  |
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area

  Scenario: Row count on the X axis in both directions
    When user right-clicks on the "y axis" area of scatter plot viewer
    Then the open menu should not list "Order..."
    When user closes the context menu
    And user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order categories: RACE" dialog should be visible
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    When user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    And user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: Row count on a categorical Y axis
    When user sets properties of scatter plot viewer:
      | X | AGE  |
      | Y | RACE |
    Then the "marker of row 2" area of scatter plot viewer should lie above the "marker of row 1" area
    When user right-clicks on the "left edge of x axis" area of scatter plot viewer
    Then the open menu should not list "Order..."
    When user closes the context menu
    And user picks "Order..." from the context menu of the "y axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 1" area of scatter plot viewer should lie above the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie above the "marker of row 10" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "marker of row 2" area of scatter plot viewer should lie above the "marker of row 1" area
    When user sets properties of scatter plot viewer:
      | X | RACE |
      | Y | AGE  |
    Then no errors should have been logged

  Scenario: The order holds while Pack and Zoom by Filter packs the categories left by the filter
    When user sets "Zoom And Filter" property of scatter plot viewer to "pack and zoom by filter"
    And user filters rows where "RACE" is one of "Black, Caucasian, Other"
    Then scatter plot viewer should not have a "marker of row 10" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "x axis max" reading of scatter plot viewer should be 2.1
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user resets the filter
    And user sets "Zoom And Filter" property of scatter plot viewer to "filter by zoom"
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And no errors should have been logged
