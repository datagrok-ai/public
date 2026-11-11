@journey @viewers @realizes:viewers.scatter-plot3d
Feature: 3D scatter plot category order
  A categorical axis offers Order > X Axis... (Y, Z alike), and the axis labels move as the dialog
  orders the categories; Default gives the column's own order back, and an order on one axis leaves
  the others alone. The WebGL frame has no readable pixels, so the order is read from the axis labels
  the viewer reports (`x labels`, `z labels`); a numerical axis offers no order. The direction is the
  Dart test `d4 | Axes | Category Order`, the markers a reversed axis moves `d4 | Axes | Custom Sort Type`.
  demog-1000, X RACE and Z SEX: rows Asian 15, Black 27, Caucasian 896, Other 62; average AGE Asian
  38.07, Caucasian 45.67, Black 46.74, Other 47.10. One journey; every scenario puts the order back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a 3d scatter plot viewer with:
      | X | RACE |
      | Z | SEX  |
    Then the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    And the "z labels" reading of 3d scatter plot viewer should be "F, M"

  Scenario: X Axis... orders the X labels by row count and by an aggregate, and Default gives them back
    When user picks "Order > X Axis..." from the context menu of the "view" area of 3d scatter plot viewer
    Then "Order categories: RACE" dialog should be visible
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "x labels" reading of 3d scatter plot viewer should be "Caucasian, Other, Black, Asian"
    And the "z labels" reading of 3d scatter plot viewer should be "F, M"
    When user selects "Aggregated value" in "Order by" input in "Order categories: RACE" dialog
    And user selects "AGE" in "Column" input in "Order categories: RACE" dialog
    Then the "x labels" reading of 3d scatter plot viewer should be "Other, Black, Caucasian, Asian"
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: An order on Z leaves X alone, and a numerical axis offers none
    When user sets "zCategoryOrder" property of 3d scatter plot viewer to '{"manual":["M","F"]}'
    Then the "z labels" reading of 3d scatter plot viewer should be "M, F"
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    When user sets "X" property of 3d scatter plot viewer to "AGE"
    And user right-clicks on the "view" area of 3d scatter plot viewer
    Then the open menu should list "Order > Z Axis..."
    And the open menu should not list "Order > X Axis..."
    When user closes the context menu
    And user sets "zCategoryOrder" property of 3d scatter plot viewer to ""
    And user sets "X" property of 3d scatter plot viewer to "RACE"
    Then the "z labels" reading of 3d scatter plot viewer should be "F, M"
    And no errors should have been logged
