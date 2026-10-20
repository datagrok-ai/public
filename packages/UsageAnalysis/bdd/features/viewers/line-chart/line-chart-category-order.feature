@journey @viewers @realizes:viewers.line-chart
Feature: Line chart category order
  A categorical X axis offers Order..., and the aggregated points move along X as the dialog orders
  the categories by another column's aggregate, in either direction; Default gives the column's own
  order back. A region drawn over an ordered category selects that category's rows (by row count
  Other, 62 rows, is the third point), also with whiskers drawn (where they meet their points is the
  Dart test `d4 | Axes | Category Order Property`).
  A numerical X axis offers no Order....
  demog-1000, AGE by RACE: X repeats, so the chart draws one averaged point per race, the rows of the
  aggregated frame being Asian, Black, Caucasian and Other (`point of row 1` to `4`). The dialog
  offers the first Y column, AGE, averaged by default: Asian 38.07, Caucasian 45.67, Black 46.74,
  Other 47.10. The dialog's modes are walked in box-plot-category-order.feature. One journey; every
  scenario puts the order back to Default.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | X | RACE |
      | Y | AGE  |
    Then the "aggregated" reading of line chart viewer should be "true"
    And the "point of row 2" area of line chart viewer should lie to the left of the "point of row 3" area

  Scenario: Aggregated value orders the points by the averaged Y column, in both directions
    When user picks "Order..." from the context menu of the "left edge of x axis" area of line chart viewer
    Then "Order categories: RACE" dialog should be visible
    When user selects "Aggregated value" in "Order by" input in "Order categories: RACE" dialog
    Then "Column" input in "Order categories: RACE" dialog should contain text "AGE"
    And the "point of row 1" area of line chart viewer should lie to the left of the "point of row 3" area
    And the "point of row 3" area of line chart viewer should lie to the left of the "point of row 2" area
    And the "point of row 2" area of line chart viewer should lie to the left of the "point of row 4" area
    When user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "point of row 4" area of line chart viewer should lie to the left of the "point of row 2" area
    And the "point of row 3" area of line chart viewer should lie to the left of the "point of row 1" area
    When user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    And user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "point of row 2" area of line chart viewer should lie to the left of the "point of row 3" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "x categories" reading of line chart viewer should be 4
    And no errors should have been logged

  Scenario: Whiskers stay with their points when the categories are reordered
    When user sets "whiskersType" property of line chart viewer to "Avg | Min, Max"
    And user picks "Order..." from the context menu of the "left edge of x axis" area of line chart viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "point of row 4" area of line chart viewer should lie to the left of the "point of row 3" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user sets "whiskersType" property of line chart viewer to "None"
    Then no errors should have been logged

  Scenario: A region over an ordered category selects that category's rows
    When user sets "xCategoryOrder" property of line chart viewer to '{"column":"RACE","aggr":"values","order":"asc"}'
    And user sets "annotationRegions" property of line chart viewer to '[{"type":"area","x":"RACE","y":"avg(AGE)","header":"Third","area":[[1.5,39],[2.5,39],[2.5,47.5],[1.5,47.5]]}]'
    Then the "point of row 4" area of line chart viewer should lie inside the "region Third" area
    When user clicks on the "bottom edge of region Third" area of line chart viewer
    Then only rows where "RACE" is "Other" should be selected
    When user clears the row selection
    And user sets "annotationRegions" property of line chart viewer to "[]"
    And user sets "xCategoryOrder" property of line chart viewer to ""
    Then the "point of row 2" area of line chart viewer should lie to the left of the "point of row 3" area
    And no errors should have been logged

  Scenario: A numerical X axis offers no Order...
    When user sets "X" property of line chart viewer to "AGE"
    And user right-clicks on the "left edge of x axis" area of line chart viewer
    Then the open menu should list "Invert X Axis"
    And the open menu should not list "Order..."
    When user closes the context menu
    And user sets "X" property of line chart viewer to "RACE"
    Then the "point of row 2" area of line chart viewer should lie to the left of the "point of row 3" area
    And no errors should have been logged
