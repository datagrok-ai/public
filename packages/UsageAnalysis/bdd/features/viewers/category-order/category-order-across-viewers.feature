@viewers @realizes:viewers.box-plot @realizes:viewers.scatter-plot @realizes:viewers.trellis-plot @realizes:viewers.filters.categorical @realizes:viewers.grid @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.scatter-plot3d
Feature: Category order across viewers
  A viewer's order stays its own until APPLY TO ALL VIEWERS makes it the column's category order:
  then every viewer on the column follows it — a scatter plot that had none, one whose own order on
  the column is dropped, the trellis plot's splits, the filter panel's list, a bar chart and a pie
  chart under the custom sort type, a 3D scatter plot — and a filter keeps the categories it was
  given by name (the rule is the Dart test `d4 | Filters | State | a reordered column keeps the
  categorical selection by name`). The viewer it came from goes back to Default, which now shows the
  column's order, and its dialog says so with a Reset that gives the column its default order back.
  The order applied is ranked over all rows, a category the filter hides included (the Dart test
  `d4 | Axes | Category Order`). The grid's Sort > Custom... sets the column's order as well, which
  the viewers without an order of their own follow, while one with its own keeps it; CANCEL after
  APPLY puts the column back as it was, but leaves an order applied elsewhere meanwhile, and RESET
  gives it its default order (its Order reads "default"); the sort they give back is the Dart test
  `d4 | Axes | Category Order Grid`.
  demog-1000: RACE rows Asian 15, Black 27, Caucasian 896, Other 62; on the scatter plot row 1 is
  Caucasian, row 2 Other and row 10 Asian. Sorted by RACE descending, the grid shows the Other rows 2
  and 133 on top, and the Caucasian rows 1 and 3 once the column's order puts Caucasian last; with every viewer
  docked the grid is too narrow for RACE, so the rows are read on USUBJID. Asian, Black and Other kept by the filter are 104 rows.
  The bars are drawn from the top, so the column's order shows bottom-up; a slice spans 0.36 degrees per row,
  so Caucasian starts at 15.12 degrees in the default order and at 37.44 by row count ascending.
  Every scenario starts from the table again.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    And user adds a scatter plot viewer with:
      | X | RACE |
      | Y | AGE  |
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area

  Scenario: APPLY TO ALL VIEWERS makes the order the column's, and the other viewers follow it
    When user double-clicks on the "header RACE" area of grid
    Then the "sort column" reading of grid should be "RACE"
    And the "cell 2 of RACE" area of grid should lie above the "cell 133 of RACE" area
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    Given user adds a trellis plot viewer with:
      | X Column Names | RACE      |
      | Viewer Type    | Histogram |
    And user adds a bar chart viewer with:
      | Split         | RACE   |
      | Bar Sort Type | custom |
    And user adds a pie chart viewer with:
      | Category      | RACE   |
      | Pie Sort Type | custom |
    And user adds a 3d scatter plot viewer with:
      | X | RACE |
    And user opens the filter panel
    Then the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And the "x label Caucasian" area of trellis plot viewer should lie to the left of the "x label Other" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "bar Caucasian" area of bar chart viewer should lie below the "bar Other" area
    And the "start angle of Caucasian" reading of pie chart viewer should be between 15.11 and 15.13
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    Then "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be disabled
    And "Order by" input in "Order categories: RACE" dialog should have value "Default"
    And "Order categories: RACE" dialog should contain text "RACE has its own order"
    And "Category Order" property of box plot viewer should be ""
    And "X Category Order" property of scatter plot viewer should be ""
    And the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "x label Other" area of trellis plot viewer should lie to the left of the "x label Caucasian" area
    And the "x label Black" area of trellis plot viewer should lie to the left of the "x label Other" area
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Other, Caucasian"
    And the "bar Other" area of bar chart viewer should lie below the "bar Caucasian" area
    And the "start angle of Caucasian" reading of pie chart viewer should be between 37.43 and 37.45
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Other, Caucasian"
    And grid should not have a "cell 2 of USUBJID" area
    And the "cell 1 of USUBJID" area of grid should lie above the "cell 3 of USUBJID" area
    When user clicks on "Reset" action link in "Order categories: RACE" dialog
    Then "Order categories: RACE" dialog should not contain text "RACE has its own order"
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "x label Caucasian" area of trellis plot viewer should lie to the left of the "x label Other" area
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And the "bar Caucasian" area of bar chart viewer should lie below the "bar Other" area
    And the "start angle of Caucasian" reading of pie chart viewer should be between 15.11 and 15.13
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    And the "cell 2 of USUBJID" area of grid should lie above the "cell 133 of USUBJID" area
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user clicks on the "checkbox Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    Then 104 rows should pass the filter
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    Then the "categories of RACE" reading of filter panel should be "Asian, Black, Other, Caucasian"
    And the "selected categories of RACE" reading of filter panel should be "Asian, Black, Other"
    And 104 rows should pass the filter
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: Sort > Custom... in the grid orders the column; a viewer with its own order keeps it
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Alphabetical" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    Then "Order categories" dialog should be visible
    When user drags the "cell 4 of RACE" area of grid in "Order categories" dialog to the "cell 1 of RACE" area
    And user clicks on OK button in "Order categories" dialog
    Then "Order categories" dialog should be absent
    And the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    And the "category Asian" area of box plot viewer should lie to the left of the "category Black" area
    And the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And no errors should have been logged

  Scenario: CANCEL in Sort > Custom... takes back its own APPLY only, and RESET gives the column its default order back
    Given the "RACE" column is the current object
    And the context panel is open
    And "Details" pane in context panel is expanded
    Then "default" action link in Details pane in context panel should be visible
    When user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user drags the "cell 4 of RACE" area of grid in "Order categories" dialog to the "cell 1 of RACE" area
    And user clicks on APPLY button in "Order categories" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    And "custom" action link in Details pane in context panel should be visible
    When user clicks on CANCEL button in "Order categories" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    And the context panel should show "RACE"
    And "default" action link in Details pane in context panel should be visible
    When user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user drags the "cell 4 of RACE" area of grid in "Order categories" dialog to the "cell 1 of RACE" area
    And user clicks on OK button in "Order categories" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    And "custom" action link in Details pane in context panel should be visible
    When user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user clicks on RESET button in "Order categories" dialog
    Then "Order categories" dialog should be absent
    And the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    And "default" action link in Details pane in context panel should be visible
    When user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user drags the "cell 4 of RACE" area of grid in "Order categories" dialog to the "cell 1 of RACE" area
    And user clicks on APPLY button in "Order categories" dialog
    And user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user clicks on CANCEL button in "Order categories" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Other" area
    And no errors should have been logged
