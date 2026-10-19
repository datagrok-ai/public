@viewers @realizes:viewers.box-plot @realizes:viewers.scatter-plot @realizes:viewers.trellis-plot @realizes:viewers.filters.categorical @realizes:viewers.grid
Feature: Category order across viewers
  A viewer's order stays its own until APPLY TO ALL VIEWERS makes it the column's category order:
  then every viewer on the column follows it — a scatter plot that had none, one whose own order on
  the column is dropped, the trellis plot's splits, the filter panel's list — and a filter keeps the
  categories it was given by name. The grid's Sort > Custom... sets the column's order as well, which
  the viewers without an order of their own follow, while one with its own keeps it; CANCEL after
  APPLY puts the column back as it was, with no order of its own (its Order reads "default").
  demog-1000: RACE rows Asian 15, Black 27, Caucasian 896, Other 62; on the scatter plot row 1 is
  Caucasian, row 2 Other and row 10 Asian. Asian, Black and Other kept by the filter are 104 rows. A
  column order cannot be taken back through the UI, so every scenario starts from the table again.

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
    Given user adds a trellis plot viewer with:
      | X Column Names | RACE      |
      | Viewer Type    | Histogram |
    And user opens the filter panel
    Then the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And the "x label Caucasian" area of trellis plot viewer should lie to the left of the "x label Other" area
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    Then "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be disabled
    And "Category Order" property of box plot viewer should be ""
    And "X Category Order" property of scatter plot viewer should be ""
    And the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "x label Other" area of trellis plot viewer should lie to the left of the "x label Caucasian" area
    And the "x label Black" area of trellis plot viewer should lie to the left of the "x label Other" area
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Other, Caucasian"
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: A filter keeps the categories it was given when the column is reordered
    Given user opens the filter panel
    Then the "selected categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user clicks on the "checkbox Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    Then 104 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Asian, Black, Other"
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "categories of RACE" reading of filter panel should be "Asian, Black, Other, Caucasian"
    And the "selected categories of RACE" reading of filter panel should be "Asian, Black, Other"
    And 104 rows should pass the filter
    And box plot viewer should not have a "category Caucasian" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Other" area
    And no errors should have been logged

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

  Scenario: CANCEL after APPLY in Sort > Custom... gives the column its default order back
    Given the "RACE" column is the current object
    And the context panel is open
    And "Details" pane in context panel is expanded
    Then "default" action link in Details pane in context panel should be visible
    When user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user drags the "cell 4 of RACE" area of grid in "Order categories" dialog to the "cell 1 of RACE" area
    And user clicks on APPLY button in "Order categories" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    When user clicks on CANCEL button in "Order categories" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    Given the "AGE" column is the current object
    And the "RACE" column is the current object
    Then the context panel should show "RACE"
    And "default" action link in Details pane in context panel should be visible
    And no errors should have been logged
