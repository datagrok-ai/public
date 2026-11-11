@viewers @realizes:viewers.scatter-plot @realizes:viewers.box-plot @realizes:viewers.line-chart @realizes:viewers.filters.categorical
Feature: Category order in everyday workflows
  What a user meets when they try things out: switching an axis to another column takes along an
  order that fits the new column (an aggregate of another column) and drops one that does not (the
  old column's names or counts), so the property and the dialog always show the order the axis has
  (the rule for every kind of order, and a renamed column keeping its order, are the Dart test
  `d4 | Axes | Category Order Property`); the dialog of the column switched away from closes, the
  box plot drops an order that does not fit its new category column, and a date shown by month can
  be ordered by an aggregate (it has no names or counts of its own). Several viewers on one column
  keep their own orders (APPLY TO ALL VIEWERS and Reset are in category-order-across-viewers.feature).
  Trying modes out in the dialog leaves the plot as the dialog shows it, Manual starts from the
  order shown, an order set elsewhere closes the dialog, and the order survives a layout round trip
  and a filter, while the removal of the column it sorts by ends it.
  demog-1000: row 1 is F / Caucasian, row 2 F / Other, row 4 M / Caucasian, row 10 M / Asian. RACE rows
  Asian 15, Black 27, Caucasian 896, Other 62; average WEIGHT Asian 70.52, Other 74.16, Caucasian
  79.49, Black 88.23; SEX rows F 553, M 447, average WEIGHT F 72.50, M 87.62. The dialog keeps the
  direction of the order it opens on. The line chart on RACE aggregates, its points of rows 1 to 4
  being Asian, Black, Caucasian and Other.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | RACE |
      | Y | AGE  |
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area

  Scenario: Switching the X column takes along only an order that fits the new column
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user picks "SEX" in the "x" column selector of scatter plot viewer
    Then "Order categories: RACE" dialog should be absent
    And "X Category Order" property of scatter plot viewer should be ""
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 4" area
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order by" input in "Order categories: SEX" dialog should have value "Default"
    When user selects "Aggregated value" in "Order by" input in "Order categories: SEX" dialog
    And user picks "WEIGHT" in the column selector of the "Order categories: SEX" dialog
    And user selects "Avg" in "Aggregation" input in "Order categories: SEX" dialog
    And user selects "Descending" in "Direction" input in "Order categories: SEX" dialog
    Then the "marker of row 4" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user clicks on "CLOSE" button in "Order categories: SEX" dialog
    And user picks "RACE" in the "x" column selector of scatter plot viewer
    Then "X Category Order" property of scatter plot viewer should contain "WEIGHT"
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order by" input in "Order categories: RACE" dialog should have value "Aggregated value"
    When user selects "Alphabetical" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user sets properties of scatter plot viewer:
      | X    | STARTED |
      | xMap | month   |
    And user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order categories: STARTED" dialog should be visible
    And "Order by" input in "Order categories: STARTED" dialog should not offer the choice "Alphabetical"
    When user selects "Aggregated value" in "Order by" input in "Order categories: STARTED" dialog
    And user picks "STARTED" in the column selector of the "Order categories: STARTED" dialog
    And user clicks on "CLOSE" button in "Order categories: STARTED" dialog
    And user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order by" input in "Order categories: STARTED" dialog should have value "Aggregated value"
    When user selects "Default" in "Order by" input in "Order categories: STARTED" dialog
    And user clicks on "CLOSE" button in "Order categories: STARTED" dialog
    Then no errors should have been logged

  Scenario: Viewers on one column keep their own orders
    Given user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    And user adds a line chart viewer with:
      | X | RACE |
      | Y | AGE  |
    And user opens the filter panel
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "point of row 3" area of line chart viewer should lie to the left of the "point of row 4" area
    When user sets "categoryOrder" property of box plot viewer to '{"column":"WEIGHT","aggr":"avg","order":"desc"}'
    Then the "category Black" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And the "x categories" reading of line chart viewer should be 4
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Black" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user picks "SEX" in the "Category 1" column selector of box plot viewer
    Then "Category Order" property of box plot viewer should be ""
    When user sets properties of line chart viewer:
      | X    | STARTED |
      | xMap | month   |
    And user right-clicks on the "left edge of x axis" area of line chart viewer
    Then the open menu should not list "Order..."
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Trying modes out leaves the plot as the dialog shows it, and the order outlives a layout, a filter and its column
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user selects "Manual" in "Order by" input in "Order categories: RACE" dialog
    And user drags the "cell 4 of RACE" area of grid in "Order categories: RACE" dialog to the "cell 1 of RACE" area
    Then the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user selects "Manual" in "Order by" input in "Order categories: RACE" dialog
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be disabled
    And the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user sets "X Category Order" property of scatter plot viewer to ""
    Then "Order categories: RACE" dialog should be absent
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user sets "xCategoryOrder" property of scatter plot viewer to '{"column":"WEIGHT","aggr":"avg","order":"asc"}'
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user saves the layout of the current table view
    And user sets "X Category Order" property of scatter plot viewer to ""
    Then the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user loads the saved layout
    Then "X Category Order" property of scatter plot viewer should contain "WEIGHT"
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user filters rows where "RACE" is one of "Other"
    Then scatter plot viewer should not have a "marker of row 1" area
    When user resets the filter
    Then the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    When user removes "WEIGHT" column
    Then "X" property of scatter plot viewer should be "RACE"
    And "X Category Order" property of scatter plot viewer should be ""
    And the "marker of row 10" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And no errors should have been logged
