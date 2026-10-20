@journey @viewers @realizes:viewers.box-plot
Feature: Box plot category order
  The X axis menu's Order... opens "Order categories: <column>", whose every change the boxes show at
  once: Default, Alphabetical, Row count, Aggregated value (a column and its aggregation) and Manual
  (the categories dragged into place), each in either direction; CLOSE keeps the order, Default
  gives the column's own order back, and APPLY TO ALL VIEWERS is off while nothing orders the
  categories. The Category Order property opens the same dialog from the context panel, an order by
  another column survives that column's rename, a Bin Color keeps each box's color through a reorder,
  and with two category levels the dialog orders the boxes inside each outer group, which keeps its place
  also when a filter leaves fewer boxes than the combinations the order ranks.
  The order of the boxes is claimed through the `category <label>` areas the box plot draws under
  them. demog-1000, AGE by RACE: Asian 15 rows, Black 27, Caucasian 896, Other 62; average WEIGHT
  Asian 70.52, Other 74.16, Caucasian 79.49, Black 88.23 (as bin colors: Other light blue #6868ff,
  Caucasian near white). With SEX first, median AGE orders the
  women Asian 44, Caucasian 47, Other 48, Black 49.5 (the men keep the alphabetical order); with RACE filtered
  to Asian and Other, four boxes stay of the eight combinations ranked.
  How each mode ranks the categories (ties, missing values, the semantic properties) is the Dart
  test `d4 | Axes | Category Order` and ApiTests `Viewers: category order`; APPLY TO ALL VIEWERS is
  in category-order/category-order-across-viewers.feature. One journey; every scenario puts the
  order back to Default.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Black" area
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area

  Scenario: Order... is on the category axis only, and opens on the column's own order
    When user right-clicks on the "y axis" area of box plot viewer
    Then the open menu should not list "Order..."
    When user closes the context menu
    And user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    Then "Order categories: RACE" dialog should be visible
    And "Order by" input in "Order categories: RACE" dialog should have value "Default"
    And "Direction" input in "Order categories: RACE" dialog should be hidden
    And "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be disabled
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then "Order categories: RACE" dialog should be absent
    And no errors should have been logged

  Scenario: Row count in both directions, then Alphabetical, previewed at once; CLOSE keeps the order
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Other" area
    And "Direction" input in "Order categories: RACE" dialog should be visible
    And "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be enabled
    When user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Asian" area
    When user selects "Alphabetical" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Asian" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    Then "Order by" input in "Order categories: RACE" dialog should have value "Alphabetical"
    And "Direction" input in "Order categories: RACE" dialog should have value "Descending"
    When user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    And user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Black" area
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be disabled
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: Aggregated value orders the boxes by another column's aggregate
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Aggregated value" in "Order by" input in "Order categories: RACE" dialog
    Then "Column" input in "Order categories: RACE" dialog should be visible
    And "Aggregation" input in "Order categories: RACE" dialog should be visible
    When user picks "WEIGHT" in the column selector of the "Order categories: RACE" dialog
    And user selects "Avg" in "Aggregation" input in "Order categories: RACE" dialog
    And user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    And the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Black" area
    When user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "category Black" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Black" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: Manual order: a category dragged to the top of the list becomes the first box
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Manual" in "Order by" input in "Order categories: RACE" dialog
    Then "Direction" input in "Order categories: RACE" dialog should be hidden
    And grid in "Order categories: RACE" dialog should be visible
    When user drags the "cell 4 of RACE" area of grid in "Order categories: RACE" dialog to the "cell 1 of RACE" area
    Then the "category Other" area of box plot viewer should lie to the left of the "category Asian" area
    And the "category Asian" area of box plot viewer should lie to the left of the "category Black" area
    And "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be enabled
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged

  Scenario: The Category Order property opens the same dialog from the context panel
    Given the context panel is open
    When user clicks on settings icon of box plot viewer
    And "X Axis" category in context panel is expanded
    Then "Category Order" property in context panel should be visible
    And "Category Order" property in context panel should have value ""
    When user clicks on button in "Category Order" property in context panel
    Then "Order categories: RACE" dialog should be visible
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And "Category Order" property in context panel should have value '{"column":"RACE","aggr":"values","property":"","order":"asc","manual":[]}'
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And "Category Order" property in context panel should have value ""
    And no errors should have been logged

  Scenario: An order by another column survives that column's rename
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Aggregated value" in "Order by" input in "Order categories: RACE" dialog
    And user picks "WEIGHT" in the column selector of the "Order categories: RACE" dialog
    And user selects "Avg" in "Aggregation" input in "Order categories: RACE" dialog
    And user selects "Ascending" in "Direction" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Black" area
    When user renames "WEIGHT" column to "BODY WEIGHT"
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Black" area
    And the "category Asian" area of box plot viewer should lie to the left of the "category Other" area
    When user renames "BODY WEIGHT" column to "WEIGHT"
    And user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "category Black" area of box plot viewer should lie to the left of the "category Caucasian" area
    And no errors should have been logged

  Scenario: Bin colors stay with their boxes when the boxes are reordered
    When user sets properties of box plot viewer:
      | binColorColumnName | WEIGHT |
      | binColorAggrType   | avg    |
    Then the "Other values" area of box plot viewer should contain the color "#6868ff"
    And the "Caucasian values" area of box plot viewer should not contain the color "#6868ff"
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "Other values" area of box plot viewer should contain the color "#6868ff"
    And the "Caucasian values" area of box plot viewer should not contain the color "#6868ff"
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    And user sets "binColorColumnName" property of box plot viewer to ""
    Then no errors should have been logged

  Scenario: With two category levels the boxes are ordered inside each outer group
    When user sets properties of box plot viewer:
      | Category 1 | SEX  |
      | Category 2 | RACE |
    Then the "category F | Black" area of box plot viewer should lie to the left of the "category F | Caucasian" area
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    Then "Order categories: RACE" dialog should be visible
    And "Order by" input in "Order categories: RACE" dialog should not offer the choice "Manual"
    When user selects "Aggregated value" in "Order by" input in "Order categories: RACE" dialog
    Then the "category F | Caucasian" area of box plot viewer should lie to the left of the "category F | Other" area
    And the "category F | Other" area of box plot viewer should lie to the left of the "category F | Black" area
    And the "category F | Black" area of box plot viewer should lie to the left of the "category M | Asian" area
    And the "category M | Black" area of box plot viewer should lie to the left of the "category M | Caucasian" area
    When user filters rows where "RACE" is one of "Asian, Other"
    Then the "category F | Asian" area of box plot viewer should lie to the left of the "category F | Other" area
    And the "category F | Other" area of box plot viewer should lie to the left of the "category M | Asian" area
    And the "category M | Asian" area of box plot viewer should lie to the left of the "category M | Other" area
    When user resets the filter
    And user selects "Alphabetical" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "category F | Other" area of box plot viewer should lie to the left of the "category F | Asian" area
    And the "category F | Asian" area of box plot viewer should lie to the left of the "category M | Other" area
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "category F | Black" area of box plot viewer should lie to the left of the "category F | Caucasian" area
    When user sets properties of box plot viewer:
      | Category 1 | RACE |
      | Category 2 |      |
    Then no errors should have been logged
