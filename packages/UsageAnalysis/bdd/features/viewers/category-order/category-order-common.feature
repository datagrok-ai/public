@viewers @realizes:viewers.box-plot @realizes:viewers.scatter-plot @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.scatter-plot3d
Feature: One category order shared by every ordering viewer
  Five viewers on RACE — box plot, scatter plot, bar chart and pie chart (both under the custom sort
  type, which shows the column's order), 3D scatter plot — all follow the order that APPLY TO ALL
  VIEWERS makes the column's, and all go back to the default one when Reset takes it away. Before
  that, an order of one viewer moves that viewer only. The line chart is left out: its points are
  the rows of its aggregated frame, which follow the column's order themselves, so no area tells a
  category from its place (its order is line-chart-category-order.feature).
  demog-1000: RACE rows Asian 15, Black 27, Caucasian 896, Other 62, so by row count descending the
  order is Caucasian, Other, Black, Asian. On the scatter plot row 1 is Caucasian, row 2 Other and
  row 10 Asian. The bars are drawn from the top under Bar Sort Order desc, so the column's order
  shows bottom-up; the pie's first slice starts at 0 degrees.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    And user adds a scatter plot viewer with:
      | X | RACE |
      | Y | AGE  |
    And user adds a bar chart viewer with:
      | Split         | RACE   |
      | Bar Sort Type | custom |
    And user adds a pie chart viewer with:
      | Category      | RACE   |
      | Pie Sort Type | custom |
    And user adds a 3d scatter plot viewer with:
      | X | RACE |
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "bar Black" area of bar chart viewer should lie below the "bar Caucasian" area
    And the "start angle of Asian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"

  Scenario: An order of one viewer moves that viewer only; APPLY TO ALL VIEWERS moves them all, and Reset brings them all back
    When user picks "Order..." from the context menu of the "x axis" area of box plot viewer
    And user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    And user selects "Descending" in "Direction" input in "Order categories: RACE" dialog
    Then the "category Other" area of box plot viewer should lie to the left of the "category Black" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "bar Black" area of bar chart viewer should lie below the "bar Caucasian" area
    And the "start angle of Asian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    When user clicks on "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog
    Then "Order categories: RACE" dialog should contain text "RACE has its own order"
    And the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "category Other" area of box plot viewer should lie to the left of the "category Black" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Asian" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    And the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 10" area
    And the "bar Caucasian" area of bar chart viewer should lie below the "bar Other" area
    And the "bar Other" area of bar chart viewer should lie below the "bar Black" area
    And the "bar Black" area of bar chart viewer should lie below the "bar Asian" area
    And the "start angle of Caucasian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "start angle of Asian" reading of pie chart viewer should be between 354.59 and 354.61
    And the "x labels" reading of 3d scatter plot viewer should be "Caucasian, Other, Black, Asian"
    When user clicks on "Reset" action link in "Order categories: RACE" dialog
    Then the "category Caucasian" area of box plot viewer should lie to the left of the "category Other" area
    And the "category Black" area of box plot viewer should lie to the left of the "category Caucasian" area
    And the "bar Black" area of bar chart viewer should lie below the "bar Caucasian" area
    And the "start angle of Asian" reading of pie chart viewer should be between -0.01 and 0.01
    And the "x labels" reading of 3d scatter plot viewer should be "Asian, Black, Caucasian, Other"
    When user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then no errors should have been logged
