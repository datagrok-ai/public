@viewers @realizes:viewers.box-plot
Feature: Box plot edge cases
  Guards for the box plot at the edges of its data and its size, on demog-1000 with a box plot of
  AGE by SEX. A table with every row deleted or every row filtered out, and a pointer moved over a
  plot shrunk below the size of its category axis, made the box plot throw on every frame or every
  mouse move; core#960 fixed them, and these scenarios hold them: no error is logged. The group
  comparison needs at least two groups: it goes when the filter leaves one SEX and comes back with
  both (help box-plot.md, "Group comparison": 2 groups give a t-test; `box_plot_group_comparison.dart`
  gates it on two categories or more).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE |
      | Category 1 | SEX |
    Then box plot viewer should be painted

  @edge-cases @realizes:core#960
  Scenario: A box plot over a table with every row deleted logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And box plot viewer should report no error
    And no errors should have been logged

  @edge-cases @realizes:core#960
  Scenario: A box plot with every row filtered out logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And box plot viewer should report no error
    And no errors should have been logged

  @edge-cases @realizes:core#960
  Scenario: The pointer over a box plot shrunk below its axes logs no error
    When user resizes box plot viewer to 170 by 150
    Then box plot viewer should not have a "x axis" area
    When user hovers over box plot viewer
    And user moves the pointer away from box plot viewer
    And user hovers over box plot viewer
    Then box plot viewer should be painted
    And no errors should have been logged

  @edge-cases
  Scenario: The group comparison goes with a single group left and comes back with two
    When user sets "Show Group Comparison" property of box plot viewer to "true"
    Then box plot viewer should have a "group comparison" area
    When user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    And box plot viewer should not have a "group comparison" area
    When user adds a categorical filter on "SEX" keeping "F,M"
    Then all rows should pass the filter
    And box plot viewer should have a "group comparison" area
    And no errors should have been logged
