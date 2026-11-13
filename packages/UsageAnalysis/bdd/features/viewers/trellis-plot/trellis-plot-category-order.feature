@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot inner viewer category order
  The inner viewer's category order is edited from its own tab of the trellis's context panel (X
  Category Order on the Scatter plot tab). It orders the categories of every cell at once and stays
  the inner viewer's own: the dialog offers no APPLY TO ALL VIEWERS, and the column keeps its order
  (the filter panel lists RACE as before). Default brings each cell's picture back.
  demog-1000, a trellis split by SEX with inner scatter plots of AGE by RACE; the cells are claimed
  through their picture signatures. Order... in an axis menu is not reachable inside a trellis: the
  cells draw no axes of their own (inner scatter plots and box plots alike). One journey; the
  scenario puts the order back to Default.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names |              |
      | Viewer Type    | Scatter plot |
    And user sets "xColumnName" inner property of trellis plot viewer to "RACE"
    And user sets "yColumnName" inner property of trellis plot viewer to "AGE"
    And user opens the filter panel
    Then the "cells" reading of trellis plot viewer should be 2

  Scenario: X Category Order of the inner scatter plots reorders every cell and stays local
    When user remembers the "cell signature F" reading of trellis plot viewer
    And user remembers the "cell signature M" reading of trellis plot viewer
    Given the context panel is open
    When user clicks on settings icon of trellis plot viewer
    And user clicks on "Scatter plot" tab in context panel
    And "X Axis" category in context panel is expanded
    Then "X Category Order" property in context panel should be visible
    When user clicks on button in "X Category Order" property in context panel
    Then "Order categories: RACE" dialog should be visible
    And "APPLY TO ALL VIEWERS" button in "Order categories: RACE" dialog should be absent
    When user selects "Row count" in "Order by" input in "Order categories: RACE" dialog
    Then the "cell signature F" reading of trellis plot viewer should not be as remembered
    And the "cell signature M" reading of trellis plot viewer should not be as remembered
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    When user selects "Default" in "Order by" input in "Order categories: RACE" dialog
    And user clicks on "CLOSE" button in "Order categories: RACE" dialog
    Then the "cell signature F" reading of trellis plot viewer should be as remembered
    And the "cell signature M" reading of trellis plot viewer should be as remembered
    And no errors should have been logged
