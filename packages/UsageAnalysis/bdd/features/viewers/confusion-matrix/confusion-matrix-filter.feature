@journey @viewers @realizes:viewers.confusion-matrix
Feature: The confusion matrix filters the way the other viewers do
  A click on a cell keeps the table on the rows of that cell, on top of the filters already there. The
  matrix re-applies it every time the table refilters, so a filter panel change keeps it, and lets go
  of it when it is closed, when its columns change or when the filter is reset (core#233). A click on
  the same cell takes the filter off, a click on another cell replaces it. On demog-1000 with SEX
  against CONTROL the "M as false" cell holds 446 rows, "F as false" 548, and 8 of the 446 are Asian.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a confusion matrix viewer with:
      | xColumnName | SEX     |
      | yColumnName | CONTROL |

  @realizes:core#233
  Scenario: Closing the matrix lets go of its filter
    When user hovers over the "cell M as false" area of confusion matrix viewer
    Then tooltip should contain text "Number of M predicted as false"
    When user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    When user clicks on close icon of confusion matrix viewer
    Then 1000 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: A second click takes the filter off, another cell replaces it
    When user adds a confusion matrix viewer with:
      | xColumnName | SEX     |
      | yColumnName | CONTROL |
    And user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    When user clicks on the "cell M as false" area of confusion matrix viewer
    Then 1000 rows should pass the filter
    When user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    When user clicks on the "cell F as false" area of confusion matrix viewer
    Then 548 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: Changing a column and resetting the filter let go of it
    When user sets "yColumnName" property of confusion matrix viewer to "DIS_POP"
    Then 1000 rows should pass the filter
    When user sets "yColumnName" property of confusion matrix viewer to "CONTROL"
    And user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    When user picks "View > Reset Filter" from the top menu
    Then 1000 rows should pass the filter
    When user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    When user clicks on settings icon of confusion matrix viewer
    And user selects "DIS_POP" in "Y" property in context panel
    Then "yColumnName" property of confusion matrix viewer should be "DIS_POP"
    And 1000 rows should pass the filter
    And confusion matrix viewer should have a "cell M as PsA" area
    And confusion matrix viewer should have a "cell PsA as PsA" area
    When user selects "CONTROL" in "Y" property in context panel
    And user clicks on the "cell M as false" area of confusion matrix viewer
    Then 446 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: The matrix's filter and the filter panel apply together, whichever comes first
    When user adds a categorical filter on "RACE" keeping "Asian"
    Then 8 rows should pass the filter
    When user picks "View > Reset Filter" from the top menu
    And user clicks on the "cell M as false" area of confusion matrix viewer
    And user adds a categorical filter on "RACE" keeping "Asian"
    Then 8 rows should pass the filter
    And no errors should have been logged
