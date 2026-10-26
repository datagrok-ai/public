@viewers @realizes:viewers.matrix-plot
Feature: Matrix plot regressions of closed issues
  core#321 (GROK-19105): switching a matrix plot's Table to a table with fewer than three numerical
  columns left stray dots on the screen, the column scroll sliders collapsed to a point. The plot now
  re-tiles to the other table's columns with every cell drawn and no scroll slider on either axis,
  and switching back gives the first table's columns and sliders again. The small table is two
  numerical columns made for the feature; the Table property is picked in the context panel.

  @realizes:core#321
  Scenario: The Table property re-tiles the matrix to a two-column table without slider dots, and back (core#321)
    Given user is logged in
    And user opens a table "two-numbers" with:
      | a | b  |
      | 1 | 10 |
      | 2 | 30 |
      | 3 | 20 |
      | 4 | 40 |
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then matrix plot viewer should be bound to table "demog-1000"
    And the "cells drawn" reading of matrix plot viewer should be 16
    And the "x scroll shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have a "x scroll slider" area
    When user remembers the "x columns" reading of matrix plot viewer
    And user clicks on settings icon of matrix plot viewer
    And user selects "two-numbers" in "Table" property in context panel
    Then matrix plot viewer should be bound to table "two-numbers"
    And the "x columns" reading of matrix plot viewer should be 2
    And the "cells drawn" reading of matrix plot viewer should be 4
    And the "blank cells" reading of matrix plot viewer should be 0
    And the "x scroll shown" reading of matrix plot viewer should be "false"
    And the "y scroll shown" reading of matrix plot viewer should be "false"
    And matrix plot viewer should not have a "x scroll slider" area
    And matrix plot viewer should not have a "y scroll slider" area
    And matrix plot viewer should report no error
    When user selects "demog-1000" in "Table" property in context panel
    Then matrix plot viewer should be bound to table "demog-1000"
    And the "x columns" reading of matrix plot viewer should be as remembered
    And the "cells drawn" reading of matrix plot viewer should be 16
    And matrix plot viewer should have a "x scroll slider" area
    And no errors should have been logged
