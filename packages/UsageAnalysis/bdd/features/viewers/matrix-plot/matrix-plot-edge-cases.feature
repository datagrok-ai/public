@viewers @edge-cases @realizes:viewers.matrix-plot
Feature: Matrix plot edge cases — no X columns, a constant column, a renamed column
  The matrix plot on data and settings its other features do not reach (a table filtered down to
  nothing or with every row deleted: `viewer-rules-edge-cases.feature`). With X Column Names cleared the plot lays out no
  cell and its render ends (the `_refresh` guard for an empty column list, `matrix_plot_core.dart`).
  A column holding one value in every row draws its
  cells without an error (the constant-column crashes of core#998 and core#1048 were in the
  histogram and density plot the cells host). A column renamed while the plot shows it is renamed
  in the plot's column list, as every viewer renames the columns its look names
  (`DataFrameViewer._renameColumnInLook`). demog-1000 opens the plot on AGE, HEIGHT, WEIGHT and
  STARTED, 16 cells.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then the "cells" reading of matrix plot viewer should be 16
    And the "cells drawn" reading of matrix plot viewer should be 16

  Scenario: With X Column Names cleared the plot draws no cell and logs no error
    When user sets "xColumnNames" property of matrix plot viewer to ""
    Then the "x columns" reading of matrix plot viewer should be 0
    And matrix plot viewer should not have a "cell 1,1" area
    And the "error" reading of matrix plot viewer should be ""
    And no errors should have been logged
    When user sets "xColumnNames" property of matrix plot viewer to "AGE, HEIGHT"
    Then the "cells drawn" reading of matrix plot viewer should be 8

  Scenario: A column holding one value in every row draws its cells without an error
    When user adds a calculated column "ONE_VALUE" with formula "5"
    And user sets properties of matrix plot viewer:
      | xColumnNames | AGE, ONE_VALUE |
      | yColumnNames | AGE, ONE_VALUE |
    Then the "cells drawn" reading of matrix plot viewer should be 4
    And the "cell viewer type of ONE_VALUE x ONE_VALUE" reading of matrix plot viewer should be "Histogram"
    And the "cell rows shown of ONE_VALUE x AGE" reading of matrix plot viewer should be 1000
    And the "error" reading of matrix plot viewer should be ""
    And no errors should have been logged

  Scenario: A column renamed while the plot shows it is renamed in the plot's columns
    Then matrix plot viewer should have a "cell HEIGHT x AGE" area
    When user renames "HEIGHT" column to "HEIGHT_CM"
    Then "xColumnNames" property of matrix plot viewer should contain "HEIGHT_CM"
    And matrix plot viewer should not have a "cell HEIGHT x AGE" area
    And the "cells drawn" reading of matrix plot viewer should be 16
    And matrix plot viewer should have a "cell HEIGHT_CM x AGE" area
    And no errors should have been logged
