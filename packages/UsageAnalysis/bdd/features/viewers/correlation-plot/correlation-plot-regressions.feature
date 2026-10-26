@viewers @realizes:viewers.correlation-plot
Feature: Correlation plot regressions of closed issues
  public#3492: hiding a column in the matrix's own Order or Hide Columns dialog (the cell menu's
  Grid > Order or Hide Columns...) took the row names away with it. Now the hidden column leaves the
  matrix (its header is gone) while the pinned name column and every row header stay; showing the
  column again brings its header back.

  Scenario: Hiding a column in Order or Hide Columns keeps the row names (public#3492)
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then the "cells" reading of correlation plot viewer should be 16
    And correlation plot viewer should have a "header WEIGHT" area
    And correlation plot viewer should have a "row header HEIGHT" area
    When user picks "Grid > Order or Hide Columns..." from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then Order or Hide Columns dialog should be visible
    When user toggles the "WEIGHT" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then Order or Hide Columns dialog should be hidden
    And correlation plot viewer should not have a "header WEIGHT" area
    And correlation plot viewer should have a "pinned band" area
    And correlation plot viewer should have a "row header HEIGHT" area
    And correlation plot viewer should have a "row header AGE" area
    And correlation plot viewer should have a "row header WEIGHT" area
    When user picks "Grid > Order or Hide Columns..." from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    And user toggles the "WEIGHT" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then correlation plot viewer should have a "header WEIGHT" area
    And correlation plot viewer should have a "row header HEIGHT" area
    And no errors should have been logged
