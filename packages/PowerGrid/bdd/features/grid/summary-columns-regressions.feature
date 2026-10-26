@viewers @realizes:viewers.grid
Feature: Grid summary columns — regressions of closed issues
  Closed issues of the summary columns this package adds to the grid, each a claim that fails if the
  bug comes back. Hiding and showing another column used to drop every summary column for good, and
  a summary column itself could not be hidden and shown again (public#1093). The column is added
  from Add > Summary Columns > Sparklines as in `summary-columns.feature`; another column is hidden
  from its header menu and brought back through Order or Hide Columns, and the summary column must
  still be drawn with its renderer; the summary column itself is unticked and ticked again in the
  same dialog, which lists it with the table columns. The Background waits for the package
  autostarts, as that feature's does, and every scenario starts on a fresh demog-1000.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"

  @realizes:public#1093
  Scenario: A summary column stays when another column is hidden and shown again (public#1093)
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    Then the "cell type of Sparklines" reading of grid should be "sparkline"
    When user picks "Hide" from the context menu of the "header AGE" area of grid
    Then grid should not have a "header AGE" area
    And the "column order" reading of grid should be "USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"
    When user picks "Order or Hide Columns..." from the context menu of the "header SEX" area of grid
    Then Order or Hide Columns dialog should be visible
    When user toggles the "AGE" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then grid should have a "header AGE" area
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"
    And the "cell type of Sparklines" reading of grid should be "sparkline"
    And no errors should have been logged

  @realizes:public#1093
  Scenario: A summary column itself is hidden and shown again in Order or Hide Columns (public#1093)
    When user picks "Add > Summary Columns > Sparklines" from the context menu of the "cell 2 of USUBJID" area of grid
    Then grid should have a "header Sparklines" area
    When user picks "Order or Hide Columns..." from the context menu of the "header SEX" area of grid
    Then Order or Hide Columns dialog should be visible
    And the "Sparklines" column should be checked in the column list of Order or Hide Columns dialog
    When user toggles the "Sparklines" column in the column list of Order or Hide Columns dialog
    Then the "Sparklines" column should not be checked in the column list of Order or Hide Columns dialog
    And grid should not have a "header Sparklines" area
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user toggles the "Sparklines" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then grid should have a "header Sparklines" area
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"
    And the "cell type of Sparklines" reading of grid should be "sparkline"
    And no errors should have been logged
