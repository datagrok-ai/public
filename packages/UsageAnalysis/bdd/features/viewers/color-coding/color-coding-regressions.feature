@viewers @realizes:viewers.grid
Feature: Colour coding regressions from closed issues
  Closed colour coding issues, each claimed where the bug showed, on demog-1000 (row 1: SEX F, RACE
  Caucasian; row 2: RACE Other; row 4: SEX M). A colouring that came back from a saved layout is
  switched off by Color Coding > Off in the header menu (public#1071). A style picked up in one table
  with Pick Up / Apply > Pick Up is applied, colouring included, to a column of another table with
  Pick Up / Apply > Apply (public#1887), here from demog-1000 to its first 100 rows (Color Coding >
  Apply Coloring after Pick Up Coloring in the other table stays disabled: a question in
  MISSING.md). The bar chart painting the category names of a colour-coded column (public#3621)
  reproduces on the current master: see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then "SEX" column should have no color coding
    And "RACE" column should have no color coding

  # public#1071
  @realizes:public#1071
  Scenario: After the layout comes back, Off in the header menu switches the colouring off (public#1071)
    When user colors "SEX" column categorically:
      | M | #3366CC |
      | F | #CC6699 |
    And user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    Then "SEX" column should have no color coding
    When user loads the saved layout
    Then "SEX" column should be color-coded categorically
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user picks "Color Coding > Off" from the context menu of the "header SEX" area of grid
    Then "SEX" column should have no color coding
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should be the same
    And the "color of cell 1 of SEX" reading of grid should be "#ffffff"
    And no errors should have been logged

  # public#1887
  @realizes:public#1887
  Scenario: A colouring picked up in one table applies to a column of another table (public#1887)
    Given user opens demog-1000 dataset keeping the first 100 rows as "demog-100"
    And user switches to the "demog-1000" table view
    When user colors "RACE" column categorically:
      | Caucasian | #3366CC |
      | Other     | #CC6699 |
    And user picks "Pick Up / Apply > Pick Up" from the context menu of the "header RACE" area of grid
    And user switches to the "demog-100" table view
    Then "RACE" column should have no color coding
    When user picks "Pick Up / Apply > Apply" from the context menu of the "header RACE" area of grid
    Then "RACE" column should be color-coded categorically
    And the categorical color of "Caucasian" in "RACE" column should be "3366CC"
    And the categorical color of "Other" in "RACE" column should be "CC6699"
    And the "color of cell 1 of RACE" reading of grid should be "#3366cc"
    And no errors should have been logged
