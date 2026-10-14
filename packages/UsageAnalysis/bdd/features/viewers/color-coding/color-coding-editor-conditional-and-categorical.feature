@journey @viewers @realizes:viewers.grid
Feature: The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors
  On demog-1000 (AGE: row 1 26, row 2 30, row 3 58, row 4 45; SEX: row 1 F, row 4 M). The
  conditional editor of AGE opens with the default four bins, where 26 and 45 sit in different bins;
  Bin... asks for a number of bins and re-bins the column, and with two bins 26 and 45 share one; Remove row takes a bin away.
  The categorical editor of SEX lists its categories, and its Type choice turns the coding Off
  (public#3413: Off after an edit really clears the column); Column Properties > Edit category
  colors... opens the same editor. Translated from the closed tickets public#3413, GROK-10142 /
  GROK-8862 (a conditional coding renders bins, not a gradient) and help visualize/viewers/grid.md
  "Color code columns". Every claim is the colour the grid's renderer resolved for a cell or the
  type the column holds.

  Not translated (MISSING.md): a category's colour picked in the editor — the colour bar that opens
  the picker has no name and is not reachable by any phrase; the conditional rows' own value and
  colour inputs, which are unnamed too.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then "AGE" column should have no color coding
    And "SEX" column should have no color coding

  Scenario: The conditional editor opens on the default bins
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    Then "AGE" column should be color-coded conditionally
    When user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    And Type input in "Color-coding: AGE" dialog should have value "Conditional"
    And "Bin..." link in "Color-coding: AGE" dialog should be visible
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    And the "color of cell 1 of AGE" and "color of cell 4 of AGE" readings of grid should differ
    And no errors should have been logged

  Scenario: Bin... re-bins the column into the number of bins asked for
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    Then "Edit" dialog should be visible
    When user enters "2" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be hidden
    And "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 4 of AGE" readings of grid should be the same
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    And no errors should have been logged

  Scenario: Remove row takes a bin away
    When user remembers the "color of cell 3 of AGE" reading of grid
    And user clicks on last Remove row button in "Color-coding: AGE" dialog
    Then the "color of cell 3 of AGE" reading of grid should not be as remembered
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: The categorical editor lists the categories, and Off in it clears the column
    When user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    Then "SEX" column should be color-coded categorically
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user picks "Color Coding > Edit..." from the context menu of the "header SEX" area of grid
    Then "Color-coding: SEX" dialog should be visible
    And "F" item in "Color-coding: SEX" dialog should be visible
    And "M" item in "Color-coding: SEX" dialog should be visible
    When user selects "Off" in Type input in "Color-coding: SEX" dialog
    Then "SEX" column should have no color coding
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should be the same
    When user clicks on CLOSE button in "Color-coding: SEX" dialog
    Then "Color-coding: SEX" dialog should be hidden
    And "SEX" column should have no color coding
    And no errors should have been logged

  Scenario: Column Properties > Edit category colors... opens the categorical editor
    When user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    And user picks "Column Properties..." from the context menu of the "header SEX" area of grid
    Then "SEX" dialog should be visible
    When user clicks on "Edit category colors..." link in "SEX" dialog
    Then "Color-coding: SEX" dialog should be visible
    And "M" item in "Color-coding: SEX" dialog should be visible
    When user clicks on CLOSE button in "Color-coding: SEX" dialog
    And user clicks on OK button in "SEX" dialog
    Then "SEX" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header SEX" area of grid
    Then "SEX" column should have no color coding
    And no errors should have been logged
