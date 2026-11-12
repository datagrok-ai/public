@journey @viewers @realizes:viewers.grid
Feature: The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors
  On demog-1000 (AGE: row 1 26, row 2 30, row 3 58, row 4 45; SEX: row 1 F, row 4 M). The
  conditional editor of AGE opens with the default four bins, where 26 and 45 sit in different bins;
  Bin... asks for a number of bins and re-bins the column, and with two bins 26 and 45 share one;
  Remove row takes a bin away. The categorical editor of SEX lists its categories, a colour picked
  from a category's bar is the colour the grid paints its cells, and its Type choice then turns the
  coding Off (public#3413: Off after an edit really clears the column); Column Properties > Edit
  category colors... opens the same editor. Translated from the closed tickets
  public#3413, GROK-10142 / GROK-8862 (a conditional coding renders bins, not a gradient) and help
  visualize/viewers/grid.md "Color code columns". Every claim is the colour the grid's renderer
  resolved for a cell or the type the column holds. Each scenario opens the editor itself and puts
  back what it changed.

  A row's pattern is typed into its "Pattern N" input: a pattern that matches more cells takes them over, a
  single value colours the cells equal to it (public#1967), and a pattern that cannot be read is marked
  invalid with an "Invalid pattern" tooltip.

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
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Bin... re-bins the column into the number of bins asked for
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    Then "Edit" dialog should be visible
    When user enters "2" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be hidden
    And "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 4 of AGE" readings of grid should be the same
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Remove row takes a bin away
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user remembers the "color of cell 3 of AGE" reading of grid
    And user clicks on last Remove row button in "Color-coding: AGE" dialog
    Then the "color of cell 3 of AGE" reading of grid should not be as remembered
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: A pattern typed into a row recolours the cells it now matches, and a bad pattern is marked
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    And user enters "4" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be hidden
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user enters ">20" into "Pattern 1" text input in "Color-coding: AGE" dialog
    Then the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should be the same
    When user enters "abc" into "Pattern 1" text input in "Color-coding: AGE" dialog
    Then "Pattern 1" text input in "Color-coding: AGE" dialog should be invalid
    When user clicks on "Pattern 2" text input in "Color-coding: AGE" dialog
    And user hovers over "Pattern 1" text input in "Color-coding: AGE" dialog
    Then tooltip should contain text "Invalid pattern"
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: A single value as a pattern colours the cells equal to it (public#1967)
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    And user enters "4" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be hidden
    And the "color of cell 1 of AGE" and "color of cell 2 of AGE" readings of grid should be the same
    When user enters "30" into "Pattern 1" text input in "Color-coding: AGE" dialog
    Then the "color of cell 1 of AGE" and "color of cell 2 of AGE" readings of grid should differ
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: The categorical editor lists the categories, a colour picked for one is the colour the grid paints it, and Off clears the column
    When user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    Then "SEX" column should be color-coded categorically
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user picks "Color Coding > Edit..." from the context menu of the "header SEX" area of grid
    Then "Color-coding: SEX" dialog should be visible
    And "F" item in "Color-coding: SEX" dialog should be visible
    And "M" item in "Color-coding: SEX" dialog should be visible
    When user clicks on color bar of "M" item in "Color-coding: SEX" dialog
    And user picks the color "#9467BD" in the color picker dialog
    Then the categorical color of "M" in "SEX" column should be "#9467BD"
    And the "color of cell 4 of SEX" reading of grid should be "#9467bd"
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
