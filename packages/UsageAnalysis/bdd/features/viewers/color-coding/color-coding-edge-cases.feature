@viewers @realizes:viewers.grid
Feature: Color coding edge cases
  Guards for colour coding where a coloring meets a column it does not fit, on demog-1000. Apply
  Coloring offers a picked-up coloring only to a column of a type it applies to: a categorical
  coloring leaves the menu item disabled on a numeric column, and applied to a selection of columns
  it colours the ones it fits and leaves the rest (`grid_popup_menu.dart`, "<type> type is not
  applicable" and the selected-columns branch). A linked coloring whose source column is renamed is
  kept, and one whose source column is removed paints nothing (core#1139). Bin... with fewer than
  one bin is the defect core#1171: the dialog takes 0, and the coloring is gone.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then "AGE" column should have no color coding
    And "RACE" column should have no color coding

  @edge-cases
  Scenario: A categorical coloring picked up cannot be applied to a numeric column
    When user picks "Color Coding > Categorical" from the context menu of the "header RACE" area of grid
    Then "RACE" column should be color-coded categorically
    When user picks "Color Coding > Pick Up Coloring" from the context menu of the "header RACE" area of grid
    And user right-clicks on the "header SEX" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Apply-Coloring" menu item should be enabled
    When user closes the context menu
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Apply-Coloring" menu item should be disabled
    When user closes the context menu
    Then "AGE" column should have no color coding
    And no errors should have been logged

  @edge-cases
  Scenario: A linear coloring applied to a selection of columns colours only the numeric ones
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    Then "AGE" column should be color-coded linearly
    When user picks "Color Coding > Pick Up Coloring" from the context menu of the "header AGE" area of grid
    And user clicks on the "header SEX" area of grid holding Control
    And user clicks on the "header HEIGHT" area of grid holding Control
    Then columns "SEX, HEIGHT" should be selected
    When user picks "Color Coding > Apply Coloring" from the context menu of the "header HEIGHT" area of grid
    Then "HEIGHT" column should be color-coded linearly
    And "SEX" column should have no color coding
    And no errors should have been logged

  @edge-cases @realizes:core#1139
  Scenario: A linked coloring keeps following its source column after the source is renamed
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000"
    And user colors "RACE" column linked to "AGE" column
    Then the coloring of "RACE" column should be linked to "AGE" column
    And the "color of cell 1 of RACE" and "color of cell 1 of AGE" readings of grid should be the same
    And the "color of cell 2 of RACE" and "color of cell 2 of AGE" readings of grid should be the same
    When user renames "AGE" column to "AGE2"
    Then "AGE2" column should be color-coded linearly
    And the "color of cell 2 of RACE" and "color of cell 2 of AGE2" readings of grid should be the same

  @edge-cases @known-failure:66 @realizes:core#1171
  Scenario: Bin... refuses a number of bins below one and keeps the coloring
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    Then "Edit" dialog should be visible
    When user enters "0" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be visible
    When user clicks on CANCEL button in "Edit" dialog
    Then "AGE" column should be color-coded conditionally
    And no errors should have been logged

  @edge-cases @realizes:core#1139
  Scenario: A linked coloring whose source column is removed paints nothing
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000"
    And user colors "RACE" column linked to "AGE" column
    Then the "color of cell 1 of RACE" and "color of cell 1 of AGE" readings of grid should be the same
    When user removes "AGE" column
    Then the "color of cell 1 of RACE" and "color of cell 1 of SEX" readings of grid should be the same
    And the "color of cell 2 of RACE" and "color of cell 2 of SEX" readings of grid should be the same
    And no errors should have been logged
