@journey @viewers @realizes:viewers.grid
Feature: Extract keeps how the columns look in the grid
  Extract opens the columns in a new table. Tags such as the format travel with the column, while a
  Back Color, a width and the rest of a grid column's settings live on the grid, and the command
  carries them over from the grid it was run from, for one column or several (core#231). What a
  column looks like is claimed by the pixels of the new grid. A Back Color is set through its API
  step; the Style section's colour picker is not what this feature is about.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  @realizes:core#231
  Scenario: Extract keeps the look of several selected columns
    When user sets the back color of the cells of "HEIGHT" column to "#F5F5A3"
    And user sets the back color of the header of "HEIGHT" column to "#904040"
    And user sets the back color of the cells of "WEIGHT" column to "#C8E6C9"
    And user drags the "column resizer HEIGHT" area of grid by 60 pixels to the right
    And user clicks on the "header HEIGHT" area of grid holding Control
    And user clicks on the "header WEIGHT" area of grid holding Control
    Then columns "HEIGHT, WEIGHT" should be selected
    Given "Actions" section in context panel is expanded
    When user clicks on "Extract" label in context panel
    Then the open table views should be exactly "demog-1000, result"
    When user switches to the "result" table view
    Then the table should have the columns "HEIGHT, WEIGHT"
    And the "cell 1 of HEIGHT" area of grid should contain the color "#F5F5A3"
    And the "header HEIGHT" area of grid should contain the color "#904040"
    And the "cell 1 of WEIGHT" area of grid should contain the color "#C8E6C9"
    And the "column width of HEIGHT" and "column width of WEIGHT" readings of grid should differ
    And no errors should have been logged

  @realizes:core#231
  Scenario: Extract keeps the look of a single column
    When user closes all views
    And user opens demog-1000 dataset
    And user sets the back color of the cells of "AGE" column to "#F5F5A3"
    And user picks "Extract" from the context menu of the "header AGE" area of grid
    Then the open table views should be exactly "demog-1000, result"
    When user switches to the "result" table view
    Then the table should have the columns "AGE"
    And the "cell 1 of AGE" area of grid should contain the color "#F5F5A3"
    And no errors should have been logged
