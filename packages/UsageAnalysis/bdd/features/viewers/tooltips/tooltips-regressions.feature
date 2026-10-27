@viewers @realizes:viewers.tooltips
Feature: Tooltip regressions from closed issues
  Closed tooltip issues, each claimed where the bug showed, on demog-1000 (spgi-100 for a viewer
  bound to another table, a three-row inline table for the Color renderer). The table's tooltip
  edited in Edit Tooltip stays what it was in a cloned view and in the original (public#1208). With
  a viewer used as the group tooltip, hovering a histogram bin leaves the grid's colour coding and
  format as they were (public#2356, public#1599). A column drawn by the Color renderer is listed in
  the row tooltip (core#826; the swatch it draws there is the renderer's, not claimed). A viewer's
  custom tooltip lists every column it was given, past the seven of the default tooltip
  (public#2297). Edit Tooltip on a viewer bound to another table lists that table's columns
  (public#373). A tooltip designed with Design custom tooltip... gives way to a column list
  confirmed in Edit Tooltip, and Tooltip > Hide then hides it (public#758, public#746). Show Column
  Names is Always on the plots, so the tooltip's columns are read by name; between two hovers the
  pointer rests above the grid.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # public#1208
  @realizes:public#1208
  Scenario: The edited table tooltip is the same in a cloned view and in the original (public#1208)
    Given user adds a scatter plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user unchecks checkbox in "Edit Tooltip" dialog
    And user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "SEX" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    And user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    And user switches to the "demog-1000" view
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  # public#2356, public#1599
  @realizes:public#2356 @realizes:public#1599
  Scenario: With a group tooltip set, hovering a histogram bin leaves the grid's colouring and format as they were (public#2356)
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000"
    And user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.00" into Custom input in Format AGE dialog
    And user clicks on OK button in Format AGE dialog
    Then Format AGE dialog should be hidden
    And the "text of cell 1 of AGE" reading of grid should be "26.00"
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user adds a scatter plot viewer
    And user adds a histogram viewer
    And user picks "Tooltip > Use as Group Tooltip" from the context menu of the "view" area of scatter plot viewer
    And user hovers over the "bin 8" area of histogram viewer
    Then scatter plot viewer in tooltip should be visible
    When user moves the pointer away from grid
    Then "AGE" column should be color-coded linearly
    And "AGE" column should have tag "format" equal to "0.00"
    And the "text of cell 1 of AGE" reading of grid should be "26.00"
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    And no errors should have been logged

  # core#826
  @realizes:core#826
  Scenario: A column drawn by the Color renderer is listed in the row tooltip (core#826)
    Given the package autostarts have completed
    And user opens a table "color-test" with:
      | id | color   | label |
      | 1  | #FF0000 | red   |
      | 2  | #00FF00 | green |
      | 3  | #0000FF | blue  |
    When user sets the semantic type of "color" column to "Color"
    Then "color" column should have semantic type "Color"
    And the "cell type of color" reading of grid should be "Color"
    When user adds a scatter plot viewer with:
      | X          | id     |
      | Y          | id     |
      | showLabels | Always |
    And user hovers over the "marker of row 2" area of scatter plot viewer
    Then the tooltip should show columns "id, color, label"
    When user moves the pointer away from grid
    Then no errors should have been logged

  # public#2297
  @realizes:public#2297
  Scenario: A viewer's custom tooltip lists every column it was given, past the default seven (public#2297)
    Given user adds a scatter plot viewer
    When user sets properties of scatter plot viewer:
      | showLabels   | Always                                                  |
      | Show Tooltip | show custom tooltip                                     |
      | Row Tooltip  | USUBJID\nAGE\nSEX\nRACE\nDIS_POP\nDEMOG\nCONTROL\nSTARTED\nSEVERITY |
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP, DEMOG, CONTROL, STARTED, SEVERITY"
    When user moves the pointer away from grid
    Then no errors should have been logged

  # public#373
  @realizes:public#373
  Scenario: Edit Tooltip on a viewer bound to another table lists that table's columns (public#373)
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    And user adds a scatter plot viewer
    When user sets "table" property of scatter plot viewer to "spgi-100"
    Then scatter plot viewer should be bound to table "spgi-100"
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    And the column list of "Edit Tooltip" dialog should start with "Id"
    And the column list of "Edit Tooltip" dialog should not list "USUBJID"
    And the column list of "Edit Tooltip" dialog should not list "AGE"
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And no errors should have been logged

  # public#758, public#746
  @realizes:public#758 @realizes:public#746
  Scenario: After Design custom tooltip, a column list picked in Edit Tooltip is the tooltip, and Hide hides it (public#758)
    Given user adds a scatter plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    And user clicks on "Design custom tooltip..." label in "Edit Tooltip" dialog
    And user clicks on "CLOSE AND APPLY" button
    Then the "demog-1000" view should be current
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then tooltip should be visible
    When user moves the pointer away from grid
    And user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user unchecks checkbox in "Edit Tooltip" dialog
    And user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "SEX" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then tooltip should be visible
    And the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    And user picks "Tooltip > Hide" from the context menu of scatter plot viewer
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then tooltip should be hidden
    And no errors should have been logged
