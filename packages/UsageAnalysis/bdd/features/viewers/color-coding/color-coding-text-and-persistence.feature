@journey @viewers @realizes:viewers.grid
Feature: A coloring applied to the text survives Extract and a saved layout, and the grid-wide All colours a boolean column
  On demog-1000. AGE gets a linear coloring from the header menu and Apply to = text in its editor,
  which sets the column's text flag (the renderer's text colour has no reading, see MISSING.md).
  Extract on the column copies it into a new table that keeps the coloring on the text
  (GROK-18365), and a saved layout applied after the coloring was switched off brings the text
  coloring back and repaints the cells (GROK-13668). Grid Color Coding > All colours the boolean
  CONTROL column too (GROK-11394); the rest of that ladder is `grid/grid-appearance.feature`'s.
  Translated from those closed tickets; Apply to on its own is `color-coding-editor-linear.feature`.
  Each scenario puts back what it changed.

  Background:
    Given user is logged in
    And the layouts named "bdd-cc-text-{run}" are deleted when the feature ends
    And user opens demog-1000 dataset
    Then "AGE" column should have no color coding

  Scenario: Extract keeps the coloring on the text
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    And user selects "text" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    And the text of "AGE" column should be color-coded
    When user picks "Extract" from the context menu of the "header AGE" area of grid
    Then the open table views should be exactly "demog-1000, result"
    When user switches to the "result" table view
    Then the table should have the columns "AGE"
    And "AGE" column should be color-coded linearly
    And the text of "AGE" column should be color-coded
    When user closes the current view
    And user switches to the "demog-1000" table view
    Then the open table views should be exactly "demog-1000"
    When user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    And user selects "background" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    And user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: A saved layout brings the text coloring back after it was switched off
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    And user selects "text" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then the text of "AGE" column should be color-coded
    When user saves the layout of the current table view as "bdd-cc-text-{run}"
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    And user selects "background" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    And user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should be the same
    When user applies the layout "bdd-cc-text-{run}" to the current table view
    Then "AGE" column should be color-coded linearly
    And the text of "AGE" column should be color-coded
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    And user selects "background" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    And user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Grid Color Coding > All colours a boolean column too
    Then "CONTROL" column should have no color coding
    And the "color of cell 1 of CONTROL" and "color of cell 3 of CONTROL" readings of grid should be the same
    When user picks "Grid Color Coding > All" from the context menu of the "cell 1 of AGE" area of grid
    Then "Color Coding" property of grid should be "All"
    And the "color of cell 1 of CONTROL" and "color of cell 3 of CONTROL" readings of grid should differ
    When user picks "Grid Color Coding > Auto" from the context menu of the "cell 1 of AGE" area of grid
    Then "Color Coding" property of grid should be "Auto"
    And the "color of cell 1 of CONTROL" and "color of cell 3 of CONTROL" readings of grid should be the same
    And no errors should have been logged
