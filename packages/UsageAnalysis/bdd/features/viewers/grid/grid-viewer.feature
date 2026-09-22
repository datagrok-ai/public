@journey @viewers @realizes:viewers.grid
Feature: A second grid as a viewer, and column tooltips
  A grid added to a view is a viewer like any other: it shows the same table, takes its own Row
  Height and Show Column Labels without touching the view's own grid, rebinds to another table, and
  closes leaving the view's grid alone. The second one is addressed as "second grid viewer" - the
  reserved "grid" always means the table view's own. The column tooltip settings live on the header
  menu as radio items, and the item that is on says so through `aria-checked`; Columns asks for the
  columns in a dialog and the row tooltip then lists them. Pick Up / Apply carries
  the main grid's look to the second one: its grid-wide Color Coding All, read back as the second
  grid's property and as the colour it resolves for a HEIGHT cell. One journey on demog-1000, with
  its first 100 rows as the second table (`grid.md` uses spgi-100; operator decision D8) - a table
  with the same columns, so the rebind is shown by the row count alone.

  The tooltip scenario walks the path of GROK-20890: with the current column's tooltip set to
  Columns and a non-empty set of columns chosen in the dialog (the dialog opens with none ticked,
  and with none the defect stayed away), a pointer resting on the header of that same column threw
  `Invalid argument (index): null` — `grid_tooltip.dart` took the Columns branch and built the row
  tooltip from a table row a header has none of. Fixed in the core on 2026-09-15 (the branch is
  gated on a data cell); the scenario carried `@known-failure` until the fix reached the stand, and
  its zero-error floor — Columns, All, OK, then the AGE header — now guards the fix.

  Not translated, and why: from `grid.md` "Column Tooltip Settings", that None shows no tooltip -
  the tooltip appears on a debounce after the pointer rests, and no signal says a hover has been
  processed, so an absence right after the hover holds whatever the setting (the menu mark of None
  is claimed instead). From `grid-ui.md` "Pick Up / Apply", the formatting and the per-column style,
  and a changed colour scheme on Average Mass - the colouring of a column is the table's, so both grids show it whatever Pick Up does;
  the grid-wide Color Coding is what the second grid only gets from Apply.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And the open tableview should have 1 grid viewer

  Scenario: A second grid shows the same table and takes its own row height
    When user adds a grid viewer
    Then grid viewer should be added to the open tableview
    And the open tableview should have 2 grid viewers
    And second grid viewer should show 1000 rows
    When user sets "Row Height" property of second grid viewer to "40"
    Then the "cell 1 of USUBJID" area of second grid viewer should be at least 34 pixels tall
    And "Row Height" property of grid should not be "40"
    When user sets "Show Column Labels" property of second grid viewer to "false"
    Then second grid viewer should not have a "header AGE" area
    And grid should have a "header AGE" area
    When user sets "Show Column Labels" property of second grid viewer to "true"
    Then second grid viewer should have a "header AGE" area
    And no errors should have been logged

  Scenario: The second grid rebinds to another table and closes without taking the view's grid
    When user opens demog-1000 dataset keeping the first 100 rows as "demog-100"
    Then the table should have 100 rows
    When user switches to the "demog-1000" table view
    And user sets "Table" property of second grid viewer to "demog-100"
    Then second grid viewer should be bound to table "demog-100"
    And second grid viewer should show 100 rows
    And grid should show 1000 rows
    When user clicks on close icon of second grid viewer
    Then the open tableview should have 1 grid viewer
    And grid should show 1000 rows
    And grid should have a "header AGE" area
    When user closes all views
    Then no errors should have been logged

  Scenario: The column tooltip menu marks the setting that is on, and Columns lists the chosen columns (GROK-20890)
    Given user opens demog-1000 dataset
    When user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "Default" menu item in context menu should be visible
    And "Form" menu item in context menu should be visible
    And "Columns" menu item in context menu should be visible
    And "None" menu item in context menu should be visible
    And "Default" menu item in context menu should be selected
    And "None" menu item in context menu should not be selected
    When user closes the context menu
    And user picks "Tooltip > Current Column > Columns" from the context menu of the "header AGE" area of grid
    Then "Select columns..." dialog should be visible
    When user clicks on All link in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be hidden
    When user hovers over the "cell 3 of AGE" area of grid
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the tooltip should show "AGE" as "58"
    When user moves the pointer away from grid
    And user hovers over the "header AGE" area of grid
    And user moves the pointer away from grid
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "Columns" menu item in context menu should be selected
    When user closes the context menu
    And user picks "Tooltip > Current Column > None" from the context menu of the "header AGE" area of grid
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "None" menu item in context menu should be selected
    And "Default" menu item in context menu should not be selected
    When user closes the context menu
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    And user clicks on "Default" menu item in context menu
    And user right-clicks on the "header AGE" area of grid
    And user hovers over "Tooltip" menu item in context menu
    And user hovers over "Current Column" menu item in context menu
    Then "Default" menu item in context menu should be selected
    And "None" menu item in context menu should not be selected
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Pick Up and Apply carry the grid's look to a second grid
    When user adds a grid viewer
    Then the open tableview should have 2 grid viewers
    When user picks "Grid Color Coding > All" from the context menu of the "cell 3 of AGE" area of grid
    And user moves the pointer away from grid
    Then "Color Coding" property of grid should be "All"
    And "Color Coding" property of second grid viewer should be "Auto"
    When user remembers the "color of cell 2 of HEIGHT" reading of grid
    Then the "color of cell 2 of HEIGHT" reading of grid should not be "#ffffff"
    And the "color of cell 2 of HEIGHT" reading of second grid viewer should not be as remembered
    When user picks "Pick Up / Apply > Pick Up" from the context menu of the "cell 3 of AGE" area of grid
    And user picks "Pick Up / Apply > Apply" from the context menu of the "cell 3 of AGE" area of second grid viewer
    Then "Color Coding" property of second grid viewer should be "All"
    And the "color of cell 2 of HEIGHT" reading of second grid viewer should be as remembered
    When user picks "Grid Color Coding > Auto" from the context menu of the "cell 3 of AGE" area of grid
    And user clicks on close icon of second grid viewer
    Then the open tableview should have 1 grid viewer
    And no errors should have been logged
