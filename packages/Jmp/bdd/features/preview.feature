@realizes:jmp.preview
Feature: A JMP file in the Browse tree previews before it is opened
  A .jmp file in the Browse tree, clicked or double-clicked, shows the package's u2 preview: a
  summary line (the JMP version that saved the file, its size in rows and columns, when it was
  saved), the table in a grid, and Open in the ribbon, which opens the table in the workspace.
  Columns that don't come across (a layout the reader doesn't know, or row states, which aren't
  imported yet) are listed in a folded "Columns not imported" section, so nothing disappears
  silently. The files are the package's own, published with it under App Data/Jmp:
  jmpreader/example1.jmp (JMP 17.2.0, 4 rows, 12 columns) and tests/unsupported-column.jmp, the
  same table with the formula column's layout byte changed to one no JMP version writes. Nothing is
  written to the server.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Jmp tree node inside browse tree is expanded

  Scenario: The preview tells what the file is, and Open opens the table
    Given no table view is open
    And Files---App-Data---Jmp---jmpreader tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Jmp---jmpreader---example1.jmp tree node inside browse tree
    Then the "example1.jmp" view should be current
    And the current view should show the text "JMP 17.2.0"
    And the current view should show the text "4 rows × 12 columns"
    And "Columns not imported" section should be absent
    When user clicks on Open button
    Then the "example1" table view should open with 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The preview lists the columns it can't import
    Given Files---App-Data---Jmp---tests tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Jmp---tests---unsupported-column.jmp tree node inside browse tree
    Then the "unsupported-column.jmp" view should be current
    And the current view should show the text "4 rows × 11 columns"
    And "Columns not imported" section should be visible
    And "Columns not imported" section should contain text "1 of 12"
    And "Columns not imported" section should be above grid
    When user expands "Columns not imported" section
    Then "Columns not imported" section should contain text "formula"
    And "Columns not imported" section should contain text "Unsupported column layout"
    And "Columns not imported" section should be above grid
    And no errors should have been logged
