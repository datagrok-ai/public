@journey
Feature: A column edited from the table view — removed, renamed, color-coded, its properties and type changed — and Ctrl+Z puts it back
  Follows help/datagrok/navigation/undo.md and the shortcuts page: Ctrl+Z reverses the last operation on the
  current table or view, Ctrl+Shift+Z (or Ctrl+Y) re-applies it, the same commands are in the Edit menu,
  undo is scoped to the view being looked at, and while a dialog is open Ctrl+Z never reaches the table.
  Claimed on a column removed from its header menu, a column renamed in Column Properties, color coding set
  from the header menu, Ctrl+Shift+R (View > Layout > Clear), a format tag, description and units set in
  Column Properties (the header tooltip shows them) and Change Type... (help/transform/change-column-type.md),
  opened from the header menu and from the Change type... link of Column Properties. F2 on a cell opens its
  column's Column Properties, not the cell editor (public#629), and Column Properties refuses a name another
  column has. Every scenario puts the column back, so the next one starts from the
  table as written here; the Change Type scenarios run on a second table, opened by the first of them.

  Left out on purpose: a grid cell edit, Shift+Delete rows (viewers/grid/grid-editing.feature) and Find and
  Replace (general/find-and-replace.feature); a viewer closed from its title bar coming back
  (pc-plot-hover-and-layout, trellis-plot-legend-and-lifecycle); Edit category colors...
  (viewers/color-coding/color-coding-editor-conditional-and-categorical.feature); the rename itself (PowerPack
  add-new-column/persistence-sources); a formula column added through Add New Column, which is not undone
  today (public#4461). The record list itself is the Dart undo_test.dart's job; the type and format matrix is
  the ddt tests of ChangeColumnsType and of the formats. Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "undo" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cy   | Oslo | 3      |

  Scenario: A column removed from its header menu comes back in its place with Ctrl+Z, and Edit > Redo and Edit > Undo do the same
    When user picks "Remove" from the context menu of the "header city" area of grid
    Then the table should not have a column "city"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And the "column order" reading of grid should be "name, city, weight"
    And the value of "city" column in row 2 should be "Rome"
    When user picks "Edit > Redo" from the top menu
    Then the table should not have a column "city"
    And the table should have 2 columns
    When user picks "Edit > Undo" from the top menu
    Then the "column order" reading of grid should be "name, city, weight"
    And no errors should have been logged

  Scenario: A column renamed in Column Properties gets its name back, and Ctrl+Y renames it again
    When user picks "Column Properties..." from the context menu of the "header city" area of grid
    And user enters "town" into "New name:" input in "city" dialog
    And user clicks on OK button in "city" dialog
    Then the table should have a column "town"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And the table should not have a column "town"
    And the value of "city" column in row 2 should be "Rome"
    When user presses Control+Y
    Then the table should have a column "town"
    And the table should not have a column "city"
    When user presses Control+Z
    Then the table should have a column "city"
    And no errors should have been logged

  Scenario: Color coding set from the header menu is taken off by Ctrl+Z and put back by redo
    When user picks "Color Coding > Linear" from the context menu of the "header weight" area of grid
    Then "weight" column should be color-coded linearly
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "weight" column should have no color coding
    When user presses Control+Shift+Z
    Then "weight" column should be color-coded linearly
    When user presses Control+Z
    Then "weight" column should have no color coding
    And no errors should have been logged

  Scenario: Ctrl+Shift+R clears the view down to the grid, Ctrl+Z brings the viewers back and redo clears again
    Given user adds a histogram viewer
    And user adds a scatter plot viewer
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Shift+R
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And grid should be visible
    When user presses Control+Z
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 scatter plot viewer
    When user presses Control+Shift+Z
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And no errors should have been logged

  # public#629
  @realizes:public#629
  Scenario: F2 opens the current column's Column Properties, not the cell editor, and a format tag set there changes the grid until Ctrl+Z (public#629)
    Then the "weight" cell of row 1 should be displayed as "1.50"
    When user clicks on the "cell 1 of weight" area of grid
    And user presses F2
    Then "weight" dialog should be visible
    And cell editor should be hidden
    And "New name:" input in "weight" dialog should have value "weight"
    And "Key 1" text input in "weight" dialog should have value "format"
    And "Value 1" text input in "weight" dialog should have value "0.00"
    When user enters "0.000" into "Value 1" text input in "weight" dialog
    And user clicks on OK button in "weight" dialog
    Then the "weight" dialog should close
    And the "weight" cell of row 1 should be displayed as "1.500"
    And the "weight" cell of row 2 should be displayed as "2.250"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the "weight" cell of row 1 should be displayed as "1.50"
    And no errors should have been logged

  Scenario: Description and units added as two tag rows show in the header tooltip
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    And user clicks on first "Add row" button in "weight" dialog
    And user enters "description" into "Key 1" text input in "weight" dialog
    And user enters "Body weight" into "Value 1" text input in "weight" dialog
    And user clicks on first "Add row" button in "weight" dialog
    Then "Key 2" text input in "weight" dialog should have value "description"
    When user enters "units" into "Key 1" text input in "weight" dialog
    And user enters "kg" into "Value 1" text input in "weight" dialog
    And user clicks on OK button in "weight" dialog
    Then the "weight" dialog should close
    And "weight" column should have units "kg"
    When user hovers over the "header weight" area of grid
    Then tooltip should contain text "Body weight"
    And tooltip should contain text "units: kg"
    When user moves the pointer away from grid
    And user hovers over the "header city" area of grid
    Then tooltip should contain text "city"
    And tooltip should not contain text "Body weight"
    When user moves the pointer away from grid
    And user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "weight" column should have units ""
    And no errors should have been logged

  Scenario: Column Properties refuses a name another column already has
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    And user enters "City" into "New name:" input in "weight" dialog
    Then "New name:" input in "weight" dialog should be invalid
    When user hovers over "New name:" input in "weight" dialog
    Then tooltip should contain text "Column 'City' already exists."
    When user enters "mass" into "New name:" input in "weight" dialog
    Then "New name:" input in "weight" dialog should be valid
    When user presses Escape
    Then "weight" dialog should be hidden
    And the table should have the columns "name, city, weight"
    And no errors should have been logged

  Scenario: While a dialog is open Ctrl+Z leaves the table alone, and after it closes Ctrl+Z works
    When user picks "Remove" from the context menu of the "header city" area of grid
    Then the table should not have a column "city"
    When user clicks on the "cell 1 of weight" area of grid
    And user presses F2
    Then "weight" dialog should be visible
    When user presses Control+Z
    Then the table should not have a column "city"
    And "weight" dialog should be visible
    When user presses Escape
    Then "weight" dialog should be hidden
    And the table should not have a column "city"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And no errors should have been logged

  Scenario: Ctrl+Z in another table's view does not undo this table's change, and switching back brings the history
    When user picks "Remove" from the context menu of the "header city" area of grid
    Then the table should not have a column "city"
    Given user opens a table "other" with:
      | id | score |
      | 1  | 10    |
      | 2  | 20    |
    When user clicks on the "cell 1 of id" area of grid
    And user presses Control+Z
    Then table "undo" should have columns "name, weight"
    And table "other" should have columns "id, score"
    When user switches to the "undo" table view
    Then the table should not have a column "city"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And table "other" should have columns "id, score"
    When user switches to the "other" table view
    And user closes the current view
    Then no errors should have been logged

  Scenario: Change Type... in the header menu shows the current type, and OK converts the column, a value that cannot be converted taking the Not convertible value
    Given user opens a table "ct" with:
      | name | code |
      | a    | 10   |
      | b    | 20   |
      | c    | n/a  |
      | d    |      |
    Then "code" column should have type "string"
    When user picks "Change Type..." from the context menu of the "header code" area of grid
    Then "Change column type" dialog should be visible
    And "Type" input in "Change column type" dialog should have value "string"
    And "Type" input in "Change column type" dialog should be disabled
    When user selects "int" in "New type" input in "Change column type" dialog
    And user enters "-1" into "Not convertible" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "int"
    And the value of "code" column in row 2 should be "20"
    And the value of "code" column in row 3 should be "-1"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "code" column should have type "string"
    And the value of "code" column in row 3 should be "n/a"
    And no errors should have been logged

  Scenario: The filter icon beside Not convertible keeps only the rows that cannot be converted
    When user picks "Change Type..." from the context menu of the "header code" area of grid
    And user selects "int" in "New type" input in "Change column type" dialog
    And user clicks on second "Filter" icon in "Change column type" dialog
    Then 1 row should pass the filter
    And all rows where "code" is "n/a" should pass the filter
    When user clicks on CANCEL button in "Change column type" dialog
    Then the "Change column type" dialog should close
    When user resets the filter
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Change type... in Column Properties opens the same dialog, and its OK converts the column
    When user picks "Column Properties..." from the context menu of the "header code" area of grid
    Then "code" dialog should be visible
    When user clicks on "Change type..." label in "code" dialog
    Then "Change column type" dialog should be visible
    When user selects "double" in "New type" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "double"
    And the value of "code" column in row 1 should be "10"
    When user presses Escape
    Then "code" dialog should be hidden
    And no errors should have been logged
