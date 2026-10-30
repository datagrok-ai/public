@journey
Feature: Ctrl+Z undoes and Ctrl+Shift+Z or Ctrl+Y redoes the table view's column and layout operations
  Follows help/datagrok/navigation/undo.md and the shortcuts page: Ctrl+Z reverses the last operation on
  the current table or view, Ctrl+Shift+Z (or Ctrl+Y) re-applies it, the same commands are in the Edit
  menu, undo is scoped to the view being looked at, and while a dialog is open Ctrl+Z never reaches the
  table. One scenario per operation the help lists as undoable that no feature claimed yet: a column
  removed from its header menu, a column renamed in Column Properties, color coding set from the header
  menu, and View > Layout > Clear. A formula column added through Edit > Add New Column is not undone
  today (MISSING.md, Findings 30 Oct 2026), so it is not claimed here. Each scenario undoes,
  redoes, and undoes again, so the next one starts from the table as written here.

  Left out on purpose: a grid cell edit, Shift+Delete rows (viewers/grid/grid-editing.feature) and Find and
  Replace (general/find-and-replace.feature) are claimed there; a viewer closed from its title bar coming
  back (pc-plot-hover-and-layout, trellis-plot-legend-and-lifecycle, with redo) likewise; changing a
  column's type is change-type.feature, tags and format edited in Column Properties are
  column-properties.feature. The record list itself (names, sizes, what is released with a closed view)
  is the Dart undo_test.dart's job. Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "undo" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cy   | Oslo | 3      |

  Scenario: A column removed from its header menu comes back in its place, and redo removes it again
    When user picks "Remove" from the context menu of the "header city" area of grid
    Then the table should not have a column "city"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And the "column order" reading of grid should be "name, city, weight"
    And the value of "city" column in row 2 should be "Rome"
    When user presses Control+Shift+Z
    Then the table should not have a column "city"
    And the table should have 2 columns
    When user presses Control+Z
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

  Scenario: View > Layout > Clear removes the viewers, Ctrl+Z brings them back and redo clears again
    Given user adds a histogram viewer
    And user adds a scatter plot viewer
    When user picks "View > Layout > Clear" from the top menu
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 scatter plot viewer
    When user presses Control+Shift+Z
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And no errors should have been logged

  Scenario: Edit > Undo and Edit > Redo do what the shortcuts do
    When user picks "Remove" from the context menu of the "header city" area of grid
    Then the table should not have a column "city"
    When user picks "Edit > Undo" from the top menu
    Then the table should have a column "city"
    When user picks "Edit > Redo" from the top menu
    Then the table should not have a column "city"
    When user picks "Edit > Undo" from the top menu
    Then the "column order" reading of grid should be "name, city, weight"
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
    And user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have a column "city"
    And table "other" should have columns "id, score"
    When user switches to the "other" table view
    And user closes the current view
    Then no errors should have been logged
