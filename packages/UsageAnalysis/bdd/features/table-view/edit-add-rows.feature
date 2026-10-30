Feature: Edit > Add Rows inserts empty rows at a position, and Ctrl+Z takes them out again
  Edit > Add Rows... opens the Add Rows dialog with the number of rows to add (Count, 1 by default) and
  their Position, which starts at the current row (help `visualize/viewers/grid.md` "Add rows and
  columns"). OK inserts that many empty rows before the row at Position (counted from 0, the row count
  appends them at the end); the existing rows keep their values and order. Edit | Add Rows is undoable
  (help `datagrok/navigation/undo.md`).

  The table is written here, so every claim names the rows that moved and the empty ones. Left out on
  purpose: the grid's own + icon (another gesture, `viewers/grid`), and large counts (the `insertRows`
  function). Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "items" with:
      | name  | qty |
      | one   | 1   |
      | two   | 2   |
      | three | 3   |

  Scenario: Count rows go in before the row at Position, the existing rows keep their values, and Ctrl+Z removes them
    When user picks "Edit > Add Rows..." from the top menu
    Then Add Rows dialog should be visible
    And Count input in Add Rows dialog should have value "1"
    When user enters "1" into Position input in Add Rows dialog
    And user enters "2" into Count input in Add Rows dialog
    And user clicks on OK button in Add Rows dialog
    Then Add Rows dialog should be hidden
    And the table should have 5 rows
    And the value of "name" column in row 1 should be "one"
    And the value of "name" column in row 2 should be ""
    And the value of "qty" column in row 2 should be ""
    And the value of "name" column in row 3 should be ""
    And the value of "name" column in row 4 should be "two"
    And the value of "qty" column in row 4 should be "2"
    And the value of "name" column in row 5 should be "three"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have 3 rows
    And the value of "name" column in row 2 should be "two"
    And the value of "name" column in row 3 should be "three"
    And no errors should have been logged

  Scenario: Position starts at the current row
    Given user makes row 3 current
    When user picks "Edit > Add Rows..." from the top menu
    Then Position input in Add Rows dialog should have value "2"
    When user clicks on OK button in Add Rows dialog
    Then the table should have 4 rows
    And the value of "name" column in row 3 should be ""
    And the value of "name" column in row 4 should be "three"
    And the value of "name" column in row 2 should be "two"
    And no errors should have been logged

  Scenario: Position equal to the row count appends the rows at the end
    When user picks "Edit > Add Rows..." from the top menu
    And user enters "3" into Position input in Add Rows dialog
    And user clicks on OK button in Add Rows dialog
    Then the table should have 4 rows
    And the value of "name" column in row 3 should be "three"
    And the value of "name" column in row 4 should be ""
    And no errors should have been logged
