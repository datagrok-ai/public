@journey
Feature: The Select menu selects rows of the current table, Edit > Add Rows inserts empty rows, and Selection to Column keeps the selection
  The table view's Select menu: Duplicates... selects the rows whose key columns repeat another row's
  (help `transform/select-duplicates.md`), Missing Values... the rows empty in the checked columns
  (`select-missing-values.md`), Random Rows... the number of rows typed (`select-random-rows.md`), and
  Selection to Column... adds a boolean column true in the selected rows, named "Selection" by default
  (`selection-to-column.md`). Edit > Add Rows... (help `visualize/viewers/grid.md` "Add rows and columns")
  inserts Count empty rows before the row at Position, which starts at the current row (counted from 0), and
  is undoable (`datagrok/navigation/undo.md`).

  One journey on a five-row table written here, so each claim names the rows selected and the ones left
  alone. Duplicates and Missing Values add to the selection, so a scenario that claims an exact selection
  clears it first. Left out on purpose: Random Rows claims the count, never which rows; Select > All is
  Ctrl+A (`viewers/grid/grid-rows.feature`); None, Invert, Extract Selected Rows, Selection to Filter and
  Filter to Column have features already (`viewers/nx/nx-regression-line`, `viewers/filter-panel/compose-with-viewers`,
  `viewers/grid/grid-extracted-rows-project`, `viewers/filter-panel/filter-panel-entry-points`); Superfluous
  rows only is the `selectDuplicates` ddt test; large counts and the append position are the `insertAt`
  ApiTests test. Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | id | name | v  |
      | 1  | a    | 10 |
      | 2  | a    |    |
      | 3  | b    | 30 |
      | 4  | b    |    |
      | 5  | c    | 50 |

  Scenario: Duplicates selects every row whose key repeats another row's
    When user picks "Select > Duplicates..." from the top menu
    Then Select Duplicates dialog should be visible
    When user toggles the "name" column in the column list of Select Duplicates dialog
    And user clicks on OK button in Select Duplicates dialog
    Then Select Duplicates dialog should be hidden
    And only rows where "name" is one of "a, b" should be selected
    And 4 rows should be selected
    And no errors should have been logged

  Scenario: Missing Values selects the rows empty in the checked column
    Given user clears the row selection
    When user picks "Select > Missing Values..." from the top menu
    Then Select Missing Values dialog should be visible
    When user toggles the "v" column in the column list of Select Missing Values dialog
    And user clicks on OK button in Select Missing Values dialog
    Then Select Missing Values dialog should be hidden
    And only rows where "id" is one of "2, 4" should be selected
    And no errors should have been logged

  Scenario: Random Rows selects as many rows as typed
    When user picks "Select > Random Rows..." from the top menu
    Then Select Random Rows dialog should be visible
    When user enters "3" into "# rows" input in Select Random Rows dialog
    And user clicks on OK button in Select Random Rows dialog
    Then Select Random Rows dialog should be hidden
    And 3 rows should be selected
    And no errors should have been logged

  Scenario: Add Rows inserts Count empty rows before the row at Position, which starts at the current row, and Ctrl+Z takes them out
    Given user makes row 3 current
    When user picks "Edit > Add Rows..." from the top menu
    Then Add Rows dialog should be visible
    And Count input in Add Rows dialog should have value "1"
    And Position input in Add Rows dialog should have value "2"
    When user enters "1" into Position input in Add Rows dialog
    And user enters "2" into Count input in Add Rows dialog
    And user clicks on OK button in Add Rows dialog
    Then Add Rows dialog should be hidden
    And the table should have 7 rows
    And the value of "name" column in row 1 should be "a"
    And the value of "name" column in row 2 should be ""
    And the value of "id" column in row 3 should be ""
    And the value of "id" column in row 4 should be "2"
    And the value of "id" column in row 7 should be "5"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should have 5 rows
    And the value of "id" column in row 2 should be "2"
    And the value of "id" column in row 3 should be "3"
    And no errors should have been logged

  Scenario: Selection to Column adds a boolean column true in the selected rows
    Given user selects rows where "name" is "b"
    When user picks "Select > Selection to Column..." from the top menu
    Then Selection to Column dialog should be visible
    And Name input in Selection to Column dialog should have value "Selection"
    When user clicks on OK button in Selection to Column dialog
    Then Selection to Column dialog should be hidden
    And a new column "Selection" should have been added
    And "Selection" column should have type "bool"
    And the value of "Selection" column in row 3 should be "true"
    And the value of "Selection" column in row 4 should be "true"
    And the value of "Selection" column in row 1 should be "false"
    And the value of "Selection" column in row 5 should be "false"
    And no errors should have been logged
