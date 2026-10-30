Feature: The Select menu selects rows and columns of the current table, and Selection to Column keeps the selection
  The table view's Select menu: Duplicates... selects the rows whose key columns repeat another row's
  (help `transform/select-duplicates.md`; Superfluous rows only leaves one row of each group unselected),
  Missing Values... the rows empty in the checked columns (`select-missing-values.md`), Random Rows... the
  number of rows typed or picked by a proportion link (`select-random-rows.md`), All every row and column,
  All Columns every column, and Selection to Column... adds a boolean column true in the selected rows,
  named "Selection" by default (`selection-to-column.md`).

  Every scenario runs on a five-row table written here, so each claim names the rows selected and the
  ones left alone. Left out on purpose: Random Rows claims the count, never which rows; None, Invert,
  Extract Selected Rows, Selection to Filter and Filter to Column have features already
  (`viewers/nx/nx-regression-line`, `viewers/filter-panel/compose-with-viewers`,
  `viewers/grid/grid-extracted-rows-project`, `viewers/filter-panel/filter-panel-entry-points`);
  Selection to Top sorts the grid, and no step reads the grid's row order (MISSING.md). The predicates
  behind the dialogs (duplicates over several keys, nulls of every type) are ddt tests. Nothing is put on
  the server.

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

  Scenario: Superfluous rows only leaves one row of each duplicate group unselected
    When user picks "Select > Duplicates..." from the top menu
    And user toggles the "name" column in the column list of Select Duplicates dialog
    And user checks "Superfluous rows only" input in Select Duplicates dialog
    And user clicks on OK button in Select Duplicates dialog
    Then only rows where "id" is one of "2, 4" should be selected
    And no errors should have been logged

  Scenario: Missing Values selects the rows empty in the checked column
    When user picks "Select > Missing Values..." from the top menu
    Then Select Missing Values dialog should be visible
    When user toggles the "v" column in the column list of Select Missing Values dialog
    And user clicks on OK button in Select Missing Values dialog
    Then Select Missing Values dialog should be hidden
    And only rows where "id" is one of "2, 4" should be selected
    And no errors should have been logged

  Scenario: Random Rows selects as many rows as typed, and a proportion link sets the count
    When user picks "Select > Random Rows..." from the top menu
    Then Select Random Rows dialog should be visible
    When user enters "3" into "# rows" input in Select Random Rows dialog
    And user clicks on OK button in Select Random Rows dialog
    Then Select Random Rows dialog should be hidden
    And 3 rows should be selected
    When user picks "Select > Random Rows..." from the top menu
    And user clicks on "100%" label in Select Random Rows dialog
    Then "# rows" input in Select Random Rows dialog should have value "5"
    When user clicks on OK button in Select Random Rows dialog
    Then all rows should be selected
    And no errors should have been logged

  Scenario: All selects every row and every column, All Columns only the columns
    When user picks "Select > All Columns" from the top menu
    Then columns "id, name, v" should be selected
    And no rows should be selected
    When user picks "Select > All" from the top menu
    Then all rows should be selected
    And columns "id, name, v" should be selected
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
