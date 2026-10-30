Feature: Data > Join Tables joins two open tables on their key columns into a new table view
  With two tables open, Data > Join Tables... opens the Join Tables dialog on them, with the column
  both tables share as the key pair (help `transform/join-tables.md`). The dialog reports per table how
  many key rows have no match ("1 / 3"), and the menu next to the number selects those rows; when no
  key row matches at all it says so and OK is disabled. The Preview result columns grid shows the join's
  rows as the settings change. OK opens the join as a new "result" table view;
  Values decides which columns of each table go into it, Join Type is passed on to the join, and a left
  join done In-place writes the right table's columns into the left table (renamed "left+right") instead
  of opening a new view.

  Every scenario runs on two small tables written here (left: ids 1-3, right: ids 2-4), so the counts
  are exact. Left out on purpose:
  - the join type matrix (inner/outer/left/right, one-to-many keys, composite keys) is the ddt unit tests
    of `joinTables`; one scenario shows that the dialog passes the type on;
  - OK "while keys are missing": the first key pair cannot be removed, and with no pair the dialog joins
    by row index (MISSING.md, question);
  - a join inside a saved project (`projects/projects-regressions-join.feature`) and the Aggregation
    tutorial (`Tutorials/bdd/features/transform/data-aggregation.feature`) already have features.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "left" with:
      | id | name |
      | 1  | a    |
      | 2  | b    |
      | 3  | c    |
    And user opens a table "right" with:
      | id | score |
      | 2  | 20    |
      | 3  | 30    |
      | 4  | 40    |

  Scenario: The dialog opens on both tables keyed by their shared column, and OK opens the inner join as a new table view
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should be visible
    And Key Columns input in Join Tables dialog should contain text "id"
    And "Join Type" input in Join Tables dialog should have value "inner"
    And grid viewer in "Preview result columns" dialog should show 2 rows
    And OK button in Join Tables dialog should be enabled
    When user clicks on OK button in Join Tables dialog
    Then Join Tables dialog should be hidden
    And the open table views should be exactly "left, right, result"
    Given user switches to the "result" table view
    Then the table should have 2 rows
    And the table should have the columns "id, name, right.id, score"
    And the value of "name" column in row 1 should be "b"
    And the value of "score" column in row 1 should be "20"
    And the value of "name" column in row 2 should be "c"
    And the value of "score" column in row 2 should be "30"
    And table "left" should have 3 rows
    And table "right" should have 3 rows
    And no errors should have been logged

  Scenario: The dialog counts the key rows without a match, and its menu selects them in that table
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should contain text "1 / 3"
    When user clicks on first "context arrow down" icon in Join Tables dialog
    And user picks "Select All Not Matching" from the open menu
    Then 1 row of table "left" should be selected
    And 0 rows of table "right" should be selected
    When user clicks on CANCEL button in Join Tables dialog
    Then the open table views should be exactly "left, right"
    And no errors should have been logged

  Scenario: Join Type is passed on to the join: outer keeps every row of both tables
    When user picks "Data > Join Tables..." from the top menu
    And user selects "outer" in "Join Type" input in Join Tables dialog
    Then grid viewer in "Preview result columns" dialog should show 4 rows
    When user clicks on OK button in Join Tables dialog
    Then the open table views should be exactly "left, right, result"
    Given user switches to the "result" table view
    Then the table should have 4 rows
    And the table should have the columns "id, name, right.id, score"
    And no errors should have been logged

  Scenario: Values decides which columns of the table go into the result
    When user picks "Data > Join Tables..." from the top menu
    And user clicks on Values input in Join Tables dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user clicks on OK button in Join Tables dialog
    Then the open table views should be exactly "left, right, result"
    Given user switches to the "result" table view
    Then the table should have 2 rows
    And the table should have the columns "id, right.id, score"
    And no errors should have been logged

  Scenario: No key row matching says so and disables OK
    Given user sets "id" column in row 1 to "7"
    And user sets "id" column in row 2 to "8"
    And user sets "id" column in row 3 to "9"
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should contain text "None of key rows are matched"
    And OK button in Join Tables dialog should be disabled
    When user clicks on CANCEL button in Join Tables dialog
    Then the open table views should be exactly "left, right"
    And no errors should have been logged

  Scenario: In-place is offered only for left and right, and a left join in place adds the columns to the left table
    When user picks "Data > Join Tables..." from the top menu
    Then "In-place" input in Join Tables dialog should be disabled
    When user selects "left" in "Join Type" input in Join Tables dialog
    Then "In-place" input in Join Tables dialog should be enabled
    When user checks "In-place" input in Join Tables dialog
    And user clicks on OK button in Join Tables dialog
    Then Join Tables dialog should be hidden
    And the open table views should be exactly "left+right, right"
    And table "left+right" should have 3 rows
    And table "left+right" should have columns "id, name, right.id, score"
    And the "score" cell of row 2 of table "left+right" should be displayed as "20"
    And table "left+right" should have missing values in "score" column
    And no errors should have been logged
