Feature: Data > Compare Tables opens the rows two tables differ in as a new color-coded "diff" table view
  Data > Compare Tables... needs two open tables (help `transform/compare-tables.md`). The Compare Tables
  dialog opens on them keyed by the first column, with the columns both tables share (same name and type) as
  Value Columns. It validates the keys (a key whose values repeat disables OK and says so). OK opens a "diff"
  table view with one row per key whose values changed, or which is in one table only, plus the "~t1 row" and
  "~t2 row" columns pointing back at the source rows; a changed value is painted orange, and a row in one table
  only in its own color. Clicking a diff row makes the matching rows current in both tables. Value Columns
  decides which columns are compared.

  Every scenario runs on two three-row tables written here (Ann is equal, Bob's city changed, Cid is only in
  the first, Dan only in the second), so the counts are exact. Left out on purpose:
  - the comparison matrix (composite keys, every column type, duplicates, columns added or removed) is the
    `TableComparison` code's job (ddt tests);
  - which color stands for "added" and which for "removed": the help does not say, and the code's comments
    contradict its colors (MISSING.md, finding, core#1231) — the feature claims only that the two differ;
  - the History button, whose states belong to every dialog with history.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | name | city | weight |
      | Ann  | Oslo | 10.5   |
      | Bob  | Rome | 10.5   |
      | Cid  | Oslo | 10.5   |

  Scenario: Compare Tables needs a second open table
    Then the top menu should list "Data > Compare Tables..." as disabled
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    Then the top menu should list "Data > Compare Tables..." as enabled
    And no errors should have been logged

  Scenario: The dialog opens keyed by name, and OK opens the rows that differ as a new table view
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    When user picks "Data > Compare Tables..." from the top menu
    Then Compare Tables dialog should be visible
    And Key Columns input in Compare Tables dialog should contain text "name"
    And OK button in Compare Tables dialog should be enabled
    When user clicks on OK button in Compare Tables dialog
    Then Compare Tables dialog should be hidden
    And the open table views should be exactly "people, people2, diff"
    Given user switches to the "diff" table view
    Then the table should have 3 rows
    And the table should have the columns "name, city, weight, ~t1 row, ~t2 row"
    And the value of "name" column in row 1 should be "Bob"
    And the value of "city" column in row 1 should be "Paris"
    And the value of "name" column in row 2 should be "Cid"
    And the value of "name" column in row 3 should be "Dan"
    And the value of "weight" column in row 3 should be "4"
    And table "people" should have 3 rows
    And table "people2" should have 3 rows
    And no errors should have been logged

  Scenario: The diff grid paints a changed value orange and the one-table rows in two other colors, and a click on a row makes it current in both tables
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    When user picks "Data > Compare Tables..." from the top menu
    And user clicks on OK button in Compare Tables dialog
    Then the open table views should be exactly "people, people2, diff"
    Given user switches to the "diff" table view
    Then the "cell 1 of city" area of grid should contain the color "#FFA500"
    And the "cell 1 of name" area of grid should not contain the color "#FFA500"
    And the "cell 2 of name" and "cell 3 of name" areas of grid should be painted in different colors
    When user clicks on the "cell 1 of name" area of grid
    Then "name" of the current row of table "people" should be "Bob"
    And "name" of the current row of table "people2" should be "Bob"
    And no errors should have been logged

  Scenario: A key whose values repeat disables OK and says why
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    When user picks "Data > Compare Tables..." from the top menu
    And user selects "city" in Key Columns input in Compare Tables dialog
    Then Compare Tables dialog should contain text "city is not unique"
    And OK button in Compare Tables dialog should be disabled
    When user clicks on CANCEL button in Compare Tables dialog
    Then the open table views should be exactly "people, people2"
    And no errors should have been logged

  Scenario: With no value columns only the rows in one table go into the diff
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    When user picks "Data > Compare Tables..." from the top menu
    And user clicks on "Clear" icon in Compare Tables dialog
    And user clicks on OK button in Compare Tables dialog
    Then the open table views should be exactly "people, people2, diff"
    Given user switches to the "diff" table view
    Then the table should have 2 rows
    And the value of "name" column in row 1 should be "Cid"
    And the value of "name" column in row 2 should be "Dan"
    And no errors should have been logged
