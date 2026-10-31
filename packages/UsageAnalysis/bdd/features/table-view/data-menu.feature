@journey
Feature: The Data menu's dialogs make a new table, new columns or new values from the open tables
  One claim per dialog of the table view's Data menu, by what OK (or the dialog's own button) visibly makes:
  Join Tables (help transform/join-tables.md) opens the join of two tables keyed by the column they share as
  a new "result" view, counts the key rows without a match and selects them from its menu, and disables OK
  when no key row matches; Compare Tables (compare-tables.md) needs a second open table, opens the rows the two
  differ in as a "diff" view with a changed value painted orange and a click on a diff row making the matching
  rows current in both tables, and refuses a key whose values repeat; Unpivot (unpivot.md) turns the merged
  columns into category and value rows of a new view; Split Column (text-to-columns.md) adds one column per
  part of a delimited text column; Batch Edit (batch-edit.md) writes a typed value into every selected row
  and follows the selection; Anonymize (anonymize-data.md) replaces the checked text and number columns in
  place, the number within the range 10^±r allows (10.5 × 10^±0.1 lies between 8.34 and 13.22); Categorize
  (categorize-data.md) bins a number column, and its SAVE COLUMN, SAVE MATRIX and APPLY COLOR CODING add a
  color-coded category column, a boolean column per category and a color coding of the source column.

  Each dialog works on small tables written here, so the counts and cells are exact; each group of
  scenarios starts by closing every view, so what one dialog left open never reaches the next. The dialogs'
  options (join types, Values, In-place, Unpivot names, Split Column's Prefix / By chars / Consecutive as
  one, Batch Edit's other row sets, Categorize without color coding, Compare without value columns) are the
  functions' tests: ddt join_test.dart, unpivot_test.dart and data_utils_test.dart (splitColumn), ApiTests
  dataframe/dataframe-join.ts and dataframe.ts ("unpivot"); an in-place join through the dialog is
  projects/projects-regressions-join.feature's. Categorize is last: its dialog has no OK or CANCEL and stays
  open. Nothing is put on the server.

  Background:
    Given user is logged in

  Scenario: Join Tables opens on both tables keyed by their shared column, and OK opens the inner join as a new table view
    Given user opens a table "left" with:
      | id | name |
      | 1  | a    |
      | 2  | b    |
      | 3  | c    |
    And user opens a table "right" with:
      | id | score |
      | 2  | 20    |
      | 3  | 30    |
      | 4  | 40    |
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should be visible
    And Key Columns input in Join Tables dialog should contain text "id"
    And "Join Type" input in Join Tables dialog should have value "inner"
    And "In-place" input in Join Tables dialog should be disabled
    And grid viewer in "Preview result columns" dialog should show 2 rows
    When user clicks on OK button in Join Tables dialog
    Then Join Tables dialog should be hidden
    And the open table views should be exactly "left, right, result"
    Given user switches to the "result" table view
    Then the table should have 2 rows
    And the table should have the columns "id, name, right.id, score"
    And the value of "name" column in row 1 should be "b"
    And the value of "score" column in row 1 should be "20"
    And no errors should have been logged

  Scenario: Join Tables counts the key rows without a match and its menu selects them, and with no key row matching it disables OK
    When user closes all views
    Given user opens a table "left" with:
      | id | name |
      | 1  | a    |
      | 2  | b    |
      | 3  | c    |
    And user opens a table "right" with:
      | id | score |
      | 2  | 20    |
      | 3  | 30    |
      | 4  | 40    |
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should contain text "1 / 3"
    When user clicks on first "context arrow down" icon in Join Tables dialog
    And user picks "Select All Not Matching" from the open menu
    Then 1 row of table "left" should be selected
    And 0 rows of table "right" should be selected
    When user clicks on CANCEL button in Join Tables dialog
    Then the open table views should be exactly "left, right"
    Given user sets "id" column in row 1 to "7"
    And user sets "id" column in row 2 to "8"
    And user sets "id" column in row 3 to "9"
    When user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should contain text "None of key rows are matched"
    And OK button in Join Tables dialog should be disabled
    When user clicks on CANCEL button in Join Tables dialog
    Then the open table views should be exactly "left, right"
    And no errors should have been logged

  Scenario: Compare Tables needs a second open table, and OK opens the rows that differ as a color-coded view linked to both tables
    When user closes all views
    Given user opens a table "people" with:
      | name | city | weight |
      | Ann  | Oslo | 10.5   |
      | Bob  | Rome | 10.5   |
      | Cid  | Oslo | 10.5   |
    Then the top menu should list "Data > Compare Tables..." as disabled
    Given user opens a table "people2" with:
      | name | city  | weight |
      | Ann  | Oslo  | 10.5   |
      | Bob  | Paris | 10.5   |
      | Dan  | Oslo  | 4      |
    Then the top menu should list "Data > Compare Tables..." as enabled
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
    And the "cell 1 of city" area of grid should contain the color "#FFA500"
    And the "cell 1 of name" area of grid should not contain the color "#FFA500"
    When user clicks on the "cell 1 of name" area of grid
    Then "name" of the current row of table "people" should be "Bob"
    And "name" of the current row of table "people2" should be "Bob"
    And no errors should have been logged

  Scenario: Compare Tables refuses a key whose values repeat: OK is disabled and the dialog says why
    When user closes all views
    Given user opens a table "people" with:
      | name | city | weight |
      | Ann  | Oslo | 10.5   |
      | Bob  | Rome | 10.5   |
      | Cid  | Oslo | 10.5   |
    And user opens a table "people2" with:
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

  # The role inputs carry no label of their own (MISSING.md): the 3rd input is student's, the 4th math's,
  # the 5th art's, after Category column name and Value column name.
  Scenario: Unpivot's Copy and Merge give one row per student and subject in a new table view
    When user closes all views
    Given user opens a table "scores" with:
      | student | math | art |
      | Ann     | 5    | 3   |
      | Bob     | 4    | 2   |
    When user picks "Data > Unpivot..." from the top menu
    Then Unpivot dialog should be visible
    And OK button in Unpivot dialog should be disabled
    When user selects "Copy" in 3rd input in Unpivot dialog
    And user selects "Merge" in 4th input in Unpivot dialog
    Then OK button in Unpivot dialog should be enabled
    When user selects "Merge" in 5th input in Unpivot dialog
    And user clicks on OK button in Unpivot dialog
    Then Unpivot dialog should be hidden
    And the open table views should be exactly "scores, result"
    Given user switches to the "result" table view
    Then the table should have 4 rows
    And the table should have the columns "student, Category, Value"
    And the value of "student" column in row 1 should be "Ann"
    And the value of "Category" column in row 1 should be "math"
    And the value of "Value" column in row 1 should be "5"
    And the value of "Category" column in row 2 should be "art"
    And the value of "Value" column in row 2 should be "3"
    And no errors should have been logged

  Scenario: Split Column adds one column per part of a delimited text column, the shorter values leaving the rest empty
    When user closes all views
    Given user opens a table "codes" with:
      | code   | n |
      | a;b    | 1 |
      | c;d;e  | 2 |
      | f      | 3 |
    And the "code" cell of row 1 is the current object
    When user picks "Data > Split Column..." from the top menu
    Then Split Column dialog should be visible
    And Column input in Split Column dialog should contain text "code"
    When user enters ";" into "Split by" input in Split Column dialog
    And user clicks on OK button in Split Column dialog
    Then Split Column dialog should be hidden
    And the table should have the columns "code, n, col1, col2, col3"
    And the value of "col1" column in row 1 should be "a"
    And the value of "col2" column in row 1 should be "b"
    And the value of "col3" column in row 1 should be ""
    And the value of "col3" column in row 2 should be "e"
    And the value of "code" column in row 2 should be "c;d;e"
    And no errors should have been logged

  Scenario: Batch Edit opens on the selected rows with their shared value, writes a typed value into them, and follows a new selection
    When user closes all views
    Given user opens a table "people" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cid  | Oslo | 3      |
    And user selects rows where "city" is "Oslo"
    When user picks "Data > Batch Edit..." from the top menu
    Then "Rows" input should have value "Selected"
    And value of "city" property should have text "Oslo"
    And value of "name" property should have text ""
    When user enters "Paris" into "city" property
    Then the "city" cell of row 1 of table "people" should be displayed as "Paris"
    And the "city" cell of row 3 of table "people" should be displayed as "Paris"
    And the "city" cell of row 2 of table "people" should be displayed as "Rome"
    And the "name" cell of row 1 of table "people" should be displayed as "Ann"
    When user selects rows of table "people" where "name" is one of "Bob"
    Then value of "name" property should have text "Bob"
    And value of "city" property should have text "Rome"
    And no errors should have been logged

  Scenario: Anonymize replaces the checked text and number columns in place and leaves the unchecked ones alone
    When user closes all views
    Given user opens a table "people" with:
      | name | city | w    | w0   |
      | Ann  | Oslo | 10.5 | 10.5 |
      | Bob  | Rome | 10.5 | 10.5 |
      | Cid  | Oslo | 10.5 | 10.5 |
    When user picks "Data > Anonymize..." from the top menu
    Then Anonymize Data dialog should be visible
    And "Number randomization factor" input in Anonymize Data dialog should have value "0.1"
    And the "city" column should not be checked in the column list of Anonymize Data dialog
    When user toggles the "city" column in the column list of Anonymize Data dialog
    And user toggles the "w" column in the column list of Anonymize Data dialog
    And user clicks on OK button in Anonymize Data dialog
    Then Anonymize Data dialog should be hidden
    And the value of "city" column in row 1 should be "city 0"
    And the value of "city" column in row 2 should be "city 1"
    And the value of "city" column in row 3 should be "city 0"
    And some value of "w" column should differ from "w0" column in the same row
    And every value of "w" column should lie between 8.34 and 13.22
    And every value of "w0" column should lie between 10.5 and 10.5
    And the value of "name" column in row 1 should be "Ann"
    And no errors should have been logged

  Scenario: Categorize bins a number column, and SAVE COLUMN, SAVE MATRIX and APPLY COLOR CODING each add what they name
    When user closes all views
    Given user opens a table "people" with:
      | name | city | v |
      | Ann  | Oslo | 1 |
      | Bob  | Rome | 2 |
      | Cid  | Oslo | 3 |
      | Dan  | Rome | 4 |
    When user picks "Data > Categorize..." from the top menu
    Then Categorize Data dialog should be visible
    When user selects "v" in "Categorize" input in Categorize Data dialog
    Then Categorize Data dialog should contain text "Matching: 0"
    When user clicks on "Bin..." link in Categorize Data dialog
    Then OK button in Edit dialog should be disabled
    When user enters "2" into "Number of bins" input in Edit dialog
    Then OK button in Edit dialog should be enabled
    When user clicks on OK button in Edit dialog
    Then Edit dialog should be hidden
    And Categorize Data dialog should contain text "Matching: 4"
    And Categorize Data dialog should contain text "Not matching: 0"
    When user clicks on "SAVE COLUMN" button in Categorize Data dialog
    Then the value of "category" column in row 1 should be "1 - 2.50"
    And the value of "category" column in row 3 should be "2.50 - 4"
    And "category" column should be color-coded categorically
    And "v" column should have no color coding
    And Categorize Data dialog should be visible
    When user clicks on "SAVE MATRIX" button in Categorize Data dialog
    Then "1 - 2.50" column should have type "bool"
    And the value of "1 - 2.50" column in row 1 should be "true"
    And the value of "1 - 2.50" column in row 3 should be "false"
    And the value of "2.50 - 4" column in row 4 should be "true"
    And the table should have 6 columns
    When user clicks on "APPLY COLOR CODING" button in Categorize Data dialog
    Then "v" column should be color-coded conditionally
    And the table should have 6 columns
    And no errors should have been logged
