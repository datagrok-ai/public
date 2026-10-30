Feature: Data > Unpivot turns the merged columns into category and value rows of a new table
  Data > Unpivot... opens the Unpivot dialog over the current table: every column gets a role (not
  included, Copy or Merge), and OK opens a new "result" table view in the tall format (help
  `transform/unpivot.md`): one row per source row and merged column, the copied columns repeated, the
  merged column's name in the category column and its cell in the value column. OK stays disabled until
  a column is set to Merge. The names of the two new columns come from the dialog.

  The table is written here (two students, two subjects), so the result's shape and cells are exact.
  Left out on purpose: the unpivot of other column types and of many columns is the `Unpivot` function's
  job (ddt tests); the help has no undo for it, and the source table is left as it was. Nothing is put on
  the server.

  Background:
    Given user is logged in
    And user opens a table "scores" with:
      | student | math | art |
      | Ann     | 5    | 3   |
      | Bob     | 4    | 2   |

  # The role inputs carry no label of their own (MISSING.md): the 3rd input is student's, the 4th math's,
  # the 5th art's, after Category column name and Value column name.
  Scenario: Copy and Merge give one row per student and subject in a new table view
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
    And the value of "student" column in row 4 should be "Bob"
    And the value of "Value" column in row 4 should be "2"
    And table "scores" should have 2 rows
    And table "scores" should have columns "student, math, art"
    And no errors should have been logged

  Scenario: The names typed in the dialog name the category and value columns, and a column left out is dropped
    When user picks "Data > Unpivot..." from the top menu
    And user enters "subject" into "Category column name" input in Unpivot dialog
    And user enters "mark" into "Value column name" input in Unpivot dialog
    And user selects "Merge" in 4th input in Unpivot dialog
    And user selects "Merge" in 5th input in Unpivot dialog
    And user clicks on OK button in Unpivot dialog
    Then the open table views should be exactly "scores, result"
    Given user switches to the "result" table view
    Then the table should have 4 rows
    And the table should have the columns "subject, mark"
    And the value of "subject" column in row 2 should be "art"
    And the value of "mark" column in row 2 should be "3"
    And no errors should have been logged
