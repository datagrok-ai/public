@journey
Feature: Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column
  Data > Compare Columns... opens the Compare Columns dialog: a table on each side, then a column on each
  side. With both columns picked it shows "Matched: N" and "Mismatched: M"; the Select link selects the
  mismatched rows of the first table, and ADD COLUMN adds a boolean "<first> == <second>" column to the first
  table, saying so in a balloon when the second column came from another table. Columns of tables with
  different row counts are not compared: the dialog says "Length mismatch" with each table's row count.
  No help page covers the dialog; the claims follow PowerPack `src/compare-columns.ts`.

  The function declares the top-menu path Data | Compare Columns..., but a table view's Data menu leaves it
  out (it takes no table input; MISSING.md, Findings): the first scenario is that known failure, and the
  others open the dialog by calling the function, the way Run a function does. ADD COLUMN over columns of two
  tables writes a formula that names the second table's column in the first table, so the column is false in
  every row (MISSING.md, Findings): the last scenario is that known failure.

  Every scenario runs on tables written here: cmp (a = 1..4, b agrees with a in rows 1 and 3), twin (four
  rows, its c agrees with cmp's a in rows 1 and 2) and short (two rows). Left out on purpose:
  - the comparison itself (`getString` equality of every row) is one rule; the scenarios show it once per
    outcome, not per column type;
  - the second table and second column selectors carry no label, so they are reached by position (MISSING.md).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "short" with:
      | x |
      | 1 |
      | 2 |
    And user opens a table "twin" with:
      | c |
      | 1 |
      | 2 |
      | 9 |
      | 9 |
    And user opens a table "cmp" with:
      | a | b |
      | 1 | 1 |
      | 2 | 3 |
      | 3 | 3 |
      | 4 | 5 |

  @known-failure:44
  Scenario: The Data menu offers Compare Columns
    Then the top menu should list "Data > Compare Columns..." as enabled
    And no errors should have been logged

  Scenario: Two columns of one table: the counts, Select picks the mismatched rows, ADD COLUMN adds the match column
    When user calls "PowerPack:_compareColumns" function
    Then "Compare Columns" dialog should be visible
    When user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "cmp" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Matched: 2"
    And "Compare Columns" dialog should contain text "Mismatched: 2"
    When user clicks on "Select" link in "Compare Columns" dialog
    Then only rows where "a" is one of "2, 4" should be selected
    When user clicks on "Add Column" button in "Compare Columns" dialog
    Then the table should have the columns "a, b, a == b"
    And the value of "a == b" column in row 1 should be "true"
    And the value of "a == b" column in row 2 should be "false"
    And the value of "a == b" column in row 3 should be "true"
    And the value of "a == b" column in row 4 should be "false"
    And "Compare Columns" dialog should be visible
    And no errors should have been logged
    When user clicks on CANCEL button in "Compare Columns" dialog
    And user removes "a == b" column
    Then the table should have the columns "a, b"

  Scenario: A column compared with itself matches in every row
    When user calls "PowerPack:_compareColumns" function
    And user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "cmp" in second choice input in "Compare Columns" dialog
    And user selects "b" in Columns input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Matched: 2"
    When user selects "b" in fourth choice input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Matched: 4"
    And "Compare Columns" dialog should contain text "Mismatched: 0"
    And no errors should have been logged
    When user clicks on CANCEL button in "Compare Columns" dialog
    Then "Compare Columns" dialog should be hidden

  Scenario: Columns of tables with different row counts are not compared
    When user calls "PowerPack:_compareColumns" function
    And user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "short" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    And user selects "x" in fourth choice input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Length mismatch"
    And "Compare Columns" dialog should contain text "cmp: 4 rows"
    And "Compare Columns" dialog should contain text "short: 2 rows"
    And "Compare Columns" dialog should not contain text "Matched"
    And "Add Column" button in "Compare Columns" dialog should be absent
    And no errors should have been logged
    When user clicks on CANCEL button in "Compare Columns" dialog
    Then "Compare Columns" dialog should be hidden

  @known-failure:109
  Scenario: A column of another table of the same length is compared, and ADD COLUMN adds the match column to the first table
    When user calls "PowerPack:_compareColumns" function
    And user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "twin" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    And user selects "c" in fourth choice input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Matched: 2"
    And "Compare Columns" dialog should contain text "Mismatched: 2"
    When user clicks on "Add Column" button in "Compare Columns" dialog
    Then an info balloon containing "Column was added to cmp" should have been shown
    And the table should have the columns "a, b, a == c"
    And the value of "a == c" column in row 1 should be "true"
    And the value of "a == c" column in row 2 should be "true"
    And the value of "a == c" column in row 3 should be "false"
    And no errors should have been logged
    When user clicks on CANCEL button in "Compare Columns" dialog
    And user removes "a == c" column
    Then the table should have the columns "a, b"
