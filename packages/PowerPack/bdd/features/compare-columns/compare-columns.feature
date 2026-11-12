@journey
Feature: Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column
  Data > Compare Columns... opens the Compare Columns dialog: a table on each side, then a column on each
  side. With both columns picked it shows "Matched: N" and "Mismatched: M"; the Select link selects the
  mismatched rows of the first table, and ADD COLUMN adds a boolean "<first> == <second>" column to the first
  table, saying so in a balloon when the second column came from another table (public#4475: that column was
  false in every row). Columns of tables with different row counts are not compared: the dialog says "Length
  mismatch". No help page covers the dialog; the claims follow PowerPack `src/compare-columns.ts`.

  Every scenario runs on tables written here: cmp (a = 1..4, b agrees with a in rows 1 and 3), twin (four
  rows, its c agrees with cmp's a in rows 1 and 2) and short (two rows). The second table and second column
  selectors carry no label, so they are reached by position. The comparison itself (`getString` equality of every
  row) is one rule: a column against itself or per column type is not claimed (MISSING.md). Nothing is put on the
  server.

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

  @realizes:core#1246
  Scenario: Two columns of one table: the counts, Select picks the mismatched rows, ADD COLUMN adds the match column
    When user picks "Data > Compare Columns..." from the top menu
    Then "Compare Columns" dialog should be visible
    When user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "cmp" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Matched: 2"
    And "Compare Columns" dialog should contain text "Mismatched: 2"
    When user clicks on "Select" link in "Compare Columns" dialog
    Then only rows where "a" is one of "2, 4" should be selected
    When user clicks on "Add Column" button in "Compare Columns" dialog
    Then a new column "a == b" should have been added
    And the value of "a == b" column in row 1 should be "true"
    And the value of "a == b" column in row 2 should be "false"
    When user clicks on CANCEL button in "Compare Columns" dialog

  Scenario: Columns of tables with different row counts are not compared
    When user picks "Data > Compare Columns..." from the top menu
    And user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "short" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    And user selects "x" in fourth choice input in "Compare Columns" dialog
    Then "Compare Columns" dialog should contain text "Length mismatch"
    And "Add Column" button in "Compare Columns" dialog should be absent
    When user clicks on CANCEL button in "Compare Columns" dialog

  @realizes:public#4475
  Scenario: A column of another table of the same length is compared, and ADD COLUMN adds the match column to the first table
    When user picks "Data > Compare Columns..." from the top menu
    And user selects "cmp" in Tables input in "Compare Columns" dialog
    And user selects "twin" in second choice input in "Compare Columns" dialog
    And user selects "a" in Columns input in "Compare Columns" dialog
    And user selects "c" in fourth choice input in "Compare Columns" dialog
    And user clicks on "Add Column" button in "Compare Columns" dialog
    Then an info balloon containing "Column was added to cmp" should have been shown
    And the value of "a == c" column in row 1 should be "true"
    And the value of "a == c" column in row 3 should be "false"
