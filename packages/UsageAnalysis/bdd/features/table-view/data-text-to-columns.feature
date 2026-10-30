Feature: Data > Split Column splits a text column into new columns by a delimiter
  Data > Split Column... (help `transform/text-to-columns.md`, "Text to columns") opens the Split Column
  dialog on the current text column. OK adds to the same table one column per part, named by the Prefix
  and a number (col1, col2, …), a shorter value leaving the rest empty. Split by is taken as one
  delimiter unless By chars splits on each of its characters, and Consecutive as one treats repeated
  delimiters as one (help, test page `text-to-columns-test.md`).

  The table is written here, so every new cell is exact; each claim also names the source column, which
  must stay as it was. Left out on purpose: type auto-detection and the delimiter matrix (letters,
  digits, spaces) are the `SplitColumn` function's job (ddt tests); one scenario per option shows the
  dialog passes it on. The help says nothing about undo, so none is claimed (MISSING.md, question).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "codes" with:
      | code   | n |
      | a;b    | 1 |
      | c;d;e  | 2 |
      | f      | 3 |
    And the "code" cell of row 1 is the current object

  Scenario: Split by a delimiter adds one column per part, the shorter values leaving the rest empty
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
    And the value of "col1" column in row 3 should be "f"
    And the value of "col2" column in row 3 should be ""
    And the value of "code" column in row 2 should be "c;d;e"
    And the table should have 3 rows
    And no errors should have been logged

  Scenario: Prefix names the new columns
    When user picks "Data > Split Column..." from the top menu
    And user enters ";" into "Split by" input in Split Column dialog
    And user enters "part" into "Prefix" input in Split Column dialog
    And user clicks on OK button in Split Column dialog
    Then the table should have the columns "code, n, part1, part2, part3"
    And the value of "part2" column in row 2 should be "d"
    And no errors should have been logged

  Scenario: By chars splits on each character of Split by
    Given user sets "code" column in row 1 to "a;b,c"
    When user picks "Data > Split Column..." from the top menu
    And user enters ";," into "Split by" input in Split Column dialog
    And user checks "By chars" input in Split Column dialog
    And user clicks on OK button in Split Column dialog
    Then the table should have the columns "code, n, col1, col2, col3"
    And the value of "col1" column in row 1 should be "a"
    And the value of "col2" column in row 1 should be "b"
    And the value of "col3" column in row 1 should be "c"
    And no errors should have been logged

  Scenario: Consecutive as one treats repeated delimiters as one
    Given user sets "code" column in row 1 to "a;;b"
    When user picks "Data > Split Column..." from the top menu
    And user enters ";" into "Split by" input in Split Column dialog
    And user checks "Consecutive as one" input in Split Column dialog
    And user clicks on OK button in Split Column dialog
    Then the value of "col1" column in row 1 should be "a"
    And the value of "col2" column in row 1 should be "b"
    And the value of "col3" column in row 2 should be "e"
    And no errors should have been logged
