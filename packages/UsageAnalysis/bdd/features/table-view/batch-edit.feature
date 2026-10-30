Feature: Data > Batch Edit shows the values a set of rows shares and writes a typed value into every row of the set
  Data > Batch Edit... docks a batch editor next to the grid (help `transform/batch-edit.md`): its Rows input
  picks the set of rows (Selected by default), and one row per column shows the value when every row of the set
  holds the same one, and nothing when they differ. A value typed into a column's row is written into that
  column in every row of the set, and the rows outside it keep theirs. The editor follows the set: a new
  selection shows that selection's values.

  Every scenario runs on a three-row table written here, so each claim names the rows changed and the one left
  alone. Left out on purpose:
  - one scenario per Rows option the help lists as a set to edit (Selected, All); the other sets (Filtered,
    Current Row, Mouse Over Row, …) differ only in the rows they hand over, which is the `RowSet` code's
    job (ddt tests);
  - every column type's editor (numbers, dates, choices) is the property grid's, not this dialog's;
  - undo: neither the help nor `datagrok/navigation/undo.md` says a batch edit can be undone (MISSING.md).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cid  | Oslo | 3      |

  Scenario: The editor opens on the selected rows, shows their shared value, and a typed value goes into every selected row
    Given user selects rows where "city" is "Oslo"
    When user picks "Data > Batch Edit..." from the top menu
    Then "Rows" input should have value "Selected"
    And value of "city" property should have text "Oslo"
    And value of "name" property should have text ""
    When user enters "Paris" into "city" property
    Then the "city" cell of row 1 of table "people" should be displayed as "Paris"
    And the "city" cell of row 3 of table "people" should be displayed as "Paris"
    And the "city" cell of row 2 of table "people" should be displayed as "Rome"
    And the "name" cell of row 1 of table "people" should be displayed as "Ann"
    And no errors should have been logged

  Scenario: The editor follows the selection
    Given user selects rows where "name" is "Ann"
    When user picks "Data > Batch Edit..." from the top menu
    Then value of "name" property should have text "Ann"
    When user selects rows of table "people" where "name" is one of "Bob"
    Then value of "name" property should have text "Bob"
    And value of "city" property should have text "Rome"
    And no errors should have been logged

  Scenario: Rows set to All writes the value into every row
    Given user selects rows where "name" is "Ann"
    When user picks "Data > Batch Edit..." from the top menu
    And user selects "All" in "Rows" input
    Then value of "city" property should have text ""
    When user enters "Paris" into "city" property
    Then the "city" cell of row 1 of table "people" should be displayed as "Paris"
    And the "city" cell of row 2 of table "people" should be displayed as "Paris"
    And the "city" cell of row 3 of table "people" should be displayed as "Paris"
    And the "name" cell of row 2 of table "people" should be displayed as "Bob"
    And no errors should have been logged
