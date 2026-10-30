Feature: Data > Categorize bins a number column into categories and saves them as a column, a matrix or a color coding
  Data > Categorize... opens the Categorize Data dialog over a column of the current table (help
  `transform/categorize-data.md`): a row per pattern with its category and the number of rows it matches, and
  the Rows / Matching / Not matching counts under them. Bin... fills the patterns with the given number of
  ranges of a number column. SAVE COLUMN adds a "category" column holding each row's category, color-coded
  when Color Coding is on; SAVE MATRIX adds a boolean column per category; APPLY COLOR CODING colors the source
  column by the categories. The dialog stays open after each button.

  Every scenario runs on a four-row table written here (v = 1..4, so two bins split it 2 + 2). Left out on
  purpose:
  - typed patterns (regular expressions on text, conditions like "> 70" on numbers): the pattern and category
    editors carry no name, label or aria-label, so no phrase reaches them (MISSING.md, needs a name in core);
    the patterns themselves are the `NumericMatcher` / RegExp code's job (ddt tests);
  - undo: neither the help nor `datagrok/navigation/undo.md` says the added columns can be undone; the
    table is the scenario's own.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | name | city | v |
      | Ann  | Oslo | 1 |
      | Bob  | Rome | 2 |
      | Cid  | Oslo | 3 |
      | Dan  | Rome | 4 |

  Scenario: Bin fills the patterns, and SAVE COLUMN adds a color-coded category column
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
    And the value of "category" column in row 2 should be "1 - 2.50"
    And the value of "category" column in row 3 should be "2.50 - 4"
    And the value of "category" column in row 4 should be "2.50 - 4"
    And "category" column should be color-coded categorically
    And "v" column should have no color coding
    And Categorize Data dialog should be visible
    And no errors should have been logged

  Scenario: With Color Coding off, SAVE COLUMN adds the category column without colors
    When user picks "Data > Categorize..." from the top menu
    And user selects "v" in "Categorize" input in Categorize Data dialog
    And user clicks on "Bin..." link in Categorize Data dialog
    And user enters "2" into "Number of bins" input in Edit dialog
    And user clicks on OK button in Edit dialog
    And user unchecks "Color Coding" input in Categorize Data dialog
    And user clicks on "SAVE COLUMN" button in Categorize Data dialog
    Then the value of "category" column in row 4 should be "2.50 - 4"
    And "category" column should have no color coding
    And no errors should have been logged

  Scenario: SAVE MATRIX adds a boolean column per category
    When user picks "Data > Categorize..." from the top menu
    And user selects "v" in "Categorize" input in Categorize Data dialog
    And user clicks on "Bin..." link in Categorize Data dialog
    And user enters "2" into "Number of bins" input in Edit dialog
    And user clicks on OK button in Edit dialog
    And user clicks on "SAVE MATRIX" button in Categorize Data dialog
    Then "1 - 2.50" column should have type "bool"
    And the value of "1 - 2.50" column in row 1 should be "true"
    And the value of "1 - 2.50" column in row 3 should be "false"
    And the value of "2.50 - 4" column in row 1 should be "false"
    And the value of "2.50 - 4" column in row 4 should be "true"
    And the table should have 5 columns
    And no errors should have been logged

  Scenario: APPLY COLOR CODING colors the source column by the categories
    When user picks "Data > Categorize..." from the top menu
    And user selects "v" in "Categorize" input in Categorize Data dialog
    And user clicks on "Bin..." link in Categorize Data dialog
    And user enters "2" into "Number of bins" input in Edit dialog
    And user clicks on OK button in Edit dialog
    And user clicks on "APPLY COLOR CODING" button in Categorize Data dialog
    Then "v" column should be color-coded conditionally
    And the table should have 3 columns
    And no errors should have been logged
