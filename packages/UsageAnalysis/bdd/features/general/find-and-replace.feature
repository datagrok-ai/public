Feature: Find and Replace replaces values in the columns the dialog names
  Ctrl+H opens Find and Replace over the current table. Where columns picks the columns by a predicate
  (public#4287, core#807): Named lists them by name in the Select columns dialog, Selected Columns takes
  the columns selected in the grid; REPLACE applies the replacement and keeps the dialog, OK applies it
  and closes the dialog, and Ctrl+Z puts the values back. Both buttons stay disabled until the predicate
  names at least one column, and Find what says what is missing (core#838); a column picked after the
  search text is typed takes part (core#837). Filter matching rows shows the rows the search matches. An
  invalid regular expression and a replacement a number column cannot hold are refused under the input
  (core#890); a replace with no rows in scope says so and keeps the earlier replace undoable (core#891).

  Every scenario runs on a small table written here, so each claim names the cells it changes and the
  ones it must leave alone. The options of the dialog (Match case, Match whole word, regular expressions,
  Selected rows only, Filtered rows only) and the kinds of columns (text, number, date) each get one
  scenario that shows the dialog passes them on. The matrix of patterns — plain text taken literally
  (core#887), the preview flagging the rows the replace changes (core#888), regex groups (core#889),
  empty searches and replacements, number queries and ranges, dates, bools and the row scopes — is the
  ddt unit tests of `findAndReplace` (`core/shared/ddt/test/data_frame/data_utils_test.dart`). The dialog
  remembers its options after OK, so each scenario sets the ones it relies on. Nothing is put on the
  server.

  Background:
    Given user is logged in
    And user opens a table "fr" with:
      | name         | code | num |
      | Apple pie    | A.1  | 5   |
      | apple        | B(2) | 50  |
      | Pineapple    | C$3  | 500 |
      | banana split | a.b  | 7   |
    And the "name" cell of row 1 is the current object

  Scenario: Named replaces only in the columns picked by name, and Ctrl+Z puts the values back
    When user presses Control+H
    Then Find and Replace dialog should be visible
    When user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user toggles the "code" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user enters "a" into "Find what" input in Find and Replace dialog
    And user enters "@" into "Replace with" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 4 should be "b@n@n@ split"
    And the value of "code" column in row 4 should be "@.b"
    And the value of "name" column in row 1 should be "Apple pie"
    And Find and Replace dialog should be visible
    When user clicks on CANCEL button in Find and Replace dialog
    And user presses Control+Z
    Then the value of "name" column in row 4 should be "banana split"
    And the value of "code" column in row 4 should be "a.b"
    And no errors should have been logged

  Scenario: A column picked after the search text is typed takes part in the search (core#837)
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user checks "Filter matching rows" input in Find and Replace dialog
    Then 2 rows should pass the filter
    When user clicks on CANCEL button in Find and Replace dialog
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario Outline: REPLACE and OK stay disabled while the predicate names no column, and Find what says what is missing (core#838)
    When user presses Control+H
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user selects "<where>" in "Where columns" input in Find and Replace dialog
    Then REPLACE button in Find and Replace dialog should be disabled
    And OK button in Find and Replace dialog should be disabled
    When user hovers over "Find what" input in Find and Replace dialog
    Then tooltip should contain text "<missing>"
    When user selects "All Columns" in "Where columns" input in Find and Replace dialog
    Then REPLACE button in Find and Replace dialog should be enabled
    And OK button in Find and Replace dialog should be enabled
    And no errors should have been logged

    Examples:
      | where            | missing             |
      | Named            | Specify columns     |
      | Tagged           | Specify tag key     |
      | Of Type          | Specify type        |
      | Selected Columns | No matching columns |

  Scenario: Selected Columns replaces only in the columns selected in the grid
    When user clicks on the "header code" area of Grid viewer holding Control
    Then columns "code" should be selected
    When user presses Control+H
    And user selects "Selected Columns" in "Where columns" input in Find and Replace dialog
    And user enters "a" into "Find what" input in Find and Replace dialog
    And user enters "@" into "Replace with" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user clicks on OK button in Find and Replace dialog
    Then Find and Replace dialog should be hidden
    And the value of "code" column in row 4 should be "@.b"
    And the value of "name" column in row 4 should be "banana split"
    And no errors should have been logged

  Scenario: The current column is the default scope and leaves the other columns alone
    When user presses Control+H
    Then "Where columns" input in Find and Replace dialog should have value "Current Column"
    When user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user enters "a" into "Find what" input in Find and Replace dialog
    And user enters "@" into "Replace with" input in Find and Replace dialog
    And user clicks on OK button in Find and Replace dialog
    Then the value of "name" column in row 4 should be "b@n@n@ split"
    And the value of "code" column in row 4 should be "a.b"

  Scenario: Without Match case the search ignores case, with it only the exact case is replaced
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 1 should be "Apple pie"
    And the value of "name" column in row 2 should be "X"
    And the value of "name" column in row 3 should be "PineX"
    When user unchecks "Match case" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 1 should be "X pie"

  Scenario: Match whole word skips the search text inside a longer word
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user unchecks "Match case" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user checks "Match whole word" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 1 should be "X pie"
    And the value of "name" column in row 2 should be "X"
    And the value of "name" column in row 3 should be "Pineapple"

  Scenario: A regular expression anchored at both ends replaces only a whole value
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user unchecks "Match case" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user checks "Use regular expressions" input in Find and Replace dialog
    And user enters "^apple$" into "Find what" input in Find and Replace dialog
    And user enters "EXACT" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 2 should be "EXACT"
    And the value of "name" column in row 1 should be "Apple pie"
    And the value of "name" column in row 3 should be "Pineapple"

  Scenario: Selected rows only replaces in the selected rows and nowhere else
    When user selects rows where "name" is "Apple pie"
    And user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user unchecks "Match case" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user checks "Selected rows only" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 1 should be "X pie"
    And the value of "name" column in row 2 should be "apple"
    And the value of "name" column in row 3 should be "Pineapple"
    And 1 row should be selected

  Scenario Outline: Filtered rows only decides whether the rows the filter hides are replaced too
    When user opens the filter panel
    And user adds a range filter on "num" from 0 to 60
    Then 3 rows should pass the filter
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user checks "Match case" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user <action> "Filtered rows only" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user clicks on OK button in Find and Replace dialog
    Then the value of "name" column in row 2 should be "X"
    And the value of "name" column in row 3 should be "<hidden row>"
    And 3 rows should pass the filter

    Examples:
      | action   | hidden row |
      | checks   | Pineapple  |
      | unchecks | PineX      |

  Scenario: Filter matching rows narrows the user's filter and gives it back on Cancel
    When user opens the filter panel
    And user adds a range filter on "num" from 0 to 60
    Then 3 rows should pass the filter
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user checks "Match case" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user checks "Filter matching rows" input in Find and Replace dialog
    Then 1 row should pass the filter
    When user clicks on CANCEL button in Find and Replace dialog
    Then 3 rows should pass the filter

  Scenario: Numerical takes the number columns, searched with a comparison
    When user presses Control+H
    And user selects "Numerical" in "Where columns" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user enters ">10" into "Find what" input in Find and Replace dialog
    And user enters "-1" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "num" column in row 2 should be "-1"
    And the value of "num" column in row 3 should be "-1"
    And the value of "num" column in row 1 should be "5"
    And the value of "num" column in row 4 should be "7"
    And the value of "code" column in row 2 should be "B(2)"

  Scenario: A replacement a number column cannot hold disables REPLACE and OK, and says why
    When user presses Control+H
    And user selects "Numerical" in "Where columns" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user enters "50" into "Find what" input in Find and Replace dialog
    And user enters "abc" into "Replace with" input in Find and Replace dialog
    Then REPLACE button in Find and Replace dialog should be disabled
    And OK button in Find and Replace dialog should be disabled
    When user hovers over "Replace with" input in Find and Replace dialog
    Then tooltip should contain text "invalid value for \"num\" (int)"
    When user enters "0" into "Replace with" input in Find and Replace dialog
    Then OK button in Find and Replace dialog should be enabled

  Scenario: A date column matches a year and takes a date as the replacement
    Given user opens a table "events" with:
      | label | day        |
      | one   | 2020-01-15 |
      | two   | 2020-06-30 |
      | three | 2021-03-02 |
    And the "day" cell of row 1 is the current object
    When user presses Control+H
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user enters "2020" into "Find what" input in Find and Replace dialog
    And user enters "1999-01-01" into "Replace with" input in Find and Replace dialog
    And user clicks on OK button in Find and Replace dialog
    Then the "day" cell of row 1 should be displayed as "1/1/1999"
    And the "day" cell of row 2 should be displayed as "1/1/1999"
    And the "day" cell of row 3 should be displayed as "3/2/2021"

  Scenario: Enter applies the replacement, and the dialog opens again with the last settings
    When user presses Control+H
    And user selects "All Columns" in "Where columns" input in Find and Replace dialog
    And user unchecks "Match whole word" input in Find and Replace dialog
    And user unchecks "Use regular expressions" input in Find and Replace dialog
    And user checks "Match case" input in Find and Replace dialog
    And user enters "banana" into "Find what" input in Find and Replace dialog
    And user enters "B" into "Replace with" input in Find and Replace dialog
    And user presses Enter in "Replace with" input in Find and Replace dialog
    Then Find and Replace dialog should be hidden
    And the value of "name" column in row 4 should be "B split"
    When user presses Control+H
    Then "Find what" input in Find and Replace dialog should have value "banana"
    And "Replace with" input in Find and Replace dialog should have value "B"
    And "Match case" input in Find and Replace dialog should be checked
    When user clicks on CANCEL button in Find and Replace dialog

  Scenario: An invalid regular expression disables REPLACE and OK, and says why (core#890)
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user checks "Use regular expressions" input in Find and Replace dialog
    And user enters "(app" into "Find what" input in Find and Replace dialog
    Then "Find what" input in Find and Replace dialog should be invalid
    And REPLACE button in Find and Replace dialog should be disabled
    And OK button in Find and Replace dialog should be disabled
    When user hovers over "Find what" input in Find and Replace dialog
    Then tooltip should contain text "Invalid regular expression"
    When user clicks on CANCEL button in Find and Replace dialog
    Then no errors should have been logged

  Scenario: A replace with no rows in scope says so and keeps the earlier replace undoable (core#891)
    When user presses Control+H
    And user selects "Named" in "Where columns" input in Find and Replace dialog
    And user clicks on Columns input in Find and Replace dialog
    And user toggles the "name" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user checks "Match case" input in Find and Replace dialog
    And user enters "apple" into "Find what" input in Find and Replace dialog
    And user enters "X" into "Replace with" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then the value of "name" column in row 2 should be "X"
    When user checks "Selected rows only" input in Find and Replace dialog
    And user enters "banana" into "Find what" input in Find and Replace dialog
    And user clicks on REPLACE button in Find and Replace dialog
    Then an info balloon containing "no rows in scope" should have been shown
    And the value of "name" column in row 4 should be "banana split"
    When user clicks on CANCEL button in Find and Replace dialog
    And user presses Control+Z
    Then the value of "name" column in row 2 should be "apple"
