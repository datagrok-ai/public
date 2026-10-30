Feature: Data > Anonymize replaces the values of the checked columns in place
  Data > Anonymize... opens the Anonymize Data dialog over the current table (help `transform/anonymize-data.md`):
  a list of the table's columns to check, and the Number randomization factor (0.1 by default). OK replaces,
  in place, the values of the checked columns only: a text column's categories become "<column> 0",
  "<column> 1", …, and a number is multiplied by 10 to a random power between -r and r. The columns left
  unchecked keep their values.

  The change cannot be undone (`datagrok/navigation/undo.md`: anonymizing data is not recorded), so every
  scenario runs on a table written here. Left out on purpose:
  - the random numbers themselves: the claim is the range the help's formula allows (10.5 × 10^±0.1 lies
    between 8.34 and 13.22) and that the values moved; the formula over other factors and types (ints, dates)
    is the API/ddt tests' job;
  - the "category" column, which maps a text column onto the categories of another open one-column table:
    its cell editor is a choice inside the dialog's grid, and the mapping is a computation (MISSING.md);
  - an invalid randomization factor logs an error (MISSING.md, finding).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | name | city | w    | w0   |
      | Ann  | Oslo | 10.5 | 10.5 |
      | Bob  | Rome | 10.5 | 10.5 |
      | Cid  | Oslo | 10.5 | 10.5 |

  Scenario: The dialog lists the columns unchecked, and OK replaces the categories of a checked text column only
    When user picks "Data > Anonymize..." from the top menu
    Then Anonymize Data dialog should be visible
    And "Number randomization factor" input in Anonymize Data dialog should have value "0.1"
    And the "city" column should not be checked in the column list of Anonymize Data dialog
    When user toggles the "city" column in the column list of Anonymize Data dialog
    And user clicks on OK button in Anonymize Data dialog
    Then Anonymize Data dialog should be hidden
    And the value of "city" column in row 1 should be "city 0"
    And the value of "city" column in row 2 should be "city 1"
    And the value of "city" column in row 3 should be "city 0"
    And the value of "name" column in row 1 should be "Ann"
    And every value of "w" column should lie between 10.5 and 10.5
    And no errors should have been logged

  Scenario: A checked number column is randomized within the factor's range, and an unchecked one is left alone
    When user picks "Data > Anonymize..." from the top menu
    And user toggles the "w" column in the column list of Anonymize Data dialog
    And user clicks on OK button in Anonymize Data dialog
    Then Anonymize Data dialog should be hidden
    And some value of "w" column should differ from "w0" column in the same row
    And every value of "w" column should lie between 8.34 and 13.22
    And every value of "w0" column should lie between 10.5 and 10.5
    And the value of "city" column in row 2 should be "Rome"
    And no errors should have been logged
