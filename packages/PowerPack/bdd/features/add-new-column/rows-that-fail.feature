@journey @realizes:powerpack.dialogs.add-new-column
Feature: Add New Column: the rows a formula fails on are left empty or filled, and their messages kept in a column
  The gear next to the Type choice (Rows that fail) opens what happens to a row the formula fails on: If
  a row fails, Leave empty or Use value with the value, and Error column with its name. Once one is set,
  the preview says how many of its rows failed and what they will get, with a Change... link that opens
  the same settings; with an error column, OK says how many rows of the table failed. The value's editor
  takes the type the formula returns, also when the gear is opened before the preview has computed it.
  On a three-row table written in the feature whose text columns hold one value that is not a date (d)
  and one that is not a number (n). Translated from the calculated-columns QA round of 13 Oct 2026
  (MISSING.md §1).

  Not translated: the value each row gets, the messages, the behaviour an edit starts from and the
  suggested error column name are the PowerPack package tests "Add new column: rows that fail" and "a
  suggested error column name follows the result name"; a renamed error column kept by a layout (core#787)
  and the type of a Use value column when every row fails (core#798) are what a layout stores and a
  function computes, the ddt tests that came with their fixes.

  Background:
    Given user is logged in
    And user opens a table "dates" with:
      | d          | n    |
      | 2024-01-05 | 12   |
      | soon       | many |
      | 2024-03-07 | 30   |

  Scenario: Rows that fail are left empty, their messages go to an error column, and the dialog says so
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "parsed" into column name input
    And user types "DateParse(${d})" into formula editor
    And user clicks on "Rows that fail" icon in "Add New Column" dialog
    Then "If a row fails" input should have value "Leave empty"
    When user checks "Error column" input
    Then "Column name" input should have value "parsed errors"
    And "Add New Column" dialog should contain text "1 of 3 preview rows failed and will be left empty."
    And "Add New Column" dialog should contain text "Messages are in \"parsed errors\"."
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And a warning balloon containing "1 of 3 rows failed and were left empty" should have been shown
    And the value of "parsed" column in row 2 should be ""
    And the value of "parsed errors" column in row 1 should be ""
    And some value of "parsed errors" column should contain "soon"
    And no errors should have been logged

  Scenario: Use value fills the rows that fail, and Change... in the preview opens the same settings
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "filled" into column name input
    And user types "ParseInt(${n})" into formula editor
    And user clicks on "Rows that fail" icon in "Add New Column" dialog
    And user selects "Use value" in "If a row fails" input
    And user enters "7" into "Value" input
    And user clicks on formula editor
    Then "Add New Column" dialog should contain text "1 of 3 preview rows failed and will be filled with 7."
    When user clicks on "Change..." link in "Add New Column" dialog
    Then "If a row fails" input should have value "Use value"
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the value of "filled" column in row 2 should be "7"
    And the value of "filled" column in row 3 should be "30"
    And no errors should have been logged
