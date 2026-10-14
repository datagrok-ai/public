@journey @realizes:powerpack.dialogs.add-new-column
Feature: Append Tables keeps the formula of a calculated column and computes it for the appended rows
  Two tables written in the feature: "first", with a calculated column dbl = ${x} * 2 made through the
  platform's API, and "second", which has only x. Data > Append Tables... with both checked opens the
  "result" table: its dbl column keeps the formula and holds the values the formula gives for the
  rows that came from "second" too (public#4068). OK with no table checked is a known failure: it
  logs a NullError (core#784). Translated from the QA round of 13 Oct 2026 on calculated columns.

  Not translated: appendAsync, appendMergeAsync, grok.data.appendTables and recalculateFormulaColumns
  keeping the formulas are ApiTests "DataFrame: Calculated columns"; this feature claims the menu
  command once.

  Background:
    Given user is logged in
    And user opens a table "first" with:
      | x |
      | 1 |
      | 2 |
    And user adds a calculated column "dbl" with formula "${x} * 2"
    And user opens a table "second" with:
      | x  |
      | 10 |

  Scenario: Both tables appended: the formula stays and covers the appended rows
    When user picks "Data > Append Tables..." from the top menu
    Then "Append Tables" dialog should be visible
    When user checks "first" checkbox in "Append Tables" dialog
    And user checks "second" checkbox in "Append Tables" dialog
    And user clicks on OK button in "Append Tables" dialog
    Then "Append Tables" dialog should be hidden
    And the open table views should be exactly "first, second, result"
    Given user switches to the "result" table view
    Then the table should have 3 rows
    And the table should have the columns "x, dbl"
    And "dbl" column should have tag "formula" equal to "${x} * 2"
    And the value of "dbl" column in row 1 should be "2"
    And the value of "dbl" column in row 3 should be "20"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  # core#784: OK with no table checked fails with a NullError
  @known-failure:47
  Scenario: OK with no table checked logs no error
    When user picks "Data > Append Tables..." from the top menu
    Then "Append Tables" dialog should be visible
    When user clicks on OK button in "Append Tables" dialog
    Then no errors should have been logged
