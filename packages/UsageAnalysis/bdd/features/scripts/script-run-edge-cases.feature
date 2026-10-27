@serial @edge-cases @realizes:views.scripts
Feature: Script run edge cases: Run... with no table to give
  Guards for the Run... dialog of a script whose table input has nothing to take. Every input is
  required unless it is nullable, and the input is validated before the function runs: "you will get
  an error if a required parameter is not specified" (help
  `datagrok/concepts/functions/func-params-annotation.md`, "Initial values and optional parameters").
  With no table open the dialog's Table input is empty, yet OK is enabled, and OK runs the script
  without a table: it throws "Cannot read properties of undefined (reading 'rowCount')". That is a
  known failure (core#1155). The same check for a scalar input (empty, out
  of its bounds, of the wrong type) is claimed once, on the query dialog
  (queries/query-parameters-edge-cases.feature), where it holds; a Run... with a table open is
  scripts-run.feature.

  The script is the feature's own, a JavaScript script (it runs in the page, see the bdd library's
  CLAUDE.md, "What never becomes a feature"), named with the run's time and deleted with its chats
  when the feature ends. Serial: it works in the Scripts view, whose search text is the account's
  own setting, like the other scripts features.

  Background:
    Given user is logged in
    And a script "BddScriptNoTable{time}" is on the server:
      """
      //language: javascript
      //input: dataframe table
      //output: int count
      count = table.rowCount;
      """
    And no table view is open
    And user opens the Scripts view

  @known-failure:41 @realizes:core#1155
  Scenario: Run... with no table open keeps OK disabled
    The defect: with no table open the Table input is empty and OK stays enabled; OK then runs the
    script without its table, which throws (core#1155). Remove the tag when the dialog keeps OK
    disabled while its required table input is empty.
    When user clears gallery search
    And user types "BddScriptNoTable{time}" into gallery search
    And user picks "Run..." from the context menu of "BddScriptNoTable{time}" link in gallery
    Then "BddScriptNoTable{time}" dialog should be visible
    And Table input in "BddScriptNoTable{time}" dialog should have value ""
    And OK button in "BddScriptNoTable{time}" dialog should be disabled
    When user clicks on CANCEL button in "BddScriptNoTable{time}" dialog
    Then the "BddScriptNoTable{time}" dialog should close
    And no errors should have been logged
    When user clears gallery search
