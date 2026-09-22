@journey @serial @realizes:views.scripts
Feature: Editing a script
  A saved script opens in the editor from the Scripts view by a double-click; Save stays disabled
  until the code changes, an appended line is saved to the server, and the editor shows it when the
  script is opened again. Translated from files/TestTrack/Scripts/edit.md and
  playwright-public/scripts/scripts-edit-debugged.test.ts, which never checked the server.

  The script is this feature's own ({time} in its name), is never run (so no container is needed),
  and is deleted with its chats at the end.

  Not translated, and why: nothing of the md is left out. The Save button's state is read from the
  class the ribbon gives it — it has no aria-disabled yet; switch to "Save button should be
  disabled" once the core names land.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And a script "BddScriptEdit{time}" is on the server:
      """
      #language: r
      #input: dataframe table
      #output: int count
      count <- nrow(table) * ncol(table)
      """
    And user opens the Scripts view

  Scenario: A double-click opens the script in the editor, with nothing to save
    When user clears gallery search
    And user types "BddScriptEdit{time}" into gallery search
    And user double-clicks on "BddScriptEdit{time}" link in gallery
    Then the "BddScriptEdit{time}" view should be current
    And code editor should contain the text "count <- nrow(table) * ncol(table)"
    And the Save button of the script view should be disabled
    And no errors should have been logged

  Scenario: An appended line is saved to the server
    When user appends "newParam = \"test\"" to code editor
    Then the Save button of the script view should be enabled
    When user saves the script
    Then an info balloon containing "Script saved." should have been shown
    And the script "BddScriptEdit{time}" on the server should contain "newParam = \"test\""
    And the Save button of the script view should be disabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Opened again, the editor shows the saved line
    When user closes the current view
    Then the "Scripts" view should be current
    When user double-clicks on "BddScriptEdit{time}" link in gallery
    Then the "BddScriptEdit{time}" view should be current
    And code editor should contain the text "newParam = \"test\""
    And no errors should have been logged
    When user closes the current view
    And user clears gallery search
