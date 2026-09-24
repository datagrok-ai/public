@journey @serial @realizes:views.scripts
Feature: Deleting a script
  The context menu of a script in the Scripts view offers Delete; the confirmation names the script,
  CANCEL keeps it, YES removes its card and the script on the server. Translated from
  files/TestTrack/Scripts/delete.md and playwright-public/scripts/scripts-delete.test.ts, whose two
  tests checked the same thing; the CANCEL case is new — the old suite never tried it.

  Not translated, and why: nothing of the md is left out. The gallery counter that read "1" with no
  card left after a delete (survey, 22 Sep) is not claimed until it is reproduced by hand.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And a script "BddScriptDelete{time}" is on the server:
      """
      #language: r
      #input: dataframe table
      #output: int count
      count <- nrow(table) * ncol(table)
      """
    And user opens the Scripts view

  Scenario: The context menu offers Delete
    When user clears gallery search
    And user types "BddScriptDelete{time}" into gallery search
    Then "BddScriptDelete{time}" link in gallery should be visible
    When user opens the context menu of "BddScriptDelete{time}" link in gallery
    Then the open menu should list "Delete"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: CANCEL keeps the script
    When user picks "Delete" from the context menu of "BddScriptDelete{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete script \"BddScriptDelete{time}\"?"
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "BddScriptDelete{time}" link in gallery should be visible
    And 1 script named "BddScriptDelete{time}" should be on the server
    And no errors should have been logged

  Scenario: YES deletes the script
    When user picks "Delete" from the context menu of "BddScriptDelete{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "BddScriptDelete{time}" link in gallery should be absent
    And 0 scripts named "BddScriptDelete{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clears gallery search
