@browse @realizes:views.browse
Feature: The shell's modes around the Browse panel
  Presentation mode puts the panels away and gives them back (its "back to design mode" link, or F7
  again — Escape does not leave it), and the Tabs toggle switches simple mode, which hides the view
  tabs. Translated from the manual
  case Browse-Modes-01 (browse_manual_tests2.md section 17, playwright-public/browse/splitmodes.test.ts,
  whose checks were best-effort `if`s).

  The toggles are the status bar's own, by their labels; the view tabs are named after their views.

  Browse-Split-01 and -02 (Split right / Split down from the view's menu) are not translated: the
  command is not in the platform any more — the view's context menu offers Close others, Rename...,
  Reset, View and Table, and no core source names a split — so the case needs product's word first.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: Presentation mode puts the panels away and its back link gives them back
    Given user opens demog-1000 dataset
    # a table view brings its Toolbox, docked as a tab over Browse
    And the toolbox pane is hidden
    And the browse panel is open
    Then the browse tree should be visible
    When user clicks on "Presentation mode" status bar toggle
    Then the browse tree should be hidden
    And status bar should be hidden
    And grid should show 1000 rows
    When user clicks on "back to design mode" link
    Then status bar should be visible
    And the browse tree should be visible
    And the "demog-1000" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: F7 enters presentation mode and leaves it
    Given user opens demog-1000 dataset
    And status bar should be visible
    When user presses F7
    Then status bar should be hidden
    And "back to design mode" link should be visible
    When user presses F7
    Then status bar should be visible
    And "back to design mode" link should be absent
    And no errors should have been logged

  # the browse panel step leaves simple mode, so the tabs start shown
  Scenario: The Tabs toggle hides the view tabs and shows them again
    Given user opens demog-1000 dataset
    Then "demog-1000" view tab should be visible
    When user clicks on "Tabs" status bar toggle
    Then "demog-1000" view tab should be hidden
    When user clicks on "Tabs" status bar toggle
    Then "demog-1000" view tab should be visible
    And the "demog-1000" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown
