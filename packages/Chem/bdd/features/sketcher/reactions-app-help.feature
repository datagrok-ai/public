@journey @realizes:chem.reactions.app @eval-b7
Feature: The Reaction Sketcher's help: its preview's buttons say what they are, and its dialogs open Chem's Reactions page
  The help pane said "Help unavailable" for Reset to example and for the Reaction dialog, and the preview's buttons
  had no tooltip. Chem ships a page on reactions (help/reactions.md). The app's view and its dialogs point the help
  pane at it, each at its own section. Chem's package tests check each dialog's link ('Reactions: help', 'Reactions:
  editor').

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  Scenario: The full-screen button names itself, and the full-screen dialog's help loads
    Given user opens Chem's Reactions app on the default example
    When user hovers over the Reactions app's full-screen button
    Then a tooltip should say "Full screen"
    When user clicks the Reactions app's full-screen button
    And user opens the help of the "Reaction" dialog
    Then help panel should show the text "The Reaction panel"
    And help panel should not show the text "Help unavailable"
    And no errors should have been logged
