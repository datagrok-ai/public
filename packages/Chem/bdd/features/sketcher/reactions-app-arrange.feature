@journey @realizes:chem.reactions.app @eval-b7
Feature: Arrange in the Reaction Sketcher lays out the drawing and leaves the Text exactly as written
  Layout once dropped four steps' molecules from the Text. Arrange changes where the molecules are drawn, not the
  reaction, so the app keeps the person's Text as written (reaction-document.ts sameReactionText). How Arrange
  wraps a route to the canvas is Crux Sketch's own scenarios. Chem's package tests ('Reactions: app') cover Clean
  Up.

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  Scenario: The default route arranged: the Text as remembered
    Given user opens Chem's Reactions app on the default example
    When user remembers the Reactions app's Text
    And user chooses "Arrange" in the Clean Up and Layout menu of the Reactions app's sketcher
    Then the Reactions app's Text should be as remembered
    And no errors should have been logged
