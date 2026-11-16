@journey @realizes:chem.reactions.app @eval-b7
Feature: The Reaction Sketcher's Text follows its drawing: a lone molecule is no error
  Drawing one molecule turned the Text into a MOL file and the preview into an error dump. Now Crux Sketch's value for a
  drawing with no arrow, a MOL block, reaches the Text as the molecule's SMILES. The status and the preview say in one
  short line that a reaction needs an arrow. An arrow drawn after the molecule makes it a reaction with no product yet.
  Clear canvas empties both the drawing and the Text, clicked as soon as the app shows its tools. What the Text holds
  for each drawing is Chem's package tests ('Reactions: app').

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  # Clear canvas clicked as the app's first open showed its tools once did nothing (the sketcher had not read the route
  # yet) and the route was then drawn, its Text as it was. First in the feature, so that it is
  # the page's first open; the sketcher's box is covered until it has read the reaction, and the click waits for that.
  Scenario: Clear canvas as the app opens empties the drawing and the Text; the preview asks for a reaction
    Given user opens Chem's Reactions app on the default example
    When user clicks on the Reactions app's Clear canvas
    Then the Reactions app's drawing should be empty
    And the Reactions app's Text should be ""
    And the Reactions app's preview should ask for a reaction
    And no errors should have been logged

  Scenario: Ethane drawn alone is no error; an arrow after it makes it a reaction
    Given user opens Chem's Reactions app on an empty reaction
    When user clicks on crux single bond tool
    And user clicks on Crux canvas 30% across and 50% down
    Then the Reactions app's Text should be "CC"
    And the Reactions app's preview should say that a reaction needs an arrow
    And no text of the Reactions app should contain "RDKit", "Crux", "M  END" or "confidence"
    When user clicks on crux reaction arrow tool
    And user drags on Crux canvas from 50% to 75% across, 50% down
    Then the Reactions app's Text should be "CC>>"
    And no errors should have been logged
