@realizes:sequencetranslator.app.oligo-toolkit
Feature: Oligo Pattern: the default example pattern and the strand switches
  The Oligo Pattern app (Browse > Apps > Peptides > Oligo Toolkit > Oligo Pattern) opens on the
  shipped <default example> pattern of the current user: the Load and Edit blocks, the pattern
  picture and a translation example for both strands; Save stays disabled until something is
  edited, and saving over the example is refused. Switching the antisense strand off hides its
  length and its example. Translated from the TestTrack case SequenceTranslator/oligo-pattern-designer,
  Block A (steps 1 and 3) and Block B step 1.

  Parked in the request document: the trash icon of Load (Block A step 2 and all of Block C: the
  icon has no name or label), the Edit strands dialog (its per-position inputs have no caption),
  and the save / reload / delete of a pattern (Blocks B and C: a saved pattern stays in the user's
  storage, and no step removes it at feature start and end). Kept without: the picture drawing a
  single strand and the sense example being 10 characters long (the example text areas have no
  names).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Peptides tree node inside browse tree is expanded
    And Apps---Peptides---Oligo-Toolkit tree node inside browse tree is expanded
    When user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree
    Then the "Oligo Pattern" view should be current
    When user hovers over "Load" heading
    Then tooltip should be hidden

  Scenario: The app opens on the default example, and saving over it is refused
    Then Author input should contain the text "(me)"
    And Pattern input should have the value "<default example>"
    And "Translation example" heading should be visible
    And "Sense strand" heading should be visible
    And "Anti sense" heading should be visible
    And Save button should be disabled
    When user enters "22" into "Sense strand length" input
    Then Save button should be enabled
    When user clicks on Save button
    Then a warning balloon containing "Cannot save default pattern" should have been shown
    And Pattern input should have the value "<default example>"
    And no errors should have been logged

  Scenario: Switching the antisense strand off hides its length and its example
    Then "Anti sense length" input should be visible
    And "Anti sense" heading should be visible
    When user enters "10" into "Sense strand length" input
    And user switches off "Anti sense strand" checkbox
    Then "Anti sense length" input should be hidden
    And "Anti sense" heading should be hidden
    And "Sense strand" heading should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged
