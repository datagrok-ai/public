@serial @realizes:sequencetranslator.app.oligo-toolkit
Feature: Oligo Pattern: the default example, the strand switches, saving, reloading and deleting a pattern
  The Oligo Pattern app (Browse > Apps > Peptides > Oligo Toolkit > Oligo Pattern) opens on the
  shipped <default example> pattern of the current user: the Load and Edit blocks, the pattern
  picture and a translation example for both strands; Save stays disabled until something is
  edited, and saving over the example is refused, as is deleting it. Switching the antisense strand
  off hides its length and its example. An edited pattern is saved to the user's storage, comes back
  when the app is reopened, and Delete removes the pattern chosen in Load, not the one named in Edit
  (GROK-16674). Translated from the TestTrack case SequenceTranslator/oligo-pattern-designer.

  A saved pattern lives in the account's shared user settings (`OligoToolkit`, written back as one
  map), so the Background removes the two patterns the scenarios save now and when the feature ends,
  and the feature is serial: another page of the account saving at the same time would overwrite it. The Edit strands dialog's per-position
  inputs are named after the position their row shows ("Sense modification 1"), and the translation
  examples after their strand and role ("Sense example output"). Kept without: the pattern picture
  being drawn (an SVG with no reading).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And no oligo pattern named "BddPattern-{time}" is in the user's settings
    And no oligo pattern named "BddDelete-{time}" is in the user's settings
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

  Scenario: The example pattern cannot be deleted
    Then Pattern input should have the value "<default example>"
    When user clicks on "Delete pattern" button
    Then a warning balloon containing "Cannot delete example pattern" should have been shown
    And "Delete pattern" dialog should be absent
    And Pattern input should have the value "<default example>"
    And no errors should have been logged

  Scenario: Switching the antisense strand off hides its length and its example
    Then "Anti sense length" input should be visible
    And "Anti sense" heading should be visible
    And "Antisense example input" text area should be visible
    When user enters "10" into "Sense strand length" input
    And user switches off "Anti sense strand" checkbox
    Then "Anti sense length" input should be hidden
    And "Anti sense" heading should be hidden
    And "Antisense example input" text area should be absent
    And "Sense strand" heading should be visible
    And "Sense example input" text area should have the value "AGCUAGCUAG"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An edited pattern is saved, and reopening the app loads it back
    When user enters "10" into "Sense strand length" input
    And user switches off "Anti sense strand" checkbox
    Then "Sense example input" text area should have the value "AGCUAGCUAG"
    When user remembers the value of "Sense example output" text area
    And user clicks on "Edit strands" button
    Then "Edit strands" dialog should be visible
    When user checks "All PTO" checkbox in "Edit strands" dialog
    And user selects "2'-Fluoro" in "Sense modification 1" choice input in "Edit strands" dialog
    And user clicks on OK button in "Edit strands" dialog
    Then the "Edit strands" dialog should close
    And "Sense example output" text area should not have the remembered value
    When user remembers the value of "Sense example output" text area
    And user enters "BddPattern-{time}" into "Pattern name" text area
    And user clicks on Save button
    Then an info balloon containing "Pattern BddPattern-{time} saved" should have been shown
    And Pattern input should offer the choice "BddPattern-{time}"
    When user closes the current view
    And user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree
    Then the "Oligo Pattern" view should be current
    When user selects "BddPattern-{time}" in Pattern input
    Then "Sense strand length" input should have the value "10"
    And "Anti sense strand" checkbox should be switched off
    And "Sense example output" text area should have the remembered value
    And no errors should have been logged

  Scenario: Delete removes the pattern chosen in Load, not the one named in Edit (GROK-16674)
    When user enters "10" into "Sense strand length" input
    And user enters "BddDelete-{time}" into "Pattern name" text area
    And user clicks on Save button
    Then an info balloon containing "Pattern BddDelete-{time} saved" should have been shown
    When user selects "BddDelete-{time}" in Pattern input
    And user enters "OtherName-{time}" into "Pattern name" text area
    And user clicks on "Delete pattern" button
    Then "Delete pattern" dialog should contain the text "Are you sure you want to delete pattern BddDelete-{time}?"
    And "Delete pattern" dialog should not contain the text "OtherName-{time}"
    When user clicks on OK button in "Delete pattern" dialog
    Then the "Delete pattern" dialog should close
    And Pattern input should not offer the choice "BddDelete-{time}"
    When user closes the current view
    And user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree
    Then the "Oligo Pattern" view should be current
    And Pattern input should not offer the choice "BddDelete-{time}"
    And no errors should have been logged
