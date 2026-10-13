@journey @spaces @realizes:views.space
Feature: What the context panel says about a space
  Clicking a space, and then a different one, and reading what the right-hand panel shows.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 9, 10 of the older
  reference and test 12 of the current spec), the part that needs no files.

  The panel is one element that swaps its content, so "contains the name" only means something if
  a second click makes it stop containing the first name — every claim here is paired that way.

  A child space is claimed in its parent's view rather than in the Browse tree: the tree builds a
  group's children when the group is opened and does not pick up one that arrives afterwards. The
  spaces are made through the API: the Create Space dialog is spaces-create's subject.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    And a space named "BDD-CP-Root" is on the server
    And a space named "BDD-CP-One" under "BDD-CP-Root" is on the server
    And a space named "BDD-CP-Two" under "BDD-CP-Root" is on the server
    And Spaces tree node inside browse tree is expanded

  Scenario: The root's view lists its two children
    When user double-clicks on Spaces---BDD-CP-Root tree node inside browse tree
    Then the "BDD-CP-Root" view should be current
    And BDD-CP-One link in gallery should be visible
    And BDD-CP-Two link in gallery should be visible

  Scenario: Selecting a space shows its details
    When user clicks on Spaces---BDD-CP-Root tree node inside browse tree
    Then context panel should be visible
    And the context panel should show "BDD-CP-Root"
    And "Details" accordion header in context panel should be visible

  # A context pane that counts its items hides itself while the count is 0 (accordion.css,
  # `.grok-prop-panel .d4-accordion-pane[d4-info="0"]`). Activity counts the space's audit log,
  # and the server writes the creation entry a moment after the space exists: on a space created
  # seconds ago the pane is in the DOM and shown only once that entry has landed.
  Scenario: The panel carries the sections a space has
    Then the following elements should be visible:
      | "Details" accordion header in context panel  |
      | "Content" accordion header in context panel  |
      | "Sharing" accordion header in context panel  |
      | "Chats" accordion header in context panel    |
    And "Activity" accordion header in context panel should be present

  Scenario: Clicking one child, then the other, switches the panel
    When user double-clicks on Spaces---BDD-CP-Root tree node inside browse tree
    Then the "BDD-CP-Root" view should be current
    When user clicks on BDD-CP-One link in gallery
    Then the context panel should show "BDD-CP-One"
    When user clicks on BDD-CP-Two link in gallery
    Then the context panel should show "BDD-CP-Two"
    And context panel should not contain text "BDD-CP-One"
    When user clicks on BDD-CP-One link in gallery
    Then the context panel should show "BDD-CP-One"
    And context panel should not contain text "BDD-CP-Two"
