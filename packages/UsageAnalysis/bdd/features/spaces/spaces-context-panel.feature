@journey @spaces @realizes:views.space
Feature: What the context panel says about a space
  Clicking a space, and then a different one, and reading what the right-hand panel shows.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 9, 10 of the older
  reference and test 12 of the current spec), the part that needs no files.

  The panel is one element that swaps its content, so "contains the name" only means something if
  a second click makes it stop containing the first name — every claim here is paired that way.

  A child space is claimed in its parent's view rather than in the Browse tree: the tree builds a
  group's children when the group is opened and does not pick up one that arrives afterwards.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-CP-Root, BDD-CP-One, BDD-CP-Two" is on the server

  Scenario: Two children to switch between
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-CP-Root" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-CP-Root" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-CP-Root tree node inside browse tree
    And user enters "BDD-CP-One" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the Create Space dialog should close
    When user picks "Create Child Space..." from the context menu of BDD-CP-Root tree node inside browse tree
    And user enters "BDD-CP-Two" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the Create Space dialog should close
    When user double-clicks on BDD-CP-Root tree node inside browse tree
    Then the "BDD-CP-Root" view should be current
    And BDD-CP-One link in space gallery should be visible
    And BDD-CP-Two link in space gallery should be visible

  Scenario: Selecting a space shows its details
    When user clicks on BDD-CP-Root tree node inside browse tree
    Then context panel should be visible
    And context panel should contain text "BDD-CP-Root"
    And "Details" accordion header in context panel should be visible

  Scenario: The panel carries the sections a space has
    Then the following elements should be visible:
      | "Details" accordion header in context panel  |
      | "Content" accordion header in context panel  |
      | "Activity" accordion header in context panel |
      | "Sharing" accordion header in context panel  |
      | "Chats" accordion header in context panel    |

  Scenario: Clicking one child, then the other, switches the panel
    When user double-clicks on BDD-CP-Root tree node inside browse tree
    Then the "BDD-CP-Root" view should be current
    When user clicks on BDD-CP-One link in space gallery
    Then context panel should contain text "BDD-CP-One"
    When user clicks on BDD-CP-Two link in space gallery
    Then context panel should contain text "BDD-CP-Two"
    And context panel should not contain text "BDD-CP-One"
    When user clicks on BDD-CP-One link in space gallery
    Then context panel should contain text "BDD-CP-One"
    And context panel should not contain text "BDD-CP-Two"
