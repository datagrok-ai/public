@journey @spaces @realizes:views.space
Feature: Creating a space
  The Spaces node of the Browse tree: what it offers, what a new space needs, and what the platform
  refuses. Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 1, 2 and 18a) and
  SPACES-general-tests.md.

  A created space reaches the server up to 8 s after OK and the tree about a second later
  (measured on dev, 2026-09-10) — the platform announces neither, so every claim here is a poll.
  A child space is claimed from the tree alone: grok.dapi.spaces.list() returns root spaces only,
  and neither projects.list() nor the parent's `children` holds it, so "how many are on the
  server" cannot be asked about a child.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Root, BDD-Dup, BDD-Parent, BDD-Child, BDD Name With Spaces" is on the server

  Scenario: A root space is created from the Spaces node
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDD-Root" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the Create Space dialog should close
    And 1 space named "BDD-Root" should be on the server
    And the browse tree should show the "BDD-Root" space

  Scenario: The space offers its actions
    When user opens the context menu of BDD-Root tree node inside browse tree
    Then the open menu should list "Share..."
    And the open menu should list "Rename..."
    And the open menu should list "Delete Space"
    And the open menu should list "Create Child Space..."
    And the open menu should list "Add to favorites"
    And the open menu should not list "Duplicate"
    When user closes the context menu

  Scenario: An empty name disables OK, and typing one enables it again
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user clears Name input in Create Space dialog
    Then OK button in Create Space dialog should be disabled
    When user enters "BDD-Dup" into Name input in Create Space dialog
    Then OK button in Create Space dialog should be enabled
    When user clicks on CANCEL button in Create Space dialog
    Then the Create Space dialog should close
    And 0 spaces named "BDD-Dup" should be on the server

  Scenario: A second root space of the same name is refused
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Dup" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Dup" should be on the server
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Dup" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then an error balloon containing "Root project with same name already exists" should have been shown
    And 1 space named "BDD-Dup" should be on the server
    And Create Space dialog should be visible
    When user clicks on CANCEL button in Create Space dialog
    Then the Create Space dialog should close

  Scenario: A child space is created under a root space
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Parent" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Parent" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-Parent tree node inside browse tree
    And user enters "BDD-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the Create Space dialog should close
    And the browse tree should show the "BDD-Child" space

  Scenario: A second child of the same name is refused
    When user picks "Create Child Space..." from the context menu of BDD-Parent tree node inside browse tree
    And user enters "BDD-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then an error balloon containing "already exists" should have been shown
    And Create Space dialog should be visible
    When user clicks on CANCEL button in Create Space dialog
    Then the Create Space dialog should close

  Scenario: A name with spaces is kept as typed
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD Name With Spaces" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD Name With Spaces" should be on the server
    And the browse tree should show the "BDD Name With Spaces" space
