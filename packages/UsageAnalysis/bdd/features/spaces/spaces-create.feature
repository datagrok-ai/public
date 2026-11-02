@journey @serial @spaces @realizes:views.space
Feature: Creating a space
  The Spaces node of the Browse tree: what it offers, what a new space needs, and what the platform
  refuses. Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 1, 2 and 18a) and
  SPACES-general-tests.md.

  A created space reaches the server up to 8 s after OK and the tree about a second later
  (measured on dev, 2026-09-10) — the platform announces neither, so every claim here is a poll.
  A child space is claimed from the tree alone: grok.dapi.spaces.list() returns root spaces only,
  and neither projects.list() nor the parent's `children` holds it, so "how many are on the
  server" cannot be asked about a child.

  A second child of the same name and a name with spaces are server rules, tested in datlas
  spaces/spaces_test.dart ('forbid duplicate subspace names', 'create space with space in name and verify files').

  Serial with spaces-delete, whose tree lost its fixture space while this feature made root spaces (twice in four
  full runs; see spaces-delete.feature).

  Edge cases: a child space dropped onto the Spaces node asks "Move <name> out?" in the "Move space"
  dialog; MOVE saves it as a root space (`_initSpacesDragDrop` in browse_panel_tree.dart). A root name
  must be unique: a child named like an existing root space is refused with the server's message in an
  error balloon, the dialog stays open, and the child keeps its parent. Those two spaces and their
  parent are made through the API and named with the run's time.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Root, BDD-Dup, BDD-Parent, BDD-Child" is on the server

  Scenario: A root space is created from the Spaces node
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDD-Root" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDD-Root" should be on the server
    And Spaces---BDD-Root tree node inside browse tree should be visible

  Scenario: The space offers its actions
    When user opens the context menu of Spaces---BDD-Root tree node inside browse tree
    Then the open menu should list "Share..."
    And the open menu should list "Rename..."
    And the open menu should list "Delete Space"
    And the open menu should list "Create Child Space..."
    And the open menu should list "Add To Favorites"
    And the open menu should not list "Duplicate"
    When user closes the context menu

  Scenario: An empty name disables OK, and typing one enables it again
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user clears Name input in Create Space dialog
    Then OK button in Create Space dialog should be disabled
    When user enters "BDD-Dup" into Name input in Create Space dialog
    Then OK button in Create Space dialog should be enabled
    When user clicks on CANCEL button in Create Space dialog
    Then the "Create Space" dialog should close
    And 0 spaces named "BDD-Dup" should be on the server

  Scenario: A second root space of the same name is refused
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Dup" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Dup" should be on the server
    And the "Create Space" dialog should close
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Dup" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then an error balloon containing "Root project with same name already exists" should have been shown
    And 1 space named "BDD-Dup" should be on the server
    And Create Space dialog should be visible
    When user clicks on CANCEL button in Create Space dialog
    Then the "Create Space" dialog should close

  Scenario: A child space is created under a root space
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Parent" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Parent" should be on the server
    And the "Create Space" dialog should close
    When user picks "Create Child Space..." from the context menu of Spaces---BDD-Parent tree node inside browse tree
    And user enters "BDD-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And BDD-Child tree node inside browse tree should be visible

  @edge-cases
  Scenario: A child space dropped on the Spaces node is moved out to the root
    Given Spaces---BDD-Parent tree node inside browse tree is expanded
    When user drags Spaces---BDD-Parent---BDD-Child tree node inside browse tree to Spaces tree node inside browse tree
    Then "Move space" dialog should be visible
    And "Move space" dialog should contain text "Move BDD-Child out?"
    When user clicks on MOVE button in "Move space" dialog
    Then the "Move space" dialog should close
    And Spaces---BDD-Child tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @edge-cases
  Scenario: A child named like a root space is refused at the root and keeps its parent
    Given no space named "BDDEdgeTaken{time}, BDDEdgeParent{time}" is on the server
    And a space named "BDDEdgeTaken{time}" is on the server
    And a space named "BDDEdgeParent{time}" is on the server
    And a space named "BDDEdgeTaken{time}" under "BDDEdgeParent{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    And Spaces---BDDEdgeParent{time} tree node inside browse tree is expanded
    Then Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible
    And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible
    When user drags Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree to Spaces tree node inside browse tree
    Then "Move space" dialog should be visible
    When user clicks on MOVE button in "Move space" dialog
    Then an error balloon containing "already exists" should have been shown
    And "Move space" dialog should be visible
    When user clicks on CANCEL button in "Move space" dialog
    Then the "Move space" dialog should close
    And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible
    And Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible
