@serial @spaces @edge-cases @realizes:views.space
Feature: Space edge cases: an emptied name, a child space moved to the root, a delete that takes a linked project with it
  Guards for the Rename dialog of a space, for a child space dragged onto the Spaces node, and for
  what deleting a space removes. The Rename dialog refuses a name with no letter in it, an emptied
  one included (`DbEntityMeta.getRenameDialog`), as the Create Space dialog does (spaces-create).

  A child space dropped onto the Spaces node asks "Move <name> out?" in the "Move space" dialog;
  MOVE saves it as a root space (`_initSpacesDragDrop` in browse_panel_tree.dart). A root name must
  be unique: a child named like an existing root space is refused with the server's message in an
  error balloon, the dialog stays open, and the child keeps its parent.

  Deleting a space removes the space and the links it holds, not a project that was only linked into
  it (core#976). Today the linked project is deleted with the space; that scenario is a known failure,
  and the tag goes when core#976 is fixed.

  The spaces and the project are this feature's own, named with the run's time, made through the API
  (the Create Space dialog is spaces-create's subject) and swept when the feature starts and ends.
  The project's name is letters and digits only: the Dashboards search misses "-" and "_".

  Serial with the other features that drag onto or list the Spaces tree.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}" is on the server
    And Spaces tree node inside browse tree is expanded

  Scenario: An emptied name disables OK in the Rename dialog
    Given a space named "BDDEdgeKeep{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of BDDEdgeKeep{time} tree node inside browse tree
    Then Rename project dialog should be visible
    And OK button in Rename project dialog should be enabled
    When user clears Name input in Rename project dialog
    Then Name input in Rename project dialog should be invalid
    And OK button in Rename project dialog should be disabled
    When user types "x" into Name input in Rename project dialog
    Then OK button in Rename project dialog should be enabled
    When user clicks on CANCEL button in Rename project dialog
    Then Rename project dialog should be hidden
    And 1 space named "BDDEdgeKeep{time}" should be on the server
    And no errors should have been logged

  Scenario: A child space dragged onto the Spaces node becomes a root space
    Given a space named "BDDEdgeParent{time}" is on the server
    And a space named "BDDEdgeOut{time}" under "BDDEdgeParent{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    And Spaces---BDDEdgeParent{time} tree node inside browse tree is expanded
    When user drags Spaces---BDDEdgeParent{time}---BDDEdgeOut{time} tree node inside browse tree to Spaces tree node inside browse tree
    Then "Move space" dialog should be visible
    And "Move space" dialog should contain text "Move BDDEdgeOut{time} out?"
    When user clicks on MOVE button in "Move space" dialog
    Then the "Move space" dialog should close
    And Spaces---BDDEdgeOut{time} tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A child space named like a root space is refused at the root and keeps its parent
    Given a space named "BDDEdgeTaken{time}" is on the server
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

  @known-failure:104 @realizes:core#976
  Scenario: Deleting a space keeps a project that was only linked into it
    core#976: the project linked into the space is deleted with the space. Remove the tag when
    core#976 is fixed.
    Given no project named "bddedgelinkproj{time}" is on the server
    And a space named "BDDEdgeLink{time}" is on the server
    And user opens demog-1000 dataset
    When user saves the current view as project "bddedgelinkproj{time}"
    And user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "bddedgelinkproj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then bddedgelinkproj{time} gallery card should be visible
    Given Spaces tree node inside browse tree is expanded
    When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDEdgeLink{time} tree node inside browse tree
    Then the "BDDEdgeLink{time}" view should be current
    And "bddedgelinkproj{time}" link in gallery should be visible
    And 1 project named "bddedgelinkproj{time}" should be on the server
    When user picks "Delete Space" from the context menu of BDDEdgeLink{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDEdgeLink{time}" should be on the server
    And 1 project named "bddedgelinkproj{time}" should be on the server
