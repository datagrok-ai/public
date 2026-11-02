@journey @serial @spaces @realizes:views.space
Feature: Deleting a space
  The confirmation the platform asks for, what a cancel preserves, and what a confirmed delete removes.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (test 5).

  How far a delete reaches (a child's sibling stays, a parent takes its children) is tested in datlas
  spaces/spaces_test.dart ('delete child space preserves parent and siblings', 'delete space with child spaces cascade').
  The space is made through the API: the Create Space dialog is spaces-create's subject.

  Serial with spaces-create: twice in four full runs, spaces-create making root spaces on the other worker at the
  moment this feature's tree listed Spaces left BDD-Del out of the tree (it was on the server), and the right-click
  landed on the Spaces root. The mechanism is not pinned down (no server push of spaces reaches the tree); the two
  do not overlap any more.

  Edge cases: deleting a space is to remove the space and the links it holds, not a project that was
  only linked into it (core#976). Today the linked project is deleted with the space; that scenario,
  after the journey's own, is a known failure, and the tag goes when core#976 is fixed. Its space and
  project are named with the run's time, the project's letters and digits only (the Dashboards search
  misses "-" and "_"), and swept when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And Spaces tree node inside browse tree is expanded
    And a space named "BDD-Del" is on the server
    And Spaces tree node inside browse tree is expanded

  Scenario: Deleting asks first
    When user picks "Delete Space" from the context menu of Spaces---BDD-Del tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "BDD-Del"

  Scenario: Cancelling the confirmation keeps the space
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 1 space named "BDD-Del" should be on the server
    And Spaces---BDD-Del tree node inside browse tree should be visible

  Scenario: Confirming removes it from the server and the tree
    When user picks "Delete Space" from the context menu of Spaces---BDD-Del tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDD-Del" should be on the server
    And Spaces---BDD-Del tree node inside browse tree should be absent

  @edge-cases @known-failure:74 @realizes:core#976
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
