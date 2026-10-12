@journey @spaces @realizes:views.space
Feature: Deleting a space
  The confirmation the platform asks for, what a cancel preserves, and what a confirmed delete removes.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (test 5).

  How far a delete reaches (a child's sibling stays, a parent takes its children) is tested in datlas
  spaces/spaces_test.dart ('delete child space preserves parent and siblings', 'delete space with child spaces cascade').
  The space is made through the API: the Create Space dialog is spaces-create's subject.

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
