@journey @spaces @realizes:views.space
Feature: Renaming a space
  The Rename dialog over a space in the Browse tree and over a child card in a space view.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (test 4).

  Every rename is claimed twice: the name the server holds, and the name the tree draws. The old
  spec's duplicate check accepted "a toast appeared OR the original name is still there", which
  passes when nothing happens at all; here the claim is that the second name never took.

  A node inside a collapsed parent is not visible, so "should be absent" about a child is only
  worth anything once the parent is expanded — every claim about a child here is preceded by that
  gesture and paired with a node that must still be there.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Ren, BDD-Ren-New, BDD-Other, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew" is on the server

  Scenario: The Rename dialog opens on the current name
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Ren" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ren" should be on the server
    When user picks "Rename..." from the context menu of BDD-Ren tree node inside browse tree
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDD-Ren"

  Scenario: Cancelling the rename keeps the old name
    When user enters "BDD-Ren-New" into Name input in Rename project dialog
    And user clicks on CANCEL button in Rename project dialog
    Then Rename project dialog should be hidden
    And 1 space named "BDD-Ren" should be on the server
    And 0 spaces named "BDD-Ren-New" should be on the server
    And the browse tree should show the "BDD-Ren" space

  Scenario: A rename reaches the server and the tree
    When user picks "Rename..." from the context menu of BDD-Ren tree node inside browse tree
    And user enters "BDD-Ren-New" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then Rename project dialog should be hidden
    And 1 space named "BDD-Ren-New" should be on the server
    And 0 spaces named "BDD-Ren" should be on the server
    And the browse tree should show the "BDD-Ren-New" space
    And the browse tree should not show the "BDD-Ren" space

  Scenario: Renaming onto an existing name is refused
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Other" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Other" should be on the server
    When user picks "Rename..." from the context menu of BDD-Ren-New tree node inside browse tree
    And user enters "BDD-Other" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then an error balloon containing "already exists" should have been shown
    And 1 space named "BDD-Ren-New" should be on the server
    And 1 space named "BDD-Other" should be on the server
    When user closes Rename project dialog
    Then Rename project dialog should be hidden

  Scenario: A child space is renamed from its card in the parent
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Ren-Parent" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ren-Parent" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-Ren-Parent tree node inside browse tree
    And user enters "BDD-Ren-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the browse tree should show the "BDD-Ren-Child" space
    When user double-clicks on BDD-Ren-Parent tree node inside browse tree
    Then the space should show the "BDD-Ren-Child" card
    When user picks "Rename..." from the context menu of BDD-Ren-Child link in space gallery
    Then Name input in Rename project dialog should have value "BDD-Ren-Child"
    When user enters "BDD-Ren-ChildNew" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the space should show the "BDD-Ren-ChildNew" card
    And the space should not show the "BDD-Ren-Child" card
    When user expands BDD-Ren-Parent tree node inside browse tree
    Then the browse tree should show the "BDD-Ren-ChildNew" space
    And the browse tree should not show the "BDD-Ren-Child" space
