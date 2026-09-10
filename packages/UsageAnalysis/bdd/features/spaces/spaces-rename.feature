@journey @spaces @realizes:views.space
Feature: Renaming a space
  The Rename dialog over a space in the Browse tree and over a child card in a space view.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (test 4).

  Every rename is claimed twice: the name the server holds, and the name the tree draws. The old
  spec's duplicate check accepted "a toast appeared OR the original name is still there", which
  passes when nothing happens at all; here the claim is that the second name never took.

  A tree node inside a collapsed parent is not visible, so "should be absent" about a child is only
  worth anything once the parent is expanded — every claim about a child here is preceded by that
  gesture and paired with a node that must still be there.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Ren, BDD-Ren-New, BDD-Other, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew" is on the server

  Scenario: The Rename dialog opens on the current name
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Ren" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ren" should be on the server
    When user picks "Rename..." from the context menu of BDD-Ren tree node
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDD-Ren"

  Scenario: Cancelling the rename keeps the old name
    When user enters "BDD-Ren-New" into Name input in Rename project dialog
    And user clicks on CANCEL button in Rename project dialog
    Then Rename project dialog should be hidden
    And 1 space named "BDD-Ren" should be on the server
    And 0 spaces named "BDD-Ren-New" should be on the server
    And BDD-Ren tree node should be visible

  Scenario: A rename reaches the server and the tree
    When user picks "Rename..." from the context menu of BDD-Ren tree node
    And user enters "BDD-Ren-New" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then Rename project dialog should be hidden
    And 1 space named "BDD-Ren-New" should be on the server
    And 0 spaces named "BDD-Ren" should be on the server
    And BDD-Ren-New tree node should be visible
    And BDD-Ren tree node should be absent

  Scenario: Renaming onto an existing name is refused
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Other" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Other" should be on the server
    When user picks "Rename..." from the context menu of BDD-Ren-New tree node
    And user enters "BDD-Other" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then an error balloon containing "already exists" should have been shown
    And 1 space named "BDD-Ren-New" should be on the server
    And 1 space named "BDD-Other" should be on the server
    When user closes Rename project dialog
    Then Rename project dialog should be hidden

  Scenario: A child space is renamed from its card in the parent
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Ren-Parent" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ren-Parent" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-Ren-Parent tree node
    And user enters "BDD-Ren-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Ren-Child tree node should be visible
    When user double-clicks on BDD-Ren-Parent tree node
    Then BDD-Ren-Child link in space gallery should be visible
    When user picks "Rename..." from the context menu of BDD-Ren-Child link in space gallery
    Then Name input in Rename project dialog should have value "BDD-Ren-Child"
    When user enters "BDD-Ren-ChildNew" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then BDD-Ren-ChildNew link in space gallery should be visible
    And BDD-Ren-Child link in space gallery should be absent
    When user expands BDD-Ren-Parent tree node
    Then BDD-Ren-ChildNew tree node should be visible
    And BDD-Ren-Child tree node should be absent
