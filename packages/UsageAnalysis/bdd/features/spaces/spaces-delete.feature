@journey @spaces @realizes:views.space
Feature: Deleting a space
  The confirmation the platform asks for, what a cancel preserves, and how far a delete reaches.
  Translated from files/TestTrack/Spaces/spaces-general.test.ts (test 5).

  A child space cannot be counted on the server (grok.dapi.spaces.list returns root spaces only),
  so the cascade is claimed from the tree — which is also where a person would look.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Del, BDD-Del-Parent, BDD-Del-Child1, BDD-Del-Child2" is on the server

  Scenario: Deleting asks first
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Del" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Del" should be on the server
    When user picks "Delete Space" from the context menu of BDD-Del tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "BDD-Del"

  Scenario: Cancelling the confirmation keeps the space
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then "Are you sure?" dialog should be hidden
    And 1 space named "BDD-Del" should be on the server
    And BDD-Del tree node inside browse tree should be visible

  Scenario: Confirming removes it from the server and the tree
    When user picks "Delete Space" from the context menu of BDD-Del tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then "Are you sure?" dialog should be hidden
    And 0 spaces named "BDD-Del" should be on the server
    And BDD-Del tree node inside browse tree should be absent

  Scenario: Deleting one child leaves its sibling
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Del-Parent" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Del-Parent" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-Del-Parent tree node inside browse tree
    And user enters "BDD-Del-Child1" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Del-Child1 tree node inside browse tree should be visible
    When user picks "Create Child Space..." from the context menu of BDD-Del-Parent tree node inside browse tree
    And user enters "BDD-Del-Child2" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Del-Child2 tree node inside browse tree should be visible
    When user picks "Delete Space" from the context menu of BDD-Del-Child1 tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then BDD-Del-Child1 tree node inside browse tree should be absent
    And BDD-Del-Child2 tree node inside browse tree should be visible

  Scenario: Deleting the parent takes the remaining child with it
    When user picks "Delete Space" from the context menu of BDD-Del-Parent tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then 0 spaces named "BDD-Del-Parent" should be on the server
    And BDD-Del-Parent tree node inside browse tree should be absent
    And BDD-Del-Child2 tree node inside browse tree should be absent
