@journey @spaces @realizes:views.space
Feature: Nested spaces and moving between them
  Three levels built the two ways the product offers — from the tree and from a card in a space's
  own view — and walked down and back up again. Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (tests 3 and 15).

  Which space is open is read from the platform (grok.shell.v.name) rather than from the URL the
  old spec matched with /\/s\//, which says a space is open but not which one.

  The part of test 15 that put a file in the parent through grok.dapi.spaces.id(...).files is left
  to the drag-and-drop feature, where a file gets into a space the way a person puts it there.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Hier-Root, BDD-Hier-Child, BDD-Hier-Grand" is on the server

  Scenario: A child is created from the tree
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Hier-Root" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Hier-Root" should be on the server
    And the "Create Space" dialog should close
    When user picks "Create Child Space..." from the context menu of BDD-Hier-Root tree node inside browse tree
    And user enters "BDD-Hier-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Hier-Child tree node inside browse tree should be visible

  Scenario: The parent's view lists the child
    When user double-clicks on BDD-Hier-Root tree node inside browse tree
    Then the "BDD-Hier-Root" view should be current
    And space gallery should be visible
    And BDD-Hier-Child link in space gallery should be visible

  Scenario: A grandchild is created from the child's card
    When user picks "Create Child Space..." from the context menu of BDD-Hier-Child link in space gallery
    And user enters "BDD-Hier-Grand" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And BDD-Hier-Grand tree node inside browse tree should be present

  Scenario: Opening the child shows the grandchild
    When user double-clicks on BDD-Hier-Child link in space gallery
    Then the "BDD-Hier-Child" view should be current
    And BDD-Hier-Grand link in space gallery should be visible
    And BDD-Hier-Child link in space gallery should be absent

  Scenario: Opening the grandchild leaves an empty space
    When user double-clicks on BDD-Hier-Grand link in space gallery
    Then the "BDD-Hier-Grand" view should be current
    And BDD-Hier-Grand link in space gallery should be absent

  Scenario: Going back up the tree finds the content again
    When user double-clicks on BDD-Hier-Root tree node inside browse tree
    Then the "BDD-Hier-Root" view should be current
    And BDD-Hier-Child link in space gallery should be visible
    When user double-clicks on BDD-Hier-Child link in space gallery
    Then the "BDD-Hier-Child" view should be current
    And BDD-Hier-Grand link in space gallery should be visible

  Scenario: Search still works after the walk
    When user enters "zzz-no-such-space" into space search
    Then BDD-Hier-Grand link in space gallery should be absent
    When user clears space search
    Then BDD-Hier-Grand link in space gallery should be visible
