@journey @spaces @realizes:views.space
Feature: A space in favorites
  Adding a space to favorites and taking it out again, claimed where the platform puts it: the
  Browse tree names the entry after its path (tree-My-stuff---Favorites---BDD-Fav). Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (test 6).

  The claim is presence in the tree rather than visibility, because it is the stronger one: a node
  that is merely inside a collapsed group would satisfy "absent" without anything having been
  removed.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Fav" is on the server

  Scenario: A space is added to favorites
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-Fav" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Fav" should be on the server
    And the "Create Space" dialog should close
    And "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be absent
    When user picks "Add to favorites" from the context menu of BDD-Fav tree node inside browse tree
    Then "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be present

  Scenario: A space is removed from favorites
    When user picks "Remove from favorites" from the context menu of BDD-Fav tree node inside browse tree
    Then "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be absent
    And "My stuff > Favorites" tree node inside browse tree should be present
    And 1 space named "BDD-Fav" should be on the server
    And BDD-Fav tree node inside browse tree should be visible
