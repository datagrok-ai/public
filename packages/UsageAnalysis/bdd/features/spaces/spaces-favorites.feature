@journey @spaces @realizes:views.space
Feature: A space in favorites
  Adding a space to favorites and taking it out again, claimed where the platform puts it: the
  Browse tree names the entry after its path (tree-My-stuff---Favorites---BDD-Fav). Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (test 6).

  The claim is presence in the tree rather than visibility, because it is the stronger one: a node
  that is merely inside a collapsed group would satisfy "absent" without anything having been
  removed. Add To Favorites is a submenu for an account that administers groups; "Only for me"
  adds the space to the user's own favorites, and picking it again takes the space out.

  For an account that administers a group, as the running one does, "Add To Favorites" is a submenu
  of "Only for me" and the groups, and its item is a check: picking it again takes the space out
  (GROK-21108 replaced "Remove from favorites" with it).

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Fav" is on the server

  Scenario: A space is added to favorites
    Given a space named "BDD-Fav" is on the server
    And Spaces tree node inside browse tree is expanded
    And "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be absent
    When user picks "Add To Favorites > Only for me" from the context menu of Spaces---BDD-Fav tree node inside browse tree
    Then "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be present

  Scenario: A space is removed from favorites
    When user picks "Add To Favorites > Only for me" from the context menu of Spaces---BDD-Fav tree node inside browse tree
    Then "My stuff > Favorites > BDD-Fav" tree node inside browse tree should be absent
    And "My stuff > Favorites" tree node inside browse tree should be present
    And 1 space named "BDD-Fav" should be on the server
    And Spaces---BDD-Fav tree node inside browse tree should be visible
