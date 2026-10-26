@users @edge-cases @realizes:views.users
Feature: User edge cases: the running account cannot disable itself from the Users view
  "You can't disable your own account" (help govern/access-control/users-and-groups.md, Disabling
  accounts): the context menu of the running account's own card in Browse > Platform > Users offers
  no Disable..., while another user's card does (users-view.feature Users-14). Nothing is changed on
  the server: users cannot be deleted, and no account is disabled, edited or made here.

  The card is found by the running account's login, which the gallery shows as the user's name; on a
  stand where the two differ the scenario cannot find the card. That the server refuses a
  self-disable through the API is not claimed: the server has no such guard today (MISSING.md).

  Background:
    Given user is logged in
    And the browse panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Users" tree node inside browse tree

  Scenario: The running account's own context menu offers no Disable...
    Then the "Users" view should be current
    When user remembers the gallery counter
    And user types "{login}" into gallery search
    Then the gallery counter should be lower than remembered
    When user opens the context menu of "{login}" link in gallery
    Then the open menu should list "Details"
    And the open menu should not list "Disable..."
    When user closes the context menu
    And user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
