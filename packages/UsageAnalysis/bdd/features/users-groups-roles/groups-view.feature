@journey @serial @groups @realizes:views.groups
Feature: The Groups view
  Browse > Platform > Groups as an administrator sees it: the list and its toolbar, the view modes,
  search, a group's context menu, context panel and chat. Translated from files/TestTrack/User
  groups/groups_manual_tests.md (Groups-01, 02, 06 to 08, 10, 16, 17, 20) and
  playwright-public/user groups/groups.test.ts.

  The group every claim is about is made by the feature and deleted at its end, with the hidden
  group its chat lives in. The groups search is fuzzy (other BDD groups of a parallel run come up
  too), and a card may be on the page before the server's result lands, so a search is claimed by
  the counter dropping below the list's first — it keeps its old number until the result lands —
  and only then by the group's link.

  Elsewhere: Groups-03 to 05, 09 and 14 are groups-lifecycle.feature, Groups-11 to 13, 15 and 18
  groups-members.feature. Groups-16 is requested by the account the feature runs as, not by a second
  signed-in user: a feature has one page; the request it sends goes with the group at the feature's
  end (Request membership threw a NullError before GROK-20906 was fixed in the core on 2026-09-15).

  Not translated: Groups-19 (favorites) — a group has no favorites entry in its context menu and no
  star in the context panel.

  Background:
    Given user is logged in
    And a group named "BDD-GV-Group-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Groups" tree node inside browse tree

  Scenario: The view opens from the Browse tree (Groups-01)
    Then the "Groups" view should be current
    And the page address should contain "/groups"
    And gallery should be visible
    And gallery counter should be visible
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The toolbar carries its controls (Groups-02)
    Then the following elements should be visible:
      | "New Group..." button                               |
      | gallery search                                      |
      | "Switch to brief view" icon inside gallery toolbar  |
      | "Switch to card view" icon inside gallery toolbar   |
      | "Switch to grid view" icon inside gallery toolbar   |
      | "Refresh" icon inside gallery toolbar               |
    And "New Group..." button should be enabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The view-mode icons switch the gallery's render mode (Groups-07)
    Then the gallery should be in brief mode
    When user clicks on "Switch to card view" icon inside gallery toolbar
    Then the gallery should be in card mode
    And "Switch to card view" icon inside gallery toolbar should be selected
    When user clicks on "Switch to grid view" icon inside gallery toolbar
    Then the gallery should be in grid mode
    And "Switch to grid view" icon inside gallery toolbar should be selected
    When user clicks on "Switch to brief view" icon inside gallery toolbar
    Then the gallery should be in brief mode
    And "Switch to brief view" icon inside gallery toolbar should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Searching by name finds the one group, and clearing brings the rest back (Groups-06)
    When user types "BDD-GV-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD-GV-Group-{time}" link in gallery should be visible
    And the page address should contain "?q=BDD-GV-Group-{time}"
    When user remembers the gallery counter
    And user clears gallery search
    Then the gallery counter should be higher than remembered
    And the page address should not contain "?q="
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A group's context menu (Groups-08)
    When user types "BDD-GV-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user opens the context menu of "BDD-GV-Group-{time}" link in gallery
    Then the open menu should list "Properties..."
    And the open menu should list "Request membership"
    And the open menu should list "Chat"
    And the open menu should list "Delete"
    When user closes the context menu
    And user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Selecting a group fills the context panel (Groups-10)
    When user types "BDD-GV-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-GV-Group-{time}" link in gallery
    Then the context panel should show "BDD-GV-Group-{time}"
    And the following elements should be visible:
      | "Actions" accordion header in context panel            |
      | "Members" accordion header in context panel            |
      | "Favorites" accordion header in context panel          |
      | "Global Permissions" accordion header in context panel |
      | "Permissions" accordion header in context panel        |
    And MANAGE button in "Members" section in context panel should be visible
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Chat opens the group's own chat (Groups-17)
    When user types "BDD-GV-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Chat" from the context menu of "BDD-GV-Group-{time}" link in gallery
    Then the "Chats" view should be current
    And chat header should contain text "BDD-GV-Group-{time}"
    When user closes the current view
    Then the "Groups" view should be current
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Request membership answers without an error (Groups-16)
    When user types "BDD-GV-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Request membership" from the context menu of "BDD-GV-Group-{time}" link in gallery
    Then an info balloon should have been shown
    And no errors should have been logged
    And no error or warning balloon should have been shown
