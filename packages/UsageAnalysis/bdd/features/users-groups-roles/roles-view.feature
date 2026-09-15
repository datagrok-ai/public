@journey @serial @roles @realizes:views.roles
Feature: The Roles view
  Browse > Platform > Roles as an administrator sees it: the list and its toolbar, the view modes,
  search, a role's context menu and context panel. Translated from files/TestTrack/User
  groups/roles_manual_tests.md (Roles-01, 02, 06 to 08, 10, 16) and
  playwright-public/user groups/roles.test.ts.

  A role is a group on the server that the Roles view lists, and the JS API cannot make one, so the
  role the claims are about is made in the New Role dialog and deleted at the end of the feature.
  The search is fuzzy, so a search for it is claimed by the counter dropping below the list's first
  before its link is read. The server count of roles by a name counts groups; the Roles view listing the role
  is what shows it is one.

  Elsewhere: Roles-03 to 05, 09 and 15 are roles-lifecycle.feature, Roles-11 to 14
  roles-assignment.feature.

  Background:
    Given user is logged in
    And no role named "BDD-RV-Role-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Roles" tree node inside browse tree

  Scenario: The view opens from the Browse tree (Roles-01)
    Then the "Roles" view should be current
    And the page address should contain "/roles"
    And gallery should be visible
    And gallery counter should be visible
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A role to look at
    When user clicks on "New Role..." button
    And user types "BDD-RV-Role-{time}" into Name input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And 1 role named "BDD-RV-Role-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The toolbar carries its controls (Roles-02)
    Then the following elements should be visible:
      | "New Role..." button                                |
      | gallery search                                      |
      | "Switch to brief view" icon inside gallery toolbar  |
      | "Switch to card view" icon inside gallery toolbar   |
      | "Switch to grid view" icon inside gallery toolbar   |
      | "Refresh" icon inside gallery toolbar               |
    And "New Role..." button should be enabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The view-mode icons switch the gallery's render mode (Roles-07)
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

  Scenario: Searching by name finds the one role, and clearing brings the rest back (Roles-06)
    When user types "BDD-RV-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD-RV-Role-{time}" link in gallery should be visible
    And the page address should contain "?q=BDD-RV-Role-{time}"
    When user remembers the gallery counter
    And user clears gallery search
    Then the gallery counter should be higher than remembered
    And the page address should not contain "?q="
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A role's context menu (Roles-08)
    When user types "BDD-RV-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user opens the context menu of "BDD-RV-Role-{time}" link in gallery
    Then the open menu should list "Properties..."
    And the open menu should list "Delete"
    When user closes the context menu
    And user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Selecting a role fills the context panel (Roles-10)
    When user types "BDD-RV-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-RV-Role-{time}" link in gallery
    Then the context panel should show "BDD-RV-Role-{time}"
    And the following elements should be visible:
      | "Actions" accordion header in context panel            |
      | "Assigned to" accordion header in context panel        |
      | "Favorites" accordion header in context panel          |
      | "Global Permissions" accordion header in context panel |
      | "Permissions" accordion header in context panel        |
      | "Sticky meta" accordion header in context panel        |
    And MANAGE button in "Assigned to" section in context panel should be visible
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
