@journey @serial @users @realizes:views.users
Feature: Managing a user
  What an administrator does to an existing user from its context menu in the Users view: the
  groups it belongs to, the roles it holds, disabling and enabling it, and marking it a favorite.
  Translated from files/TestTrack/User groups/users_manual_tests.md (Users-18 to Users-21) and
  playwright-public/user groups/users.test.ts.

  The user is the bddmanaged fixture user, shared with the other @serial features that change a
  user, which never run at the same time. The group and the role it joins are made by the feature
  and deleted at its end, and its memberships go with them; the fixture step puts it back to active
  and out of the favorites when an earlier run died between the halves.

  The context menu of a disabled user offers Enable at once: the gallery card that stayed stale
  until a reload in the older suite follows the change now, and the claim that Disable... is gone
  holds it to that.

  Every membership is claimed on the server, which is the fact, and a new one also in the pane that
  shows it, which is what the manual case reads. A removal is claimed on the server only: the pane
  empties and refills while it reloads, and an absence read off it would prove nothing.

  Background:
    Given user is logged in
    And a user "bddmanaged" is on the server
    And a group named "BDD-UM-Group-{time}" is on the server
    And no role named "BDD-UM-Role-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree

  Scenario: A role to assign
    When user clicks on "Platform > Roles" tree node inside browse tree
    Then the "Roles" view should be current
    When user clicks on "New Role..." button
    And user types "BDD-UM-Role-{time}" into Name input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And 1 role named "BDD-UM-Role-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Groups... adds the user to a group (Users-18)
    When user clicks on "Platform > Users" tree node inside browse tree
    Then the "Users" view should be current
    When user remembers the gallery counter
    And user types "bddmanaged" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Groups..." from the context menu of "bddmanaged" link in gallery
    Then "bddmanaged groups" dialog should be visible
    When user types "BDD-UM-Group-{time}" into membership search
    And user clicks on add button of "BDD-UM-Group-{time}" membership candidate
    Then "BDD-UM-Group-{time}" membership row should be visible
    When user clicks on SAVE button in "bddmanaged groups" dialog
    Then the "bddmanaged groups" dialog should close
    And "bddmanaged" should be a member of "BDD-UM-Group-{time}" on the server
    When user clicks on "bddmanaged" link in gallery
    Then the context panel should show "bddmanaged"
    When user expands "Member of" section in context panel
    Then "Member of" section in context panel should contain text "BDD-UM-Group-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Groups... takes the user out of the group again (Users-18)
    When user picks "Groups..." from the context menu of "bddmanaged" link in gallery
    Then "BDD-UM-Group-{time}" membership row should be visible
    When user clicks on remove button of "BDD-UM-Group-{time}" membership row
    Then "BDD-UM-Group-{time}" membership row should be absent
    When user clicks on SAVE button in "bddmanaged groups" dialog
    Then the "bddmanaged groups" dialog should close
    And "bddmanaged" should not be a member of "BDD-UM-Group-{time}" on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Roles... gives the user a role (Users-19)
    When user picks "Roles..." from the context menu of "bddmanaged" link in gallery
    Then "bddmanaged roles" dialog should be visible
    When user types "BDD-UM-Role-{time}" into membership search
    And user clicks on add button of "BDD-UM-Role-{time}" membership candidate
    Then "BDD-UM-Role-{time}" membership row should be visible
    When user clicks on SAVE button in "bddmanaged roles" dialog
    Then the "bddmanaged roles" dialog should close
    And "bddmanaged" should be a member of "BDD-UM-Role-{time}" on the server
    When user clicks on "bddmanaged" link in gallery
    And user expands "Roles" section in context panel
    Then "Roles" section in context panel should contain text "BDD-UM-Role-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Roles... takes the role away again (Users-19)
    When user picks "Roles..." from the context menu of "bddmanaged" link in gallery
    Then "BDD-UM-Role-{time}" membership row should be visible
    When user clicks on remove button of "BDD-UM-Role-{time}" membership row
    Then "BDD-UM-Role-{time}" membership row should be absent
    When user clicks on SAVE button in "bddmanaged roles" dialog
    Then the "bddmanaged roles" dialog should close
    And "bddmanaged" should not be a member of "BDD-UM-Role-{time}" on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Disable... disables the user (Users-20)
    Then the user "bddmanaged" should be active on the server
    When user picks "Disable..." from the context menu of "bddmanaged" link in gallery
    Then "Disable user" dialog should be visible
    When user clicks on DISABLE button in "Disable user" dialog
    Then the "Disable user" dialog should close
    And the user "bddmanaged" should be disabled on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A disabled user can be enabled again (Users-20)
    When user opens the context menu of "bddmanaged" link in gallery
    Then the open menu should list "Enable"
    And the open menu should not list "Disable..."
    When user picks "Enable" from the open menu
    Then "Enable user" dialog should be visible
    When user clicks on ENABLE button in "Enable user" dialog
    Then the "Enable user" dialog should close
    And the user "bddmanaged" should be active on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A user is added to favorites and removed again (Users-21)
    Then "My stuff > Favorites > bddmanaged" tree node inside browse tree should be absent
    When user picks "Add to favorites" from the context menu of "bddmanaged" link in gallery
    Then "My stuff > Favorites > bddmanaged" tree node inside browse tree should be present
    When user opens the context menu of "bddmanaged" link in gallery
    Then the open menu should list "Remove from favorites"
    When user picks "Remove from favorites" from the open menu
    Then "My stuff > Favorites > bddmanaged" tree node inside browse tree should be absent
    When user opens the context menu of "bddmanaged" link in gallery
    Then the open menu should list "Add to favorites"
    When user closes the context menu
    And user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
