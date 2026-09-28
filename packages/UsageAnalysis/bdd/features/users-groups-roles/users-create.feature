@journey @users @realizes:views.users
Feature: Creating users
  A user and a service user made through the New menu of the Users view, and found again. Translated
  from files/TestTrack/User groups/users_manual_tests.md (Users-05, Users-07) and
  playwright-public/user groups/users.test.ts.

  Every run makes two users named by the time it ran, and they stay: a user can never be deleted.
  On a fresh CI database that is two users per run; on a shared stand it adds up, which is the
  price of checking the creation itself rather than a dialog that was cancelled.

  OK in the new user dialog does not save anyone: it opens the new user's profile, and the user
  exists only once the profile is saved (cmdUsersAddNew in users_browser.dart). The feature claims
  both halves: the profile is open and no such user exists, then one does after Save.

  A service user is saved by OK, which then shows its API token dialog (UserMeta.session) — the only
  way this dialog opens, and so the claim that what was made is a service account: the JS API does
  not expose the service flag.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Users" tree node inside browse tree

  Scenario: A new user is saved from the profile that OK opens (Users-05)
    Then the "Users" view should be current
    When user clicks on New button
    And user picks "User..." from the open menu
    And user types "bdd{time}@datagrok.ai" into Email input in "Create new user" dialog
    And user types "bdd{time}" into Login input in "Create new user" dialog
    And user types "BDD" into "First Name" input in "Create new user" dialog
    And user types "User {time}" into "Last Name" input in "Create new user" dialog
    Then OK button in "Create new user" dialog should be enabled
    When user clicks on OK button in "Create new user" dialog
    Then the "Create new user" dialog should close
    And the "BDD User {time}" view should be current
    And 0 users with login "bdd{time}" should be on the server
    When user clicks on Save button
    Then 1 user with login "bdd{time}" should be on the server
    And the user "bdd{time}" should be active on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The new user is found in the Users view (Users-05)
    When user switches to the "Users" view
    And user remembers the gallery counter
    And user types "bdd{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD User {time}" link in gallery should be visible
    When user remembers the gallery counter
    And user clears gallery search
    Then the gallery counter should be higher than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A service user is saved by OK and shown its token (Users-07)
    When user clicks on New button
    And user picks "Service User..." from the open menu
    Then "Create new service user" dialog should be visible
    And OK button in "Create new service user" dialog should be disabled
    When user types "bdd-svc{time}" into Login input in "Create new service user" dialog
    Then OK button in "Create new service user" dialog should be enabled
    When user clicks on OK button in "Create new service user" dialog
    Then the "Create new service user" dialog should close
    And "API token" dialog should be visible
    And "API token" input in "API token" dialog should be visible
    And 1 user with login "bdd-svc{time}" should be on the server
    When user clicks on CLOSE button in "API token" dialog
    Then the "API token" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The service user is found in the Users view (Users-07)
    When user remembers the gallery counter
    And user types "bdd-svc{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "bdd-svc{time}" link in gallery should be visible
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
