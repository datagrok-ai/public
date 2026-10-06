@journey @users @realizes:views.users
Feature: Creating users
  The New menu of the Users view: the new user dialog and the profile its OK opens, the new service
  user dialog, and a user found in the view. Translated from files/TestTrack/User groups/
  users_manual_tests.md (Users-05, Users-07) and playwright-public/user groups/users.test.ts.

  A user can never be deleted, so nothing here saves one: a saved user per run would stay on the
  stand for good. OK in the new user dialog does not save anyone: it opens the new user's profile,
  and the user exists only once the profile is saved (cmdUsersAddNew in users_browser.dart). The
  feature claims the profile opens with no such user on the server, and that closing it unsaved
  leaves none. The user the Users view finds is a fixture made once per stand, "bddcreated". The
  new service user dialog is claimed up to its OK, which a login enables; OK would save the service
  user and show its API token, so the dialog is cancelled.

  Not translated, and why: saving a user and a service user, and the API token dialog — each run
  would leave a user behind.

  Background:
    Given user is logged in
    And a user "bddcreated" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Users" tree node inside browse tree

  Scenario: OK in the new user dialog opens the new user's profile and saves nobody (Users-05)
    Then the "Users" view should be current
    When user clicks on New button
    And user picks "User..." from the open menu
    And user types "bddunsaved@datagrok.ai" into Email input in "Create new user" dialog
    And user types "bddunsaved" into Login input in "Create new user" dialog
    And user types "BDD" into "First Name" input in "Create new user" dialog
    And user types "Unsaved" into "Last Name" input in "Create new user" dialog
    Then OK button in "Create new user" dialog should be enabled
    When user clicks on OK button in "Create new user" dialog
    Then the "Create new user" dialog should close
    And the "BDD Unsaved" view should be current
    And 0 users with login "bddunsaved" should be on the server
    When user closes "BDD Unsaved" view
    Then the "Users" view should be current
    And 0 users with login "bddunsaved" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A user is found in the Users view (Users-05)
    When user remembers the gallery counter
    And user types "bddcreated" into gallery search
    Then the gallery counter should be lower than remembered
    And "bddcreated" link in gallery should be visible
    When user remembers the gallery counter
    And user clears gallery search
    Then the gallery counter should be higher than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The new service user dialog asks for a login before OK (Users-07)
    When user clicks on New button
    And user picks "Service User..." from the open menu
    Then "Create new service user" dialog should be visible
    And OK button in "Create new service user" dialog should be disabled
    When user types "bdd-svc-unsaved" into Login input in "Create new service user" dialog
    Then OK button in "Create new service user" dialog should be enabled
    When user clicks on CANCEL button in "Create new service user" dialog
    Then the "Create new service user" dialog should close
    And 0 users with login "bdd-svc-unsaved" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown
