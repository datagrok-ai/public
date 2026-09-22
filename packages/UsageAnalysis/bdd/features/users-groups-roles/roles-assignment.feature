@journey @serial @roles @realizes:views.roles
Feature: Who holds a role, and what it grants
  The Assigned to pane of a role and the editor behind its MANAGE button — assigning a user, letting
  it assign the role on, taking the role away — and the role's Global Permissions. Translated from
  files/TestTrack/User groups/roles_manual_tests.md (Roles-11 to Roles-14) and
  playwright-public/user groups/roles.test.ts.

  The role is made in the New Role dialog (the JS API cannot make a role) and deleted at the end of
  the feature; the user assigned is the bddmanaged fixture user, whose assignment goes with the role.
  An assignment is claimed on the server, and an addition also in the pane; a removal on the server
  only, since the pane empties while it reloads. A global permission has no reading in the JS API,
  so it is claimed in the MANAGE dialog, reopened after SAVE: the pane cannot tell the role's own
  grant from the ones it lists of other groups (GROK-20902). The cleanup revokes the grant before it
  deletes the role.

  Background:
    Given user is logged in
    And a user "bddmanaged" is on the server
    And no role named "BDD-RA-Role-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Roles" tree node inside browse tree

  Scenario: A role to assign
    Then the "Roles" view should be current
    When user remembers the gallery counter
    And user clicks on "New Role..." button
    And user types "BDD-RA-Role-{time}" into Name input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And 1 role named "BDD-RA-Role-{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be higher than remembered
    When user remembers the gallery counter
    And user types "BDD-RA-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-RA-Role-{time}" link in gallery
    Then the context panel should show "BDD-RA-Role-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: MANAGE assigns the role to a user (Roles-11)
    When user clicks on MANAGE button in "Assigned to" section in context panel
    Then "BDD-RA-Role-{time} members" dialog should be visible
    When user types "bddmanaged" into membership search
    And user clicks on add button of "bddmanaged" membership candidate
    Then "bddmanaged" membership row should be visible
    And checkbox label of "bddmanaged" membership row should have text "Can assign"
    When user clicks on SAVE button in "BDD-RA-Role-{time} members" dialog
    Then the "BDD-RA-Role-{time} members" dialog should close
    And "bddmanaged" should be a plain member of "BDD-RA-Role-{time}" on the server
    And "Assigned to" section in context panel should contain text "bddmanaged"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Can assign is saved for the assignee (Roles-13)
    When user clicks on MANAGE button in "Assigned to" section in context panel
    Then checkbox of "bddmanaged" membership row should be unchecked
    When user checks checkbox of "bddmanaged" membership row
    And user clicks on SAVE button in "BDD-RA-Role-{time} members" dialog
    Then the "BDD-RA-Role-{time} members" dialog should close
    And "bddmanaged" should be an admin member of "BDD-RA-Role-{time}" on the server
    When user clicks on MANAGE button in "Assigned to" section in context panel
    Then checkbox of "bddmanaged" membership row should be checked
    When user clicks on CANCEL button in "BDD-RA-Role-{time} members" dialog
    Then the "BDD-RA-Role-{time} members" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Removing the assignment takes the role away (Roles-12)
    When user clicks on MANAGE button in "Assigned to" section in context panel
    Then "bddmanaged" membership row should be visible
    When user clicks on remove button of "bddmanaged" membership row
    Then "bddmanaged" membership row should be absent
    When user clicks on SAVE button in "BDD-RA-Role-{time} members" dialog
    Then the "BDD-RA-Role-{time} members" dialog should close
    And "bddmanaged" should not be a member of "BDD-RA-Role-{time}" on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # GROK-20902, fixed 2026-09-21 in privileges_service.dart: the global permissions of a group
  # were every group's ("via Administrators", "via All users"); they are the group's own grants.
  Scenario: A new role has no global permissions in its pane (Roles-14)
    When user expands "Global Permissions" section in context panel
    Then "Global Permissions" section in context panel should contain text "No global permissions"

  Scenario: Global Permissions grants the role a permission (Roles-14)
    When user expands "Global Permissions" section in context panel
    And user clicks on MANAGE button in "Global Permissions" section in context panel
    Then "BDD-RA-Role-{time}: Global Permissions" dialog should be visible
    When user expands "Browse" tree node in "BDD-RA-Role-{time}: Global Permissions" dialog
    Then "Browse > Browse Apps" tree node in "BDD-RA-Role-{time}: Global Permissions" dialog should be unchecked
    When user checks "Browse > Browse Apps" tree node in "BDD-RA-Role-{time}: Global Permissions" dialog
    And user clicks on SAVE button in "BDD-RA-Role-{time}: Global Permissions" dialog
    Then the "BDD-RA-Role-{time}: Global Permissions" dialog should close
    When user clicks on MANAGE button in "Global Permissions" section in context panel
    And user expands "Browse" tree node in "BDD-RA-Role-{time}: Global Permissions" dialog
    Then "Browse > Browse Apps" tree node in "BDD-RA-Role-{time}: Global Permissions" dialog should be checked
    When user clicks on CANCEL button in "BDD-RA-Role-{time}: Global Permissions" dialog
    Then the "BDD-RA-Role-{time}: Global Permissions" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # GROK-20904, fixed 2026-09-22 by db_up/20260922_0_permissions_group_cascade.sql: deleting a role
  # that held a global permission violated permissions_user_group_id_fkey and the role stayed; the
  # key now cascades, so a group's grants go with it. This scenario needs the grant of the one before
  # it. The feature's own cleanup still revokes the role's global permissions before it deletes the role.
  Scenario: A role that grants a permission can still be deleted (Roles-14, Roles-15)
    When user types "BDD-RA-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Delete" from the context menu of "BDD-RA-Role-{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 roles named "BDD-RA-Role-{time}" should be on the server
    And no errors should have been logged
