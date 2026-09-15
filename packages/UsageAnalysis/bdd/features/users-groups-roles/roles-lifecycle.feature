@journey @serial @roles @realizes:views.roles
Feature: A role from creation to deletion
  The New Role... dialog, a role created with it, renamed through Properties... and deleted.
  Translated from files/TestTrack/User groups/roles_manual_tests.md (Roles-03, 04, 05, 09, 15) and
  playwright-public/user groups/roles.test.ts.

  Every change is claimed on the server as well as in the list, and the list through a search for
  its name, whose counter must drop below the list's first before the link is read (the search is
  fuzzy, and a card can be on the page before the result lands). The roles are removed at the
  feature's end whatever it got to.

  The Name field of the New Role dialog starts as "New Role"; Roles-05 is claimed on the field
  cleared, where it is marked invalid and OK is disabled.

  Background:
    Given user is logged in
    And no role named "BDD-RL-Role-{time}, BDD-RL-Renamed-{time}, BDD-RL-Cancelled-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Roles" tree node inside browse tree

  Scenario: The New Role dialog, cancelled (Roles-03)
    Then the "Roles" view should be current
    When user remembers the gallery counter
    And user clicks on "New Role..." button
    Then "Create New Role" dialog should be visible
    And the following elements should be visible:
      | Name input in "Create New Role" dialog        |
      | Description input in "Create New Role" dialog |
    When user types "BDD-RL-Cancelled-{time}" into Name input in "Create New Role" dialog
    And user clicks on CANCEL button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The New Role dialog refuses an empty name (Roles-05)
    When user clicks on "New Role..." button
    Then Name input in "Create New Role" dialog should have value "New Role"
    And OK button in "Create New Role" dialog should be enabled
    When user clears Name input in "Create New Role" dialog
    Then Name input in "Create New Role" dialog should be invalid
    And OK button in "Create New Role" dialog should be disabled
    When user clicks on CANCEL button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A role is created (Roles-04)
    When user clicks on "New Role..." button
    And user types "BDD-RL-Role-{time}" into Name input in "Create New Role" dialog
    And user types "created by the BDD suite" into Description input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And 1 role named "BDD-RL-Role-{time}" should be on the server
    And 0 roles named "BDD-RL-Cancelled-{time}" should be on the server
    When user types "BDD-RL-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD-RL-Role-{time}" link in gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Properties... renames the role and rewrites its description (Roles-09)
    When user picks "Properties..." from the context menu of "BDD-RL-Role-{time}" link in gallery
    Then "BDD-RL-Role-{time} Properties" dialog should be visible
    And Name input in "BDD-RL-Role-{time} Properties" dialog should have value "BDD-RL-Role-{time}"
    And Description input in "BDD-RL-Role-{time} Properties" dialog should have value "created by the BDD suite"
    When user types "BDD-RL-Renamed-{time}" into Name input in "BDD-RL-Role-{time} Properties" dialog
    And user types "renamed by the BDD suite" into Description input in "BDD-RL-Role-{time} Properties" dialog
    And user clicks on OK button in "BDD-RL-Role-{time} Properties" dialog
    Then the "BDD-RL-Role-{time} Properties" dialog should close
    And 1 role named "BDD-RL-Renamed-{time}" should be on the server
    And 0 roles named "BDD-RL-Role-{time}" should be on the server
    When user types "BDD-RL-Renamed-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Properties..." from the context menu of "BDD-RL-Renamed-{time}" link in gallery
    Then "BDD-RL-Renamed-{time} Properties" dialog should be visible
    And Description input in "BDD-RL-Renamed-{time} Properties" dialog should have value "renamed by the BDD suite"
    When user clicks on CANCEL button in "BDD-RL-Renamed-{time} Properties" dialog
    Then the "BDD-RL-Renamed-{time} Properties" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete removes the role after confirmation (Roles-15)
    When user picks "Delete" from the context menu of "BDD-RL-Renamed-{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 roles named "BDD-RL-Renamed-{time}" should be on the server
    When user clears gallery search
    And user types "BDD-RL-Renamed-{time}" into gallery search
    Then "BDD-RL-Renamed-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
