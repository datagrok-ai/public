@journey @serial @groups @realizes:views.groups
Feature: A group's members
  The Members pane of a group and the editor behind its MANAGE button: adding a user, making it an
  admin, nesting a group, removing a member. Translated from files/TestTrack/User
  groups/groups_manual_tests.md (Groups-11, 12, 13, 15, 18) and playwright-public/user groups/groups.test.ts.

  The group and the group nested into it are made by the feature and deleted at its end; the user
  added is the bddmanaged fixture user, whose membership goes with the group. Each change is claimed on the server;
  an addition also in the pane. A removal is claimed on the server only: the pane empties and
  refills while it reloads, and an absence read off it would prove nothing.

  Every user has a personal security group, named by the login; the Groups view does not list those,
  so the user's own group is claimed on the server and its login found nowhere in the gallery.

  Background:
    Given user is logged in
    And a user "bddmanaged" is on the server
    And a group named "BDD-GM-Group-{time}" is on the server
    And a group named "BDD-GM-Child-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Groups" tree node inside browse tree

  Scenario: MANAGE adds a user to the group (Groups-11)
    When user remembers the gallery counter
    And user types "BDD-GM-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-GM-Group-{time}" link in gallery
    Then the context panel should show "BDD-GM-Group-{time}"
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GM-Group-{time} members" dialog should be visible
    When user types "bddmanaged" into membership search
    And user clicks on add button of "bddmanaged" membership candidate
    Then "bddmanaged" membership row should be visible
    When user clicks on SAVE button in "BDD-GM-Group-{time} members" dialog
    Then the "BDD-GM-Group-{time} members" dialog should close
    And "bddmanaged" should be a plain member of "BDD-GM-Group-{time}" on the server
    And "Members" section in context panel should contain text "bddmanaged"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Admin box makes the member an admin (Groups-13)
    When user clicks on MANAGE button in "Members" section in context panel
    Then checkbox label of "bddmanaged" membership row should have text "Admin"
    And checkbox of "bddmanaged" membership row should be unchecked
    When user checks checkbox of "bddmanaged" membership row
    And user clicks on SAVE button in "BDD-GM-Group-{time} members" dialog
    Then the "BDD-GM-Group-{time} members" dialog should close
    And "bddmanaged" should be an admin member of "BDD-GM-Group-{time}" on the server
    When user clicks on MANAGE button in "Members" section in context panel
    Then checkbox of "bddmanaged" membership row should be checked
    When user clicks on CANCEL button in "BDD-GM-Group-{time} members" dialog
    Then the "BDD-GM-Group-{time} members" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A group added as a member is nested (Groups-15)
    When user clicks on MANAGE button in "Members" section in context panel
    And user types "BDD-GM-Child-{time}" into membership search
    And user clicks on add button of "BDD-GM-Child-{time}" membership candidate
    Then "BDD-GM-Child-{time}" membership row should be visible
    When user clicks on SAVE button in "BDD-GM-Group-{time} members" dialog
    Then the "BDD-GM-Group-{time} members" dialog should close
    And "BDD-GM-Child-{time}" should be a member of "BDD-GM-Group-{time}" on the server
    And "Members" section in context panel should contain text "BDD-GM-Child-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Removing a member leaves the others (Groups-12)
    When user clicks on MANAGE button in "Members" section in context panel
    Then "bddmanaged" membership row should be visible
    When user clicks on remove button of "bddmanaged" membership row
    Then "bddmanaged" membership row should be absent
    When user clicks on SAVE button in "BDD-GM-Group-{time} members" dialog
    Then the "BDD-GM-Group-{time} members" dialog should close
    And "bddmanaged" should not be a member of "BDD-GM-Group-{time}" on the server
    And "BDD-GM-Child-{time}" should be a member of "BDD-GM-Group-{time}" on the server
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A user's personal group is on the server but not in the Groups view (Groups-18)
    Then the user "bddmanaged" should have a personal group on the server
    When user types "bddmanaged" into gallery search
    Then gallery counter should have text "0"
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
