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
  so the user's login is found nowhere in the gallery. That the personal group is made is tested in
  datlas users_groups/users_test.dart 'Save user'.

  Edge cases: CANCEL in MANAGE leaves the members as they were. The account that creates a group
  becomes its only admin (`groups_repository.dart` save), and "if user is the only admin in group,
  they can't leave the group or revoke their admin privileges" (help
  govern/access-control/users-and-groups.md): unticking its Admin box or removing itself keeps the
  MANAGE dialog open with the server's message in an error balloon, and nothing changes. These run on
  a group made through the API before each of them (the Background of their Rule), so the running
  account is its only admin and bddmanaged is not a member. The running account is named by its login
  in the membership editor, which lists a user by name: on a stand where the two differ those
  scenarios cannot find its row. A group added as a member of its own member is refused as a circular
  reference: two groups of the scenario's own, made through the API, one nested into the other through
  MANAGE first. Every group these scenarios make is named with the run's time and deleted at the end.

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

  Scenario: A user's personal group is not in the Groups view (Groups-18)
    When user types "bddmanaged" into gallery search
    Then gallery counter should have text "0"
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  @edge-cases
  Scenario: A group added as a member of its own member is refused as a circular reference
    Given a group named "BDD-GM-Outer-{time}" is on the server
    And a group named "BDD-GM-Inner-{time}" is on the server
    When user clears gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user remembers the gallery counter
    And user types "BDD-GM-Outer-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-GM-Outer-{time}" link in gallery
    Then the context panel should show "BDD-GM-Outer-{time}"
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GM-Outer-{time} members" dialog should be visible
    When user types "BDD-GM-Inner-{time}" into membership search
    And user clicks on add button of "BDD-GM-Inner-{time}" membership candidate
    Then "BDD-GM-Inner-{time}" membership row should be visible
    When user clicks on SAVE button in "BDD-GM-Outer-{time} members" dialog
    Then the "BDD-GM-Outer-{time} members" dialog should close
    And "BDD-GM-Inner-{time}" should be a member of "BDD-GM-Outer-{time}" on the server
    When user clears gallery search
    And user types "BDD-GM-Inner-{time}" into gallery search
    And user clicks on "BDD-GM-Inner-{time}" link in gallery
    Then the context panel should show "BDD-GM-Inner-{time}"
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GM-Inner-{time} members" dialog should be visible
    When user types "BDD-GM-Outer-{time}" into membership search
    And user clicks on add button of "BDD-GM-Outer-{time}" membership candidate
    Then "BDD-GM-Outer-{time}" membership row should be visible
    When user clicks on SAVE button in "BDD-GM-Inner-{time} members" dialog
    Then an error balloon containing "Circular reference detected" should have been shown
    And "BDD-GM-Outer-{time}" should not be a member of "BDD-GM-Inner-{time}" on the server
    And "BDD-GM-Inner-{time}" should be a member of "BDD-GM-Outer-{time}" on the server
    When user clicks on CANCEL button in "BDD-GM-Inner-{time} members" dialog
    Then the "BDD-GM-Inner-{time} members" dialog should close
    When user clears gallery search

  Rule: MANAGE on a group of the scenario's own, whose only admin is the running account
    Background:
      Given a group named "BDD-GM-Edge-{time}" is on the server
      When user clears gallery search
      And user clicks on "Refresh" icon inside gallery toolbar
      And user remembers the gallery counter
      And user types "BDD-GM-Edge-{time}" into gallery search
      Then the gallery counter should be lower than remembered
      When user clicks on "BDD-GM-Edge-{time}" link in gallery
      Then the context panel should show "BDD-GM-Edge-{time}"

    @edge-cases
    Scenario: CANCEL in MANAGE leaves the members as they were
      When user clicks on MANAGE button in "Members" section in context panel
      Then "BDD-GM-Edge-{time} members" dialog should be visible
      When user types "bddmanaged" into membership search
      And user clicks on add button of "bddmanaged" membership candidate
      Then "bddmanaged" membership row should be visible
      When user clicks on CANCEL button in "BDD-GM-Edge-{time} members" dialog
      Then the "BDD-GM-Edge-{time} members" dialog should close
      And "bddmanaged" should not be a member of "BDD-GM-Edge-{time}" on the server
      When user clears gallery search
      Then no errors should have been logged
      And no error or warning balloon should have been shown

    @edge-cases
    Scenario Outline: The only admin <change> is refused and stays the admin
      Then "{login}" should be an admin member of "BDD-GM-Edge-{time}" on the server
      When user clicks on MANAGE button in "Members" section in context panel
      Then "BDD-GM-Edge-{time} members" dialog should be visible
      When user <gesture> of "{login}" membership row
      And user clicks on SAVE button in "BDD-GM-Edge-{time} members" dialog
      Then an error balloon containing "should have at least one admin user" should have been shown
      And "BDD-GM-Edge-{time} members" dialog should be visible
      And "{login}" should be an admin member of "BDD-GM-Edge-{time}" on the server
      When user clicks on CANCEL button in "BDD-GM-Edge-{time} members" dialog
      Then the "BDD-GM-Edge-{time} members" dialog should close
      When user clears gallery search

      Examples:
        | change                        | gesture                      |
        | unticking their Admin box     | unchecks checkbox            |
        | removing themselves           | clicks on remove button      |
