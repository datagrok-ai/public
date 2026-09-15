@journey @serial @groups @realizes:views.groups
Feature: A group from creation to deletion
  The New Group... dialog, a group created with it, renamed through Properties... and deleted.
  Translated from files/TestTrack/User groups/groups_manual_tests.md (Groups-03, 04, 05, 09, 14) and
  playwright-public/user groups/groups.test.ts.

  Every step that changes a group is claimed on the server as well as in the list, and the list
  through a search for its name, whose counter must drop below the list's first before the link is
  read (the search is fuzzy, and a card can be on the page before the result lands). The groups are removed at the feature's end whatever it got to.

  Groups-05 in the manual case says the dialog has no name validation. It has: the Name field starts
  as "New Group", and cleared it is marked invalid and OK is disabled (grok_group_meta.dart
  propertiesDialog). The feature claims what the dialog does.

  Background:
    Given user is logged in
    And no group named "BDD-GL-Group-{time}, BDD-GL-Renamed-{time}, BDD-GL-Cancelled-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Groups" tree node inside browse tree

  Scenario: The New Group dialog, cancelled (Groups-03)
    Then the "Groups" view should be current
    When user remembers the gallery counter
    And user clicks on "New Group..." button
    Then "Create New Group" dialog should be visible
    And the following elements should be visible:
      | Name input in "Create New Group" dialog        |
      | Description input in "Create New Group" dialog |
    When user types "BDD-GL-Cancelled-{time}" into Name input in "Create New Group" dialog
    And user clicks on CANCEL button in "Create New Group" dialog
    Then the "Create New Group" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The New Group dialog refuses an empty name (Groups-05)
    When user clicks on "New Group..." button
    Then Name input in "Create New Group" dialog should have value "New Group"
    And OK button in "Create New Group" dialog should be enabled
    When user clears Name input in "Create New Group" dialog
    Then Name input in "Create New Group" dialog should be invalid
    And OK button in "Create New Group" dialog should be disabled
    When user clicks on CANCEL button in "Create New Group" dialog
    Then the "Create New Group" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A group is created (Groups-04)
    When user clicks on "New Group..." button
    And user types "BDD-GL-Group-{time}" into Name input in "Create New Group" dialog
    And user types "created by the BDD suite" into Description input in "Create New Group" dialog
    And user clicks on OK button in "Create New Group" dialog
    Then the "Create New Group" dialog should close
    And 1 group named "BDD-GL-Group-{time}" should be on the server
    And 0 groups named "BDD-GL-Cancelled-{time}" should be on the server
    When user types "BDD-GL-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD-GL-Group-{time}" link in gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Properties... renames the group and rewrites its description (Groups-09)
    When user picks "Properties..." from the context menu of "BDD-GL-Group-{time}" link in gallery
    Then "BDD-GL-Group-{time} Properties" dialog should be visible
    And Name input in "BDD-GL-Group-{time} Properties" dialog should have value "BDD-GL-Group-{time}"
    And Description input in "BDD-GL-Group-{time} Properties" dialog should have value "created by the BDD suite"
    When user types "BDD-GL-Renamed-{time}" into Name input in "BDD-GL-Group-{time} Properties" dialog
    And user types "renamed by the BDD suite" into Description input in "BDD-GL-Group-{time} Properties" dialog
    And user clicks on OK button in "BDD-GL-Group-{time} Properties" dialog
    Then the "BDD-GL-Group-{time} Properties" dialog should close
    And 1 group named "BDD-GL-Renamed-{time}" should be on the server
    And 0 groups named "BDD-GL-Group-{time}" should be on the server
    When user types "BDD-GL-Renamed-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Properties..." from the context menu of "BDD-GL-Renamed-{time}" link in gallery
    Then "BDD-GL-Renamed-{time} Properties" dialog should be visible
    And Description input in "BDD-GL-Renamed-{time} Properties" dialog should have value "renamed by the BDD suite"
    When user clicks on CANCEL button in "BDD-GL-Renamed-{time} Properties" dialog
    Then the "BDD-GL-Renamed-{time} Properties" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete removes the group after confirmation (Groups-14)
    When user picks "Delete" from the context menu of "BDD-GL-Renamed-{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 groups named "BDD-GL-Renamed-{time}" should be on the server
    When user clears gallery search
    And user types "BDD-GL-Renamed-{time}" into gallery search
    Then "BDD-GL-Renamed-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
