@serial @groups @edge-cases @realizes:views.groups
Feature: Group edge cases: the only admin, a group nested into its own member, a taken name, cancelled edits
  Guards for a group's editors under Browse > Platform > Groups. The account that creates a group
  becomes its only admin (`groups_repository.dart` save), and "if user is the only admin in group,
  they can't leave the group or revoke their admin privileges" (help
  govern/access-control/users-and-groups.md): MANAGE keeps its dialog open with the server's message
  in an error balloon (`grok_group_meta.dart` `_editMembers`) and nothing changes. A group added as a
  member of its own member is refused as a circular reference, and a second group made in the New
  Group dialog under a name already made there is refused (`groups_repository.dart`: "A group or role
  with the same name already exists"; the guard compares the internal names, which the dialog derives
  from the name typed, so a group named the same through the API is not one: MISSING.md). The Properties dialog refuses a name of
  spaces only ("Name must be at least 1 character", `grok_group_meta.dart` `propertiesDialog`).
  CANCEL in MANAGE and in Properties leaves the group as it was on the server.

  The groups are made through the API and deleted at feature end. The running account is named by
  its login in the membership editor, which lists a user by name: on a stand where the running
  account's name differs from its login the two only-admin scenarios cannot find its row. That the
  server refuses a group with no admin whatever route saves it is a datlas matter (MISSING.md).

  Background:
    Given user is logged in
    And a group named "BDD-GE-Group-{time}" is on the server
    And a group named "BDD-GE-Child-{time}" is on the server
    And no group named "BDD-GE-Twice-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Groups" tree node inside browse tree
    And user remembers the gallery counter
    And user types "BDD-GE-Group-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-GE-Group-{time}" link in gallery
    Then the context panel should show "BDD-GE-Group-{time}"

  Scenario: The only admin unticking their Admin box is refused and stays the admin
    Then "{login}" should be an admin member of "BDD-GE-Group-{time}" on the server
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GE-Group-{time} members" dialog should be visible
    And checkbox of "{login}" membership row should be checked
    When user unchecks checkbox of "{login}" membership row
    And user clicks on SAVE button in "BDD-GE-Group-{time} members" dialog
    Then an error balloon containing "should have at least one admin user" should have been shown
    And "BDD-GE-Group-{time} members" dialog should be visible
    And "{login}" should be an admin member of "BDD-GE-Group-{time}" on the server
    When user clicks on CANCEL button in "BDD-GE-Group-{time} members" dialog
    Then the "BDD-GE-Group-{time} members" dialog should close

  Scenario: The only admin removing themselves is refused and stays the admin
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GE-Group-{time} members" dialog should be visible
    When user clicks on remove button of "{login}" membership row
    Then "{login}" membership row should be absent
    When user clicks on SAVE button in "BDD-GE-Group-{time} members" dialog
    Then an error balloon containing "should have at least one admin user" should have been shown
    And "BDD-GE-Group-{time} members" dialog should be visible
    And "{login}" should be an admin member of "BDD-GE-Group-{time}" on the server
    When user clicks on CANCEL button in "BDD-GE-Group-{time} members" dialog
    Then the "BDD-GE-Group-{time} members" dialog should close

  Scenario: A group added as a member of its own member is refused as a circular reference
    When user clicks on MANAGE button in "Members" section in context panel
    And user types "BDD-GE-Child-{time}" into membership search
    And user clicks on add button of "BDD-GE-Child-{time}" membership candidate
    And user clicks on SAVE button in "BDD-GE-Group-{time} members" dialog
    Then the "BDD-GE-Group-{time} members" dialog should close
    And "BDD-GE-Child-{time}" should be a member of "BDD-GE-Group-{time}" on the server
    When user remembers the gallery counter
    And user clears gallery search
    And user types "BDD-GE-Child-{time}" into gallery search
    And user clicks on "BDD-GE-Child-{time}" link in gallery
    Then the context panel should show "BDD-GE-Child-{time}"
    When user clicks on MANAGE button in "Members" section in context panel
    Then "BDD-GE-Child-{time} members" dialog should be visible
    When user types "BDD-GE-Group-{time}" into membership search
    And user clicks on add button of "BDD-GE-Group-{time}" membership candidate
    Then "BDD-GE-Group-{time}" membership row should be visible
    When user clicks on SAVE button in "BDD-GE-Child-{time} members" dialog
    Then an error balloon containing "Circular reference detected" should have been shown
    And "BDD-GE-Group-{time}" should not be a member of "BDD-GE-Child-{time}" on the server
    And "BDD-GE-Child-{time}" should be a member of "BDD-GE-Group-{time}" on the server
    When user clicks on CANCEL button in "BDD-GE-Child-{time} members" dialog
    Then the "BDD-GE-Child-{time} members" dialog should close

  Scenario: CANCEL in MANAGE leaves the members as they were
    When user clicks on MANAGE button in "Members" section in context panel
    And user types "BDD-GE-Child-{time}" into membership search
    And user clicks on add button of "BDD-GE-Child-{time}" membership candidate
    Then "BDD-GE-Child-{time}" membership row should be visible
    When user clicks on CANCEL button in "BDD-GE-Group-{time} members" dialog
    Then the "BDD-GE-Group-{time} members" dialog should close
    And "BDD-GE-Child-{time}" should not be a member of "BDD-GE-Group-{time}" on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A second group made in the New Group dialog under the same name is refused
    When user clicks on "New Group..." button
    And user types "BDD-GE-Twice-{time}" into Name input in "Create New Group" dialog
    And user clicks on OK button in "Create New Group" dialog
    Then the "Create New Group" dialog should close
    And 1 group named "BDD-GE-Twice-{time}" should be on the server
    When user clicks on "New Group..." button
    And user types "BDD-GE-Twice-{time}" into Name input in "Create New Group" dialog
    And user clicks on OK button in "Create New Group" dialog
    Then an error balloon containing "A group or role with the same name already exists" should have been shown
    And "Create New Group" dialog should be visible
    And 1 group named "BDD-GE-Twice-{time}" should be on the server
    When user clicks on CANCEL button in "Create New Group" dialog
    Then the "Create New Group" dialog should close

  Scenario: Properties refuses a name of spaces only
    When user picks "Properties..." from the context menu of "BDD-GE-Group-{time}" link in gallery
    Then "BDD-GE-Group-{time} Properties" dialog should be visible
    When user types "   " into Name input in "BDD-GE-Group-{time} Properties" dialog
    Then Name input in "BDD-GE-Group-{time} Properties" dialog should be invalid
    And OK button in "BDD-GE-Group-{time} Properties" dialog should be disabled
    When user clicks on CANCEL button in "BDD-GE-Group-{time} Properties" dialog
    Then the "BDD-GE-Group-{time} Properties" dialog should close
    And 1 group named "BDD-GE-Group-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown
