@serial @roles @edge-cases @realizes:views.roles
Feature: Role edge cases: a name a group already has, a name of spaces only, a cancelled grant
  Guards for the New Role dialog and a role's Global Permissions. Groups and roles share one set of
  names: a role named in the New Role dialog after a group made in the New Group dialog is refused
  ("A group or role with the same name already exists", `groups_repository.dart` save) and the
  dialog stays open with the message in an error balloon (`grok_group_meta.dart` `propertiesDialog`).
  The dialog refuses a name of spaces only ("Name must be at least 1 character", the same
  `propertiesDialog`). CANCEL in the Global Permissions dialog grants nothing: the box ticked before
  it is unticked when the dialog is opened again.

  The group and the role are made in their dialogs (the JS API cannot make a role) and deleted at
  feature end with any global permission they hold. A global permission has no reading in the JS
  API, so it is claimed in the dialog reopened, as roles-assignment.feature does. The Global
  Permissions pane after a role holding a grant is deleted through the API (core#561) is not here:
  see MISSING.md.

  Background:
    Given user is logged in
    And no group named "BDD-RE-Group-{time}, BDD-RE-Role-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree

  Scenario: A new role under the name of a group is refused
    When user clicks on "Platform > Groups" tree node inside browse tree
    Then the "Groups" view should be current
    When user clicks on "New Group..." button
    And user types "BDD-RE-Group-{time}" into Name input in "Create New Group" dialog
    And user clicks on OK button in "Create New Group" dialog
    Then the "Create New Group" dialog should close
    And 1 group named "BDD-RE-Group-{time}" should be on the server
    When user clicks on "Platform > Roles" tree node inside browse tree
    Then the "Roles" view should be current
    When user clicks on "New Role..." button
    And user types "BDD-RE-Group-{time}" into Name input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then an error balloon containing "A group or role with the same name already exists" should have been shown
    And "Create New Role" dialog should be visible
    And 1 group named "BDD-RE-Group-{time}" should be on the server
    When user clicks on CANCEL button in "Create New Role" dialog
    Then the "Create New Role" dialog should close

  Scenario: The New Role dialog refuses a name of spaces only
    When user clicks on "Platform > Roles" tree node inside browse tree
    Then the "Roles" view should be current
    When user clicks on "New Role..." button
    And user types "   " into Name input in "Create New Role" dialog
    Then Name input in "Create New Role" dialog should be invalid
    And OK button in "Create New Role" dialog should be disabled
    When user clicks on CANCEL button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: CANCEL in Global Permissions grants nothing
    When user clicks on "Platform > Roles" tree node inside browse tree
    Then the "Roles" view should be current
    When user remembers the gallery counter
    And user clicks on "New Role..." button
    And user types "BDD-RE-Role-{time}" into Name input in "Create New Role" dialog
    And user clicks on OK button in "Create New Role" dialog
    Then the "Create New Role" dialog should close
    And 1 role named "BDD-RE-Role-{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be higher than remembered
    When user remembers the gallery counter
    And user types "BDD-RE-Role-{time}" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "BDD-RE-Role-{time}" link in gallery
    Then the context panel should show "BDD-RE-Role-{time}"
    When user expands "Global Permissions" section in context panel
    And user clicks on MANAGE button in "Global Permissions" section in context panel
    Then "BDD-RE-Role-{time}: Global Permissions" dialog should be visible
    When user expands "Browse" tree node in "BDD-RE-Role-{time}: Global Permissions" dialog
    And user checks "Browse > Browse Apps" tree node in "BDD-RE-Role-{time}: Global Permissions" dialog
    Then "Browse > Browse Apps" tree node in "BDD-RE-Role-{time}: Global Permissions" dialog should be checked
    When user clicks on CANCEL button in "BDD-RE-Role-{time}: Global Permissions" dialog
    Then the "BDD-RE-Role-{time}: Global Permissions" dialog should close
    When user clicks on MANAGE button in "Global Permissions" section in context panel
    And user expands "Browse" tree node in "BDD-RE-Role-{time}: Global Permissions" dialog
    Then "Browse > Browse Apps" tree node in "BDD-RE-Role-{time}: Global Permissions" dialog should be unchecked
    When user clicks on CANCEL button in "BDD-RE-Role-{time}: Global Permissions" dialog
    Then the "BDD-RE-Role-{time}: Global Permissions" dialog should close
    And "Global Permissions" section in context panel should contain text "No global permissions"
    And no errors should have been logged
    And no error or warning balloon should have been shown
