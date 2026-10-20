@realizes:sharing.share-dialog
Feature: A table saved in a project, shared from its Sharing pane
  A table the Save project dialog uploads is listed under its project in My stuff > My dashboards. Its
  own context menu there offers no Share...; clicked, it fills the context panel, whose Sharing pane
  names the owner and opens the Share dialog, which shares it as any entity is shared: the dialog starts
  at View and use, and the grant it makes shows when the dialog opens again. Translated from the owner
  side of the TestTrack case Sharing/share-table-permissions (playwright-public/Sharing
  share-table-permissions).

  The table is the first 20 rows of demog, named with the run's time and saved through the ribbon's Save
  in a project of the same name: a table saved through the JS API has no friendly name, and the tree
  lists it without a label. The project, with its table and view, the grant with them, is removed when
  the feature starts and ends. Neither the project's Share dialog the save opens (closed with Escape) nor
  the table's share sends a notification.

  Not translated, and why: what the second account may then do with the table (read its data, not delete
  or share it further) and the per-privilege matrix of the PermissionsView are the server's answers, an API
  test (ApiTests, Dapi: permissions); the dialog has no Advanced editor link any more.

  Background:
    Given user is logged in
    And no project named "BDD-ShareTable-{time}" is on the server
    And user opens demog dataset keeping the first 20 rows as "BDD-ShareTable-{time}"
    When user opens the Save project dialog from the ribbon
    And user types "BDD-ShareTable-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDD-ShareTable-{time}" should be on the server
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > My dashboards" tree node inside browse tree is expanded
    And "My stuff > My dashboards > BDD-ShareTable-{time}" tree node inside browse tree is expanded

  Scenario: Share... on the table shares it with another user
    Given the context panel is open
    When user clicks on "My stuff > My dashboards > BDD-ShareTable-{time} > BDD-ShareTable-{time}" tree node inside browse tree
    Then the context panel should show "BDD-ShareTable-{time}"
    When user expands Sharing section in context panel
    Then Sharing section in context panel should contain text "You are the owner"
    When user clicks on "Share..." button in Sharing section in context panel
    Then "Share BDD-ShareTable-{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-ShareTable-{time}" dialog
    And user unchecks "Send notifications" input in "Share BDD-ShareTable-{time}" dialog
    And user clicks on OK button in "Share BDD-ShareTable-{time}" dialog
    Then the "Share BDD-ShareTable-{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "Share..." button in Sharing section in context panel
    Then "Share BDD-ShareTable-{time}" dialog should be visible
    And the access level of the sharing user in "Share BDD-ShareTable-{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDD-ShareTable-{time}" dialog
    Then the "Share BDD-ShareTable-{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown
