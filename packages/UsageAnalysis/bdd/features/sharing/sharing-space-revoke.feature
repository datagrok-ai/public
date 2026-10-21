@journey @serial @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me @realizes:views.space
Feature: A space, shared and taken back
  A space the owner shares and then takes back from its Share dialog: the second account sees it
  under My stuff > Shared with me while it is shared, and not once the share is gone. Translated from
  playwright-public/Sharing share-files-and-spaces-permissions (blocks C and D, whose second-account
  steps were API reads); what the second account may do inside a shared space is sharing-recipient.

  The spaces are this feature's own root spaces, made through the API and deleted when the feature
  starts and ends, their grants with them: ShareRevokeSpace1, shared and taken back, and
  ShareRevokeSpace2, shared too and never taken back, so that the second account still has the
  owner's node under Shared with me when ShareRevokeSpace1's absence there is claimed. The second
  account is the library's sharing user; it signs in on the feature's page, and every share is made
  with Send notifications off.

  Not translated, and why: that a revoked space takes its content's access with it is the server's
  answer (ApiTests, Dapi: permissions, "a space shared for viewing"); that the second account may not
  edit, delete or share the space itself is an API test still to write (sharing/MISSING.md).

  Serial: it signs in as the second account on its page.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a space named "ShareRevokeSpace1" is on the server
    And a space named "ShareRevokeSpace2" is on the server
    And the context panel is open
    And the browse panel is open
    When user refreshes the browse tree
    Given Spaces tree node inside browse tree is expanded

  Scenario: The owner shares both spaces
    When user picks "Share..." from the context menu of Spaces---ShareRevokeSpace2 tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share ShareRevokeSpace2" dialog
    And user unchecks "Send notifications" input in "Share ShareRevokeSpace2" dialog
    And user clicks on OK button in "Share ShareRevokeSpace2" dialog
    Then the "Share ShareRevokeSpace2" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user picks "Share..." from the context menu of Spaces---ShareRevokeSpace1 tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share ShareRevokeSpace1" dialog
    And user unchecks "Send notifications" input in "Share ShareRevokeSpace1" dialog
    And user clicks on OK button in "Share ShareRevokeSpace1" dialog
    Then the "Share ShareRevokeSpace1" dialog should close
    And an info balloon containing "Shared" should have been shown
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account sees both spaces under Shared with me
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareRevokeSpace1" shared by the running account should be visible
    And "ShareRevokeSpace2" shared by the running account should be visible

  Scenario: The owner takes the share of one space back
    When user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given Spaces tree node inside browse tree is expanded
    When user clicks on Spaces---ShareRevokeSpace1 tree node inside browse tree
    Then the context panel should show "ShareRevokeSpace1"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of Spaces---ShareRevokeSpace1 tree node inside browse tree
    Then "Share ShareRevokeSpace1" dialog should be visible
    And the access level of the sharing user in "Share ShareRevokeSpace1" dialog should be "View and use"
    When user removes the sharing user from "Share ShareRevokeSpace1" dialog
    And user clicks on OK button in "Share ShareRevokeSpace1" dialog
    Then the "Share ShareRevokeSpace1" dialog should close
    When user clicks on Spaces---ShareRevokeSpace1 tree node inside browse tree
    Then the context panel should show "ShareRevokeSpace1"
    And the sharing pane should not list the sharing user
    And no errors should have been logged

  Scenario: The second account no longer has the space
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareRevokeSpace2" shared by the running account should be visible
    And "ShareRevokeSpace1" shared by the running account should be absent
    When user signs in as themselves again
    Then the running account should be signed in
