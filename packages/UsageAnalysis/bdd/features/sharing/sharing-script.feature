@journey @serial @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me
Feature: A script, shared and taken back
  A script shared through its Share dialog, seen from the other side: the second account finds it
  under My stuff > Shared with me and runs it, and loses it when the owner takes the share back.
  Translated from the recipient side of playwright-public/Sharing share-script-permissions (blocks D
  and F), whose second-account steps were API calls; the owner's Share dialog on a script is
  scripts/scripts-browser, and the dialog's own controls and the Sharing pane before a share are
  sharing-dialog (on a query); the Advanced editor of the old block C is gone from the dialog
  (GROK-20322).

  The scripts are this feature's own, made through the API and deleted when the feature starts and
  ends, their grants with them: ShareScript1, a JavaScript script with no inputs whose output is the
  string "shared-ok", and ShareScript2, shared too and never taken back, so that the second account
  still has the owner's node under Shared with me when ShareScript1's absence there is claimed. The second account is the library's sharing user; it signs in on the feature's page,
  and the share is made with Send notifications off.

  Not translated, and why: that the second account may not edit, delete or share the script further,
  and that View and use carries Execute, are the server's answers (ApiTests, Dapi: permissions).

  Serial: it signs in as the second account on its page.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And no script named "ShareScript1" is on the server
    And no script named "ShareScript2" is on the server
    And a script "ShareScript1" is on the server:
      """
      //language: javascript
      //output: string res
      res = 'shared-ok';
      """
    And a script "ShareScript2" is on the server:
      """
      //language: javascript
      //output: string res
      res = 'anchor';
      """
    And the context panel is open
    And the browse panel is open

  Scenario: The owner shares both scripts to view and use
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user enters "ShareScript2" into gallery search
    And user picks "Share..." from the context menu of "ShareScript2" gallery card
    And user picks the sharing user in "User, group, or email" input in "Share ShareScript2" dialog
    And user unchecks "Send notifications" input in "Share ShareScript2" dialog
    And user clicks on OK button in "Share ShareScript2" dialog
    Then the "Share ShareScript2" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user enters "ShareScript1" into gallery search
    And user picks "Share..." from the context menu of "ShareScript1" gallery card
    Then "Share ShareScript1" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share ShareScript1" dialog
    And user unchecks "Send notifications" input in "Share ShareScript1" dialog
    And user clicks on OK button in "Share ShareScript1" dialog
    Then the "Share ShareScript1" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "ShareScript1" gallery card
    Then the context panel should show "ShareScript1"
    And the sharing pane should list the sharing user
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account finds the script under Shared with me and runs it
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareScript1" shared by the running account should be visible
    When user picks "Run" from the context menu of "ShareScript1" shared by the running account
    Then the script results should show "res" as "shared-ok"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner takes the share back
    When user signs in as themselves again
    Then the running account should be signed in
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user enters "ShareScript1" into gallery search
    And user picks "Share..." from the context menu of "ShareScript1" gallery card
    Then "Share ShareScript1" dialog should be visible
    And the access level of the sharing user in "Share ShareScript1" dialog should be "View and use"
    When user removes the sharing user from "Share ShareScript1" dialog
    And user clicks on OK button in "Share ShareScript1" dialog
    Then the "Share ShareScript1" dialog should close
    When user clicks on "ShareScript1" gallery card
    Then the context panel should show "ShareScript1"
    And the sharing pane should not list the sharing user
    When user clears gallery search
    Then gallery search should have value ""

  Scenario: The second account no longer has the script
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareScript2" shared by the running account should be visible
    And "ShareScript1" shared by the running account should be absent
    When user signs in as themselves again
    Then the running account should be signed in
