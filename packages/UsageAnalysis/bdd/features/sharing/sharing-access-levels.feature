@journey @serial @realizes:sharing.share-dialog
Feature: Sharing at Full access, and with All users
  The Share dialog grants more than its default: a level checked in its list of levels before OK is
  the level the grant gets (a click on a level's name only selects the row; its check box picks it), and All users is shared with as a user is. Translated from
  playwright-public/Sharing sharing-all-users-and-owner (block A) and the access-level part of
  share-query-permissions, whose checks were API reads; sharing-dialog opens the list and cancels.

  The queries are this feature's own, ShareLevel1 and ShareLevel2 on System:Datagrok, made through
  the API and deleted when the feature starts and ends, their grants with them. The All users grant on
  ShareLevel2 is removed through the dialog before the end (`user removes "All users" from …`), and the
  Sharing pane is checked to no longer list it. The second account is the library's sharing user; every
  share is made with Send notifications off.

  Not translated, and why: what Full access lets the second account do, and that an All users grant
  reaches every account, are the server's answers (ApiTests, Dapi: permissions); that the All users
  grant leaves with its revoke needs a step the dialog lacks, and that the owner keeps every
  privilege with no grant left (sharing-all-users-and-owner C, sharing-edge-cases 2) is an API test
  still to write (sharing/MISSING.md); the Share with everyone permission gate of sharing-edge-cases 1
  read only a server flag.

  Serial: the Queries view's search text is the account's own setting.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a query "ShareLevel1" on "System:Datagrok" reads "select 1 as x"
    And a query "ShareLevel2" on "System:Datagrok" reads "select 2 as y"
    And the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareLevel" into gallery search
    Then "ShareLevel1" gallery card should be visible
    And "ShareLevel2" gallery card should be visible

  Scenario: A level checked in the list is the level the grant gets
    When user picks "Share..." from the context menu of "ShareLevel1" gallery card
    Then "Share ShareLevel1" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share ShareLevel1" dialog
    And user clicks on share access selector
    And user checks "Full access" tree node inside privilege tree
    And user presses Escape
    Then privilege tree should be hidden
    And share access selector should contain text "Full access"
    When user unchecks "Send notifications" input in "Share ShareLevel1" dialog
    And user clicks on OK button in "Share ShareLevel1" dialog
    Then the "Share ShareLevel1" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "ShareLevel2" gallery card
    Then the context panel should show "ShareLevel2"
    When user clicks on "ShareLevel1" gallery card
    Then the context panel should show "ShareLevel1"
    And the sharing pane should list the sharing user
    And Sharing section in context panel should contain text "has full access"
    When user picks "Share..." from the context menu of "ShareLevel1" gallery card
    Then "Share ShareLevel1" dialog should be visible
    And the access level of the sharing user in "Share ShareLevel1" dialog should be "Full access"
    When user clicks on CANCEL button in "Share ShareLevel1" dialog
    Then the "Share ShareLevel1" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: All users is shared with as a user is
    When user clicks on "ShareLevel1" gallery card
    Then the context panel should show "ShareLevel1"
    When user clicks on "ShareLevel2" gallery card
    Then the context panel should show "ShareLevel2"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should not contain text "All users"
    When user picks "Share..." from the context menu of "ShareLevel2" gallery card
    Then "Share ShareLevel2" dialog should be visible
    When user types "All" into "User, group, or email" input in "Share ShareLevel2" dialog
    And user clicks on "All users" table row in "Share ShareLevel2" dialog
    And user unchecks "Send notifications" input in "Share ShareLevel2" dialog
    And user clicks on OK button in "Share ShareLevel2" dialog
    Then the "Share ShareLevel2" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "ShareLevel1" gallery card
    Then the context panel should show "ShareLevel1"
    When user clicks on "ShareLevel2" gallery card
    Then the context panel should show "ShareLevel2"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should contain text "All users"
    When user picks "Share..." from the context menu of "ShareLevel2" gallery card
    Then "Share ShareLevel2" dialog should be visible
    When user removes "All users" from "Share ShareLevel2" dialog
    And user clicks on OK button in "Share ShareLevel2" dialog
    Then the "Share ShareLevel2" dialog should close
    When user clicks on "ShareLevel1" gallery card
    And user clicks on "ShareLevel2" gallery card
    Then the context panel should show "ShareLevel2"
    And Sharing section in context panel should not contain text "All users"
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown
