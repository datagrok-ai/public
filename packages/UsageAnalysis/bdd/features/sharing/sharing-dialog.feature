@journey @serial @realizes:sharing.share-dialog
Feature: The Share dialog and the Sharing pane
  What a user sees when about to share an entity: the Sharing pane of the context panel, the Share
  dialog it opens, and the Resulting permissions dialog. Nothing is shared here; the share itself
  and what the other user then gets are sharing-recipient. Translated from the TestTrack case
  Sharing/sharing-dialog.

  The entity is a query of this feature's own, ShareDialog1 on System:Datagrok, made through the API
  and deleted when the feature starts and ends; it is found in Browse > Platform > Functions >
  Queries, whose search the feature clears at the end. A second query, ShareDialog2, is there only to
  be made current in between: the context panel rebuilds the Sharing pane only when the current
  object changes, so each scenario that cancels makes ShareDialog2 current and then ShareDialog1
  again before it reads the pane. The second account is the library's sharing user. Nothing is
  granted: every dialog here is cancelled.

  Send notifications is checked on its own only while the account setting
  sendNotificationOnShareCheckbox is on (permissions_browser.dart takes the default from it; unchecking
  in the dialog does not write it back): the scenario claims the setting first.

  Not claimed: that the running account's own grant row is not shown in the dialog (the dialog keeps
  it as a hidden row). The suggestions are reached as the "All users" row of the typeahead's table, a
  group every stand has. Which privileges each level holds for each kind of entity, and that the
  server refuses the rest, is an API test (ApiTests, Dapi: permissions).

  Serial: the Queries view's search text is the account's own setting, shared with every feature
  that searches it.

  Background:
    Given user is logged in
    And a query "ShareDialog1" on "System:Datagrok" reads "select 1 as x"
    And a query "ShareDialog2" on "System:Datagrok" reads "select 2 as y"
    And the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareDialog" into gallery search
    Then "ShareDialog1" gallery card should be visible
    And "ShareDialog2" gallery card should be visible

  Scenario: The Sharing pane names the owner and opens the dialog
    When user clicks on "ShareDialog1" gallery card
    Then the context panel should show "ShareDialog1"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    And "Share..." button in Sharing section in context panel should be visible
    And "Calculate resulting permissions for this entity" button in Sharing section in context panel should be visible
    When user clicks on "Share..." button in Sharing section in context panel
    Then "Share ShareDialog1" dialog should be visible
    When user clicks on CANCEL button in "Share ShareDialog1" dialog
    Then the "Share ShareDialog1" dialog should close
    When user clicks on "ShareDialog2" gallery card
    Then the context panel should show "ShareDialog2"
    When user clicks on "ShareDialog1" gallery card
    Then the context panel should show "ShareDialog1"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    And no errors should have been logged

  Scenario: The dialog asks who and how much, and offers a message once someone is picked
    When user picks "Share..." from the context menu of "ShareDialog1" gallery card
    Then "Share ShareDialog1" dialog should be visible
    And "User, group, or email" input in "Share ShareDialog1" dialog should be visible
    And share access selector should contain text "View and use"
    And OK button in "Share ShareDialog1" dialog should be visible
    And CANCEL button in "Share ShareDialog1" dialog should be visible
    And "Send notifications" input in "Share ShareDialog1" dialog should be hidden
    And share message should be hidden
    And "Remove" icon in "Share ShareDialog1" dialog should be hidden
    When user picks the sharing user in "User, group, or email" input in "Share ShareDialog1" dialog
    Then "Remove" icon in "Share ShareDialog1" dialog should be visible
    And "Send notifications" input in "Share ShareDialog1" dialog should be visible
    And share message should be visible
    And the "sendNotificationOnShareCheckbox" shell setting should be true
    And "Send notifications" input in "Share ShareDialog1" dialog should be checked
    And "Share ShareDialog1" dialog should contain text "(Will be sent to 1 user)"
    When user clicks on CANCEL button in "Share ShareDialog1" dialog
    Then the "Share ShareDialog1" dialog should close
    When user clicks on "ShareDialog2" gallery card
    Then the context panel should show "ShareDialog2"
    When user clicks on "ShareDialog1" gallery card
    Then the context panel should show "ShareDialog1"
    And the sharing pane should not list the sharing user
    And no errors should have been logged

  Scenario: Escape closes the list of levels and leaves the dialog open, and CANCEL closes it while the suggestions are open (GROK-20758)
    When user picks "Share..." from the context menu of "ShareDialog1" gallery card
    Then "Share ShareDialog1" dialog should be visible
    When user clicks on share access selector
    Then the following elements should be visible:
      | "Full access" tree node inside privilege tree  |
      | "View and use" tree node inside privilege tree |
      | "Edit" tree node inside privilege tree         |
      | "Delete" tree node inside privilege tree       |
      | "Share" tree node inside privilege tree        |
    And "View and use" tree node inside privilege tree should be checked
    When user presses Escape
    Then privilege tree should be hidden
    And "Share ShareDialog1" dialog should be visible
    And "All users" table row in "Share ShareDialog1" dialog should be hidden
    When user types "All" into "User, group, or email" input in "Share ShareDialog1" dialog
    Then "All users" table row in "Share ShareDialog1" dialog should be visible
    When user clicks on CANCEL button in "Share ShareDialog1" dialog
    Then the "Share ShareDialog1" dialog should close
    And no errors should have been logged

  Scenario: Your own account is not offered (GROK-20755)
    When user picks "Share..." from the context menu of "ShareDialog1" gallery card
    Then "Share ShareDialog1" dialog should be visible
    When user types the current user's login into "User, group, or email" input in "Share ShareDialog1" dialog
    Then the user suggestions should not offer the current user
    When user clicks on CANCEL button in "Share ShareDialog1" dialog
    Then the "Share ShareDialog1" dialog should close
    And no errors should have been logged

  Scenario: Resulting permissions shows what a user gets
    When user clicks on "ShareDialog1" gallery card
    Then the context panel should show "ShareDialog1"
    And the sharing pane should not list the sharing user
    When user clicks on "Calculate resulting permissions for this entity" button in Sharing section in context panel
    Then "Resulting permissions" dialog should be visible
    And "Resulting permissions" dialog should contain text "ShareDialog1"
    And "User, group or role" input in "Resulting permissions" dialog should be visible
    And CANCEL button in "Resulting permissions" dialog should be visible
    And "Resulting permissions" dialog should not contain text "No access"
    When user picks the sharing user in "User, group or role" input in "Resulting permissions" dialog
    Then "Resulting permissions" dialog should contain text "No access"
    When user clicks on CANCEL button in "Resulting permissions" dialog
    Then the "Resulting permissions" dialog should close
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown
