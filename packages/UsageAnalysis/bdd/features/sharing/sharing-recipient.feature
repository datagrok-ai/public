@journey @serial @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me @realizes:views.queries @realizes:views.space
Feature: What the other user gets, and losing it
  A query shared once through the Share dialog, seen from the other side: the second account finds
  it under My stuff > Shared with me and runs it, cannot put a query shared with it into a space
  shared with it for viewing only, and loses the query when the owner takes the share back.
  Translated from the TestTrack case Sharing/sharing-recipient.

  The entities are this feature's own, made through the API and deleted when the feature starts and
  ends, their grants with them: the queries ShareRecipient1 (select 1 as x) and ShareRecipient2
  (select 2 as y) on System:Datagrok, and the root space ShareSpace1, which holds one file,
  ShareSpaceStay1.csv (the project's fixtures/cars-small.csv): the second account's view of the space
  shows it before the query's absence there is claimed. The second account's
  starting state, ShareRecipient2 and ShareSpace1 shared with it at View and use, is made through the
  Share dialog in the first scenario: there is no step that grants through the API. The second account
  is the library's sharing user; it signs in on the feature's page, and every share is made with Send
  notifications off. The Datagrok connection is not shared with it by this feature.

  That the second account may not edit, rename, delete or share the query further is the server's
  answer and an API test (ApiTests, Dapi: permissions); the feature claims the refusal the Share
  dialog shows.

  Adding to a space takes Edit on it, for a Link as for a Move (core#693: a Link went through on View
  until 2026-10-12). The Link's claim is the refusal and then the outcome: after Escape, which closes
  the dialog if the refusal keeps it open, the space's view lists its file and not the query. The
  balloon is also the barrier: the absence is read only once the server has answered the Link. The
  Move entity dialog opens on Link for the second account, so the Move attempt picks Move first.

  After a refused Move the Browse tree takes no second drop of the same query onto the same space
  until it is refreshed (a drop onto another space still opens the dialog); the refusal is remembered
  for that pair, as designed, so the Link attempt starts from a refreshed tree.

  Serial: it signs in as the second account on its page, and the Queries view's search text is the
  account's own setting (cleared at the end of the owner's last visit).

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a query "ShareRecipient1" on "System:Datagrok" reads "select 1 as x"
    And a query "ShareRecipient2" on "System:Datagrok" reads "select 2 as y"
    And a space named "ShareSpace1" is on the server
    And the "fixtures/cars-small.csv" file of the project is in the space "ShareSpace1" as "ShareSpaceStay1.csv"
    And the context panel is open
    And the browse panel is open

  Scenario: The second account's starting state: ShareRecipient2 and ShareSpace1 shared with it
    Given Spaces tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of Spaces---ShareSpace1 tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share ShareSpace1" dialog
    And user unchecks "Send notifications" input in "Share ShareSpace1" dialog
    And user clicks on OK button in "Share ShareSpace1" dialog
    Then the "Share ShareSpace1" dialog should close
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareRecipient2" into gallery search
    And user picks "Share..." from the context menu of "ShareRecipient2" gallery card
    And user picks the sharing user in "User, group, or email" input in "Share ShareRecipient2" dialog
    And user unchecks "Send notifications" input in "Share ShareRecipient2" dialog
    And user clicks on OK button in "Share ShareRecipient2" dialog
    Then the "Share ShareRecipient2" dialog should close
    And no error or warning balloon should have been shown

  Scenario: The owner shares the query to view and use (GROK-19235)
    When user enters "ShareRecipient1" into gallery search
    And user picks "Share..." from the context menu of "ShareRecipient1" gallery card
    Then "Share ShareRecipient1" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share ShareRecipient1" dialog
    And user unchecks "Send notifications" input in "Share ShareRecipient1" dialog
    And user clicks on OK button in "Share ShareRecipient1" dialog
    Then the "Share ShareRecipient1" dialog should close
    And an info balloon containing "Shared" should have been shown
    And no error or warning balloon should have been shown
    When user clicks on "ShareRecipient1" gallery card
    Then the context panel should show "ShareRecipient1"
    And the sharing pane should list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    And Sharing section in context panel should contain text "On ShareRecipient1"
    And Sharing section in context panel should contain text "can view and use"

  Scenario: The second account finds the query under Shared with me and runs it (GROK-19401, GROK-18525, GROK-18210)
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    And user double-clicks on "ShareRecipient1" shared by the running account
    Then the "ShareRecipient1" view should be current
    And the table should have 1 row
    And the table should have the columns "x"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account cannot share the query further
    Given the browse panel is open
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    And user picks "Share..." from the context menu of "ShareRecipient1" shared by the running account
    Then "Share ShareRecipient1" dialog should be visible
    When user picks the running account in "User, group, or email" input in "Share ShareRecipient1" dialog
    And user unchecks "Send notifications" input in "Share ShareRecipient1" dialog
    And user clicks on OK button in "Share ShareRecipient1" dialog
    Then an error balloon containing "You don't have a permission to share this object" should have been shown
    And "Share ShareRecipient1" dialog should be visible
    When user clicks on CANCEL button in "Share ShareRecipient1" dialog
    Then the "Share ShareRecipient1" dialog should close

  Scenario: The second account cannot Move a query into a space shared for viewing
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    And user drags "ShareRecipient2" shared by the running account to "ShareSpace1" shared by the running account
    Then "Move entity" dialog should be visible
    And "Move entity" dialog should contain text "to ShareSpace1 project"
    And "Move entity" dialog should contain text "ShareRecipient2"
    And choice input in "Move entity" dialog should have value "Link"
    When user selects "Move" in choice input in "Move entity" dialog
    Then choice input in "Move entity" dialog should have value "Move"
    When user clicks on YES button in "Move entity" dialog
    Then an error balloon containing "insufficient privileges to edit the project" should have been shown
    And "Move entity" dialog should be visible
    When user clicks on CANCEL button in "Move entity" dialog
    Then the "Move entity" dialog should close

  Scenario: The second account asks for a Link into the same space
    When user refreshes the browse tree
    And user expands "." shared by the running account
    And user drags "ShareRecipient2" shared by the running account to "ShareSpace1" shared by the running account
    Then "Move entity" dialog should be visible
    And choice input in "Move entity" dialog should have value "Link"

  @realizes:core#693
  Scenario: The Link is refused, and the space does not get the query (core#693)
    When user clicks on YES button in "Move entity" dialog
    Then an error balloon containing "insufficient privileges to edit the project" should have been shown
    When user presses Escape
    And user double-clicks on "ShareSpace1" shared by the running account
    Then the "ShareSpace1" view should be current
    And "ShareSpaceStay1.csv" link in gallery should be visible
    And "ShareRecipient2" link in gallery should be absent

  Scenario: The owner takes the share back
    When user signs in as themselves again
    Then the running account should be signed in
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareRecipient1" into gallery search
    And user picks "Share..." from the context menu of "ShareRecipient1" gallery card
    Then "Share ShareRecipient1" dialog should be visible
    And the access level of the sharing user in "Share ShareRecipient1" dialog should be "View and use"
    When user removes the sharing user from "Share ShareRecipient1" dialog
    And user clicks on OK button in "Share ShareRecipient1" dialog
    Then the "Share ShareRecipient1" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "ShareRecipient1" gallery card
    Then the context panel should show "ShareRecipient1"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    When user clears gallery search
    Then gallery search should have value ""

  Scenario: The second account no longer has the query
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareRecipient2" shared by the running account should be visible
    And "ShareRecipient1" shared by the running account should be absent
    When user signs in as themselves again
    Then the running account should be signed in
