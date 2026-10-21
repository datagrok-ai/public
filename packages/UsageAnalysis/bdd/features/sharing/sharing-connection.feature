@journey @serial @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me
Feature: A connection, shared and taken back
  A database connection shared from its Browse node, seen from the other side: the second account
  finds it under My stuff > Shared with me, and loses it when the owner takes the share back.
  Translated from playwright-public/Sharing share-connection-permissions (blocks D and F, whose
  second-account steps were API calls) and the deleted TestTrack case
  Sharing/share-connection-permissions; the owner's Share dialog from the Connections gallery is
  connections/connections-browser, and the dialog's own controls and the Sharing pane before a share
  are sharing-dialog (on a query); the Advanced editor of the old block C is gone from the dialog
  (GROK-20322).

  The connections are this feature's own, Postgres connections the library points at the stand's
  test Northwind, made through the API and deleted when the feature starts and ends, their grants
  with them: ShareConn1, shared and taken back, and ShareConn2, shared too and never taken back, so
  that the second account still has the owner's node under Shared with me when ShareConn1's absence
  there is claimed. The second account is the library's sharing user; it signs in on the feature's
  page, and every share is made with Send notifications off.

  Known failure: core#977 — right after a View and use share the Sharing pane names the grant "has
  special permissions", not "can view and use"; its scenario claims only that label. Not translated, and why: what View and use and Full access let the second
  account do with the connection (query, schema, files, no edit) is the server's answer (ApiTests,
  Dapi: permissions); browsing the schema as the second account needs the database's credentials,
  which a feature's own connection does not carry; the "will also be shared" notice of a query's
  Share dialog, and the TestTrack case's claim that sharing a
  connection shares nothing along with it, are in sharing/MISSING.md.

  Serial: it signs in as the second account on its page.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a "Postgres" connection named "ShareConn1" is on the server
    And a "Postgres" connection named "ShareConn2" is on the server
    And the context panel is open
    And the browse panel is open

  Scenario: The owner shares both connections from their Browse nodes
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of Databases---Postgres---ShareConn2 tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share ShareConn2" dialog
    And user unchecks "Send notifications" input in "Share ShareConn2" dialog
    And user clicks on OK button in "Share ShareConn2" dialog
    Then the "Share ShareConn2" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user picks "Share..." from the context menu of Databases---Postgres---ShareConn1 tree node inside browse tree
    Then "Share ShareConn1" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share ShareConn1" dialog
    And user unchecks "Send notifications" input in "Share ShareConn1" dialog
    And user clicks on OK button in "Share ShareConn1" dialog
    Then the "Share ShareConn1" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree
    Then the context panel should show "ShareConn1"
    And the sharing pane should list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @known-failure:66 @realizes:core#977
  Scenario: The Sharing pane names the level of the connection's grant (core#977)
    When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree
    Then the context panel should show "ShareConn1"
    And the sharing pane should list the sharing user
    And Sharing section in context panel should contain text "can view and use"

  Scenario: The second account finds the connection under Shared with me
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareConn1" shared by the running account should be visible
    And "ShareConn2" shared by the running account should be visible
    And no errors should have been logged

  Scenario: The owner takes the share back
    When user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of Databases---Postgres---ShareConn1 tree node inside browse tree
    Then "Share ShareConn1" dialog should be visible
    And the access level of the sharing user in "Share ShareConn1" dialog should be "View and use"
    When user removes the sharing user from "Share ShareConn1" dialog
    And user clicks on OK button in "Share ShareConn1" dialog
    Then the "Share ShareConn1" dialog should close
    When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree
    Then the context panel should show "ShareConn1"
    And the sharing pane should not list the sharing user
    And no errors should have been logged

  Scenario: The second account no longer has the connection
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareConn2" shared by the running account should be visible
    And "ShareConn1" shared by the running account should be absent
    When user signs in as themselves again
    Then the running account should be signed in
