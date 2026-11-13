@journey @realizes:sharing.share-dialog
Feature: A connection shared from its Browse node, and what its Sharing pane says
  A database connection shared from its Browse node, and the Sharing pane of the context panel after
  the share. Translated from playwright-public/Sharing share-connection-permissions and the deleted
  TestTrack case Sharing/share-connection-permissions; the owner's Share dialog from the Connections
  gallery is connections/connections-browser, and the dialog's own controls and the Sharing pane
  before a share are sharing-dialog (on a query); the Advanced editor of the old block C is gone from
  the dialog (GROK-20322).

  The connection is this feature's own, ShareConn1, a Postgres connection the library points at the
  stand's test Northwind, made through the API and deleted when the feature starts and ends, its
  grants with it. The share goes to the library's sharing user, with Send notifications off.

  core#977 (fixed): right after a View and use share the Sharing pane names the grant "can view and
  use" (it used to say "has special permissions"); the last scenario guards it.

  Not translated, and why: what the second account gets — the connection to query, its schema and
  files, no edit, its loss after the revoke, and none of the connection's queries along with it — is
  the server's answer (ApiTests, Dapi: permissions: "a database connection", "a database connection
  shared on its own", "a revoke"); the "will also be shared" notice of a query's Share dialog is in
  sharing/MISSING.md.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a "Postgres" connection named "ShareConn1" is on the server
    And the context panel is open
    And the browse panel is open

  Scenario: The owner shares the connection from its Browse node
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
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

  @realizes:core#977
  Scenario: The Sharing pane names the level of the connection's grant (core#977)
    When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree
    Then the context panel should show "ShareConn1"
    And the sharing pane should list the sharing user
    And Sharing section in context panel should contain text "can view and use"
