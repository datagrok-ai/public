@connections
Feature: Deleting a database connection
  Delete... asks "Are you sure?" and removes the connection for everyone on DELETE; CANCEL keeps it.
  Translated from TestTrack Connections/delete.md (playwright-public connections/05-delete.test.ts).
  The subjects are Postgres connections saved through the API without credentials; whatever a
  scenario leaves is deleted at feature end and checked gone.

  "Deleting a connection also deletes a query" (help `access/databases/databases.md`, the caution
  under Query Editor): a query saved under a connection in the query editor goes with it.

  Not translated, and why: the md's "Browse > Platform > Connections" — the Platform section has no
  Connections node any more; the list of connections is the provider's "Browse connections" view,
  and deleting from its gallery is the third scenario. The md's YES button is DELETE now.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario: DELETE removes the connection from the server and the tree, and a query saved under it goes too
    Given a "Postgres" connection named "BDD-Conn-Delete-{run}" is on the server
    And no query named "BDDConnCascadeQ{time}" is on the server
    And a "Postgres" connection named "BDDConnCascade{time}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Delete..." from the context menu of Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete connection \"BDD-Conn-Delete-{run}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 connections named "BDD-Conn-Delete-{run}" should be on the server
    And Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user picks "New Query..." from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDConnCascadeQ{time}" into Name input
    And user replaces the code of code editor with "select * from products"
    And user clicks on Save button
    Then 1 query named "BDDConnCascadeQ{time}" should be on the server
    When user closes all views
    And user picks "Delete..." from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 connections named "BDDConnCascade{time}" should be on the server
    And 0 queries named "BDDConnCascadeQ{time}" should be on the server
    And Databases---Postgres---BDDConnCascade{time} tree node inside browse tree should be absent
    And no errors should have been logged

  Scenario: CANCEL keeps the connection
    Given a "Postgres" connection named "BDD-Conn-Keep-{run}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Delete..." from the context menu of Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 1 connection named "BDD-Conn-Keep-{run}" should be on the server
    And Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: A connection is deleted from the connections gallery
    Given a "Postgres" connection named "BDD-Conn-Gallery-{run}" is on the server
    When user picks "Browse connections" from the context menu of Databases---Postgres tree node inside browse tree
    Then the "Postgres" view should be current
    When user types "BDD-Conn-Gallery-{run}" into gallery search
    Then "BDD-Conn-Gallery-{run}" link in gallery should be visible
    When user picks "Delete..." from the context menu of "BDD-Conn-Gallery-{run}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 connections named "BDD-Conn-Gallery-{run}" should be on the server
    And "BDD-Conn-Gallery-{run}" link in gallery should be absent
    And no error or warning balloon should have been shown
