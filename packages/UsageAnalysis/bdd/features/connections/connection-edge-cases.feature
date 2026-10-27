@connections @edge-cases
Feature: Connection edge cases: names the dialogs refuse, CANCEL, and a delete that takes the queries
  Guards for the connection dialogs and for the delete of a connection that has a query. The New
  connection and Edit Connection dialogs refuse a name of spaces ("Can't be empty", `Validators.notEmpty`
  on the name, `data_connection_meta.dart`); the Rename dialog refuses a name with no letter in it
  ("Name must contain at least one letter", `DbEntityMeta.getRenameDialog`). CANCEL saves nothing, a
  changed name included. "Deleting a connection also deletes a query" (help
  `access/databases/databases.md`, the caution under Query Editor): the query saved under a
  connection goes with it.

  Every connection and query is the feature's own: Postgres connections to db.datagrok.ai:54322/northwind
  saved through the API without credentials (nothing here connects to the database), named with the
  run's time and deleted, with their queries, when the feature ends. A connection with the same name as
  another (core#546's class) is an ApiTests case, and what TEST answers for a wrong host or password is
  the database server's reply, not the dialog's: see MISSING.md.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario: A name of spaces keeps OK disabled, and CANCEL saves nothing
    Given no connection named "BDDConnSpaces{time}" is on the server
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    When user enters "   " into Name input in "Add new connection" dialog
    And user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be disabled
    When user enters "BDDConnSpaces{time}" into Name input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be enabled
    When user clicks on CANCEL button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 0 connections named "BDDConnSpaces{time}" should be on the server
    And no errors should have been logged

  Scenario: CANCEL in the Edit dialog keeps the name, and a name of spaces is refused there too
    Given a "Postgres" connection named "BDDConnEditEdge{time}" is on the server
    And no connection named "BDDConnEditedEdge{time}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Edit..." from the context menu of Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    When user enters "   " into Name input in "Edit Connection" dialog
    Then OK button in "Edit Connection" dialog should be disabled
    When user enters "BDDConnEditedEdge{time}" into Name input in "Edit Connection" dialog
    Then OK button in "Edit Connection" dialog should be enabled
    When user clicks on CANCEL button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 1 connection named "BDDConnEditEdge{time}" should be on the server
    And 0 connections named "BDDConnEditedEdge{time}" should be on the server
    And Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Rename... refuses a name with no letter
    Given a "Postgres" connection named "BDDConnRenameEdge{time}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree
    Then "Rename dataconnection" dialog should be visible
    When user enters "2026" into Name input in "Rename dataconnection" dialog
    Then OK button in "Rename dataconnection" dialog should be disabled
    When user enters "BDDConnRenamedEdge{time}" into Name input in "Rename dataconnection" dialog
    Then OK button in "Rename dataconnection" dialog should be enabled
    When user clicks on CANCEL button in "Rename dataconnection" dialog
    Then the "Rename dataconnection" dialog should close
    And 1 connection named "BDDConnRenameEdge{time}" should be on the server
    And no errors should have been logged

  Scenario: Deleting a connection deletes the query saved under it
    Given no query named "BDDConnCascadeQ{time}" is on the server
    And a "Postgres" connection named "BDDConnCascade{time}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
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
