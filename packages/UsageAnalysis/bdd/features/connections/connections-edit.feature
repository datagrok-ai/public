@connections @journey
Feature: Editing a database connection
  A connection is renamed from its Edit dialog and from Rename..., and copied by Clone....
  Translated from TestTrack Connections/edit.md (playwright-public connections/03-edit.test.ts).
  The subject is a Postgres connection to db.datagrok.ai:54322/northwind saved through the API
  without credentials, BDD-Conn-Edit-{run}; it and every name it takes are deleted at feature end
  and checked gone.

  Not translated, and why: "Test connection" with the right and with wrong credentials, and "repeat
  for Oracle, MariaDB, MySQL, MS SQL" — what a test answers is the database server's reply to a
  login, not the dialog's (see the bdd library's CLAUDE.md, "What never becomes a feature"), and the
  rename is the same dialog for every provider. The old spec's 1.5 s sleep after Save went with it.

  Edge cases: on two more connections of the same kind, named with the run's time and deleted at
  feature end, the Edit Connection dialog refuses a name of spaces ("Can't be empty") and its CANCEL
  keeps the name, a changed one included; the Rename dialog refuses a name with no letter ("Name must
  contain at least one letter", `DbEntityMeta.getRenameDialog`). A clone that carries the credentials
  (core#1176) is in `connections-add.feature`.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And no connection named "BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}" is on the server
    And a "Postgres" connection named "BDD-Conn-Edit-{run}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded

  Scenario: The Edit dialog shows what the connection was saved with
    When user picks "Edit..." from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    And Name input in "Edit Connection" dialog should have value "BDD-Conn-Edit-{run}"
    And Server input in "Edit Connection" dialog should have value "db.datagrok.ai"
    And Db input in "Edit Connection" dialog should have value "northwind"
    When user clicks on CANCEL button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And no errors should have been logged

  Scenario: Renaming through the Edit dialog
    When user picks "Edit..." from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    When user enters "BDD-Conn-Edited-{run}" into Name input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 1 connection named "BDD-Conn-Edited-{run}" should be on the server
    And 0 connections named "BDD-Conn-Edit-{run}" should be on the server
    And Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree should be visible
    And no error or warning balloon should have been shown

  Scenario: Renaming through Rename...
    When user picks "Rename..." from the context menu of Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree
    Then "Rename dataconnection" dialog should be visible
    When user enters "BDD-Conn-Renamed-{run}" into Name input in "Rename dataconnection" dialog
    And user clicks on OK button in "Rename dataconnection" dialog
    Then the "Rename dataconnection" dialog should close
    And 1 connection named "BDD-Conn-Renamed-{run}" should be on the server
    And 0 connections named "BDD-Conn-Edited-{run}" should be on the server
    And Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree should be visible
    And no error or warning balloon should have been shown

  Scenario: Clone... offers a copy under another name and saves nothing when cancelled
    When user picks "Clone..." from the context menu of Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    And Name input in "Edit Connection" dialog should have value "Copy of BDD-Conn-Renamed-{run}"
    And Server input in "Edit Connection" dialog should have value "db.datagrok.ai"
    When user clicks on CANCEL button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 0 connections named "Copy of BDD-Conn-Renamed-{run}" should be on the server
    And no errors should have been logged

  @edge-cases
  Scenario: The Edit dialog refuses a name of spaces and keeps the name on CANCEL
    Given a "Postgres" connection named "BDDConnEditEdge{time}" is on the server
    And no connection named "BDDConnEditedEdge{time}" is on the server
    When user refreshes the browse tree
    Given Databases tree node inside browse tree is expanded
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

  @edge-cases
  Scenario: Rename... refuses a name with no letter and keeps the name on CANCEL
    Given a "Postgres" connection named "BDDConnRenameEdge{time}" is on the server
    And no connection named "BDDConnRenamedEdge{time}" is on the server
    When user refreshes the browse tree
    Given Databases tree node inside browse tree is expanded
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
