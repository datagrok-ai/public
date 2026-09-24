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
