@connections
Feature: Adding a database connection
  A provider's New connection dialog asks for the fields that provider needs and keeps OK disabled
  until the connection has a name; OK saves it under the provider. Translated from TestTrack
  Connections/adding.md (playwright-public connections/01-adding.test.ts); the dialog fields of
  every provider were probed on dev.

  Every connection a scenario makes is named BDD-Conn-…-{run} and deleted at feature end, with any
  chat on it; the server is checked for it before and after.

  Not translated, and why: what TEST answers — whether the connector reaches the database and logs
  in is the database server's answer, not the dialog's (see the bdd library's CLAUDE.md, "What never
  becomes a feature"); the scenarios claim that TEST is offered. The md's "check the connection in
  the Databases list" as a separate step — the saved connection is claimed both on the server and as
  its tree node here; the second connection of the md repeats the first with other names, so it is
  one scenario, not two.

  The fields are claimed on Postgres and on MS SQL, a provider behind the tree's "Show more"; Oracle,
  MySQL and MariaDB ask for the same fields and are not repeated.

  Edge cases: a name of spaces keeps OK disabled ("Can't be empty", `Validators.notEmpty` on the name,
  `data_connection_meta.dart`), and CANCEL with a name typed saves nothing. A clone of the stand's
  NorthwindTest (the DBTests package's connection), saved from Clone... under a name with the run's
  time, is expected to carry the credentials and list its schemas (core#1176: today the clone says
  "Unavailable", a known failure). That a rename in the clone's Edit dialog keeps the credentials is
  claimed once the clone lists its schemas.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And the hidden providers of the Databases tree are shown

  Scenario Outline: The New connection dialog of <provider> asks for its fields
    When user picks "New connection..." from the context menu of Databases---<node> tree node inside browse tree
    Then "Add new connection" dialog should be visible
    And the following elements should be visible:
      | Name input in "Add new connection" dialog     |
      | Server input in "Add new connection" dialog   |
      | Db input in "Add new connection" dialog       |
      | Login input in "Add new connection" dialog    |
      | Password input in "Add new connection" dialog |
      | TEST button in "Add new connection" dialog    |
    And OK button in "Add new connection" dialog should be disabled
    When user clicks on CANCEL button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And no errors should have been logged

    Examples:
      | provider | node     |
      | Postgres | Postgres |
      | MS SQL   | MS-SQL   |

  Scenario: A connection string replaces the server fields
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    When user selects "Connection string" in Configure input in "Add new connection" dialog
    Then Conn-String input in "Add new connection" dialog should be visible
    And Server input in "Add new connection" dialog should be hidden
    When user clicks on CANCEL button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And no errors should have been logged

  Scenario: A named connection is saved under its provider
    Given no connection named "BDD-Conn-Add-{run}" is on the server
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    And OK button in "Add new connection" dialog should be disabled
    When user enters "BDD-Conn-Add-{run}" into Name input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be enabled
    When user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    And user enters "54322" into Port input in "Add new connection" dialog
    And user enters "northwind" into Db input in "Add new connection" dialog
    And user enters "datagrok" into Login input in "Add new connection" dialog
    And user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "BDD-Conn-Add-{run}" should be on the server
    And the "BDD-Conn-Add-{run}" connection on the server should have the data source "Postgres"
    Given Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @edge-cases
  Scenario: A name of spaces keeps OK disabled, and CANCEL with a name typed saves nothing
    Given no connection named "BDDConnSpaces{time}" is on the server
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    When user enters "BDDConnSpaces{time}" into Name input in "Add new connection" dialog
    And user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be enabled
    When user enters "   " into Name input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be disabled
    When user enters "BDDConnSpaces{time}" into Name input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be enabled
    When user clicks on CANCEL button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 0 connections named "BDDConnSpaces{time}" should be on the server
    And no errors should have been logged

  @edge-cases @known-failure:115 @realizes:core#1176
  Scenario: A clone keeps the credentials and lists its schemas
    Today a clone of NorthwindTest says "Unavailable" under its node: the copy does not carry the
    credentials (core#1176). Remove the tag when the clone lists its schemas.
    Given no connection named "BDDConnClone{time}" is on the server
    And the stand has a reachable "PostgresTest" connection
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Clone..." from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree
    Then "Edit Connection" dialog should be visible
    When user enters "BDDConnClone{time}" into Name input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 1 connection named "BDDConnClone{time}" should be on the server
    Given Databases---Postgres---BDDConnClone{time} tree node inside browse tree is expanded
    Then Databases---Postgres---BDDConnClone{time}---Schemas tree node inside browse tree should become visible within 60 seconds
