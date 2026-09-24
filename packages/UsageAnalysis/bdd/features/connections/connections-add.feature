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

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

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
      | Oracle   | Oracle   |
      | MySQL    | MySQL    |
      | MariaDB  | MariaDB  |

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
