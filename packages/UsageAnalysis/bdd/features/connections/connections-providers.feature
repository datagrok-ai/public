@connections @full-stand
Feature: Editing connections of the other providers
  The edit.md case repeated for Oracle, MariaDB, MySQL and MS SQL: a connection of each is renamed
  through Rename..., and without credentials its Test connection ends on an error balloon that
  names it. The connections are saved through the API with the coordinates of that provider's
  Samples Northwind and no credentials, and deleted at feature end — checked gone. These providers
  exist only on a full stand.

  The renamed connection is named by its tree label: the tree keeps the old name on a relabelled
  node until it is refreshed (switch to the path name once the core renames the node).

  Not translated, and why: "set the right login and password — test OK" for these providers — the
  test databases' passwords are not in the environment the suite reads (only the Postgres ones
  are: DG_PG_*); the Postgres flow claims it in connections-edit.feature.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario Outline: A <provider> connection is renamed, and its test without credentials fails
    Given no connection named "BDD-Conn-<short>-{run}, BDD-Conn-<short>-Renamed-{run}" is on the server
    And a "<provider>" connection named "BDD-Conn-<short>-{run}" is on the server
    And Databases---<node> tree node inside browse tree is expanded
    When user right-clicks on Databases---<node>---BDD-Conn-<short>-{run} tree node inside browse tree
    And user picks "Rename..." from the open menu
    Then "Rename dataconnection" dialog should be visible
    When user enters "BDD-Conn-<short>-Renamed-{run}" into Name input in "Rename dataconnection" dialog
    And user clicks on OK button in "Rename dataconnection" dialog
    Then the "Rename dataconnection" dialog should close
    And 1 connection named "BDD-Conn-<short>-Renamed-{run}" should be on the server
    And "BDD-Conn-<short>-Renamed-{run}" tree node inside browse tree should be visible
    Given user watches the task bar
    When user right-clicks on "BDD-Conn-<short>-Renamed-{run}" tree node inside browse tree
    And user picks "Test connection" from the open menu
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an error balloon containing "BDD-Conn-<short>-Renamed-{run}"

    Examples:
      | provider | node    | short   |
      | Oracle   | Oracle  | Oracle  |
      | MariaDB  | MariaDB | MariaDB |
      | MySQL    | MySQL   | MySQL   |
      | MS SQL   | MS-SQL  | MSSQL   |
