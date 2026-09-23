@connections
Feature: Adding a database connection
  A provider's New connection dialog asks for the fields that provider needs and keeps OK disabled
  until the connection has a name; TEST reports on a balloon what the connection can reach, and OK
  saves it under the provider. Translated from TestTrack Connections/adding.md (playwright-public
  connections/01-adding.test.ts); the dialog fields of every provider were probed on dev.

  Every connection a scenario makes is named BDD-Conn-…-{run} and deleted at feature end, with any
  chat on it; the server is checked for it before and after. The providers other than Postgres
  live only on a full stand (@full-stand), and a TEST that connects needs the Northwind password
  from DG_PG_PASSWORD (@needs-credentials) — typed from the environment, never printed. Without it
  the same flow still claims the save and a TEST that fails for the missing password, which the old
  spec never checked: it accepted any balloon after TEST.

  Not translated, and why: the md's "check the connection in the Databases list" as a separate
  step — the saved connection is claimed both on the server and as its tree node here; the second
  connection of the md repeats the first with other names, so it is one scenario, not two.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario Outline: The New connection dialog of <provider> asks for its fields
    When user right-clicks on Databases---<node> tree node inside browse tree
    And user picks "New connection..." from the open menu
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

    @full-stand
    Examples:
      | provider | node    |
      | MS SQL   | MS-SQL  |
      | Oracle   | Oracle  |
      | MySQL    | MySQL   |
      | MariaDB  | MariaDB |

  Scenario: A connection string replaces the server fields
    When user right-clicks on Databases---Postgres tree node inside browse tree
    And user picks "New connection..." from the open menu
    Then "Add new connection" dialog should be visible
    When user selects "Connection string" in Configure input in "Add new connection" dialog
    Then Conn-String input in "Add new connection" dialog should be visible
    And Server input in "Add new connection" dialog should be hidden
    When user clicks on CANCEL button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And no errors should have been logged

  Scenario: A connection without a password is saved, and its TEST fails for the password
    Given no connection named "BDD-Conn-Add-{run}" is on the server
    When user right-clicks on Databases---Postgres tree node inside browse tree
    And user picks "New connection..." from the open menu
    Then "Add new connection" dialog should be visible
    And OK button in "Add new connection" dialog should be disabled
    When user enters "BDD-Conn-Add-{run}" into Name input in "Add new connection" dialog
    Then OK button in "Add new connection" dialog should be enabled
    When user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    And user enters "54322" into Port input in "Add new connection" dialog
    And user enters "northwind" into Db input in "Add new connection" dialog
    And user enters "datagrok" into Login input in "Add new connection" dialog
    Given user watches the task bar
    When user clicks on TEST button in "Add new connection" dialog
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an error balloon containing "failed to connect"
    When user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "BDD-Conn-Add-{run}" should be on the server
    Given Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree should be visible

  @needs-credentials
  Scenario: A connection with the right credentials passes its TEST and is saved
    Given no connection named "BDD-Conn-Add-Ok-{run}" is on the server
    When user right-clicks on Databases---Postgres tree node inside browse tree
    And user picks "New connection..." from the open menu
    Then "Add new connection" dialog should be visible
    When user enters "BDD-Conn-Add-Ok-{run}" into Name input in "Add new connection" dialog
    And user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    And user enters "54322" into Port input in "Add new connection" dialog
    And user enters "northwind" into Db input in "Add new connection" dialog
    And user enters "datagrok" into Login input in "Add new connection" dialog
    And user enters the DG_PG_PASSWORD secret into Password input in "Add new connection" dialog
    Given user watches the task bar
    When user clicks on TEST button in "Add new connection" dialog
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an info balloon containing "connected successfully"
    When user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "BDD-Conn-Add-Ok-{run}" should be on the server
    Given Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---BDD-Conn-Add-Ok-{run} tree node inside browse tree should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged
