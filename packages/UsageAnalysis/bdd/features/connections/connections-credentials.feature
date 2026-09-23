@connections @needs-credentials
Feature: Connections that log in with real credentials
  The two cases whose claim is a successful login: an existing connection's Edit dialog passes TEST
  once it holds the right login and password (edit.md, "set the right login/password — test OK"),
  and the external provider's connection is created from the dialog against the writable test
  database (external-provider.md, its first step). The login of the test server is "datagrok", the
  same literal the old specs default to; only the passwords are secrets (DG_PG_PASSWORD,
  DG_PG_EXT_PASSWORD), typed from the environment and never
  printed; a run without them leaves this feature out (--grep-invert @needs-credentials). The
  connections are named BDD-Conn-…-{run} and deleted at feature end — checked gone.

  Not translated, and why: nothing of the two cases beyond what the other connections features claim
  without secrets.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario: The right credentials pass the Edit dialog's TEST
    Given a "Postgres" connection named "BDD-Conn-Creds-{run}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree
    And user picks "Edit..." from the open menu
    Then "Edit Connection" dialog should be visible
    When user enters "datagrok" into Login input in "Edit Connection" dialog
    And user enters the DG_PG_PASSWORD secret into Password input in "Edit Connection" dialog
    Given user watches the task bar
    When user clicks on TEST button in "Edit Connection" dialog
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an info balloon containing "connected successfully"
    When user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    Given user watches the task bar
    When user right-clicks on Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree
    And user picks "Test connection" from the open menu
    Then the connection test should have ended on an info balloon containing "connected successfully"

  Scenario: The external provider's connection is created from the dialog and passes its test
    Given no connection named "BDD-Conn-Ext-{run}" is on the server
    When user right-clicks on Databases---Postgres tree node inside browse tree
    And user picks "New connection..." from the open menu
    Then "Add new connection" dialog should be visible
    When user enters "BDD-Conn-Ext-{run}" into Name input in "Add new connection" dialog
    And user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    And user enters "54327" into Port input in "Add new connection" dialog
    And user enters "test" into Db input in "Add new connection" dialog
    And user enters the DG_PG_EXT_LOGIN secret into Login input in "Add new connection" dialog
    And user enters the DG_PG_EXT_PASSWORD secret into Password input in "Add new connection" dialog
    Given user watches the task bar
    When user clicks on TEST button in "Add new connection" dialog
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an info balloon containing "connected successfully"
    When user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "BDD-Conn-Ext-{run}" should be on the server
