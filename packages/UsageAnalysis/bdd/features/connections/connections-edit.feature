@connections @journey
Feature: Editing a database connection
  A connection is renamed from its Edit dialog and from Rename..., copied by Clone..., and "Test
  connection" reports the credentials it was saved with. Translated from TestTrack
  Connections/edit.md (playwright-public connections/03-edit.test.ts). The subject is a Postgres
  connection to db.datagrok.ai:54322/northwind saved through the API without credentials,
  BDD-Conn-Edit-{run}; it and every name it takes are deleted at feature end and checked gone.

  "Repeat for Oracle, MariaDB, MySQL, MS SQL" is connections-providers.feature: an outline over
  connections of those providers, each renamed and tested without credentials. "Set the right login
  and password — test OK" is connections-credentials.feature (@needs-credentials): a journey is one
  test, so a scenario that needs secrets lives where a run without them can leave it out.

  A renamed connection is named in the tree by its label: the tree relabels the node at once but
  keeps the old name on it until the tree is refreshed (a core naming gap, not something a person
  sees) — switch to the Databases---Postgres---<new name> tree node once the core renames the node.

  Not translated, and why: the old spec's 1.5 s sleep after Save — it covered a Test connection that
  answered from the credentials cached before the save. The positive test here runs the Edit
  dialog's own TEST on what the dialog holds, which is what the md asks; whether the context-menu
  test after a save reads stale credentials is a question for the platform, not a wait for the test.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And no connection named "BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}" is on the server
    And a "Postgres" connection named "BDD-Conn-Edit-{run}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded

  Scenario: The Edit dialog shows what the connection was saved with
    When user right-clicks on Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree
    And user picks "Edit..." from the open menu
    Then "Edit Connection" dialog should be visible
    And Name input in "Edit Connection" dialog should have value "BDD-Conn-Edit-{run}"
    And Server input in "Edit Connection" dialog should have value "db.datagrok.ai"
    And Db input in "Edit Connection" dialog should have value "northwind"
    When user clicks on CANCEL button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And no errors should have been logged

  Scenario: Renaming through the Edit dialog
    When user right-clicks on Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree
    And user picks "Edit..." from the open menu
    Then "Edit Connection" dialog should be visible
    When user enters "BDD-Conn-Edited-{run}" into Name input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 1 connection named "BDD-Conn-Edited-{run}" should be on the server
    And 0 connections named "BDD-Conn-Edit-{run}" should be on the server
    And "BDD-Conn-Edited-{run}" tree node inside browse tree should be visible
    And no error or warning balloon should have been shown

  Scenario: Renaming through Rename...
    When user right-clicks on "BDD-Conn-Edited-{run}" tree node inside browse tree
    And user picks "Rename..." from the open menu
    Then "Rename dataconnection" dialog should be visible
    When user enters "BDD-Conn-Renamed-{run}" into Name input in "Rename dataconnection" dialog
    And user clicks on OK button in "Rename dataconnection" dialog
    Then the "Rename dataconnection" dialog should close
    And 1 connection named "BDD-Conn-Renamed-{run}" should be on the server
    And 0 connections named "BDD-Conn-Edited-{run}" should be on the server
    And "BDD-Conn-Renamed-{run}" tree node inside browse tree should be visible
    And no error or warning balloon should have been shown

  Scenario: Clone... offers a copy under another name and saves nothing when cancelled
    When user right-clicks on "BDD-Conn-Renamed-{run}" tree node inside browse tree
    And user picks "Clone..." from the open menu
    Then "Edit Connection" dialog should be visible
    And Name input in "Edit Connection" dialog should have value "Copy of BDD-Conn-Renamed-{run}"
    And Server input in "Edit Connection" dialog should have value "db.datagrok.ai"
    When user clicks on CANCEL button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    And 0 connections named "Copy of BDD-Conn-Renamed-{run}" should be on the server
    And no errors should have been logged

  # @full-stand: a stand whose connector answers a refused login. On a local stand the connector
  # waits out its 180 s socket timeout and says nothing at all, so the claim cannot be made there.
  @slow @full-stand
  Scenario: Wrong credentials make Test connection fail
    When user right-clicks on "BDD-Conn-Renamed-{run}" tree node inside browse tree
    And user picks "Edit..." from the open menu
    Then "Edit Connection" dialog should be visible
    When user enters "bdd_nobody" into Login input in "Edit Connection" dialog
    And user enters "wrong" into Password input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    Given user watches the task bar
    When user right-clicks on "BDD-Conn-Renamed-{run}" tree node inside browse tree
    And user picks "Test connection" from the open menu
    Then the task bar should have shown "Testing"
    And the connection test should have ended on an error balloon containing "password authentication failed"
