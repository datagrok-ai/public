@tutorials @serial @realizes:tutorials.data-connectors
Feature: The Data Connectors tutorial
  Walks Data access > Data Connectors from its card to the end: a Postgres connection added from the
  Browse tree's context menu and filled in, a query created on it from the same menu and named, and the
  query run. Each step is claimed as ticked and as done — the dialog, every field, the connection and
  the query view, and the rows the query brings.
  Translated from playwright-tests/e2e/tutorials/data-connectors.test.ts. The connection's parameters
  are the ones the tutorial shows the learner (a public demo database).

  The tutorial needs Grok Connect (its prerequisite), and the query runs on db.datagrok.ai, outside the
  stand: where the stand cannot reach it, the walk ends
  there (skipped, not failed). The connection and the query have fixed names learners share, so only the
  running user's own are removed, before the walk and after it.
  Fixed in the tutorial for this translation: the counter was 12 for 11 steps; it told the learner to
  click "Add connection...", and the menu item is "New connection...".
  Fails at "no errors should have been logged" while GROK-20982 is open: the Activity count of the new
  connection and query (/log/count) answers with an error body that the client parses as a number.
  Walked end to end on a local stand with the Grok Connect prerequisite lifted, since that stand
  reports no service health.

  Serial, for the fixed names and because a finished tutorial writes its completion record into the
  account's settings, which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the user's own connection "Starbucks" and query "Get Starbucks US" are removed now and at feature end
    And the "Data Connectors" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Data Connectors tutorial
    # the tutorial will not start without the connector service (its prerequisite)
    Given the stand runs the "Grok Connect" service
    When user starts the "Data Connectors" tutorial
    Then the tutorial progress should be 1 of 11
    Given the tutorial step "Create a connection to Postgres server" should not be done yet
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then the tutorial step "Create a connection to Postgres server" should be done
    And "Add new connection" dialog should be visible
    When user enters "Starbucks" into "Name" input in "Add new connection" dialog
    Then the tutorial step "Set \"Name\" to \"Starbucks\"" should be done
    When user enters "db.datagrok.ai" into "Server" input in "Add new connection" dialog
    Then the tutorial step "Set \"Server\" to \"db.datagrok.ai\"" should be done
    When user enters "54324" into "Port" input in "Add new connection" dialog
    Then the tutorial step "Set \"Port\" to \"54324\"" should be done
    When user enters "starbucks" into "Db" input in "Add new connection" dialog
    Then the tutorial step "Set \"Db\" to \"starbucks\"" should be done
    When user enters "datagrok" into "Login" input in "Add new connection" dialog
    Then the tutorial step "Set \"Login\" to \"datagrok\"" should be done
    When user enters "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF" into "Password" input in "Add new connection" dialog
    Then the tutorial step "Set \"Password\" to \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\"" should be done
    When user clicks on OK button in "Add new connection" dialog
    Then the tutorial step "Click \"OK\"" should be done
    And 1 connection named "Starbucks" should be on the server

    Given the tutorial step "Create a data query to the \"Starbucks\" data connection" should not be done yet
    When user picks "New Query..." from the context menu of Databases---Postgres---Starbucks tree node inside browse tree
    Then the tutorial step "Create a data query to the \"Starbucks\" data connection" should be done
    Given the tutorial step "Set \"Name\" to \"Get Starbucks US\"" should not be done yet
    When user enters "Get Starbucks US" into "Name" input
    Then the tutorial step "Set \"Name\" to \"Get Starbucks US\"" should be done

    # the query runs on the demo database outside the stand
    Given the stand can reach the database of the "Starbucks" connection
    When user puts "select * from starbucks_us" on the first line of code editor
    And user clicks on play icon
    Then the tutorial step "Add \"select * from starbucks_us\" to the editor and hit \"Play\"" should be done
    And the "Data Connectors" tutorial should be completed
    And the tutorial should have listed 11 steps
    And the tutorial progress should be 11 of 11
    And no hint should be shown
    And no errors should have been logged
