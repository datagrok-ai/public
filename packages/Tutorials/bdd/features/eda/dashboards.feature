@tutorials @serial @realizes:tutorials.dashboards
Feature: The Dashboards tutorial
  Walks Exploratory data analysis > Dashboards from its card to the end: a parameterised query written
  on the Starbucks connection, run from Browse for one state, a bar chart added, the result saved as a
  project and closed, reopened from Dashboards and refreshed for another state. Each step is claimed as
  ticked and as done — the connection and the query on the server, the 645 New York stores, the bar
  chart, the saved project, the 84 Louisiana stores after the refresh.
  Translated from playwright-tests/e2e/tutorials/dashboards.test.ts.

  The tutorial needs Grok Connect (its prerequisite), and its query runs on db.datagrok.ai, outside the
  stand, so the walk ends there where the stand cannot reach it (skipped, not failed). The connection,
  the query and the project have fixed names learners share: only the running user's own are removed,
  before the walk and after it.
  Fixed in the tutorial for this translation: it told the learner to click "Add connection...", and the
  menu item is "New connection..."; the Dashboards step looked its tree row up before Close all rebuilt
  the tree, and listened on the caption only.
  Walked end to end on a local stand with the Grok Connect prerequisite lifted, since that stand
  reports no service health.

  Serial, for the fixed names and because a finished tutorial writes its completion record into the
  account's settings, which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the user's own project "Coffee sales dashboard" is removed now and at feature end
    And the user's own connection "Starbucks" and query "Stores in @state" are removed now and at feature end
    And the "Dashboards" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Dashboards tutorial
    # the tutorial will not start without the connector service (its prerequisite)
    Given the stand runs the "Grok Connect" service
    When user starts the "Dashboards" tutorial
    Then the tutorial progress should be 1 of 27
    Given the tutorial step "Create a connection to Postgres server" should not be done yet
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then the tutorial step "Create a connection to Postgres server" should be done
    When user enters "Starbucks" into "Name" input in "Add new connection" dialog
    Then the tutorial step "Set \"Name\" to \"Starbucks\"" should be done
    When user enters "db.datagrok.ai" into "Server" input in "Add new connection" dialog
    And user enters "54324" into "Port" input in "Add new connection" dialog
    And user enters "starbucks" into "Db" input in "Add new connection" dialog
    And user enters "datagrok" into "Login" input in "Add new connection" dialog
    And user enters "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF" into "Password" input in "Add new connection" dialog
    Then the tutorial step "Set \"Password\" to \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\"" should be done
    When user clicks on OK button in "Add new connection" dialog
    Then the tutorial step "Click \"OK\"" should be done
    Given the tutorial step "Create a data query to the \"Starbucks\" data connection" should not be done yet
    When user picks "New Query..." from the context menu of Databases---Postgres---Starbucks tree node inside browse tree
    Then the tutorial step "Create a data query to the \"Starbucks\" data connection" should be done
    Given the tutorial step "Set \"Name\" to \"Stores in @state\"" should not be done yet
    When user enters "Stores in @state" into "Name" input
    Then the tutorial step "Set \"Name\" to \"Stores in @state\"" should be done
    When user puts "select * from starbucks_us where state = @state;" on the first line of code editor
    Then the tutorial step "Add \"select * from starbucks_us where state = @state;\" to the editor" should be done
    When user puts "--input: string state" on the first line of code editor
    Then the tutorial step "Add \"--input: string state\" as the first line of the query" should be done
    When user clicks on SAVE button
    Then the tutorial step "Save the query" should be done
    And 1 query named "Stores in @state" should be on the server
    Given the tutorial step "Find Browse on the sidebar and click" should not be done yet
    When user clicks on browse tab
    Then the tutorial step "Find Browse on the sidebar and click" should be done
    Given the tutorial step "Find the created query in the browse view, right-click it and hit Run" should not be done yet
    When user expands Databases---Postgres---Starbucks tree node inside browse tree
    And user picks "Run" from the context menu of Databases---Postgres---Starbucks---Stores-in-@state tree node inside browse tree
    Then the tutorial step "Find the created query in the browse view, right-click it and hit Run" should be done
    When user enters "NY" into "State" input in "Stores in @state" dialog
    Then the tutorial step "Set state to \"NY\"" should be done

    # the query runs on the demo database outside the stand
    Given the stand can reach the database of the "Starbucks" connection
    When user clicks on OK button in "Stores in @state" dialog
    Then the tutorial step "Click \"OK\" to run the query" should be done
    And the table should have 645 rows

    Given the tutorial step "Open bar chart" should not be done yet
    When user clicks on bar-chart icon in toolbox
    Then the tutorial step "Open bar chart" should be done
    Given the tutorial step "Save a project" should not be done yet
    When user clicks on SAVE button
    Then the tutorial step "Save a project" should be done
    When user enters "Coffee sales dashboard" into "Name" text input in "Save project" dialog
    Then the tutorial step "Set the project name to \"Coffee sales dashboard\"" should be done
    When user clicks on OK button in "Save project" dialog
    Then the tutorial step "Click \"OK\"" should be done 2 times
    When user clicks on CANCEL button in "Share Coffee sales dashboard" dialog
    Then the tutorial step "Skip the sharing step" should be done
    And 1 project named "Coffee sales dashboard" should be on the server
    Given the tutorial step "Close the project" should not be done yet
    When user picks "Close All" from the context menu of browse tab
    Then the tutorial step "Close the project" should be done
    Given the tutorial step "Open browse and click on Dashboards" should not be done yet
    When user clicks on Dashboards tree node inside browse tree
    Then the tutorial step "Open browse and click on Dashboards" should be done
    When user double-clicks on "Coffee sales dashboard" gallery card
    Then the tutorial step "Find and open your project" should be done
    Given the tutorial step "Set State to LA" should not be done yet
    When user enters "LA" into "State" input
    Then the tutorial step "Set State to LA" should be done
    When user clicks on REFRESH button
    Then the tutorial step "Click REFRESH button" should be done
    # the dashboard re-runs its query for the other state
    And the table should have 84 rows

    And the "Dashboards" tutorial should be completed
    And the tutorial should have listed 27 steps
    And the tutorial progress should be 27 of 27
    And no hint should be shown
    And no errors should have been logged
