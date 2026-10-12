@dev-only @journey @serial @realizes:views.projects
Feature: A project of the NorthwindTest query PostgresAll: saved and reopened with its creation script
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the first scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  Test 1 of the md: the saved query PostgresAll is run by a double-click under NorthwindTest (830
  rows), saved through the ribbon's Save dialog with Data sync on and its creation script calling
  PostgresAll, reopened from the Dashboards gallery (830 rows, re-run by data sync), and the Save
  dialog of the reopened project still shows Data sync on and the call; Delete Project removes it.
  Translated from the TestTrack case Projects/projects-lifecycle-db; Test 2 (the orders table) is
  projects-northwind-lifecycle-db-table. Nothing of Test 1 is parked; the System:Datagrok version
  for other stands is parked (see the request document: its row counts differ per server and need
  a remembered row count).

  The project name carries the run's time (letters and digits only) and the project is removed
  (with its table and view) at the start and at the end. It is serial: the Dashboards search is
  shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDNwLifeDbQuery{time}" is on the server

  Scenario: The saved query runs from the tree into a table view
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the current view should be a TableView view
    And the "PostgresAll" view should be current
    And the table should have 830 rows
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Saved with Data sync, the creation script calls the query
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be checked
    When user enters "BDDNwLifeDbQuery{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "PostgresAll" project table in "Save project" dialog
    Then "PostgresAll" project table in "Save project" dialog should contain text ":PostgresAll()"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDNwLifeDbQuery{time}" uploaded' should have been shown
    And 1 project named "BDDNwLifeDbQuery{time}" should be on the server
    And the "PostgresAll" table of the "BDDNwLifeDbQuery{time}" project should be saved with data sync
    And "Share BDDNwLifeDbQuery{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDNwLifeDbQuery{time}" dialog
    Then the "Share BDDNwLifeDbQuery{time}" dialog should close
    And no errors should have been logged

  Scenario: The query project reopens from Dashboards by re-running the query
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbQuery{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeDbQuery{time} gallery card
    Then the current view should be a TableView view
    And the "PostgresAll" view should be current
    And the table should have 830 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And "Creation script" button in "PostgresAll" project table in "Save project" dialog should be visible
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be checked
    When user clicks on "Creation script" button in "PostgresAll" project table in "Save project" dialog
    Then "PostgresAll" project table in "Save project" dialog should contain text ":PostgresAll()"
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: Delete Project removes the query project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbQuery{time}" into gallery search
    Then BDDNwLifeDbQuery{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDNwLifeDbQuery{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDNwLifeDbQuery{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDNwLifeDbQuery{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDNwLifeDbQuery{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
