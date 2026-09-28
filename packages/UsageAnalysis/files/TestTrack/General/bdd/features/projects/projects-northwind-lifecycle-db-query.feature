@dev-only @journey @serial @realizes:views.projects
Feature: A project of the NorthwindTest query PostgresAll: saved, reopened, shared and renamed
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the first scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  Test 1 of the md, as written: the saved query PostgresAll is run by a double-click under
  NorthwindTest (830 rows), saved through the ribbon's Save dialog with Data sync on and its
  creation script calling PostgresAll, reopened from the Dashboards gallery, shared through the
  tile's Share... with the second account, which opens it too, renamed and reopened under its new
  name, and removed through Delete Project. Every reopen is claimed by the md's 830 rows and by the
  data-sync mark (the query was re-run, not a snapshot loaded). Translated from the TestTrack case
  Projects/projects-lifecycle-db; Test 2 (the orders table) is
  projects-northwind-lifecycle-db-table, and the System:Datagrok version of both is
  packages/UsageAnalysis/bdd/features/projects/projects-lifecycle-db.feature.

  The md's setup asks for "a second user who can access the NorthwindTest connection". The
  connection is not shared with the second account on dev, and the project's Share dialog does not
  share it (without a grant the second account's open fails with "Connection not found"), so the
  feature grants the second account "View and use" on the connection for the run and revokes it at
  the end, reading the server back both times. PostgresAll and the connection belong to the Dbtests
  package and are otherwise never changed.

  The two project names carry the run's time and are removed (with their tables and views) at the
  start and at the end. @serial: the Dashboards search is shared with every feature that saves a
  project.

  Not translated, and why: Logout and signing in with the second user's credentials — the platform's
  Logout ends every session of the account, which all workers of a run share, so the second account
  is entered through its own session ("user signs in as the sharing user"); Close All from the left
  sidebar's context menu is done through the shell (closing views is not the claim).

  Background:
    Given user is logged in
    And the browse panel is open
    And the sharing user may use the "Dbtests:PostgresTest" connection until the feature ends
    And no project named "BDDNwLifeDbQuery{time}" is on the server
    And no project named "BDDNwLifeDbQueryRenamed{time}" is on the server

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
    Then creation script text in "PostgresAll" project table in "Save project" dialog should be visible
    And creation script text in "PostgresAll" project table in "Save project" dialog should contain text ":PostgresAll()"
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

  Scenario: The query project is shared with the second account
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbQuery{time}" into gallery search
    And user picks "Share..." from the context menu of BDDNwLifeDbQuery{time} gallery card
    Then "Share BDDNwLifeDbQuery{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDNwLifeDbQuery{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDNwLifeDbQuery{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDNwLifeDbQuery{time}" dialog
    Then the "Share BDDNwLifeDbQuery{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDNwLifeDbQuery{time} gallery card
    Then the context panel should show "BDDNwLifeDbQuery{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of BDDNwLifeDbQuery{time} gallery card
    Then the access level of the sharing user in "Share BDDNwLifeDbQuery{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDNwLifeDbQuery{time}" dialog
    Then the "Share BDDNwLifeDbQuery{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account opens the shared query project
    Given user signs in as the sharing user
    And the browse panel is open
    And the sharing user may use the "Dbtests:PostgresTest" connection until the feature ends
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

  Scenario: Renamed, the query project opens under its new name
    Given user signs in as themselves again
    And the browse panel is open
    And the sharing user may use the "Dbtests:PostgresTest" connection until the feature ends
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbQuery{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDNwLifeDbQuery{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDNwLifeDbQueryRenamed{time}" into Name input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDNwLifeDbQueryRenamed{time}" should be on the server
    And 0 projects named "BDDNwLifeDbQuery{time}" should be on the server
    When user enters "BDDNwLifeDbQueryRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeDbQueryRenamed{time} gallery card
    Then the current view should be a TableView view
    And the table should have 830 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the query project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbQueryRenamed{time}" into gallery search
    Then BDDNwLifeDbQueryRenamed{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDNwLifeDbQueryRenamed{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDNwLifeDbQueryRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDNwLifeDbQueryRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDNwLifeDbQueryRenamed{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
