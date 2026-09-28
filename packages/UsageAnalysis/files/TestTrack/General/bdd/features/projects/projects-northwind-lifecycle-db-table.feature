@dev-only @journey @serial @realizes:views.projects
Feature: A project of the NorthwindTest orders table: saved, reopened and shared
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the first scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  Test 2 of the md, as written: NorthwindTest > Schemas > public > orders is opened with Get All
  (830 rows), saved through the ribbon's Save dialog with Data sync on and the creation script
  DbQuery(Dbtests:PostgresTest, "public.orders", …), reopened from the Dashboards gallery, shared
  through the tile's Share... with the second account, which opens it too, and removed through
  Delete Project. Every reopen is claimed by the md's 830 rows and by the data-sync mark (the table
  was re-read from the database, not a snapshot loaded). Translated from the TestTrack case
  Projects/projects-lifecycle-db; Test 1 (the PostgresAll query) is
  projects-northwind-lifecycle-db-query, and the System:Datagrok version of both is
  packages/UsageAnalysis/bdd/features/projects/projects-lifecycle-db.feature.

  The md's setup asks for "a second user who can access the NorthwindTest connection". The
  connection is not shared with the second account on dev, and the project's Share dialog does not
  share it, so the feature grants the second account "View and use" on the connection for the run
  and revokes it at the end, reading the server back both times.

  The project name carries the run's time; the project (with its table and view) is removed at the
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
    And no project named "BDDNwLifeDbTable{time}" is on the server

  Scenario: Get All opens the database table into a table view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the current view should be a TableView view
    And the "orders" view should be current
    And the table should have 830 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Saved with Data sync, the creation script reads the table from the database
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "orders" project table in "Save project" dialog should be checked
    When user enters "BDDNwLifeDbTable{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "orders" project table in "Save project" dialog
    Then creation script text in "orders" project table in "Save project" dialog should be visible
    And creation script text in "orders" project table in "Save project" dialog should contain text "DbQuery(Dbtests:PostgresTest, \"public.orders\""
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDNwLifeDbTable{time}" uploaded' should have been shown
    And 1 project named "BDDNwLifeDbTable{time}" should be on the server
    And the "orders" table of the "BDDNwLifeDbTable{time}" project should be saved with data sync
    And "Share BDDNwLifeDbTable{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDNwLifeDbTable{time}" dialog
    Then the "Share BDDNwLifeDbTable{time}" dialog should close
    And no errors should have been logged

  Scenario: The table project reopens from Dashboards by re-reading the table
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbTable{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeDbTable{time} gallery card
    Then the current view should be a TableView view
    And the "orders" view should be current
    And the table should have 830 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table project is shared with the second account
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbTable{time}" into gallery search
    And user picks "Share..." from the context menu of BDDNwLifeDbTable{time} gallery card
    Then "Share BDDNwLifeDbTable{time}" dialog should be visible
    And "Share BDDNwLifeDbTable{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDNwLifeDbTable{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDNwLifeDbTable{time}" dialog
    Then the "Share BDDNwLifeDbTable{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDNwLifeDbTable{time} gallery card
    Then the context panel should show "BDDNwLifeDbTable{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of BDDNwLifeDbTable{time} gallery card
    Then the access level of the sharing user in "Share BDDNwLifeDbTable{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDNwLifeDbTable{time}" dialog
    Then the "Share BDDNwLifeDbTable{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account opens the shared table project
    Given user signs in as the sharing user
    And the browse panel is open
    And the sharing user may use the "Dbtests:PostgresTest" connection until the feature ends
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbTable{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeDbTable{time} gallery card
    Then the current view should be a TableView view
    And the "orders" view should be current
    And the table should have 830 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the table project
    Given user signs in as themselves again
    And the browse panel is open
    And the sharing user may use the "Dbtests:PostgresTest" connection until the feature ends
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeDbTable{time}" into gallery search
    Then BDDNwLifeDbTable{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDNwLifeDbTable{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDNwLifeDbTable{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDNwLifeDbTable{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDNwLifeDbTable{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
