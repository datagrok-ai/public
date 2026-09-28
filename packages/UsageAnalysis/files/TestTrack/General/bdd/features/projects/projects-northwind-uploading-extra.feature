@dev-only @serial @realizes:views.projects @realizes:file.menu.save.tables-as-project
Feature: Projects uploaded from a Get Top 100 result of the NorthwindTest orders table
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and every row fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  The md's "Project from the Get Top 100 query result" as written, on the Northwind orders it names:
  Get Top 100 on NorthwindTest > Schemas > public > orders opens 100 of its 830 rows, its creation
  script reads the table with limit = 100, and the result is saved with Data sync on
  (UI_Top100_Sync) and off (UI_Top100_NoSync). Each reopens from its Dashboards card: with Data sync
  on the query is re-run, with it off the stored result opens — both with 100 rows. Translated from
  the TestTrack case Projects/uploading-ui; the md's other cases (an SDF file, the scratchpad's
  upload dialog) and the System:Datagrok version of this one are
  packages/UsageAnalysis/bdd/features/projects/projects-uploading-extra.feature, and "Project from
  two local files" is not translated (the operating system's file picker with files of the local
  machine, which a CI agent does not have).

  A reopen is claimed to be a reopen: the table open before the save is marked in memory, Close All
  is claimed to leave no table, and the table that comes back must be a frame without the mark.

  Every project (with its table and view) is named with the run's time and removed before its row
  starts and when the feature ends. The Save dialog's preview logs "Unable to find element in cloned
  iframe" (GROK-18606, known noise), so a save is claimed to log nothing else; the reopen is claimed
  to log nothing at all. It is @serial: the Dashboards search and the uploads are shared with every
  feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario Outline: A Get Top 100 result of the NorthwindTest orders is saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "Get Top 100" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 100 rows
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on "Creation script" button in "orders" project table in "Save project" dialog
    Then "orders" project table in "Save project" dialog should contain text "limit = 100"
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "orders" project table in "Save project" dialog
    Then "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "orders" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on <Project> gallery card
    Then table "orders" should have been <how> with 100 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                  | switch   | script  | saved          | how                   |
      | ON   | BDDNwUpXTopSync{time}    | checks   | visible | with data sync | reloaded by data sync |
      | OFF  | BDDNwUpXTopNoSync{time}  | unchecks | hidden  | as a snapshot  | loaded as a snapshot  |
