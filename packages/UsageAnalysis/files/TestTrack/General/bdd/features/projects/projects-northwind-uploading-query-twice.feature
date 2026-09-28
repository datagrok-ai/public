@dev-only @serial @realizes:views.projects @realizes:data.menu.link-tables @realizes:file.menu.save.tables-as-project
Feature: The NorthwindTest query PostgresAll run twice, linked and saved with Data sync on and off
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the first scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  TestCase2 of the md's source matrix (TestCase2Sync and TestCase2NoSync): the query PostgresAll
  (830 rows) is run twice under NorthwindTest — the second run opens as "PostgresAll (2)" — and the
  two are linked selection-to-filter by customerid through Data > Link Tables...; the first two orders
  of the first table (VINET, TOMSP) leave the md's 11 rows in the second. Both tables are saved in one
  project through the ribbon's Save dialog, once with Data sync on for both and once with it off; the
  project reopens from its Dashboards card with both tables, the Save dialog of the reopened project
  shows each table's Creation script only when it was saved with Data sync, and the link still
  filters. Translated from the TestTrack case Projects/uploading; the other NorthwindTest rows are
  projects-northwind-uploading, and the System:Datagrok version of this row is in
  packages/UsageAnalysis/bdd/features/projects/projects-uploading.feature.

  The md says to double-click the query. Both runs here are the query's context menu > Run in the
  Browse tree (it runs a query without parameters at once, no dialog): a table opened by a second
  double-click on the same query gets no Data sync switch in the Save dialog (its row holds the switch
  hidden), and the md's "set Data sync for both tables" needs the switch on both.

  The link is checked with the first pair before the save and with the third and fourth orders
  (HANAR, VICTE) after the reopen, so a filter the project merely restored cannot pass for a link that
  no longer works. A reopen is claimed to be a reopen: the tables open before the save are marked in
  memory, Close All is claimed to leave no table, and every table that comes back must be a frame
  without the mark. Before a reopened table's Creation script is claimed shown or hidden, its Data
  sync switch is claimed on or off. The query belongs to the Dbtests package and is never changed.

  Every project (with its tables and views) is named with the run's time and is removed before its
  row starts and again when the feature ends. The Save dialog's preview logs "Unable to find element
  in cloned iframe" (GROK-18606, known noise), so the save is claimed to log nothing else; the reopen
  is claimed to log nothing at all. It is @serial: the Dashboards search and the uploads are shared
  with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open

  Scenario Outline: PostgresAll run twice, linked by customerid and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the "PostgresAll" view should be current
    And table "PostgresAll" should have 830 rows
    When user clicks on browse tab
    And user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the "PostgresAll (2)" view should be current
    And table "PostgresAll (2)" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "PostgresAll" in link tables first table
    And user selects "PostgresAll (2)" in link tables second table
    And user selects "customerid" in link tables first key
    And user selects "customerid" in link tables second key
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "PostgresAll -> PostgresAll"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on "PostgresAll" view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "PostgresAll" should be selected
    When user clicks on "PostgresAll (2)" view
    Then 11 rows of table "PostgresAll (2)" should pass the filter
    And status bar should contain text "Filtered: 11"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # the dialog fills its table rows after it shows; a switch clicked before its row is drawn is missed
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be checked
    And Data sync switch in "PostgresAll (2)" project table in "Save project" dialog should be checked
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "PostgresAll" project table in "Save project" dialog
    And user <switch> Data sync switch in "PostgresAll (2)" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "PostgresAll" table of the "<Project>" project should be saved <saved>
    And the "PostgresAll (2)" table of the "<Project>" project should be saved <saved>
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
    Then table "PostgresAll" should have been <how> with 830 rows
    And table "PostgresAll (2)" should have been <how> with 830 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "PostgresAll" project table in "Save project" dialog should be <script>
    And Data sync switch in "PostgresAll (2)" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "PostgresAll (2)" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on "PostgresAll" view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "PostgresAll" should be selected
    When user clicks on "PostgresAll (2)" view
    Then 24 rows of table "PostgresAll (2)" should pass the filter
    And status bar should contain text "Filtered: 24"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                  | switch   | saved          | state | script  | how                   |
      | ON   | BDDNwUpload2Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDNwUpload2NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
