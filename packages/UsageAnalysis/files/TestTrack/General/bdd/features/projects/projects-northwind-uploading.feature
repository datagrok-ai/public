@dev-only @serial @realizes:views.projects @realizes:data.menu.link-tables @realizes:data.menu.aggregate-rows @realizes:file.menu.save.tables-as-project @realizes:viewers.pivot-viewer
Feature: Projects uploaded from the NorthwindTest query and table, with Data sync on and off
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and every row fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  The rows of the md's source matrix that read NorthwindTest (TestCase3, 6, 8), as written. Two
  tables opened from Browse are linked selection-to-filter through Data > Link Tables..., the link is checked in the
  status bar of the second table, and both are saved in one project through the ribbon's Save
  dialog — once with Data sync on for both tables and once with it off. The project reopens from its
  Dashboards card with both tables (rebuilt by data sync, or loaded from the snapshot) and a link
  that still filters; the Save dialog of the reopened project shows each table's Creation script
  only when it was saved with Data sync. Then an aggregation from Data > Aggregate Rows... over the
  NorthwindTest orders table (GROK-19103) is saved and reopened inside its project. Every row is one
  Scenario Outline example that makes and removes its own fixtures. Translated from the TestTrack
  case Projects/uploading; the System:Datagrok version of these rows, and the rows that read no
  database (TestCase1, 4, 5, 7), are packages/UsageAnalysis/bdd/features/projects/projects-uploading.feature.

  The md's rows here:
  - TestCase2 (PostgresAll run twice) is projects-northwind-uploading-query-twice.
  - TestCase3 and TestCase6: customers.csv (91 rows) from Files > Demo > northwind or from a Space,
    linked to PostgresAll by the customer key (CustomerID / customerid); the first two customers
    (ALFKI, ANATR) leave the md's 10 orders. After the reopen the next two (ANTON, AROUT) leave 20.
    The Space is made per row through the JS API, with a copy of customers.csv in its Files; the
    Create Space dialog and the drag with Copy of the md's setup are claimed in
    projects-lifecycle-spaces. The Files row makes the Space too and does not open it.
  - TestCase8: NorthwindTest > Schemas > public > orders opened with Get All (830 rows) and
    aggregated by customerid with count(orderid): one row per customer, 89. The pivot's values are
    not compared against the table (the pivot table reports only the rows its grid shows, and this
    one has 89): its readings and the published table's rows and columns are claimed. After the
    reopen the Dashboards panel lists New Dashboard and the project, with both tables under the
    project and no separate project for the aggregation.

  The link is checked with the first pair before the save and with the second pair after the
  reopen, so a filter the project merely restored cannot pass for a link that no longer works. A
  reopen is claimed to be a reopen: the tables open before the save are marked in memory, Close All
  is claimed to leave no table, and every table that comes back must be a frame without the mark.
  Before a reopened table's Creation script is claimed shown or hidden, its Data sync switch is
  claimed on or off. The query and the table belong to the Dbtests package and are never changed.

  Every project (with its tables and views) and every space is named with the run's time and is
  removed before its row starts and again when the feature ends. The Save dialog's preview logs
  "Unable to find element in cloned iframe" (GROK-18606, known noise), so the save is claimed to log
  nothing else; the reopen is claimed to log nothing at all. It is @serial: the Dashboards search
  and the uploads are shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And user clears the saved pivot table parameters

  Scenario Outline: customers.csv from <Source> and the query PostgresAll, linked and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And a space "<Space>" holds copies of the northwind files "customers.csv"
    When user clicks on "Refresh" icon inside browse toolbar
    And Spaces tree node inside browse tree is expanded
    And Spaces---<Space> tree node inside browse tree is expanded
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks <customers node> tree node inside browse tree
    Then the "customers" view should be current
    And table "customers" should have 91 rows
    When user clicks on browse tab
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the "PostgresAll" view should be current
    And table "PostgresAll" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "customers" in link tables first table
    And user selects "PostgresAll" in link tables second table
    And user selects "CustomerID" in link tables first key
    And user selects "customerid" in link tables second key
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "customers -> PostgresAll"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on "PostgresAll" view
    Then 10 rows of table "PostgresAll" should pass the filter
    And status bar should contain text "Filtered: 10"
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # the dialog fills its table rows after it shows; a switch clicked before its row is drawn is missed
    And Data sync switch in "customers" project table in "Save project" dialog should be checked
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be checked
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "customers" project table in "Save project" dialog
    And user <switch> Data sync switch in "PostgresAll" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "customers" table of the "<Project>" project should be saved <saved>
    And the "PostgresAll" table of the "<Project>" project should be saved <saved>
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
    Then table "customers" should have been <how> with 91 rows
    And table "PostgresAll" should have been <how> with 830 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "PostgresAll" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "PostgresAll" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on "PostgresAll" view
    Then 20 rows of table "PostgresAll" should pass the filter
    And status bar should contain text "Filtered: 20"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Source                   | Sync | Project                  | Space                     | customers node                                     | switch   | saved          | state | script  | how                   |
      | Files > Demo > northwind | ON   | BDDNwUpload3Sync{time}   | BDDNwUploadS3Sync{time}   | Files---Demo---northwind---customers.csv           | checks   | with data sync | on    | visible | reloaded by data sync |
      | Files > Demo > northwind | OFF  | BDDNwUpload3NoSync{time} | BDDNwUploadS3NoSync{time} | Files---Demo---northwind---customers.csv           | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
      | a space                  | ON   | BDDNwUpload6Sync{time}   | BDDNwUploadS6Sync{time}   | Spaces---BDDNwUploadS6Sync{time}---customers.csv   | checks   | with data sync | on    | visible | reloaded by data sync |
      | a space                  | OFF  | BDDNwUpload6NoSync{time} | BDDNwUploadS6NoSync{time} | Spaces---BDDNwUploadS6NoSync{time}---customers.csv | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |

  Scenario Outline: The NorthwindTest orders table aggregated with Aggregate Rows stays in its project, saved with Data sync <Sync> (GROK-19103)
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    When user picks "Data > Aggregate Rows..." from the top menu
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip shipname" area of pivot table viewer
    And user adds "customerid" to the "group by" row of first pivot table viewer
    And user clicks on the "remove pivot chip customerid" area of pivot table viewer
    And user picks "Aggregation > count" from the context menu of the "aggregate chip avg(orderid)" area of pivot table viewer
    And user closes the context menu
    Then the "group by" reading of pivot table viewer should be "customerid"
    And the "pivot" reading of pivot table viewer should be ""
    And the "aggregate" reading of pivot table viewer should be "count(orderid)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And table "orders aggregation" should have 89 rows
    And table "orders aggregation" should have columns "customerid, count(orderid)"
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # the dialog fills its table rows after it shows; a switch clicked before its row is drawn is missed
    And Data sync switch in "orders" project table in "Save project" dialog should be checked
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be checked
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "orders" project table in "Save project" dialog
    And user <switch> Data sync switch in "orders aggregation" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "orders" table of the "<Project>" project should be saved <saved>
    And the "orders aggregation" table of the "<Project>" project should be saved <saved>
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
    Then table "orders" should have been <how> with 830 rows
    And table "orders aggregation" should have been <how> with 89 rows
    And the table views "orders, orders aggregation" should be open
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders aggregation" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    Given the dashboards panel of the left sidebar is open
    Then <Project> tree node inside browse tree should be visible
    Given <Project> tree node inside browse tree is expanded
    Then <Project>---orders-aggregation tree node inside browse tree should be visible
    And <Project>---orders tree node inside browse tree should be visible
    And New-Dashboard tree node inside browse tree should be visible
    And there should be 2 visible dashboards project nodes
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                  | switch   | saved          | state | script  | how                   |
      | ON   | BDDNwUpload8Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDNwUpload8NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
