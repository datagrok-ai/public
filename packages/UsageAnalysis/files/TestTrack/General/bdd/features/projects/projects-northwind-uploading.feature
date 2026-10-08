@dev-only @serial @realizes:views.projects @realizes:data.menu.link-tables @realizes:file.menu.save.tables-as-project
Feature: Projects uploaded from the NorthwindTest query PostgresAll, linked and saved with Data sync on and off
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) and its query PostgresAll exist only
  on dev.datagrok.ai; on any other stand the tree has no such node and every row fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  TestCase2 (PostgresAll run twice) and TestCase3 (customers.csv from Browse > Files > Demo >
  northwind and PostgresAll) of the TestTrack case Projects/uploading: the two tables are linked
  selection-to-filter on the customer key through Data > Link Tables..., the link is checked in the
  status bar of the second table, and both are saved in one project through the ribbon's Save dialog —
  once with Data sync on for both tables and once with it off. The project reopens from its Dashboards
  card with both tables (re-run by data sync, or loaded from the snapshot) and a link that still
  filters; the reopened project's Save dialog shows each table's Creation script only when it was
  saved with Data sync (each table's switch claimed first, so a dialog that has not filled its rows
  cannot pass for "hidden"). PostgresAll is run with Run from its context menu, as the md says; its
  830 rows are the orders. After the reopen the link is checked with rows 3 and 4 of the first table
  (24 orders for TestCase2, as the md does; 20 for TestCase3, the md's 1 and 2 giving 10 again): a
  filter the project merely restored cannot pass for a link that still works.

  Kept without (see the request document): the proof that a reopened table is not a frame that
  stayed open (here the reopen is claimed by each table view's rows and its data-sync mark), and the
  console after the save (the Save dialog's preview logs "Unable to find element in cloned iframe"
  for views of several tables, GROK-18606, won't fix). Parked: TestCase6 (customers.csv from a
  space) and TestCase8 (GROK-19103, Aggregate Rows of the orders table, whose row count has to be
  remembered and whose Dashboards panel nodes have to be counted).

  The project (with its tables and views) is named with the run's time and removed before its row
  starts and again when the feature ends. It is serial: the Dashboards search and the uploads are
  shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open

  Scenario Outline: <First view> and <Second view>, linked and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user <open first>
    Then the "<First view>" view should be current
    And table "<First view>" should have <First rows> rows
    When user clicks on browse tab
    And user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the "<Second view>" view should be current
    And table "<Second view>" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user sets the tables of the Link Tables dialog to "<First view>" and "<Second view>"
    And user sets key columns 1 of the Link Tables dialog to "<First key>" and "customerid"
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    And user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on <First view> view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "<First view>" should be selected
    When user clicks on <Second view> view
    Then <Filtered> rows of table "<Second view>" should pass the filter
    And status bar should contain text "Filtered: <Filtered>"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "<First view>" project table in "Save project" dialog
    And user <switch> Data sync switch in "<Second view>" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And the "<First view>" table of the "<Project>" project should be saved <saved>
    And the "<Second view>" table of the "<Project>" project should be saved <saved>
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on <Project> gallery card
    Then the "<First view>" table view should open with <First rows> rows
    And the "<Second view>" table view should open with 830 rows
    Given user switches to the "<First view>" table view
    Then the table should have been <how>
    Given user switches to the "<Second view>" table view
    Then the table should have been <how>
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "<First view>" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "<First view>" project table in "Save project" dialog should be <script>
    And Data sync switch in "<Second view>" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "<Second view>" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user clicks on <First view> view
    And user clicks on the "row header <Row A>" area of grid
    And user clicks on the "row header <Row B>" area of grid holding Shift
    Then 2 rows of table "<First view>" should be selected
    When user clicks on <Second view> view
    Then <Filtered after> rows of table "<Second view>" should pass the filter
    And status bar should contain text "Filtered: <Filtered after>"

    Examples:
      | open first                                                                                                            | First view  | First rows | First key  | Second view     | Filtered | Row A | Row B | Filtered after | Sync | Project              | switch   | saved          | state | script  | how                   |
      | picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree | PostgresAll | 830        | customerid | PostgresAll (2) | 11       | 3     | 4     | 24             | ON   | BDDNwUp2Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree | PostgresAll | 830        | customerid | PostgresAll (2) | 11       | 3     | 4     | 24             | OFF  | BDDNwUp2NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
      | double-clicks Files---Demo---northwind---customers.csv tree node inside browse tree                                  | customers   | 91         | CustomerID | PostgresAll     | 10       | 3     | 4     | 20             | ON   | BDDNwUp3Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | double-clicks Files---Demo---northwind---customers.csv tree node inside browse tree                                  | customers   | 91         | CustomerID | PostgresAll     | 10       | 3     | 4     | 20             | OFF  | BDDNwUp3NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
