@serial @realizes:views.projects @realizes:data.menu.link-tables @realizes:file.menu.save.tables-as-project
Feature: Projects uploaded from different sources, with Data sync on and off
  Scenario 1, TestCase1: two tables opened from Browse > Files > Demo > northwind are linked
  selection-to-filter through Data > Link Tables..., the link is checked in the status bar of the
  second table, and both are saved in one project through the ribbon's Save dialog — once with Data
  sync on for both tables and once with it off. The project reopens from its Dashboards card with
  both tables (rebuilt by data sync, or loaded from the snapshot) and a link that still filters; the
  Save dialog of the reopened project shows each table's Creation script only when it was saved
  with Data sync. Scenario 2, TestCase7: orders.csv with a pivot table (Group by ShipCountry, Pivot
  ShipVia, count(OrderID)) published with ADD is saved with Data sync on and off, and reopens with
  both tables. Translated from the TestTrack case Projects/uploading; TestCase4 and 5 (files from a
  space) are projects-uploading-space.

  Fixture numbers: customers.csv has 91 rows and orders.csv 830; the first two customers (ALFKI,
  ANATR) have 6 and 4 orders, the next two (ANTON, AROUT) 7 and 13; orders has 40 ship countries.
  The link is checked with the first pair before the save and with the second pair after the
  reopen, so a filter the project merely restored cannot pass for a link that no longer works.
  Before a reopened table's Creation script is claimed shown or hidden, its Data sync switch is
  claimed on or off, so a dialog that has not filled its rows cannot pass for "hidden". The pivot's
  configuration is set through its properties: it prepares the table the case is about.

  Parked (see the request document): TestCase2, 3 and 6 (the NorthwindTest query PostgresAll,
  dev only; the Link Tables phrases are not in the library the General/bdd project sees) and
  Scenario 3, TestCase8 (GROK-19103: Aggregate Rows of the NorthwindTest orders table, whose row
  count must be remembered), with their System:Datagrok versions; the proof that a reopened table is
  not a frame that stayed open. Here the reopen is claimed by each table view's rows and its
  data-sync mark. After the reopen the link is checked with customers' rows 3 and 4 (20 orders), not
  the md table's 1 and 2: a filter the project merely restored would pass 10 whatever the link does
  (the md's own note does this for TestCase2).

  The project (with its tables and views) is named with the run's time and removed before its row
  starts and again when the feature ends. The console errors are not claimed: the Save dialog's
  preview logs "Unable to find element in cloned iframe" for these views (GROK-18606, won't fix),
  and a check that lets only that message through is requested in the request document. It is
  serial: the Dashboards search and the uploads are shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open

  Scenario Outline: Two northwind files from Files > Demo > northwind, linked and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And table "customers" should have 91 rows
    When user clicks on browse tab
    And user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user sets the tables of the Link Tables dialog to "customers" and "orders"
    And user sets key columns 1 of the Link Tables dialog to "CustomerID" and "CustomerID"
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "customers -> orders"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 10 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 10"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "customers" project table in "Save project" dialog
    And user <switch> Data sync switch in "orders" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "customers" table of the "<Project>" project should be saved <saved>
    And the "orders" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then the "customers" table view should open with 91 rows
    And the "orders" table view should open with 830 rows
    Given user switches to the "customers" table view
    Then the table should have been <how>
    Given user switches to the "orders" table view
    Then the table should have been <how>
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 20 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 20"
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | saved          | state | script  | how               |
      | ON   | BDDUpload1Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDUpload1NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: orders.csv with a pivot table added to the workspace, saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And user clears the saved pivot table parameters
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user sets properties of pivot table viewer:
      | Group By Column Names  | ShipCountry |
      | Pivot Column Names     | ShipVia     |
      | Aggregate Agg Types    | count       |
      | Aggregate Column Names | OrderID     |
    Then the "group by" reading of pivot table viewer should be "ShipCountry"
    And the "pivot" reading of pivot table viewer should be "ShipVia"
    And the "aggregate" reading of pivot table viewer should be "count(OrderID)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And the table should have 40 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "orders" project table in "Save project" dialog
    And user <switch> Data sync switch in "orders aggregation" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And the "orders" table of the "<Project>" project should be saved <saved>
    And the "orders aggregation" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on <Project> gallery card
    Then the "orders" table view should open with 830 rows
    And the "orders aggregation" table view should open with 40 rows
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders aggregation" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

    Examples:
      | Sync | Project                | switch   | saved          | state | script  |
      | ON   | BDDUpload7Sync{time}   | checks   | with data sync | on    | visible |
      | OFF  | BDDUpload7NoSync{time} | unchecks | as a snapshot  | off   | hidden  |
