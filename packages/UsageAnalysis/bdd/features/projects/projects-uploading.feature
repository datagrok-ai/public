@serial @realizes:views.projects @realizes:data.menu.link-tables @realizes:data.menu.aggregate-rows @realizes:file.menu.save.tables-as-project @realizes:viewers.pivot-viewer
Feature: Projects uploaded from different sources, with Data sync on and off
  Two tables opened from Browse are linked selection-to-filter through Data > Link Tables..., the
  link is checked in the status bar of the second table, and both are saved in one project through
  the ribbon's Save dialog — once with Data sync on for both tables and once with it off. The
  project reopens from its Dashboards card with both tables (rebuilt by data sync, or loaded from
  the snapshot) and a link that still filters; the Save dialog of the reopened project shows each
  table's Creation script only when it was saved with Data sync. Then a pivot table published into
  the workspace with ADD, and an aggregation from Data > Aggregate Rows... over a database table
  (GROK-19103), are saved and reopened inside their project. Every source is one Scenario Outline
  row that makes and removes its own fixtures, so no row depends on another. Translated from the
  TestTrack case Projects/uploading.

  The md's source matrix, and what stands in for it:
  - TestCase1 (two files from Files > Demo > northwind), TestCase4 (both from a Space), TestCase5
    (customers from a Space, orders from Files): as written, the files opened by double-clicking
    their Browse tree nodes (a space lists its files under its node). The Space is made per row
    through the JS API, with copies of the northwind files in its Files; the Create Space dialog
    and the drag with Copy of the md's setup are claimed in projects-lifecycle-spaces.
  - TestCase2 (the NorthwindTest query PostgresAll twice): NorthwindTest exists only on dev, so this
    row opens System:Datagrok's public.entity_types twice with Get All and links the two by name;
    two rows selected in the first leave exactly two in the second, three leave three. The table's
    row count differs per server: it is remembered after the first open and compared on reopen.
  - TestCase3 and TestCase6 (customers.csv from Files > Demo > northwind or from a Space, linked to
    the PostgresAll query by the customer key): the query is the row's own, saved on System:Datagrok
    through the JS API and run by double-clicking it under the connection. It selects literal
    customer IDs (ALFKI once, ANATR twice, ANTON three times, AROUT four times; no table is read), so
    customers rows 1-2 leave 3 of its 10 rows and rows 3-4 leave 7. The Link Tables dialog lists
    the new link with the table names cut to 15 characters ("customers -> BDDUploadQ3NoSy..."),
    so the listed link is claimed by that prefix. Every row makes a Space holding
    customers.csv; the Files rows do not open it.
  - TestCase7 (orders.csv with a pivot): as written. The pivot's aggregated values are not compared
    against the table: the pivot table reports only the rows its grid shows (ten), and this pivot
    has 40 — Demo's orders.csv reads 40 distinct ShipCountry values, because a shifted address
    field puts postal codes into the column in some rows. Its readings (group by, pivot, aggregate)
    and the published table's rows and columns are claimed instead.
  - TestCase8 (an aggregation of a database table, GROK-19103): over public.entity_types, by
    is_package_entity with count(id), instead of the NorthwindTest orders. After the reopen the
    Dashboards panel lists New Dashboard and the project, with both tables under the project and no
    separate project for the aggregation.

  Fixture numbers: customers.csv has 91 rows and orders.csv 830; the first two customers (ALFKI,
  ANATR) have 6 and 4 orders, the next two (ANTON, AROUT) 7 and 13. The link is checked with the
  first pair before the save and with the second pair after the reopen, so a filter the project
  merely restored cannot pass for a link that no longer works.

  A reopen is claimed to be a reopen: the tables open before the save are marked in memory, Close
  All is claimed to leave no table, and every table that comes back must be a frame without the
  mark. Before a reopened table's Creation script is claimed shown or hidden, its Data sync switch
  is claimed on or off, so a dialog that has not filled its rows cannot pass for "hidden". The
  scenarios with a derived table claim both table views open after the reopen, not only the tables.

  Every project (with its tables and views), every space and every query is named with the run's
  time and is removed before its row starts and again when the feature ends. The Save dialog's preview logs
  "Unable to find element in cloned iframe" (GROK-18606, known noise), so the save is claimed to
  log nothing else; the reopen is claimed to log nothing at all. @serial: the Dashboards search and
  the uploads are shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And user clears the saved pivot table parameters

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
    When user selects "customers" in link tables first table
    And user selects "orders" in link tables second table
    And user selects "CustomerID" in link tables first key
    And user selects "CustomerID" in link tables second key
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
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
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
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "customers" should have been <how> with 91 rows
    And table "orders" should have been <how> with 830 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 20 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 20"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | saved          | state | script  | how               |
      | ON   | BDDUpload1Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDUpload1NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: customers.csv from a space and orders.csv from <Source>, linked and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And a space "<Space>" holds copies of the northwind files "customers.csv, orders.csv"
    When user clicks on "Refresh" icon inside browse toolbar
    And Spaces tree node inside browse tree is expanded
    And Spaces---<Space> tree node inside browse tree is expanded
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Spaces---<Space>---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And table "customers" should have 91 rows
    When user clicks on browse tab
    And user double-clicks <orders node> tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "customers" in link tables first table
    And user selects "orders" in link tables second table
    And user selects "CustomerID" in link tables first key
    And user selects "CustomerID" in link tables second key
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
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
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
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "customers" should have been <how> with 91 rows
    And table "orders" should have been <how> with 830 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 20 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 20"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Source                   | Sync | Project                | Space                   | orders node                                    | switch   | saved          | state | script  | how               |
      | the same space           | ON   | BDDUpload4Sync{time}   | BDDUploadS4Sync{time}   | Spaces---BDDUploadS4Sync{time}---orders.csv    | checks   | with data sync | on    | visible | reloaded by data sync |
      | the same space           | OFF  | BDDUpload4NoSync{time} | BDDUploadS4NoSync{time} | Spaces---BDDUploadS4NoSync{time}---orders.csv  | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |
      | Files > Demo > northwind | ON   | BDDUpload5Sync{time}   | BDDUploadS5Sync{time}   | Files---Demo---northwind---orders.csv          | checks   | with data sync | on    | visible | reloaded by data sync |
      | Files > Demo > northwind | OFF  | BDDUpload5NoSync{time} | BDDUploadS5NoSync{time} | Files---Demo---northwind---orders.csv          | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: A database table opened twice with Get All, linked by name and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    When user remembers the row count of table "entity_types" as "entity types"
    And user clicks on browse tab
    And user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types (2)" view should be current
    And table "entity_types (2)" should be open
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "entity_types" in link tables first table
    And user selects "entity_types (2)" in link tables second table
    And user selects "name" in link tables first key
    And user selects "name" in link tables second key
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "entity_types -> entity_types"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on "entity_types" view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "entity_types" should be selected
    When user clicks on "entity_types (2)" view
    Then 2 rows of table "entity_types (2)" should pass the filter
    And status bar should contain text "Filtered: 2"
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "entity_types" project table in "Save project" dialog
    And user <switch> Data sync switch in "entity_types (2)" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "entity_types" table of the "<Project>" project should be saved <saved>
    And the "entity_types (2)" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "entity_types" should have been <how> with the "entity types" row count
    And table "entity_types (2)" should have been <how> with the "entity types" row count
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "entity_types" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "entity_types" project table in "Save project" dialog should be <script>
    And Data sync switch in "entity_types (2)" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "entity_types (2)" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on "entity_types" view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 5" area of grid holding Shift
    Then 3 rows of table "entity_types" should be selected
    When user clicks on "entity_types (2)" view
    Then 3 rows of table "entity_types (2)" should pass the filter
    And status bar should contain text "Filtered: 3"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | saved          | state | script  | how               |
      | ON   | BDDUpload2Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDUpload2NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: customers.csv from <Source> and a database query of customer IDs, linked and saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And no query named "<Query>" is on the server
    And a query "<Query>" on the Datagrok connection is:
      """
      select customerid from (values ('ALFKI'), ('ANATR'), ('ANATR'), ('ANTON'), ('ANTON'), ('ANTON'), ('AROUT'), ('AROUT'), ('AROUT'), ('AROUT')) as t(customerid)
      """
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
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---<Query> tree node inside browse tree
    Then the "<Query>" view should be current
    And table "<Query>" should have 10 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "customers" in link tables first table
    And user selects "<Query>" in link tables second table
    And user selects "CustomerID" in link tables first key
    And user selects "customerid" in link tables second key
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "customers -> <Query shown>"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on "<Query>" view
    Then 3 rows of table "<Query>" should pass the filter
    And status bar should contain text "Filtered: 3"
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "customers" project table in "Save project" dialog
    And user <switch> Data sync switch in "<Query>" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "customers" table of the "<Project>" project should be saved <saved>
    And the "<Query>" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "customers" should have been <how> with 91 rows
    And table "<Query>" should have been <how> with 10 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "<Query>" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "<Query>" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on "<Query>" view
    Then 7 rows of table "<Query>" should pass the filter
    And status bar should contain text "Filtered: 7"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Source                   | Sync | Project                | Query                   | Query shown     | Space                   | customers node                                   | switch   | saved          | state | script  | how               |
      | Files > Demo > northwind | ON   | BDDUpload3Sync{time}   | BDDUploadQ3Sync{time}   | BDDUploadQ3Sync | BDDUploadS3Sync{time}   | Files---Demo---northwind---customers.csv         | checks   | with data sync | on    | visible | reloaded by data sync |
      | Files > Demo > northwind | OFF  | BDDUpload3NoSync{time} | BDDUploadQ3NoSync{time} | BDDUploadQ3NoSy | BDDUploadS3NoSync{time} | Files---Demo---northwind---customers.csv         | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |
      | a space                  | ON   | BDDUpload6Sync{time}   | BDDUploadQ6Sync{time}   | BDDUploadQ6Sync | BDDUploadS6Sync{time}   | Spaces---BDDUploadS6Sync{time}---customers.csv   | checks   | with data sync | on    | visible | reloaded by data sync |
      | a space                  | OFF  | BDDUpload6NoSync{time} | BDDUploadQ6NoSync{time} | BDDUploadQ6NoSy | BDDUploadS6NoSync{time} | Spaces---BDDUploadS6NoSync{time}---customers.csv | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: orders.csv with a pivot table added to the workspace, saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip ShipName" area of pivot table viewer
    And user adds "ShipCountry" to the "group by" row of pivot table viewer
    And user clicks on the "remove pivot chip CustomerID" area of pivot table viewer
    And user adds "ShipVia" to the "pivot" row of pivot table viewer
    And user picks "Aggregation > count" from the context menu of the "aggregate chip avg(OrderID)" area of pivot table viewer
    And user closes the context menu
    Then the "group by" reading of pivot table viewer should be "ShipCountry"
    And the "pivot" reading of pivot table viewer should be "ShipVia"
    And the "aggregate" reading of pivot table viewer should be "count(OrderID)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And table "orders aggregation" should have 40 rows
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
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
    And user double-clicks on <Project> gallery card
    Then table "orders" should have been <how> with 830 rows
    And table "orders aggregation" should have been <how> with 40 rows
    And table "orders aggregation" should have columns "ShipCountry, 1 count(OrderID), 2 count(OrderID), 3 count(OrderID)"
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
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | saved          | state | script  | how               |
      | ON   | BDDUpload7Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDUpload7NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |

  Scenario Outline: A database table aggregated with Aggregate Rows stays in its project, saved with Data sync <Sync> (GROK-19103)
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    When user remembers the row count of table "entity_types" as "entity types"
    And user picks "Data > Aggregate Rows..." from the top menu
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip id" area of pivot table viewer
    And user adds "is_package_entity" to the "group by" row of first pivot table viewer
    And user adds "id" to the "aggregate" row of first pivot table viewer
    Then the "group by" reading of pivot table viewer should be "is_package_entity"
    When user picks "Aggregation > count" from the context menu of the "aggregate chip values(id)" area of pivot table viewer
    And user closes the context menu
    Then the "aggregate" reading of pivot table viewer should be "count(id)"
    And the aggregated values of pivot table viewer should match "count(id)" grouped by "is_package_entity"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "entity_types aggregation" view should be current
    And table "entity_types aggregation" should have 2 rows
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "entity_types" project table in "Save project" dialog
    And user <switch> Data sync switch in "entity_types aggregation" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "entity_types" table of the "<Project>" project should be saved <saved>
    And the "entity_types aggregation" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "entity_types" should have been <how> with the "entity types" row count
    And table "entity_types aggregation" should have been <how> with 2 rows
    And the table views "entity_types, entity_types aggregation" should be open
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "entity_types" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "entity_types" project table in "Save project" dialog should be <script>
    And Data sync switch in "entity_types aggregation" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "entity_types aggregation" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors but the project preview's should have been logged
    Given the dashboards panel of the left sidebar is open
    Then <Project> tree node inside browse tree should be visible
    Given <Project> tree node inside browse tree is expanded
    Then <Project>---entity-types-aggregation tree node inside browse tree should be visible
    And <Project>---entity-types tree node inside browse tree should be visible
    And New-Dashboard tree node inside browse tree should be visible
    And there should be 2 visible dashboards project nodes
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | saved          | state | script  | how               |
      | ON   | BDDUpload8Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | OFF  | BDDUpload8NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot |
