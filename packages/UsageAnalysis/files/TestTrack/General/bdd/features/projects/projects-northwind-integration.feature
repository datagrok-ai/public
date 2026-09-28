@dev-only @journey @serial @realizes:views.projects @realizes:viewers.pivot-viewer @realizes:data.menu.aggregate-rows @realizes:data.menu.join-tables @realizes:views.space @realizes:views.databases
Feature: One project of tables from nine sources, with the NorthwindTest table and query
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the sources scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  The md as written: nine tables are opened in one workspace, each the way a user gets it — orders.csv
  from Browse > Files > Demo > northwind (830 rows), a copy of customers.csv dragged into a Space the
  feature creates (91 rows), NorthwindTest's products table through Get All (77 rows), the
  NorthwindTest query PostgresAll run from the Browse tree (830 rows), a JavaScript script run from
  the Scripts view (demog, 5850 rows), a pivot of orders added from the toolbox and published with
  ADD, an aggregation of demog from Data > Aggregate Rows..., a join of orders and customers from
  Data > Join Tables..., and a clone of products from the table's Context Panel > Actions > Clone.
  The Dashboards panel of the left sidebar lists all nine under New Dashboard and no other project;
  the ribbon's Save saves them in one project, every table with its own Data sync switch on, and the
  server holds all nine in that project; reopened from the Dashboards gallery, the card's Content
  pane lists the nine and the project opens the nine tables, each rebuilt by data sync with the rows
  it had. Translated from the TestTrack case Projects/complex-integration; its System:Datagrok
  version is packages/UsageAnalysis/bdd/features/projects/projects-integration.feature. The md gives
  no row count for products: 77 is the count NorthwindTest's products holds.

  Fixture substituted: the script is made through the JS API rather than the script editor ("a
  script … is on the server"; projects-lifecycle-script makes its script through the UI). The Space
  is created and filled through the UI.

  The project, the Space and the script are named with the run's time and removed (the project with
  its tables and views, the Space with its file) at the start and at the end; Delete Project removes
  the project through the UI at the end. PostgresAll and the connection belong to the Dbtests
  package and are never changed. @serial: the Dashboards search is shared with every feature that
  saves a project.

  Not translated, and why: the md's cleanup of the Space and the script through their context menus
  (Delete Space, Delete) — the feature-end cleanup removes them and reads the server back, and the
  delete dialogs are the Spaces and Scripts features' subject; a clean console around the save (the
  publish preview logs "Unable to find element in cloned iframe", known noise with no ticket).

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDNwIntegr{time}" is on the server
    And no space named "BDDNwIntegr{time}" is on the server
    And a script "BDDNwIntegrScript{time}" is on the server:
      """
      //language: javascript
      //output: dataframe df
      df = await grok.data.getDemoTable('demog.csv');
      """
    And user clears the saved pivot table parameters

  Scenario: A Space gets a copy of customers.csv
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDDNwIntegr{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDNwIntegr{time}" should be on the server
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---northwind tree node inside browse tree
    Then customers.csv link in gallery should be visible
    When user drags customers.csv link in gallery to Spaces---BDDNwIntegr{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given Spaces---BDDNwIntegr{time} tree node inside browse tree is expanded
    Then Spaces---BDDNwIntegr{time}---customers.csv tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Two files, a database table, a query and a script open into the workspace
    Given Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And the table should have 830 rows
    Given the browse panel is open
    When user double-clicks on Spaces---BDDNwIntegr{time}---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And the table should have 91 rows
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree
    Then the "products" view should be current
    And table "products" should have 77 rows
    Given the browse panel is open
    When user double-clicks on Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree
    Then the "PostgresAll" view should be current
    And table "PostgresAll" should have 830 rows
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    And user clears gallery search
    And user types "BDDNwIntegrScript{time}" into gallery search
    And user picks "Run..." from the context menu of "BDDNwIntegrScript{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A pivot, an aggregation, a join and a clone are derived from them
    When user switches to the "orders" table view
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip ShipName" area of pivot table viewer
    And user adds "ShipCountry" to the "group by" row of first pivot table viewer
    And user clicks on the "remove pivot chip CustomerID" area of pivot table viewer
    And user picks "Aggregation > count" from the context menu of the "aggregate chip avg(OrderID)" area of pivot table viewer
    And user closes the context menu
    Then the "group by" reading of pivot table viewer should be "ShipCountry"
    And the "pivot" reading of pivot table viewer should be ""
    And the "aggregate" reading of pivot table viewer should be "count(OrderID)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    # one row per ShipCountry of orders.csv
    And table "orders aggregation" should have 40 rows
    And table "orders aggregation" should have columns "ShipCountry, count(OrderID)"
    When user switches to the "demog" table view
    And user picks "Data > Aggregate Rows..." from the top menu
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    And user adds "RACE" to the "group by" row of first pivot table viewer
    And user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    And the "pivot" reading of pivot table viewer should be ""
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "RACE"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    And table "demog aggregation" should have 4 rows
    And table "demog aggregation" should have columns "RACE, avg(AGE)"
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "orders" in join left table selector
    And user selects "customers" in join right table selector
    And user picks column "CustomerID" in join left key selector
    And user picks column "CustomerID" in join right key selector
    Then "Join Type" input in "Join Tables" dialog should have value "inner"
    When user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have 830 rows
    When user switches to the "products" table view
    Given the context panel is open
    When user clicks on status bar table name
    Then the context panel should show "products"
    When user expands "Actions" pane in context panel
    And user clicks on Clone action in context panel
    Then the "products (2)" view should be current
    And table "products (2)" should have 77 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Dashboards panel lists the nine tables under New Dashboard only
    Given the dashboards panel of the left sidebar is open
    Then the following elements should be visible:
      | New-Dashboard---orders tree node inside browse tree                |
      | New-Dashboard---customers tree node inside browse tree             |
      | New-Dashboard---products tree node inside browse tree          |
      | New-Dashboard---PostgresAll tree node inside browse tree      |
      | New-Dashboard---demog tree node inside browse tree                 |
      | New-Dashboard---orders-aggregation tree node inside browse tree    |
      | New-Dashboard---demog-aggregation tree node inside browse tree     |
      | New-Dashboard---result tree node inside browse tree                |
      | New-Dashboard---products-(2) tree node inside browse tree      |
    And there should be 1 visible dashboards project node
    And no errors should have been logged

  Scenario: Saved with Data sync for all nine tables
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # the dialog fills its table rows after it shows; the last row's Creation script is the barrier
    And "Creation script" button in "products (2)" project table in "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "orders" project table in "Save project" dialog             |
      | Data sync switch in "customers" project table in "Save project" dialog          |
      | Data sync switch in "products" project table in "Save project" dialog       |
      | Data sync switch in "PostgresAll" project table in "Save project" dialog   |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "orders aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog  |
      | Data sync switch in "result" project table in "Save project" dialog             |
      | Data sync switch in "products (2)" project table in "Save project" dialog   |
    When user enters "BDDNwIntegr{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDNwIntegr{time}\" uploaded" should have been shown
    And 1 project named "BDDNwIntegr{time}" should be on the server
    # the join is saved in this project, not in a separate one (GROK-19103)
    And the "BDDNwIntegr{time}" project on the server should hold the tables "orders, customers, products, PostgresAll, demog, orders aggregation, demog aggregation, result, products (2)"
    And the "orders" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "customers" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "products" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "PostgresAll" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "demog" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "orders aggregation" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "demog aggregation" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "result" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And the "products (2)" table of the "BDDNwIntegr{time}" project should be saved with data sync
    And "Share BDDNwIntegr{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDNwIntegr{time}" dialog
    Then the "Share BDDNwIntegr{time}" dialog should close

  Scenario: Reopened from Dashboards, all nine tables come back with their rows
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwIntegr{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given the context panel is open
    When user clicks on BDDNwIntegr{time} gallery card
    Then the context panel should show "BDDNwIntegr{time}"
    When user expands "Content" pane in context panel
    Then the following elements should be visible:
      | BDDNwIntegr{time}---orders tree node in context panel             |
      | BDDNwIntegr{time}---customers tree node in context panel          |
      | BDDNwIntegr{time}---products tree node in context panel       |
      | BDDNwIntegr{time}---PostgresAll tree node in context panel   |
      | BDDNwIntegr{time}---demog tree node in context panel              |
      | BDDNwIntegr{time}---orders-aggregation tree node in context panel |
      | BDDNwIntegr{time}---demog-aggregation tree node in context panel  |
      | BDDNwIntegr{time}---result tree node in context panel             |
      | BDDNwIntegr{time}---products-(2) tree node in context panel   |
    When user double-clicks on BDDNwIntegr{time} gallery card
    Then the table views "orders, customers, products, PostgresAll, demog, orders aggregation, demog aggregation, result, products (2)" should be open
    And table "orders" should have been reloaded by data sync with 830 rows
    And table "customers" should have been reloaded by data sync with 91 rows
    And table "products" should have been reloaded by data sync with 77 rows
    And table "PostgresAll" should have been reloaded by data sync with 830 rows
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "orders aggregation" should have been reloaded by data sync with 40 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "result" should have been reloaded by data sync with 830 rows
    And table "products (2)" should have been reloaded by data sync with 77 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwIntegr{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDNwIntegr{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDNwIntegr{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDNwIntegr{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDNwIntegr{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
