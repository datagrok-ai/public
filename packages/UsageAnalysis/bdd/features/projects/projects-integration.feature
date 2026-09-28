@journey @serial @realizes:views.projects @realizes:viewers.pivot-viewer @realizes:data.menu.aggregate-rows @realizes:data.menu.join-tables @realizes:views.space @realizes:views.databases
Feature: One project of tables from nine sources
  Nine tables are opened in one workspace, each the way a user gets it: a file from Browse > Files >
  Demo > northwind (orders.csv), a file from a Space (a copy of customers.csv dragged into a Space
  the feature creates), a database table through Get All, a saved query run from the Browse tree, a
  JavaScript script run from the Scripts view, a pivot of orders added from the toolbox and
  published with ADD, an aggregation of the script's demog from Data > Aggregate Rows..., a join of
  orders and customers from Data > Join Tables..., and a clone of the database table from the
  table's Context Panel > Actions > Clone. The Dashboards panel of the left sidebar lists all nine
  under New Dashboard and no other project; the ribbon's Save saves them in one project, every table
  with its own Data sync switch on, and the server holds all nine in that project; reopened from
  the Dashboards gallery, the card's Content pane
  lists the nine and the project opens the nine tables, each rebuilt by data sync with the rows it
  had. Translated from the TestTrack case Projects/complex-integration.

  Fixtures substituted: NorthwindTest (a dev-only connection) is replaced by System:Datagrok, which
  every user may read and query — the table products by public.entity_types (Get All), the query
  PostgresAll by a query of the feature's own, "BDDIntegrQ{time}" (select id, name from
  public.entity_types), saved through the JS API; their row counts are the ones read at the first
  open, since a table of the Datagrok database differs between servers. The script is made through
  the JS API rather than the script editor ("a script … is on the server"; projects-lifecycle-script
  makes its script through the UI). The Space is created and filled through the UI.

  Everything is named with the run's time: the project, the Space, the query, the script. The
  project (with its tables and views), the Space (with its file), the query and the script are
  removed at the start and at the end; Delete Project removes the project through the UI at the end.
  It is @serial: the Dashboards search is shared with every feature that saves a project.

  Not translated, and why: the md's cleanup of the Space and the script through their context menus
  (Delete Space, Delete) — the feature-end cleanup removes them and reads the server back, and the
  delete dialogs are the Spaces and Scripts features' subject; a clean console around the save (the
  publish preview logs "Unable to find element in cloned iframe", known noise with no ticket). Covered
  elsewhere: the derived tables shared, renamed and reopened by a second account in
  projects-lifecycle-derived; a table added to a saved project by a drag in projects-augment.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDIntegr{time}" is on the server
    And no space named "BDDIntegr{time}" is on the server
    And no query named "BDDIntegrQ{time}" is on the server
    And a query "BDDIntegrQ{time}" on the Datagrok connection is:
      """
      select id, name from public.entity_types
      """
    And a script "BDDIntegrScript{time}" is on the server:
      """
      //language: javascript
      //output: dataframe df
      df = await grok.data.getDemoTable('demog.csv');
      """
    And user clears the saved pivot table parameters

  Scenario: A Space gets a copy of customers.csv
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDDIntegr{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDIntegr{time}" should be on the server
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---northwind tree node inside browse tree
    Then customers.csv link in gallery should be visible
    When user drags customers.csv link in gallery to Spaces---BDDIntegr{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given Spaces---BDDIntegr{time} tree node inside browse tree is expanded
    Then Spaces---BDDIntegr{time}---customers.csv tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Two files, a database table, a query and a script open into the workspace
    Given Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And the table should have 830 rows
    Given the browse panel is open
    When user double-clicks on Spaces---BDDIntegr{time}---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And the table should have 91 rows
    Given the browse panel is open
    # the tree keeps the children it listed before the query was saved; a refresh lists it
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    When user remembers the row count of table "entity_types" as "entity types"
    Given the browse panel is open
    When user double-clicks on Databases---Postgres---Datagrok---BDDIntegrQ{time} tree node inside browse tree
    Then the "BDDIntegrQ{time}" view should be current
    When user remembers the row count of table "BDDIntegrQ{time}" as "integration query"
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    And user clears gallery search
    And user types "BDDIntegrScript{time}" into gallery search
    And user picks "Run..." from the context menu of "BDDIntegrScript{time}" link in gallery
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
    When user switches to the "entity_types" table view
    Given the context panel is open
    When user clicks on status bar table name
    Then the context panel should show "entity_types"
    When user expands "Actions" pane in context panel
    And user clicks on Clone action in context panel
    Then the "entity_types (2)" view should be current
    And table "entity_types (2)" should have the "entity types" row count
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Dashboards panel lists the nine tables under New Dashboard only
    Given the dashboards panel of the left sidebar is open
    Then the following elements should be visible:
      | New-Dashboard---orders tree node inside browse tree                |
      | New-Dashboard---customers tree node inside browse tree             |
      | New-Dashboard---entity_types tree node inside browse tree          |
      | New-Dashboard---BDDIntegrQ{time} tree node inside browse tree      |
      | New-Dashboard---demog tree node inside browse tree                 |
      | New-Dashboard---orders-aggregation tree node inside browse tree    |
      | New-Dashboard---demog-aggregation tree node inside browse tree     |
      | New-Dashboard---result tree node inside browse tree                |
      | New-Dashboard---entity_types-(2) tree node inside browse tree      |
    And there should be 1 visible dashboards project node
    And no errors should have been logged

  Scenario: Saved with Data sync for all nine tables
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # the dialog fills its table rows after it shows; the last row's Creation script is the barrier
    And "Creation script" button in "entity_types (2)" project table in "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "orders" project table in "Save project" dialog             |
      | Data sync switch in "customers" project table in "Save project" dialog          |
      | Data sync switch in "entity_types" project table in "Save project" dialog       |
      | Data sync switch in "BDDIntegrQ{time}" project table in "Save project" dialog   |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "orders aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog  |
      | Data sync switch in "result" project table in "Save project" dialog             |
      | Data sync switch in "entity_types (2)" project table in "Save project" dialog   |
    When user enters "BDDIntegr{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDIntegr{time}\" uploaded" should have been shown
    And 1 project named "BDDIntegr{time}" should be on the server
    # the join is saved in this project, not in a separate one (GROK-19103)
    And the "BDDIntegr{time}" project on the server should hold the tables "orders, customers, entity_types, BDDIntegrQ{time}, demog, orders aggregation, demog aggregation, result, entity_types (2)"
    And the "orders" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "customers" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "entity_types" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "BDDIntegrQ{time}" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "demog" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "orders aggregation" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "demog aggregation" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "result" table of the "BDDIntegr{time}" project should be saved with data sync
    And the "entity_types (2)" table of the "BDDIntegr{time}" project should be saved with data sync
    And "Share BDDIntegr{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDIntegr{time}" dialog
    Then the "Share BDDIntegr{time}" dialog should close

  Scenario: Reopened from Dashboards, all nine tables come back with their rows
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDIntegr{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given the context panel is open
    When user clicks on BDDIntegr{time} gallery card
    Then the context panel should show "BDDIntegr{time}"
    When user expands "Content" pane in context panel
    Then the following elements should be visible:
      | BDDIntegr{time}---orders tree node in context panel             |
      | BDDIntegr{time}---customers tree node in context panel          |
      | BDDIntegr{time}---entity_types tree node in context panel       |
      | BDDIntegr{time}---BDDIntegrQ{time} tree node in context panel   |
      | BDDIntegr{time}---demog tree node in context panel              |
      | BDDIntegr{time}---orders-aggregation tree node in context panel |
      | BDDIntegr{time}---demog-aggregation tree node in context panel  |
      | BDDIntegr{time}---result tree node in context panel             |
      | BDDIntegr{time}---entity_types-(2) tree node in context panel   |
    When user double-clicks on BDDIntegr{time} gallery card
    Then the table views "orders, customers, entity_types, BDDIntegrQ{time}, demog, orders aggregation, demog aggregation, result, entity_types (2)" should be open
    And table "orders" should have been reloaded by data sync with 830 rows
    And table "customers" should have been reloaded by data sync with 91 rows
    And table "entity_types" should have been reloaded by data sync with the "entity types" row count
    And table "BDDIntegrQ{time}" should have been reloaded by data sync with the "integration query" row count
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "orders aggregation" should have been reloaded by data sync with 40 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "result" should have been reloaded by data sync with 830 rows
    And table "entity_types (2)" should have been reloaded by data sync with the "entity types" row count
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDIntegr{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDIntegr{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDIntegr{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDIntegr{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDIntegr{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
