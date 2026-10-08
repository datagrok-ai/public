@journey @serial @realizes:views.projects @realizes:data.menu.join-tables @realizes:views.space
Feature: A saved project when what it is built from is renamed or moved
  Three parts, each with its own project. Part 1 (GROK-19212): orders.csv and customers.csv from
  Browse > Files > Demo > northwind, a pivot of orders published with ADD and an inner join of the
  two made with Data > Join Tables... are saved with Data sync; the two source tables are renamed
  from their view tabs (Table > Rename...), the project is saved again and reopened from its
  Dashboards card. Part 2 (GROK-21026): a query made with New SQL Query... on the Datagrok
  database's entity_types, its result and entity_types itself (Get All) are joined on id and saved;
  the query is renamed in the Browse tree and the project reopened. Part 3: a project made of a
  query's result and a script's result is moved into a space from its card, the query and the script
  are dragged into the space from the Browse tree and the Scripts view, and the project opens from
  the space. Translated from the TestTrack case Projects/complex.

  The pivot's configuration (Group by ShipCountry, count(OrderID)) is set through the viewer's
  properties: it prepares the table the part is about. Row counts of the northwind files are the
  md's (orders 830, customers 91); the pivot and the join were measured once (40 and 830 rows). The
  Datagrok database's tables differ per server, so its row counts are not claimed literally: the
  query's result is claimed to have rows.

  Kept without one claim each, restored once the library can remember a row count (see the request
  document): Part 2's "entity_types opens with the same row count as the query" and the reopened
  entity_types and query result having those counts.

  Known failures: GROK-19212 (reopened): the project whose pivot and join sources were renamed
  does not reopen whole (on localhost 1.28.0 no view of it opened within two minutes; the md reports
  the pivot and the join missing and "2 tasks" left in the task bar); GROK-21026: the join
  over the renamed query does not come back ("Could not resolve table", or "Opening project" left
  in the task bar on dev). Each is one scenario claiming what the md expects.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The three projects, both queries (under both names of the renamed one), the script and the
  space are removed when the feature starts and ends. It is serial: the Dashboards search and the
  uploads are shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And user clears the saved pivot table parameters
    And no project named "BDDComplexRename{time}" is on the server
    And no project named "BDDComplexJoin{time}" is on the server
    And no project named "BDDComplexMove{time}" is on the server
    And no query named "BDDComplexQuery{time}, BDDComplexQueryRenamed{time}, BDDComplexMoveQuery{time}" is on the server
    And no script named "BDDComplexScript{time}" is on the server
    And no space named "BDDComplexSpace{time}" is on the server

  Scenario: Part 1 — two northwind files, a pivot and a join are saved with Data sync
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    When user clicks on browse tab
    And user double-clicks Files---Demo---northwind---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And table "customers" should have 91 rows
    When user clicks on orders view
    Then the "orders" view should be current
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user sets "Group By Column Names" property of pivot table viewer to "ShipCountry"
    And user sets properties of pivot table viewer:
      | Aggregate Agg Types    | count   |
      | Aggregate Column Names | OrderID |
    Then the "group by" reading of pivot table viewer should be "ShipCountry"
    And the "aggregate" reading of pivot table viewer should be "count(OrderID)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And table "orders aggregation" should have 40 rows
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "orders" in Tables input in "Join Tables" dialog
    And user selects "customers" in second choice input in "Join Tables" dialog
    # the preview of the result opens once both tables are chosen, and takes the focus
    Then "Preview result columns" dialog should be visible
    When user selects "CustomerID" in first column input in "Join Tables" dialog
    And user selects "CustomerID" in second column input in "Join Tables" dialog
    And user selects "inner" in "Join Type" input in "Join Tables" dialog
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the current view should be a TableView view
    And the table should have 830 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "orders, customers, orders aggregation, result" with data sync
    When user enters "BDDComplexRename{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDComplexRename{time}" uploaded' should have been shown
    And "Share BDDComplexRename{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplexRename{time}" dialog
    Then the "Share BDDComplexRename{time}" dialog should close

  Scenario: Part 1 — the source tables are renamed from their view tabs and saved again
    When user picks "Table > Rename..." from the context menu of orders view
    Then "Rename table" dialog should be visible
    When user enters "ordersR" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And ordersR view should be visible
    When user picks "Table > Rename..." from the context menu of customers view
    Then "Rename table" dialog should be visible
    When user enters "customersR" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And customersR view should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "ordersR, customersR, orders aggregation, result" with data sync
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDComplexRename{time}" uploaded' should have been shown

  @known-failure @realizes:GROK-19212
  Scenario: Part 1 — reopened, the renamed sources, the pivot and the join come back
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRename{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRename{time} gallery card
    Then the "ordersR" table view should open with 830 rows
    And the "customersR" table view should open with 91 rows
    And the "orders aggregation" table view should open with 40 rows
    And the "result" table view should open with 830 rows
    And no error or warning balloon should have been shown
    And "Data loading error" dialog should be absent
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "ordersR, customersR, orders aggregation, result" with data sync
    And "Creation script" button in "ordersR" project table in "Save project" dialog should be visible
    And "Creation script" button in "customersR" project table in "Save project" dialog should be visible
    And "Creation script" button in "orders aggregation" project table in "Save project" dialog should be visible
    And "Creation script" button in "result" project table in "Save project" dialog should be visible
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: Part 2 — a query is made on entity_types with New SQL Query...
    When user presses Escape
    And user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDComplexQuery{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDComplexQuery{time}" should be on the server
    When user closes the current view
    Then no errors should have been logged

  Scenario: Part 2 — the query's result and entity_types are joined and saved
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexQuery{time} tree node inside browse tree
    Then the "BDDComplexQuery{time}" view should be current
    And the "rows" reading of grid should be at least 1
    When user clicks on browse tab
    Given Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    And the "rows" reading of grid should be at least 1
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "entity_types" in Tables input in "Join Tables" dialog
    And user selects "BDDComplexQuery{time}" in second choice input in "Join Tables" dialog
    # the preview of the result opens once both tables are chosen, and takes the focus
    Then "Preview result columns" dialog should be visible
    When user selects "id" in first column input in "Join Tables" dialog
    And user selects "id" in second column input in "Join Tables" dialog
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the current view should be a TableView view
    And the "rows" reading of grid should be at least 1
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "BDDComplexQuery{time}, entity_types, result" with data sync
    When user enters "BDDComplexJoin{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDComplexJoin{time}" uploaded' should have been shown
    And "Share BDDComplexJoin{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplexJoin{time}" dialog
    Then the "Share BDDComplexJoin{time}" dialog should close

  Scenario: Part 2 — the query is renamed in the Browse tree
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    # the double-click that ran the query unfolded its node over the runs it logged
    When user collapses Databases---Postgres---Datagrok---BDDComplexQuery{time} tree node inside browse tree
    And user picks "Rename..." from the context menu of Databases---Postgres---Datagrok---BDDComplexQuery{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "BDDComplexQueryRenamed{time}" into Name input in "Rename dataquery" dialog
    And user clicks on OK button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDComplexQueryRenamed{time}" should be on the server
    And 0 queries named "BDDComplexQuery{time}" should be on the server

  @known-failure @realizes:GROK-21026
  Scenario: Part 2 — reopened, the join over the renamed query comes back
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexJoin{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexJoin{time} gallery card
    Then table "result" should be open
    Given user switches to the "result" table view
    Then the "rows" reading of grid should be at least 1
    Given user switches to the "entity_types" table view
    Then the "rows" reading of grid should be at least 1
    Given user switches to the "BDDComplexQuery{time}" table view
    Then the "rows" reading of grid should be at least 1
    And "Data loading error" dialog should be absent

  Scenario: Part 3 — a space, a query and a script are made
    When user presses Escape
    And user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDComplexSpace{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDComplexSpace{time}" should be on the server
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDComplexMoveQuery{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDComplexMoveQuery{time}" should be on the server
    When user closes the current view
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    When user replaces the code of code editor with "//name: BDDComplexScript{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = await grok.data.getDemoTable('demog.csv');" to code editor
    And user saves the script
    Then 1 script named "BDDComplexScript{time}" should be on the server
    And no errors should have been logged

  Scenario: Part 3 — the query's and the script's results are saved as a project
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexMoveQuery{time} tree node inside browse tree
    Then the "BDDComplexMoveQuery{time}" view should be current
    And the "rows" reading of grid should be at least 1
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDComplexScript{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Run..." from the context menu of "BDDComplexScript{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "BDDComplexMoveQuery{time}, demog" with data sync
    When user enters "BDDComplexMove{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDComplexMove{time}" uploaded' should have been shown
    And "Share BDDComplexMove{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplexMove{time}" dialog
    Then the "Share BDDComplexMove{time}" dialog should close

  Scenario: Part 3 — the project, the query and the script are moved into the space
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexMove{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Move to Space..." from the context menu of BDDComplexMove{time} gallery card
    Then "Move to space" dialog should be visible
    When user selects "BDDComplexSpace{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    Given Spaces tree node inside browse tree is expanded
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user collapses Databases---Postgres---Datagrok---BDDComplexMoveQuery{time} tree node inside browse tree
    And user drags Databases---Postgres---Datagrok---BDDComplexMoveQuery{time} tree node inside browse tree to BDDComplexSpace{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDComplexScript{time}" into gallery search
    And user drags "BDDComplexScript{time}" link in gallery to BDDComplexSpace{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user clears gallery search
    And user double-clicks on BDDComplexSpace{time} tree node inside browse tree
    Then the "BDDComplexSpace{time}" view should be current
    And "BDDComplexMove{time}" link in gallery should be visible
    And "BDDComplexMoveQuery{time}" link in gallery should be visible
    And "BDDComplexScript{time}" link in gallery should be visible

  Scenario: Part 3 — the moved project opens from the space with its data
    When user double-clicks on "BDDComplexMove{time}" link in gallery
    Then the "demog" table view should open with 5850 rows
    Given user switches to the "BDDComplexMoveQuery{time}" table view
    Then the "rows" reading of grid should be at least 1
    And "Data loading error" dialog should be absent
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "BDDComplexMoveQuery{time}, demog" with data sync
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And "Creation script" button in "BDDComplexMoveQuery{time}" project table in "Save project" dialog should be visible
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no error or warning balloon should have been shown

  Scenario: Cleanup through the UI — the projects, the queries, the script and the space
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplex" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDComplexRename{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user picks "Delete Project" from the context menu of BDDComplexJoin{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDComplexRename{time}" should be on the server
    And 0 projects named "BDDComplexJoin{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDComplexSpace{time} tree node inside browse tree
    Then the "BDDComplexSpace{time}" view should be current
    When user picks "Delete Project" from the context menu of "BDDComplexMove{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDComplexMove{time}" should be on the server
    When user picks "Delete" from the context menu of "BDDComplexMoveQuery{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDDComplexMoveQuery{time}" should be on the server
    When user picks "Delete" from the context menu of "BDDComplexScript{time}" link in gallery
    And user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 scripts named "BDDComplexScript{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user collapses Databases---Postgres---Datagrok---BDDComplexQueryRenamed{time} tree node inside browse tree
    And user picks "Delete" from the context menu of Databases---Postgres---Datagrok---BDDComplexQueryRenamed{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDDComplexQueryRenamed{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    When user picks "Delete Space" from the context menu of BDDComplexSpace{time} tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete space "BDDComplexSpace{time}"?'
    And "Are you sure?" dialog should contain text "This will delete space and its related data"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 spaces named "BDDComplexSpace{time}" should be on the server
    And no errors should have been logged
