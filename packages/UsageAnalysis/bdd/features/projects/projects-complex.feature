@journey @serial @realizes:views.projects @realizes:data.menu.join-tables @realizes:views.space
Feature: A saved project when what it is built from is renamed or moved
  Three parts, each with its own project. Part 1 (core#188): orders.csv and customers.csv from
  Browse > Files > Demo > northwind, a pivot of orders published with ADD and an inner join of the
  two made with Data > Join Tables... are saved with Data sync; the two source tables are renamed
  from their view tabs (Table > Rename...), the project is saved again and reopened from its
  Dashboards card. Part 2 (core#628): a query on the Datagrok database (put on the server through the API)
  over its entity_types, its result and entity_types itself (Get All) are joined on id and saved;
  the query is renamed in the Browse tree and the project reopened. Part 3: a project made of a
  query's result and a script's result is moved into a space from its card, the query and the script
  are dragged into the space from the Browse tree and the Scripts view, and the project opens from
  the space. Translated from the TestTrack case Projects/complex.

  The pivot is configured through its own rows and chips: ShipCountry added to Group by and the
  default ShipName chip removed, the default avg(OrderID) chip turned into count(OrderID) from its
  menu. Row counts of the northwind files are the md's (orders 830, customers 91); the pivot and the
  join were measured once (40 and 830 rows). The Datagrok database's tables differ per server, so its
  row counts are not claimed literally: the query's result is claimed to have rows.

  Kept without one claim each, restored once the library can remember a row count (see the request
  document): Part 2's "entity_types opens with the same row count as the query" and the reopened
  entity_types and query result having those counts.

  core#188 (open): the project whose pivot and join sources were renamed did not
  reopen whole (localhost 1.28.0 in September: no view within two minutes; the md reports the pivot
  and the join missing). On localhost f98337ae1b (2026-10-12, two runs, also with the feature as
  committed) it reopens with all four tables, so its scenario is a plain one now. Known failure
  core#628: the join over the renamed query does not come back ("Could not resolve table", or
  "Opening project" left in the task bar on dev); one scenario claims what the md expects. The Part 1
  project is then shared with the second account (the library's sharing user), who opens it with its
  pivot and join — the recipient's open of derived tables of projects-lifecycle-derived.md step 8,
  which is parked as a whole.

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The three projects, both queries (under both names of the renamed one), the script and the
  space are removed when the feature starts and ends.

  Serial (`@serial`): Part 3 drags a query and a script onto a space in the Browse tree, whose
  Spaces and Datagrok queries every feature of the account adds to and removes from; a node moved by
  another feature during the drag takes the drop (seen in parallel runs, 2026-10-19).

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
    And a query "BDDComplexQuery{time}" on "System:Datagrok" reads "select * from public.entity_types"
    And a query "BDDComplexMoveQuery{time}" on "System:Datagrok" reads "select * from public.entity_types"
    And a script "BDDComplexScript{time}" is on the server:
      """
      //language: javascript
      //output: dataframe df
      df = await grok.data.getDemoTable('demog.csv');
      """
    And a space named "BDDComplexSpace{time}" is on the server

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
    # a new pivot of orders groups by ShipName, pivots on CustomerID and aggregates avg(OrderID)
    And the "group by" reading of pivot table viewer should be "ShipName"
    When user adds "ShipCountry" to the "group by" row of pivot table viewer
    And user clicks on the "remove group by chip ShipName" area of pivot table viewer
    And user picks "Aggregation > count" from the context menu of the "aggregate chip avg(OrderID)" area of pivot table viewer
    And user closes the context menu
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

  @realizes:core#188
  Scenario: Part 1 — reopened, the renamed sources, the pivot and the join come back
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRename{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRename{time} gallery card and the project opens
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

  Scenario: Part 1 — shared, the project's derived tables open for the second account
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRename{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDComplexRename{time} gallery card
    Then "Share BDDComplexRename{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDComplexRename{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDComplexRename{time}" dialog
    And user clicks on OK button in "Share BDDComplexRename{time}" dialog
    Then the "Share BDDComplexRename{time}" dialog should close
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRename{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRename{time} gallery card and the project opens
    Then the "ordersR" table view should open with 830 rows
    And the "customersR" table view should open with 91 rows
    And the "orders aggregation" table view should open with 40 rows
    And the "result" table view should open with 830 rows
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: Part 2 — the query's result and entity_types are joined and saved
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexQuery{time} tree node inside browse tree and the query it runs completes
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

  @known-failure:253 @realizes:core#628
  Scenario: Part 2 — reopened, the join over the renamed query comes back
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexJoin{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexJoin{time} gallery card
    Then table "result" should be open
    When user clicks on result view
    Then the "rows" reading of grid should be at least 1
    When user clicks on entity_types view
    Then the "rows" reading of grid should be at least 1
    When user clicks on BDDComplexQuery{time} view
    Then the "rows" reading of grid should be at least 1
    And "Data loading error" dialog should be absent

  Scenario: Part 3 — the query's and the script's results are saved as a project
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexMoveQuery{time} tree node inside browse tree and the query it runs completes
    Then the "BDDComplexMoveQuery{time}" view should be current
    And the "rows" reading of grid should be at least 1
    Given the browse panel is open
    When user opens the address "/scripts"
    Then the "Scripts" view should be current
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
    Given the browse panel is open
    When user opens the address "/scripts"
    Then the "Scripts" view should be current
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
    When user clicks on BDDComplexMoveQuery{time} view
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
