@journey @serial @realizes:views.projects @realizes:views.space @realizes:sharing.share-dialog @realizes:data.menu.join-tables @realizes:data.menu.aggregate-rows
Feature: A project whose tables are renamed, grown by drags, and whose query and script move into a Space
  The TestTrack case Projects/complex (with its manual companion complex-ui). A project is built from
  four sources — orders.csv from Browse > Files > Demo > northwind, the result of a saved query run
  from the Browse tree, a database table through Get All, and the demog table of a JavaScript script
  run from the Scripts view — and three tables derived from the script's demog: a pivot published
  with ADD, an aggregation from Data > Aggregate Rows..., and a join of demog with the pivot from
  Data > Join Tables.... It is saved through the ribbon's Save dialog with Data sync, and the server
  holds the join in the same project (GROK-19103: the join once went into a separate, broken
  project). Every table but demog and its pivot (see below) is renamed from its view tab
  (right-click > Table > Rename...), and the project is saved as a copy with Data sync; the Dashboards panel then names the copy, and the copy
  reopens from the Dashboards gallery with the new names, the rows and Data sync still on
  (GROK-19212: a renamed table lost its Data sync on reopen). A file of a Space, the result of a
  second query and the database table, opened again, are dragged in the Dashboards panel onto the
  open copy's node and saved with the node's SAVE; the copy reopens with them. The query is then
  renamed in the Browse tree and the script in its editor, and both are moved into the Space (the
  Move entity dialog's Move; complex-ui step 10, the Space half only — an entity cannot be moved
  into a file share). Both projects still open with every table — the pivot, the aggregation and
  the join over the renamed and moved script included — rebuilt by data sync and no error dialog,
  and their Save dialogs show Data sync on and a creation script for every table saved with it.
  Last, the original is shared with the second account with Full access (the Share dialog's
  privilege tree), which opens it and saves the original: the server holds one project of that
  name, changed after the recipient's OK, and no copy of it.

  Fixtures substituted: NorthwindTest (a dev-only connection) is replaced by System:Datagrok, which
  every user may read and query — the md's complexQuery over orders becomes "BDDComplexQ{time}"
  (select id, name from public.entity_types) and the second query "BDDComplexQB{time}", both saved
  through the JS API with the script, and the database tables customers and products become
  public.entity_types (Get All); their row count is the one read at the first open, since a table of
  the Datagrok database differs between servers. The md's pivot of orders, aggregation of customers
  and join of customers with complexQuery become a pivot (RACE by SEX, avg(AGE)), an aggregation
  (count(USUBJID) by DIS_POP) and a join on RACE, all over the script's demog. The Space holds a copy
  of customers.csv rather than of demog.csv (a second "demog" in the workspace would take a
  generated name). The tables are dragged onto the copy rather than the original (the md's step 12),
  and the second drag is of a second query, since the query already in the project keeps its name.

  Not translated here, and why, and where it is claimed: the copy without Data sync and its re-save
  with Data sync (steps 13-14) — projects-data-sync (complex-save-copy) for one table; here the copy
  is saved with Data sync, and step 13's claim on the node's name is read after that save. The
  project renamed and moved into a Space and opened from it (steps 17, 20, 22's open from the
  Space) — projects-lifecycle-* and projects-move; here the moved query and script are what the
  reopened projects depend on. The view-only recipient's Save dialog and the shares of a Space
  table and an unshared script (steps 24, 26, GROK-18345, GROK-19403) — projects-lifecycle-files,
  projects-lifecycle-spaces and projects-lifecycle-script; the broken script (steps 28-29,
  GROK-19728) — projects-lifecycle-script. complex-ui step 4's Link / Clone / Move / Copy choice:
  a table dropped onto a project node opens Move entity with the project and the table and a YES,
  no such choice. Logout and signing in with the second user's credentials — the platform's Logout
  ends every session of the account, which all workers of a run share, so the second account is
  entered through its own session. A clean console wherever the Save dialog was opened (its
  preview clones the views into an iframe and logs "Unable to find element in cloned iframe", known
  noise with no ticket).

  Known failure, GROK-19212 (reproduced 2026-09-28 on localhost through both view-tab paths, the
  tab's Table > Rename... and the tab click's Context Panel > Actions > Rename..., 2 and 3 times):
  a project saved after the source of a derived table was renamed does not reopen the derived
  table. orders.csv and its pivot (Toolbox > Pivot table > ADD, "orders aggregation", 90 rows) are
  saved with Data sync, orders is renamed to ordersR and the project saved again; reopened from the
  Dashboards gallery, ordersR is re-read by data sync with 830 rows, but the pivot never arrives,
  no table view opens, the status bar keeps "2 tasks" and the console logs 'Could not resolve
  table "orders"' (the pivot's creation script still names the old table); no dialog. A join does
  the same (reproduced 2 of 2 for either input): the result of the saved query BDDComplexQB{time}
  joined with entity_types (Get All) on id, the query's table renamed before the first save; on
  reopen both inputs are re-read, the join never arrives, "2 tasks" stays and the console logs
  'Could not resolve table "BDDComplexQB…"'; no dialog. Renaming the join's own table reopens fine.
  Each case is a normal scenario that reopens the project and waits for the tables that do load,
  followed by a @known-failure scenario that claims the derived table with its rows and the views.
  For the same reason demog and demog aggregation, the sources of the pivot and the join, keep
  their names in the renamed copy.

  The md's join of a database table with the query result is not built here: once the query is
  renamed, a project holding such a join does not open (GROK-21026), which
  projects-regressions-query-rename-join claims. Observed on localhost and dev (2026-09-25), left
  out of the claims until the operator decides: the tables dragged onto the project node arrive
  with Data sync off, so only their rows are claimed.

  Everything is named with the run's time and removed at the start and the end: the four projects
  (with their tables and views), a stray "Copy of" the original, both queries (the first under both
  names), the script under both names (the editor's save removes the one it saved), and the Space
  last, after what was moved into it. It is @serial: the Dashboards search is shared with every
  feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDComplex{time}" is on the server
    And no project named "BDDComplexRenamed{time}" is on the server
    And no project named "Copy of BDDComplex{time}" is on the server
    And no project named "BDDComplexPivot{time}" is on the server
    And no project named "BDDComplexJoin{time}" is on the server
    And no query named "BDDComplexQ{time}, BDDComplexQRen{time}, BDDComplexQB{time}" is on the server
    And no script named "BDDComplexScriptRen{time}" is on the server
    And a query "BDDComplexQ{time}" on the Datagrok connection is:
      """
      select id, name from public.entity_types
      """
    And a query "BDDComplexQB{time}" on the Datagrok connection is:
      """
      select id, friendly_name from public.entity_types
      """
    And a script "BDDComplexScript{time}" is on the server:
      """
      //language: javascript
      //output: dataframe df
      df = await grok.data.getDemoTable('demog.csv');
      """
    And no space named "BDDComplex{time}" is on the server
    And user clears the saved pivot table parameters

  Scenario: A file, a query, a database table and a script open into the workspace
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And the table should have 830 rows
    Given the browse panel is open
    # the tree keeps the children it listed before the query was saved; a refresh lists it
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexQ{time} tree node inside browse tree
    Then the "BDDComplexQ{time}" view should be current
    When user remembers the row count of table "BDDComplexQ{time}" as "complex query"
    Given the browse panel is open
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    And table "entity_types" should have the "complex query" row count
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    And user clears gallery search
    And user types "BDDComplexScript{time}" into gallery search
    And user picks "Run..." from the context menu of "BDDComplexScript{time}" link in gallery
    Then the "demog" view should be current
    And the table should have 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # the md's join of a database table with the query result is not built: GROK-21026, claimed in
  # projects-regressions-query-rename-join
  Scenario: A pivot, an aggregation and a join are derived from the script's table (GROK-19103)
    When user switches to the "demog" table view
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    And user adds "RACE" to the "group by" row of first pivot table viewer
    And user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    And user adds "SEX" to the "pivot" row of first pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    And the "pivot" reading of pivot table viewer should be "SEX"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "RACE" pivoted on "SEX"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    And table "demog aggregation" should have 4 rows
    And table "demog aggregation" should have columns "RACE, F avg(AGE), M avg(AGE)"
    When user switches to the "demog" table view
    And user picks "Data > Aggregate Rows..." from the top menu
    Then the open tableview should have 2 pivot table viewers
    When user clicks on the "remove pivot chip SEVERITY" area of second pivot table viewer
    And user picks "Column > USUBJID" from the context menu of the "aggregate chip avg(AGE)" area of second pivot table viewer
    And user closes the context menu
    And user picks "Aggregation > count" from the context menu of the "aggregate chip first(USUBJID)" area of second pivot table viewer
    And user closes the context menu
    Then the "group by" reading of second pivot table viewer should be "DIS_POP"
    And the "pivot" reading of second pivot table viewer should be ""
    And the "aggregate" reading of second pivot table viewer should be "count(USUBJID)"
    And the aggregated values of second pivot table viewer should match "count(USUBJID)" grouped by "DIS_POP"
    When user clicks on the "add to workspace" area of second pivot table viewer
    Then the "demog aggregation (2)" view should be current
    And table "demog aggregation (2)" should have 6 rows
    And table "demog aggregation (2)" should have columns "DIS_POP, count(USUBJID)"
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "demog" in join left table selector
    And user selects "demog aggregation" in join right table selector
    And user picks column "RACE" in join left key selector
    And user picks column "RACE" in join right key selector
    Then "Join Type" input in "Join Tables" dialog should have value "inner"
    When user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have 5850 rows
    Given the dashboards panel of the left sidebar is open
    Then New-Dashboard---result tree node inside browse tree should be visible
    And there should be 1 visible dashboards project node
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The seven tables are saved in one project with Data sync
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "orders" project table in "Save project" dialog             |
      | Data sync switch in "BDDComplexQ{time}" project table in "Save project" dialog  |
      | Data sync switch in "entity_types" project table in "Save project" dialog       |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation (2)" project table in "Save project" dialog  |
      | Data sync switch in "result" project table in "Save project" dialog             |
    When user enters "BDDComplex{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDComplex{time}\" uploaded" should have been shown
    And 1 project named "BDDComplex{time}" should be on the server
    And the "BDDComplex{time}" project on the server should hold the tables "orders, BDDComplexQ{time}, entity_types, demog, demog aggregation, demog aggregation (2), result"
    And the "orders" table of the "BDDComplex{time}" project should be saved with data sync
    And the "BDDComplexQ{time}" table of the "BDDComplex{time}" project should be saved with data sync
    And the "entity_types" table of the "BDDComplex{time}" project should be saved with data sync
    And the "demog" table of the "BDDComplex{time}" project should be saved with data sync
    And the "demog aggregation" table of the "BDDComplex{time}" project should be saved with data sync
    And the "demog aggregation (2)" table of the "BDDComplex{time}" project should be saved with data sync
    And the "result" table of the "BDDComplex{time}" project should be saved with data sync
    And "Share BDDComplex{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplex{time}" dialog
    Then the "Share BDDComplex{time}" dialog should close

  # demog and demog aggregation, the sources of the pivot and the join, keep their names: renamed, the
  # project reopens without the tables derived from them (GROK-19212, claimed by the last two scenarios)
  Scenario: The file, query, database, aggregation and join tables are renamed from their view tabs (GROK-19212)
    Given simple mode is off
    When user picks "Table > Rename..." from the context menu of "orders" view
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "orders"
    When user enters "ordersR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "ordersR" should be open
    When user picks "Table > Rename..." from the context menu of "BDDComplexQ{time}" view
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "BDDComplexQ{time}"
    When user enters "complexQueryR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "complexQueryR" should be open
    When user picks "Table > Rename..." from the context menu of "entity_types" view
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "entity_types"
    When user enters "entity_typesR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "entity_typesR" should be open
    When user picks "Table > Rename..." from the context menu of "demog aggregation (2)" view
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "demog aggregation (2)"
    When user enters "demog aggregation (2)R" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "demog aggregation (2)R" should be open
    When user picks "Table > Rename..." from the context menu of "result" view
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "result"
    When user enters "resultR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "resultR" should be open
    And the table views "ordersR, complexQueryR, entity_typesR, demog, demog aggregation, demog aggregation (2)R, resultR" should be open
    And no errors should have been logged

  Scenario: A copy saved with Data sync keeps the new names, and the Dashboards panel names the copy (GROK-19212)
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then the following elements should be checked:
      | Data sync switch in "ordersR" project table in "Save project" dialog             |
      | Data sync switch in "complexQueryR" project table in "Save project" dialog   |
      | Data sync switch in "entity_typesR" project table in "Save project" dialog       |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation (2)R" project table in "Save project" dialog  |
      | Data sync switch in "resultR" project table in "Save project" dialog             |
    When user enters "BDDComplexRenamed{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDComplexRenamed{time}\" uploaded" should have been shown
    And 1 project named "BDDComplexRenamed{time}" should be on the server
    And the "BDDComplexRenamed{time}" project on the server should hold the tables "ordersR, complexQueryR, entity_typesR, demog, demog aggregation, demog aggregation (2)R, resultR"
    And the "ordersR" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "complexQueryR" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "entity_typesR" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "demog" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "demog aggregation" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "demog aggregation (2)R" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "resultR" table of the "BDDComplexRenamed{time}" project should be saved with data sync
    And the "BDDComplex{time}" project on the server should hold the tables "orders, BDDComplexQ{time}, entity_types, demog, demog aggregation, demog aggregation (2), result"
    Given the dashboards panel of the left sidebar is open
    Then BDDComplexRenamed{time} tree node inside browse tree should be visible
    And BDDComplex{time} tree node inside browse tree should be absent
    And New-Dashboard tree node inside browse tree should be visible
    And there should be 2 visible dashboards project nodes

  Scenario: A Space gets a copy of customers.csv
    Given the browse panel is open
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDDComplex{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDComplex{time}" should be on the server
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---northwind tree node inside browse tree
    Then customers.csv link in gallery should be visible
    When user drags customers.csv link in gallery to Spaces---BDDComplex{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given Spaces---BDDComplex{time} tree node inside browse tree is expanded
    Then Spaces---BDDComplex{time}---customers.csv tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Reopened, the copy's renamed tables keep their names, rows and Data sync (GROK-19212)
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRenamed{time} gallery card
    Then the table views "ordersR, complexQueryR, entity_typesR, demog, demog aggregation, demog aggregation (2)R, resultR" should be open
    And table "ordersR" should have been reloaded by data sync with 830 rows
    And table "complexQueryR" should have been reloaded by data sync with the "complex query" row count
    And table "entity_typesR" should have been reloaded by data sync with the "complex query" row count
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)R" should have been reloaded by data sync with 6 rows
    And table "resultR" should have been reloaded by data sync with 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be visible:
      | "Creation script" button in "ordersR" project table in "Save project" dialog             |
      | "Creation script" button in "complexQueryR" project table in "Save project" dialog   |
      | "Creation script" button in "entity_typesR" project table in "Save project" dialog       |
      | "Creation script" button in "demog" project table in "Save project" dialog              |
      | "Creation script" button in "demog aggregation" project table in "Save project" dialog |
      | "Creation script" button in "demog aggregation (2)R" project table in "Save project" dialog  |
      | "Creation script" button in "resultR" project table in "Save project" dialog             |
    And the following elements should be checked:
      | Data sync switch in "ordersR" project table in "Save project" dialog             |
      | Data sync switch in "complexQueryR" project table in "Save project" dialog   |
      | Data sync switch in "entity_typesR" project table in "Save project" dialog       |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation (2)R" project table in "Save project" dialog  |
      | Data sync switch in "resultR" project table in "Save project" dialog             |
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: A Space file, a query's result and a database table are dragged onto the open copy
    Given the browse panel is open
    And Spaces---BDDComplex{time} tree node inside browse tree is expanded
    When user double-clicks on Spaces---BDDComplex{time}---customers.csv tree node inside browse tree
    Then the "customers" view should be current
    And the table should have 91 rows
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexQB{time} tree node inside browse tree
    Then the "BDDComplexQB{time}" view should be current
    And table "BDDComplexQB{time}" should have the "complex query" row count
    Given the browse panel is open
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    And table "entity_types" should have the "complex query" row count
    Given the dashboards panel of the left sidebar is open
    Then the following elements should be visible:
      | New-Dashboard---customers tree node inside browse tree          |
      | New-Dashboard---BDDComplexQB{time} tree node inside browse tree |
      | New-Dashboard---entity_types tree node inside browse tree       |
    When user collapses BDDComplexRenamed{time} tree node inside browse tree
    And user drags New-Dashboard---customers tree node inside browse tree to BDDComplexRenamed{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDComplexRenamed{time} project"
    And "customers" project table in Move entity dialog should be visible
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    # the panel rebuilds its tree after a move; the next drag waits for it
    And New-Dashboard---customers tree node inside browse tree should be absent
    When user collapses BDDComplexRenamed{time} tree node inside browse tree
    And user drags New-Dashboard---BDDComplexQB{time} tree node inside browse tree to BDDComplexRenamed{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And "BDDComplexQB{time}" project table in Move entity dialog should be visible
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    # the panel rebuilds its tree after a move; the next drag waits for it
    And New-Dashboard---BDDComplexQB{time} tree node inside browse tree should be absent
    When user collapses BDDComplexRenamed{time} tree node inside browse tree
    And user drags New-Dashboard---entity_types tree node inside browse tree to BDDComplexRenamed{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And "entity_types" project table in Move entity dialog should be visible
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user expands BDDComplexRenamed{time} tree node inside browse tree
    Then the following elements should be visible:
      | BDDComplexRenamed{time}---customers tree node inside browse tree          |
      | BDDComplexRenamed{time}---BDDComplexQB{time} tree node inside browse tree |
      | BDDComplexRenamed{time}---entity_types tree node inside browse tree       |
    And New-Dashboard---entity_types tree node inside browse tree should be absent
    When user clicks on Save button in BDDComplexRenamed{time} tree node inside browse tree
    Then "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "BDDComplexRenamed{time}" project on the server should hold the tables "ordersR, complexQueryR, entity_typesR, demog, demog aggregation, demog aggregation (2)R, resultR, customers, BDDComplexQB{time}, entity_types"
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRenamed{time} gallery card
    Then the table views "customers, BDDComplexQB{time}, entity_types" should be open
    And table "customers" should have 91 rows
    And table "BDDComplexQB{time}" should have the "complex query" row count
    And table "entity_types" should have the "complex query" row count
    And no error or warning balloon should have been shown

  Scenario: The query and the script are renamed and moved into the Space
    When user closes all views
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of Databases---Postgres---Datagrok---BDDComplexQ{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "BDDComplexQRen{time}" into Name input in "Rename dataquery" dialog
    And user clicks on OK button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDComplexQRen{time}" should be on the server
    And 0 queries named "BDDComplexQ{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user drags Databases---Postgres---Datagrok---BDDComplexQRen{time} tree node inside browse tree to Spaces---BDDComplex{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDComplex{time} project"
    And Move entity dialog should contain text "BDDComplexQRen{time}"
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    And the "BDDComplex{time}" space should hold the query "BDDComplexQRen{time}" on the server
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDComplexScript{time}" into gallery search
    And user picks "Edit..." from the context menu of "BDDComplexScript{time}" link in gallery
    Then the "BDDComplexScript{time}" view should be current
    When user replaces the code of code editor with "//name: BDDComplexScriptRen{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = await grok.data.getDemoTable('demog.csv');" to code editor
    And user saves the script
    Then 1 script named "BDDComplexScriptRen{time}" should be on the server
    And 0 scripts named "BDDComplexScript{time}" should be on the server
    Given user opens the Scripts view
    And the browse panel is open
    When user clears gallery search
    And user types "BDDComplexScriptRen{time}" into gallery search
    And user drags "BDDComplexScriptRen{time}" link in gallery to Spaces---BDDComplex{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDComplex{time} project"
    And Move entity dialog should contain text "BDDComplexScriptRen{time}"
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    And the "BDDComplex{time}" space should hold the script "BDDComplexScriptRen{time}" on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: After the moves, the copy's derived tables over the moved query and script are rebuilt
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexRenamed{time} gallery card
    Then the table views "ordersR, complexQueryR, entity_typesR, demog, demog aggregation, demog aggregation (2)R, resultR, customers, BDDComplexQB{time}, entity_types" should be open
    And table "ordersR" should have been reloaded by data sync with 830 rows
    And table "complexQueryR" should have been reloaded by data sync with the "complex query" row count
    And table "entity_typesR" should have been reloaded by data sync with the "complex query" row count
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)R" should have been reloaded by data sync with 6 rows
    And table "resultR" should have been reloaded by data sync with 5850 rows
    # tables dragged onto a project node arrive with Data sync off (switched on separately, not a defect),
    # so only their rows are claimed
    And table "customers" should have 91 rows
    And table "BDDComplexQB{time}" should have the "complex query" row count
    And table "entity_types" should have the "complex query" row count
    And "Data loading error" dialog should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be visible:
      | "Creation script" button in "ordersR" project table in "Save project" dialog             |
      | "Creation script" button in "complexQueryR" project table in "Save project" dialog   |
      | "Creation script" button in "entity_typesR" project table in "Save project" dialog       |
      | "Creation script" button in "demog" project table in "Save project" dialog              |
      | "Creation script" button in "demog aggregation" project table in "Save project" dialog |
      | "Creation script" button in "demog aggregation (2)R" project table in "Save project" dialog  |
      | "Creation script" button in "resultR" project table in "Save project" dialog             |
    And the following elements should be checked:
      | Data sync switch in "ordersR" project table in "Save project" dialog             |
      | Data sync switch in "complexQueryR" project table in "Save project" dialog   |
      | Data sync switch in "entity_typesR" project table in "Save project" dialog       |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation (2)R" project table in "Save project" dialog  |
      | Data sync switch in "resultR" project table in "Save project" dialog             |
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: After the moves, the original opens every table with its rows, Data sync and creation script
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplex{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplex{time} gallery card
    Then the table views "orders, BDDComplexQ{time}, entity_types, demog, demog aggregation, demog aggregation (2), result" should be open
    And table "orders" should have been reloaded by data sync with 830 rows
    And table "BDDComplexQ{time}" should have been reloaded by data sync with the "complex query" row count
    And table "entity_types" should have been reloaded by data sync with the "complex query" row count
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)" should have been reloaded by data sync with 6 rows
    And table "result" should have been reloaded by data sync with 5850 rows
    And "Data loading error" dialog should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be visible:
      | "Creation script" button in "orders" project table in "Save project" dialog             |
      | "Creation script" button in "BDDComplexQ{time}" project table in "Save project" dialog  |
      | "Creation script" button in "entity_types" project table in "Save project" dialog       |
      | "Creation script" button in "demog" project table in "Save project" dialog              |
      | "Creation script" button in "demog aggregation" project table in "Save project" dialog |
      | "Creation script" button in "demog aggregation (2)" project table in "Save project" dialog  |
      | "Creation script" button in "result" project table in "Save project" dialog             |
    And the following elements should be checked:
      | Data sync switch in "orders" project table in "Save project" dialog             |
      | Data sync switch in "BDDComplexQ{time}" project table in "Save project" dialog  |
      | Data sync switch in "entity_types" project table in "Save project" dialog       |
      | Data sync switch in "demog" project table in "Save project" dialog              |
      | Data sync switch in "demog aggregation" project table in "Save project" dialog |
      | Data sync switch in "demog aggregation (2)" project table in "Save project" dialog  |
      | Data sync switch in "result" project table in "Save project" dialog             |
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: Shared with Full access, the recipient saves the original project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplex{time}" into gallery search
    And user picks "Share..." from the context menu of BDDComplex{time} gallery card
    Then "Share BDDComplex{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDComplex{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDComplex{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on share access selector
    Then Full-access tree node inside privilege tree should be unchecked
    When user checks Full-access tree node inside privilege tree
    And user clicks outside the privilege tree
    Then share access selector should contain text "Full access"
    When user clicks on OK button in "Share BDDComplex{time}" dialog
    Then the "Share BDDComplex{time}" dialog should close
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplex{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplex{time} gallery card
    Then the table views "orders, BDDComplexQ{time}, entity_types, demog, demog aggregation, demog aggregation (2), result" should be open
    And table "orders" should have been reloaded by data sync with 830 rows
    And table "BDDComplexQ{time}" should have been reloaded by data sync with the "complex query" row count
    And table "entity_types" should have been reloaded by data sync with the "complex query" row count
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)" should have been reloaded by data sync with 6 rows
    And table "result" should have been reloaded by data sync with 5850 rows
    And "Data loading error" dialog should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user remembers when the "BDDComplex{time}" project was saved
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be enabled
    And "Save original project" radio choice in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDComplex{time}\" uploaded" should have been shown
    And 1 project named "BDDComplex{time}" should be on the server
    And 0 projects named "Copy of BDDComplex{time}" should be on the server
    And the "BDDComplex{time}" project should have been saved again since remembered
    When user closes all views
    Given user signs in as themselves again

  Scenario: A pivot's source renamed from its view tab is saved, and the project is reopened (GROK-19212)
    When user closes all views
    Given simple mode is off
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And table "orders aggregation" should have 90 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "orders" project table in "Save project" dialog             |
      | Data sync switch in "orders aggregation" project table in "Save project" dialog |
    When user enters "BDDComplexPivot{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDComplexPivot{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplexPivot{time}" dialog
    Then the "Share BDDComplexPivot{time}" dialog should close
    When user picks "Table > Rename..." from the context menu of "orders" view
    Then "Rename table" dialog should be visible
    When user enters "ordersR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "ordersR" should be open
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "BDDComplexPivot{time}" project on the server should hold the tables "ordersR, orders aggregation"
    And the "ordersR" table of the "BDDComplexPivot{time}" project should be saved with data sync
    And the "orders aggregation" table of the "BDDComplexPivot{time}" project should be saved with data sync
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexPivot{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexPivot{time} gallery card
    Then table "ordersR" should have been reloaded by data sync with 830 rows

  @known-failure @realizes:GROK-19212
  Scenario: The reopened project holds the pivot of the renamed source (GROK-19212)
    Then table "orders aggregation" should have been reloaded by data sync with 90 rows
    And the table views "ordersR, orders aggregation" should be open

  Scenario: A join whose input was renamed from its view tab is saved, and the project is reopened (GROK-19212)
    When user closes all views
    Given simple mode is off
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDComplexQB{time} tree node inside browse tree
    Then the "BDDComplexQB{time}" view should be current
    And table "BDDComplexQB{time}" should have the "complex query" row count
    Given the browse panel is open
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "BDDComplexQB{time}" in join left table selector
    And user selects "entity_types" in join right table selector
    And user picks column "id" in join left key selector
    And user picks column "id" in join right key selector
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have the "complex query" row count
    When user picks "Table > Rename..." from the context menu of "BDDComplexQB{time}" view
    Then "Rename table" dialog should be visible
    When user enters "complexQueryBR" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "complexQueryBR" should be open
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "complexQueryBR" project table in "Save project" dialog |
      | Data sync switch in "entity_types" project table in "Save project" dialog   |
      | Data sync switch in "result" project table in "Save project" dialog         |
    When user enters "BDDComplexJoin{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDComplexJoin{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDComplexJoin{time}" dialog
    Then the "Share BDDComplexJoin{time}" dialog should close
    And the "BDDComplexJoin{time}" project on the server should hold the tables "complexQueryBR, entity_types, result"
    And the "result" table of the "BDDComplexJoin{time}" project should be saved with data sync
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDComplexJoin{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDComplexJoin{time} gallery card
    Then table "complexQueryBR" should have been reloaded by data sync with the "complex query" row count
    And table "entity_types" should have been reloaded by data sync with the "complex query" row count

  @known-failure @realizes:GROK-19212
  Scenario: The reopened project holds the join of the renamed input (GROK-19212)
    Then table "result" should have been reloaded by data sync with the "complex query" row count
    And the table views "complexQueryBR, entity_types, result" should be open
