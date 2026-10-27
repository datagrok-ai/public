@realizes:views.projects
Feature: Data sync edge cases: the query behind a dynamic dashboard deleted, recreated, or without a column
  A query on the Datagrok database (made through the API: writing a query is the database features'
  subject) is run from the Browse tree, a bar chart split by its "name" column is added, and the
  result is saved through the ribbon's Save dialog with Data sync on. What the dashboard does when
  the query changes is stated in `help/datagrok/concepts/project/dashboard.md`, "When the source
  changes":
  - the query deleted: "The dashboard fails to open with an error". Reopened, the dashboard shows
    the Data loading error dialog naming the query, and CLOSE PROJECT leaves no table open;
  - "A new query with the same name on the same connection repairs it": the query deleted and a new
    one saved under its name, the dashboard reopens with the rows the new query returns, rebuilt by
    data sync;
  - "A column renamed or removed: viewers bound to it show a placeholder with the missing column
    name. The dashboard opens": a new query of the name without the "name" column, and the bar
    chart says that the column "name" does not exist;
  - a column a calculated column reads removed from the query: the dashboard opens, a warning balloon
    names the calculated column and the missing column (`openProject` in shell_project.dart,
    "Some calculated columns may not work"), and the calculated column is empty;
  - the file a dashboard was saved from with Data sync deleted from its space: reopened, the
    dashboard shows the Data loading error dialog, as for a deleted query.

  Each reopen follows a page reload, so the dashboard is read from the server and not from the
  session that saved it. The query and the dashboards carry the run's time (letters and digits only:
  the Dashboards search misses "-" and "_") and are removed when a scenario starts and when the
  feature ends.

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And no project named "BDDEdgeSyncDel{time}" is on the server
    And no project named "BDDEdgeSyncRep{time}" is on the server
    And no project named "BDDEdgeSyncCol{time}" is on the server
    And no project named "BDDEdgeSyncCalc{time}" is on the server
    And no project named "BDDEdgeSyncFile{time}" is on the server

  @edge-cases
  Scenario: With its query deleted, the dashboard does not open and names the query in the Data loading error
    Given a query "BDDEdgeSyncQDel{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
    And the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDEdgeSyncQDel{time} tree node inside browse tree and the query it runs completes
    Then the "BDDEdgeSyncQDel{time}" view should be current
    And the table should have 3 rows
    Given user adds a bar chart viewer with:
      | splitColumnName | name |
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDEdgeSyncQDel{time}" project table in "Save project" dialog should be checked
    When user enters "BDDEdgeSyncDel{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSyncDel{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeSyncDel{time}" dialog
    Then the "Share BDDEdgeSyncDel{time}" dialog should close
    And the "BDDEdgeSyncQDel{time}" table of the "BDDEdgeSyncDel{time}" project should be saved with data sync
    When user closes all views
    Given no query named "BDDEdgeSyncQDel{time}" is on the server
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncDel{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncDel{time} gallery card
    Then "Data loading error" dialog should become visible within 60 seconds
    And "Data loading error" dialog should contain text "could not load some of its data"
    And "Data loading error" dialog should contain text "BDDEdgeSyncQDel{time}"
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close
    And no table should be left in the workspace

  @edge-cases
  Scenario: A new query saved under the deleted query's name repairs the dashboard
    Given a query "BDDEdgeSyncQRep{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
    And the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDEdgeSyncQRep{time} tree node inside browse tree and the query it runs completes
    Then the "BDDEdgeSyncQRep{time}" view should be current
    And the table should have 3 rows
    Given user adds a bar chart viewer with:
      | splitColumnName | name |
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDEdgeSyncQRep{time}" project table in "Save project" dialog should be checked
    When user enters "BDDEdgeSyncRep{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSyncRep{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeSyncRep{time}" dialog
    Then the "Share BDDEdgeSyncRep{time}" dialog should close
    And the "BDDEdgeSyncQRep{time}" table of the "BDDEdgeSyncRep{time}" project should be saved with data sync
    When user closes all views
    Given no query named "BDDEdgeSyncQRep{time}" is on the server
    And a query "BDDEdgeSyncQRep{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project')"
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncRep{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncRep{time} gallery card and the project opens
    Then the "BDDEdgeSyncQRep{time}" view should be current
    And the table should have been reloaded by data sync
    And the table should have 2 rows
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown

  @edge-cases
  Scenario: Without the column its bar chart is split by, the dashboard opens and the chart names the missing column
    Given a query "BDDEdgeSyncQCol{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
    And the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDEdgeSyncQCol{time} tree node inside browse tree and the query it runs completes
    Then the "BDDEdgeSyncQCol{time}" view should be current
    And the table should have 3 rows
    Given user adds a bar chart viewer with:
      | splitColumnName | name |
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDEdgeSyncQCol{time}" project table in "Save project" dialog should be checked
    When user enters "BDDEdgeSyncCol{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSyncCol{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeSyncCol{time}" dialog
    Then the "Share BDDEdgeSyncCol{time}" dialog should close
    And the "BDDEdgeSyncQCol{time}" table of the "BDDEdgeSyncCol{time}" project should be saved with data sync
    When user closes all views
    Given a query "BDDEdgeSyncQCol{time}" on "System:Datagrok" reads "select id from entity_types where name in ('Script', 'Project', 'User')"
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncCol{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncCol{time} gallery card and the project opens
    Then the "BDDEdgeSyncQCol{time}" view should be current
    And the table should have been reloaded by data sync
    And the table should have 3 rows
    And the table should not have a column "name"
    And bar chart viewer should contain text 'Column "name" does not exist'
    And "Data loading error" dialog should be absent

  @edge-cases
  Scenario: Without a column its calculated column reads, the dashboard opens with a warning naming it
    Given a query "BDDEdgeSyncQCalc{time}" on "System:Datagrok" reads "select name, length(name) as len from entity_types where name in ('Script', 'Project', 'User')"
    And the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDEdgeSyncQCalc{time} tree node inside browse tree and the query it runs completes
    Then the "BDDEdgeSyncQCalc{time}" view should be current
    And the table should have 3 rows
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user enters "LenPlus" into Name input in "Add New Column" dialog
    And user replaces the code of code editor in "Add New Column" dialog with "${len} + 1"
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "LenPlus" column should have no missing values
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDEdgeSyncQCalc{time}" project table in "Save project" dialog should be checked
    When user enters "BDDEdgeSyncCalc{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSyncCalc{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeSyncCalc{time}" dialog
    Then the "Share BDDEdgeSyncCalc{time}" dialog should close
    When user closes all views
    Given a query "BDDEdgeSyncQCalc{time}" on "System:Datagrok" reads "select name from entity_types where name in ('Script', 'Project', 'User')"
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncCalc{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncCalc{time} gallery card and the project opens
    Then the "BDDEdgeSyncQCalc{time}" view should be current
    And the table should have 3 rows
    And the table should not have a column "len"
    And an error or warning balloon matching "Some calculated columns may not work.*LenPlus.*missing columns: len" should have been shown
    And "LenPlus" column should have missing values
    And "Data loading error" dialog should be absent

  @edge-cases
  Scenario: With its file deleted, the dashboard does not open and shows the Data loading error
    Given a space named "BDDEdgeSyncSpace{time}" is on the server
    And the "fixtures/cars-small.csv" file of the project is in the space "BDDEdgeSyncSpace{time}" as "syncedcars{time}.csv"
    And the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDEdgeSyncSpace{time} tree node inside browse tree
    Then the "BDDEdgeSyncSpace{time}" view should be current
    When user double-clicks on "syncedcars{time}.csv" link in gallery
    Then the "syncedcars{time}" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "syncedcars{time}" project table in "Save project" dialog should be checked
    When user enters "BDDEdgeSyncFile{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSyncFile{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeSyncFile{time}" dialog
    Then the "Share BDDEdgeSyncFile{time}" dialog should close
    When user closes all views
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDEdgeSyncSpace{time} tree node inside browse tree
    And user picks "Delete..." from the context menu of "syncedcars{time}.csv" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "syncedcars{time}.csv" link in gallery should be absent
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncFile{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncFile{time} gallery card
    Then "Data loading error" dialog should become visible within 60 seconds
    And "Data loading error" dialog should contain text "could not load some of its data"
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close
