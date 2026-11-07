@journey @realizes:views.projects
Feature: A project saved with and without Data sync
  demog.csv opened from Browse > Files > Demo is saved through the ribbon's Save dialog with Data
  sync on, then saved as a copy with Data sync off, and the copy is saved again with Data sync on.
  What the server keeps is claimed for every save (the table's data-sync flag and creation script).
  Translated from the TestTrack case Projects/complex-save-copy.

  How a reopened project brings the table back — rebuilt by its creation script, or loaded from the
  uploaded snapshot — is ApiTests `Dapi: projects: data sync` (src/dapi/projects.ts).

  Both projects are named with the run's time and removed, with their tables and views, when the
  feature ends.

  Not translated, and why: opening the projects from the Dashboards tiles (the tile gallery is
  Browse's subject) and the Share dialog the first save
  shows, which is cancelled (a copy of a saved project shows none).

  Edge cases: a table opened with Open local file has no creation script, so the Save dialog offers it
  no Data sync switch and no Creation script button (`help/datagrok/concepts/project/dashboard.md`,
  "Data sync"). What a dashboard saved with Data sync does when its source changes ("When the source
  changes" there), each reopen following a page reload so the dashboard is read from the server: its
  query deleted, it does not open and the Data loading error dialog names the query; a new query saved
  under the deleted query's name repairs it, rebuilt by data sync; its query replaced by one of the
  same name without the column its bar chart is split by, it opens, the chart saying the column does
  not exist; replaced by one without the column a calculated column reads, it opens with a warning
  balloon naming both, the calculated column empty; the file it was saved from deleted from its space,
  it shows the Data loading error, reopened without a reload. These queries are made through the API
  on System:Datagrok. Their projects are named with the run's time in letters and digits only (the Dashboards search misses "-"
  and "_") and removed, with their tables and views, when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDD-Sync-Orig-{time}" is on the server
    And no project named "BDD-Sync-Copy-{time}" is on the server

  Scenario: A file opened from Browse is saved with Data sync on
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 5850 rows
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDD-Sync-Orig-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDD-Sync-Orig-{time}" should have been shown
    And 1 project named "BDD-Sync-Orig-{time}" should be on the server
    And the "demog" table of the "BDD-Sync-Orig-{time}" project should be saved with data sync
    And "Share BDD-Sync-Orig-{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Sync-Orig-{time}" dialog
    Then the "Share BDD-Sync-Orig-{time}" dialog should close
    And no errors should have been logged

  Scenario: A copy saved with Data sync off leaves the original synced
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDD-Sync-Orig-{time}"
    And choice input in "demog" project table in "Save project" dialog should have value "Clone"
    When user enters "BDD-Sync-Copy-{time}" into Name text input in "Save project" dialog
    And user unchecks Data sync switch in "demog" project table in "Save project" dialog
    Then Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    And "Creation script" button in "demog" project table in "Save project" dialog should be hidden
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Sync-Copy-{time}" should be on the server
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved as a snapshot
    And the "demog" table of the "BDD-Sync-Orig-{time}" project should be saved with data sync
    And no errors should have been logged

  Scenario: The copy's Save dialog keeps Data sync off
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be unchecked
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  Scenario: Turning Data sync on for the copy saves its creation script
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user checks Data sync switch in "demog" project table in "Save project" dialog
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text "System:DemoFiles/demog.csv"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "BDD-Sync-Copy-{time}" project should be saved with data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @edge-cases
  Scenario: A table opened from a local file is offered no Data sync
    Given the context panel is open
    And simple mode is off
    And the browse panel is open
    When user closes all views
    When user uploads "fixtures/browse-import.csv" through "Open local file" icon inside browse toolbar
    Then the "browse-import" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "browse-import" project table in "Save project" dialog should be visible
    And Data sync switch in "browse-import" project table in "Save project" dialog should be hidden
    And "Creation script" button in "browse-import" project table in "Save project" dialog should be hidden
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged
    When user closes all views

  @edge-cases
  Scenario: With its query deleted, the dashboard does not open and names the query in the Data loading error
    Given no project named "BDDEdgeSyncDel{time}" is on the server
    And a query "BDDEdgeSyncQDel{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
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
    Given no project named "BDDEdgeSyncRep{time}" is on the server
    And a query "BDDEdgeSyncQRep{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
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
    When user closes all views

  @edge-cases
  Scenario: With its query replaced by one of the same name without the bar chart's split column, the dashboard opens and the chart names that column
    Given no project named "BDDEdgeSyncCol{time}" is on the server
    And a query "BDDEdgeSyncQCol{time}" on "System:Datagrok" reads "select name, id from entity_types where name in ('Script', 'Project', 'User')"
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
    When user closes all views

  @edge-cases
  Scenario: With its query replaced by one of the same name without a column its calculated column reads, the dashboard opens with a warning naming it
    Given no project named "BDDEdgeSyncCalc{time}" is on the server
    And a query "BDDEdgeSyncQCalc{time}" on "System:Datagrok" reads "select name, length(name) as len from entity_types where name in ('Script', 'Project', 'User')"
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
    When user closes all views

  @edge-cases
  Scenario: With its file deleted, the dashboard shows the Data loading error
    Given no project named "BDDEdgeSyncFile{time}" is on the server
    And a space named "BDDEdgeSyncSpace{time}" is on the server
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
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSyncFile{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSyncFile{time} gallery card
    Then "Data loading error" dialog should become visible within 60 seconds
    And "Data loading error" dialog should contain text "could not load some of its data"
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close
