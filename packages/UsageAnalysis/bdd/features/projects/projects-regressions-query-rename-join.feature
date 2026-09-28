@journey @serial @realizes:views.projects @realizes:data.menu.join-tables
Feature: Projects regressions: a join over a query result after the query is renamed
  GROK-21026: a project holding a join over a query's result does not open after the query is
  renamed, when the query's name has no spaces (so its result table takes the same name). The
  operator's repro (projects-suspected-defects-repro.md, item 6), step by step through the UI: New
  SQL Query... on Browse > Databases > Postgres > Datagrok > Schemas > public > entity_types, saved
  under a name of letters and digits; the query run by a double-click in the tree; entity_types
  opened with Get All; Data > Join Tables... of entity_types with the query's table on id; the three
  tables saved through the ribbon's Save dialog with Data sync; everything closed; the query renamed
  in the tree (right-click > Rename...); everything closed; the project reopened from the Dashboards
  gallery.

  What the server does on the rename, claimed as it is: the creation scripts of the query's table
  and of the join ("result") are rewritten to name the query's new name, while the query's table in
  the project keeps the old one. That mismatch is the root of the defect: the reopen re-reads the
  query's table and entity_types, but the join asks for a table by the query's new name, and
  "result" never arrives (the console logs 'Could not resolve table "<new name>"'). On core
  3e53f6462e the open never completes: no table view opens, no dialog, and "Opening project" stays
  in the task bar (an older core showed a "Data loading error" dialog instead). The scenario that
  reopens therefore waits only on the part that does load (the query's table and entity_types,
  re-read by data sync with their rows), and the known-failure scenario after it claims the defect
  itself: "result" open with its rows and the three table views. The last scenario reloads the
  page: the open that never completes would otherwise leave its "Opening project" entry in the task
  bar, and a later feature on the same page that claims "the task bar should have finished" that
  entry would wait on it.

  Why a feature of its own rather than a scenario of projects-complex: the defect needs only a query,
  a database table and their join, and its reopen leaves the page with a hung open, which a long
  journey would have to recover from before its later scenarios.

  Fixtures: the query reads public.entity_types of System:Datagrok (readable by every user); its row
  count differs between servers, so it is remembered at the first run and every later count is
  claimed against it (the join on id keeps every row). Close All from the left sidebar's context
  menu is done through the shell (closing views is not the claim). The query (under both names) and
  the project with its tables and views are named with the run's time and removed at the start and
  at the end. @serial: the Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDRenJoinProj{time}" is on the server
    And no query named "BDDRenJoinQ{time}, BDDRenJoinQR{time}" is on the server

  Scenario: A join of entity_types with a query's result is saved with Data sync, and the query is renamed in the tree
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    And code editor should hold the code "select * from public.entity_types"
    When user enters "BDDRenJoinQ{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDRenJoinQ{time}" should be on the server
    When user closes the current view
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree
    Then the "BDDRenJoinQ{time}" view should be current
    When user remembers the row count of the table as "rename join query"
    Given the browse panel is open
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    And table "entity_types" should have the "rename join query" row count
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "entity_types" in join left table selector
    And user selects "BDDRenJoinQ{time}" in join right table selector
    And user picks column "id" in join left key selector
    And user picks column "id" in join right key selector
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have the "rename join query" row count
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the following elements should be checked:
      | Data sync switch in "BDDRenJoinQ{time}" project table in "Save project" dialog |
      | Data sync switch in "entity_types" project table in "Save project" dialog      |
      | Data sync switch in "result" project table in "Save project" dialog            |
    When user enters "BDDRenJoinProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRenJoinProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRenJoinProj{time}" dialog
    Then the "Share BDDRenJoinProj{time}" dialog should close
    And the "BDDRenJoinProj{time}" project on the server should hold the tables "BDDRenJoinQ{time}, entity_types, result"
    And the "result" table of the "BDDRenJoinProj{time}" project should be saved with data sync
    And the creation script of the "result" table of the "BDDRenJoinProj{time}" project on the server should contain 'JoinTables("entity_types", "BDDRenJoinQ{time}"'
    And no errors but the project preview's should have been logged
    When user closes all views
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    # the double-click that ran the query also unfolded its node over the runs it logged, and a right-click
    # aimed at the unfolded node can land on one of those rows (a function call's menu)
    When user collapses Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree
    And user picks "Rename..." from the context menu of Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "BDDRenJoinQR{time}" into Name input in "Rename dataquery" dialog
    And user clicks on OK button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDRenJoinQR{time}" should be on the server
    And 0 queries named "BDDRenJoinQ{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: After the rename, the join's creation script names the new query while the table keeps the old name
    Then the creation script of the "BDDRenJoinQ{time}" table of the "BDDRenJoinProj{time}" project on the server should contain ':BDDRenJoinQR{time}()'
    And the creation script of the "result" table of the "BDDRenJoinProj{time}" project on the server should contain 'JoinTables("entity_types", "BDDRenJoinQR{time}"'
    And the "BDDRenJoinProj{time}" project on the server should hold the tables "BDDRenJoinQ{time}, entity_types, result"

  Scenario: Reopened from the Dashboards gallery, the query's table and entity_types are re-read
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRenJoinProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRenJoinProj{time} gallery card
    # the open never completes on GROK-21026, so the barrier is the two tables that do load, not the task bar
    Then table "BDDRenJoinQ{time}" should have been reloaded by data sync with the "rename join query" row count
    And table "entity_types" should have been reloaded by data sync with the "rename join query" row count

  @known-failure @realizes:GROK-21026
  Scenario: The reopened project holds the join over the renamed query's result (GROK-21026)
    Then table "result" should have been reloaded by data sync with the "rename join query" row count
    And the table views "BDDRenJoinQ{time}, entity_types, result" should be open

  Scenario: The page is reloaded, dropping the open that never completed
    When user reloads the page
    Then the table view tabs should read ""
