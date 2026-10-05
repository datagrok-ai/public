@dev-only @journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project of the user's own NorthwindTest query survives renaming the query
  Dev only: NorthwindTest (the Postgres connection Dbtests:PostgresTest of the Dbtests package,
  listed as NorthwindTest under Browse > Databases > Postgres) exists only on dev.datagrok.ai; on
  any other stand the tree has no such node and the first scenario fails. Run it with
  DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev, and with the budgets of a slow stand:
  BDD_EXPECT_TIMEOUT=30000 and `run --timeout 300000` (on dev every "no project named" cleanup
  lists all the stand's projects, 35-65 s each, and the Dashboards gallery answers a search in up to
  30 s, past the default 120 s per test and afterAll hook and 15 s per check).

  As written in the md: a query is made in the Browse tree with New SQL Query... on NorthwindTest >
  Schemas > public > orders (select * from public.orders), saved under its own name, run from the
  tree (830 rows), and its result saved as a project with Data sync on. The project alone is shared
  with the second account (notifications off). Then the query is renamed (github-3550): the owner
  opens the project again and gets the 830 rows re-read from the database. The project goes through
  Delete Project and the query through its Delete. Translated from the TestTrack case
  Projects/projects-lifecycle-query.

  The query's SQL is then changed to "limit 100" in its editor and, after the page is reloaded as
  the md does, the owner's reopen shows 100 rows.

  Parked (see the request document): the recipient's opens (before the rename, and after the
  rename and the SQL change). Signed in as the sharing user on dev (1.28.0), the double-click on the
  card leaves "Opening project" in the task bar for over two minutes with the Dashboards view
  current, and no table and no message come — a suspected defect, described outside the
  repository. Kept without one claim, restored with the phrase: that the Save dialog lists only the
  query's table (the project-table rows cannot be counted apart from the dialog's other rows with
  library phrases). The System:Datagrok version for other stands is parked (row counts differ per
  server).

  The query and the project are named with the run's time (letters and digits only) and removed,
  the query under both its names, at the start and at the end. It is serial: the Dashboards search
  is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDNwLifeQProj{time}" is on the server
    And no query named "BDDNwLifeQ{time}, BDDNwLifeQRenamed{time}" is on the server

  Scenario: New SQL Query... on a table makes the query, saved under the feature's name
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the current view should be a DataQueryView view
    And Name input should have value "orders"
    And code editor should hold the code "select * from public.orders"
    When user enters "BDDNwLifeQ{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDNwLifeQ{time}" should be on the server
    When user closes the current view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query runs from the tree into a table view
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the "BDDNwLifeQ{time}" view should be current
    And the table should have 830 rows
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query result is saved as a project with Data sync
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And "BDDNwLifeQ{time}" project table in "Save project" dialog should be visible
    And Data sync switch in "BDDNwLifeQ{time}" project table in "Save project" dialog should be checked
    When user enters "BDDNwLifeQProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDNwLifeQProj{time}" should be on the server
    And the "BDDNwLifeQ{time}" table of the "BDDNwLifeQProj{time}" project should be saved with data sync
    And "Share BDDNwLifeQProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDNwLifeQProj{time}" dialog
    Then the "Share BDDNwLifeQProj{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Only the project is shared with the second account, without notifications
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeQProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDNwLifeQProj{time} gallery card
    Then "Share BDDNwLifeQProj{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDNwLifeQProj{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDNwLifeQProj{time}" dialog
    And user clicks on OK button in "Share BDDNwLifeQProj{time}" dialog
    Then the "Share BDDNwLifeQProj{time}" dialog should close
    When user clicks on BDDNwLifeQProj{time} gallery card
    Then the context panel should show "BDDNwLifeQProj{time}"
    And the sharing pane should list the sharing user

  Scenario: The query is renamed in the Browse tree
    When user closes all views
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the double-click that ran the query also unfolded its node over the runs it logged, and a right-click
    # aimed at the unfolded node can land on one of those rows (a function call's menu)
    When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree
    And user picks "Rename..." from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "BDDNwLifeQRenamed{time}" into Name input in "Rename dataquery" dialog
    And user clicks on OK button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDNwLifeQRenamed{time}" should be on the server
    And 0 queries named "BDDNwLifeQ{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner opens the project after the query was renamed (github-3550)
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeQProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeQProj{time} gallery card
    Then the current view should be a TableView view
    And the table should have 830 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user closes all views

  Scenario: With the query's SQL changed, the reopened project shows the new result
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree
    And user picks "Edit..." from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user replaces the code of code editor with "select * from public.orders limit 100"
    And user clicks on Save button
    Then code editor should hold the code "select * from public.orders limit 100"
    When user closes the current view
    And user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeQProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDNwLifeQProj{time} gallery card
    Then the current view should be a TableView view
    And the table should have 100 rows
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the project, and Delete the renamed query
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDNwLifeQProj{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDNwLifeQProj{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDNwLifeQProj{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDNwLifeQProj{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the double-click that ran the query also unfolded its node over the runs it logged, and a right-click
    # aimed at the unfolded node can land on one of those rows (a function call's menu)
    When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree
    When user picks "Delete" from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete query \"BDDNwLifeQRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDDNwLifeQRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    # the connection's children listed again, then read
    Given Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree should be visible
    And Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
