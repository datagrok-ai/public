@journey @serial @realizes:views.projects
Feature: A project built from the user's own query survives renaming the query
  A query is made in the Browse tree with New SQL Query... on a table of the platform's own database
  (Browse > Databases > Postgres > Datagrok), saved, run from the tree, and its result saved as a
  project with Data sync on. The project is shared with the second account through the Dashboards
  tile's Share..., with notifications off; the second account opens it. Then the query is renamed
  (github-3550): the owner and the second account both open the project again and get the same
  rows re-read from the database. At last the project is renamed and reopened, deleted through
  Delete Project, and the query through its Delete. Translated from the TestTrack case
  Projects/projects-lifecycle-query.

  Fixtures substituted: NorthwindTest (a dev-only connection) is replaced by System:Datagrok, which
  every user may read and query, and its orders table by public.entity_types. The md's 830 rows
  become the count read at the first run, since a table of the Datagrok database differs between
  servers; every reopen is claimed against it and by the data-sync mark (re-read, not a snapshot).

  The query and the project are named with the run's time and removed, under both their names, at
  the start and at the end; the grant the Share dialog adds lives on the project and the query and
  goes with them. @serial: the Dashboards search is shared with every feature that saves a project.

  Not translated, and why: Logout and signing in with the second user's credentials — the platform's
  Logout ends every session of the account, which all workers of a run share, so the second account
  is entered through its own session ("user signs in as the sharing user"); Close All from the left
  sidebar's context menu is done through the shell (closing views is not the claim); "Load more"
  under the connection is not needed — the query's name sorts among the first nodes listed.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeQueryProj{time}" is on the server
    And no project named "BDDLifeQueryProjRenamed{time}" is on the server
    And no query named "BDDLifeQuery{time}" is on the server
    And no query named "BDDLifeQueryRenamed{time}" is on the server

  Scenario: New SQL Query... on a table makes the query, saved under the feature's name
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    And Name input should have value "entity_types"
    And code editor should hold the code "select * from public.entity_types"
    When user enters "BDDLifeQuery{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDLifeQuery{time}" should be on the server
    When user closes the current view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query runs from the tree into a table view
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDLifeQuery{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the "BDDLifeQuery{time}" view should be current
    When user remembers the row count of the table as "entity types query"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query result is saved as a project with Data sync
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And there should be 1 visible saved table row
    And "BDDLifeQuery{time}" project table in "Save project" dialog should be visible
    And Data sync switch in "BDDLifeQuery{time}" project table in "Save project" dialog should be checked
    When user enters "BDDLifeQueryProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDLifeQueryProj{time}" should be on the server
    And the "BDDLifeQuery{time}" table of the "BDDLifeQueryProj{time}" project should be saved with data sync
    And "Share BDDLifeQueryProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeQueryProj{time}" dialog
    Then the "Share BDDLifeQueryProj{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Only the project is shared with the second account, without notifications
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDLifeQueryProj{time} gallery card
    Then "Share BDDLifeQueryProj{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDLifeQueryProj{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeQueryProj{time}" dialog
    Then share access selector should contain text "View and use"
    And "Send notifications" input in "Share BDDLifeQueryProj{time}" dialog should be checked
    When user unchecks "Send notifications" input in "Share BDDLifeQueryProj{time}" dialog
    Then "Send notifications" input in "Share BDDLifeQueryProj{time}" dialog should be unchecked
    When user clicks on OK button in "Share BDDLifeQueryProj{time}" dialog
    Then the "Share BDDLifeQueryProj{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDLifeQueryProj{time} gallery card
    Then the context panel should show "BDDLifeQueryProj{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of BDDLifeQueryProj{time} gallery card
    Then the access level of the sharing user in "Share BDDLifeQueryProj{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDLifeQueryProj{time}" dialog
    Then the "Share BDDLifeQueryProj{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account opens the shared project
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeQueryProj{time} gallery card
    Then the current view should be a TableView view
    And the table should have the "entity types query" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query is renamed in the Browse tree
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    # the double-click that ran the query also unfolded its node over the runs it logged, and a right-click
    # aimed at the unfolded node can land on one of those rows (a function call's menu)
    When user collapses Databases---Postgres---Datagrok---BDDLifeQuery{time} tree node inside browse tree
    And user picks "Rename..." from the context menu of Databases---Postgres---Datagrok---BDDLifeQuery{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "BDDLifeQueryRenamed{time}" into Name input in "Rename dataquery" dialog
    And user clicks on OK button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDLifeQueryRenamed{time}" should be on the server
    And 0 queries named "BDDLifeQuery{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner opens the project after the query was renamed (github-3550)
    # a fresh page: the one that renamed the query holds the query entity under its old name
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeQueryProj{time} gallery card
    Then the current view should be a TableView view
    And the table should have the "entity types query" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user closes all views

  Scenario: The second account opens the project after the query was renamed (github-3550)
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeQueryProj{time} gallery card
    Then the current view should be a TableView view
    And the table should have the "entity types query" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Renamed, the project opens under its new name
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProj{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDLifeQueryProj{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDLifeQueryProjRenamed{time}" into Name input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeQueryProjRenamed{time}" should be on the server
    And 0 projects named "BDDLifeQueryProj{time}" should be on the server
    When user enters "BDDLifeQueryProjRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeQueryProjRenamed{time} gallery card
    Then the current view should be a TableView view
    And the table should have the "entity types query" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the project, and Delete the renamed query
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeQueryProjRenamed{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDLifeQueryProjRenamed{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDLifeQueryProjRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeQueryProjRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    # the double-click that ran the query also unfolded its node over the runs it logged, and a right-click
    # aimed at the unfolded node can land on one of those rows (a function call's menu)
    When user collapses Databases---Postgres---Datagrok---BDDLifeQueryRenamed{time} tree node inside browse tree
    When user picks "Delete" from the context menu of Databases---Postgres---Datagrok---BDDLifeQueryRenamed{time} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete query \"BDDLifeQueryRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDDLifeQueryRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    # the connection's children listed again, then read
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    Then Databases---Postgres---Datagrok---Schemas tree node inside browse tree should be visible
    And Databases---Postgres---Datagrok---BDDLifeQueryRenamed{time} tree node inside browse tree should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
