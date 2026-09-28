@journey @serial @realizes:views.projects
Feature: Projects built from a database: saved, reopened, shared and renamed
  Two projects are built from the platform's own database (Browse > Databases > Postgres >
  Datagrok): one from a saved query run from the tree, one from the table public.entity_types
  opened with Get All. Each is saved through the ribbon's Save dialog with Data sync on and its
  creation script shown, reopened from the Dashboards gallery, shared through the tile's Share...
  with the second account, which opens it too, and removed through Delete Project; the query-based
  one is also renamed and reopened under its new name. Every reopen is claimed by the row count of
  the first open and by the data-sync mark: the table was re-read from the database, not loaded
  from a snapshot. Translated from the TestTrack case Projects/projects-lifecycle-db.

  Fixtures substituted: NorthwindTest (a dev-only connection) is replaced by System:Datagrok, which
  every user may read and query; the md's PostgresAll query by a query of the feature's own,
  "BDDLifeDbQ{time}" over public.entity_types, saved through the JS API (the query editor is not
  this case's subject; projects-lifecycle-query makes its query through the UI); orders by
  public.entity_types. The md's 830 rows become the count read at the first open, since a table of
  the Datagrok database differs between servers.

  Everything is named with the run's time: the query, the three project names. The projects (with
  their tables and views) and the query are removed at the start and at the end; the grant the
  Share dialog adds lives on the project and the query and goes with them (the share adds nothing
  to the System:Datagrok connection, checked on localhost 2026-09-25). @serial: the Dashboards
  search is shared with every feature that saves a project.

  Not translated, and why: Logout and signing in with the second user's credentials — the platform's
  Logout ends every session of the account, which all workers of a run share, so the second account
  is entered through its own session ("user signs in as the sharing user"); Close All from the left
  sidebar's context menu is done through the shell (closing views is not the claim).

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeDbQuery{time}" is on the server
    And no project named "BDDLifeDbQueryRenamed{time}" is on the server
    And no project named "BDDLifeDbTable{time}" is on the server
    And no query named "BDDLifeDbQ{time}" is on the server
    And a query "BDDLifeDbQ{time}" on the Datagrok connection is:
      """
      select * from public.entity_types
      """

  Scenario: The saved query runs from the tree into a table view
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDLifeDbQ{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the "BDDLifeDbQ{time}" view should be current
    When user remembers the row count of the table as "entity types"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Saved with Data sync, the creation script calls the query
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDLifeDbQ{time}" project table in "Save project" dialog should be checked
    When user enters "BDDLifeDbQuery{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "BDDLifeDbQ{time}" project table in "Save project" dialog
    Then creation script text in "BDDLifeDbQ{time}" project table in "Save project" dialog should be visible
    And creation script text in "BDDLifeDbQ{time}" project table in "Save project" dialog should contain text ":BDDLifeDbQ{time}()"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeDbQuery{time}" uploaded' should have been shown
    And 1 project named "BDDLifeDbQuery{time}" should be on the server
    And the "BDDLifeDbQ{time}" table of the "BDDLifeDbQuery{time}" project should be saved with data sync
    And "Share BDDLifeDbQuery{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeDbQuery{time}" dialog
    Then the "Share BDDLifeDbQuery{time}" dialog should close
    And no errors should have been logged

  Scenario: The query project reopens from Dashboards by re-running the query
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbQuery{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeDbQuery{time} gallery card
    Then the current view should be a TableView view
    And the "BDDLifeDbQ{time}" view should be current
    And the table should have the "entity types" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The query project is shared with the second account
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbQuery{time}" into gallery search
    And user picks "Share..." from the context menu of BDDLifeDbQuery{time} gallery card
    Then "Share BDDLifeDbQuery{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDLifeDbQuery{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeDbQuery{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDLifeDbQuery{time}" dialog
    Then the "Share BDDLifeDbQuery{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDLifeDbQuery{time} gallery card
    Then the context panel should show "BDDLifeDbQuery{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of BDDLifeDbQuery{time} gallery card
    Then the access level of the sharing user in "Share BDDLifeDbQuery{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDLifeDbQuery{time}" dialog
    Then the "Share BDDLifeDbQuery{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account opens the shared query project
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbQuery{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeDbQuery{time} gallery card
    Then the current view should be a TableView view
    And the "BDDLifeDbQ{time}" view should be current
    And the table should have the "entity types" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Renamed, the query project opens under its new name
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbQuery{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDLifeDbQuery{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDLifeDbQueryRenamed{time}" into Name input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeDbQueryRenamed{time}" should be on the server
    And 0 projects named "BDDLifeDbQuery{time}" should be on the server
    When user enters "BDDLifeDbQueryRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeDbQueryRenamed{time} gallery card
    Then the current view should be a TableView view
    And the table should have the "entity types" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the query project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbQueryRenamed{time}" into gallery search
    Then BDDLifeDbQueryRenamed{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDLifeDbQueryRenamed{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDLifeDbQueryRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeDbQueryRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDLifeDbQueryRenamed{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Get All opens the database table into a table view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get All" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a TableView view
    And the "entity_types" view should be current
    And the table should have the "entity types" row count
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Saved with Data sync, the creation script reads the table from the database
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And Data sync switch in "entity_types" project table in "Save project" dialog should be checked
    When user enters "BDDLifeDbTable{time}" into Name text input in "Save project" dialog
    And user clicks on "Creation script" button in "entity_types" project table in "Save project" dialog
    Then creation script text in "entity_types" project table in "Save project" dialog should be visible
    And creation script text in "entity_types" project table in "Save project" dialog should contain text "DbQuery(System:Datagrok, \"public.entity_types\""
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeDbTable{time}" uploaded' should have been shown
    And 1 project named "BDDLifeDbTable{time}" should be on the server
    And the "entity_types" table of the "BDDLifeDbTable{time}" project should be saved with data sync
    And "Share BDDLifeDbTable{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeDbTable{time}" dialog
    Then the "Share BDDLifeDbTable{time}" dialog should close
    And no errors should have been logged

  Scenario: The table project reopens from Dashboards by re-reading the table
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbTable{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeDbTable{time} gallery card
    Then the current view should be a TableView view
    And the "entity_types" view should be current
    And the table should have the "entity types" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table project is shared with the second account
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbTable{time}" into gallery search
    And user picks "Share..." from the context menu of BDDLifeDbTable{time} gallery card
    Then "Share BDDLifeDbTable{time}" dialog should be visible
    And "Share BDDLifeDbTable{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeDbTable{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDLifeDbTable{time}" dialog
    Then the "Share BDDLifeDbTable{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDLifeDbTable{time} gallery card
    Then the context panel should show "BDDLifeDbTable{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of BDDLifeDbTable{time} gallery card
    Then the access level of the sharing user in "Share BDDLifeDbTable{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDLifeDbTable{time}" dialog
    Then the "Share BDDLifeDbTable{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account opens the shared table project
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbTable{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeDbTable{time} gallery card
    Then the current view should be a TableView view
    And the "entity_types" view should be current
    And the table should have the "entity types" row count
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the table project
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeDbTable{time}" into gallery search
    Then BDDLifeDbTable{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDLifeDbTable{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDLifeDbTable{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeDbTable{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDLifeDbTable{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
