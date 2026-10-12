@journey @serial @realizes:views.projects
Feature: Projects regressions: projects of a parameterized query
  Regression guards for bugs about projects built on a parameterized query: GROK-20006 (a
  changed parameter is not saved to the project or its copy), core#463 (the project does not open
  after a parameter was removed from its query) and GROK-20174 (a misleading "related entities were
  deleted" error after a parameter was added). The reports use Northwind and Dbtests queries; here
  each scenario makes its own query on System:Datagrok over public.entity_types (the platform's
  entity type names, the same on every server): one row per type name the parameter asks for.

  The query is prepared through the JS API (it sets the scene); everything the regression is about
  goes through the UI: the query is run from Browse > Databases > Postgres > Datagrok, its
  parameter changed in Toolbox > Source and refreshed, projects saved through the ribbon's Save
  dialog and reopened from the Dashboards gallery, and the query's SQL changed in its editor (Edit...
  from the Browse tree). Each claim is what the reopened project shows: its one row, read in the
  grid. What the saved creation script holds (the value of the parameter) is an API test (request
  document, section 5, E1, E2).

  GROK-20174 (done) asked only for an accurate error message: a project whose query gained a
  parameter does not open, and its Data loading error names the parameters ("2 required input
  parameters, 1 passed" on localhost 1.28.0, 2026-10-12) instead of "related entities were deleted";
  that is what its scenario claims. core#463 (open): a project whose query lost a parameter does
  not open either ("1 input parameters, 2 passed"); the reopen is a known failure claiming the
  project opens with its row, set up by the scenario before it. GROK-20006's second case
  (parameters carried by a copied link, then saved) is not translated: no scenario here opens the
  project through a link. The feature is a journey, so the known failure never runs without its
  setup.

  The console is not claimed after a save: the Save dialog's preview logs "Unable to find element in
  cloned iframe" for some views (GROK-18606, won't fix); after each reopen it is. The query and every
  project are named with the run's time and removed when each scenario starts and when the feature
  ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-20006
  Scenario: A parameter changed in Toolbox > Source is saved to the project and to its copy
    Given no project named "BDDRegParamProj{time}" is on the server
    And no project named "BDDRegParamCopy{time}" is on the server
    And a query "BDDRegParams{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegParams{time} tree node inside browse tree
    Then the "BDDRegParams{time}" view should be current
    And the "rows" reading of grid should be 1
    And the "text of cell 1 of name" reading of grid should be "Project"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamProj{time}" dialog
    Then the "Share BDDRegParamProj{time}" dialog should close
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegParamProj{time} gallery card
    Then the "BDDRegParams{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Project"
    Given the toolbox pane is shown
    When user enters "Script" into "Type Name" input in Source pane in toolbox
    And user clicks on REFRESH button in Source pane in toolbox
    Then the "text of cell 1 of name" reading of grid should be "Script"
    And the "rows" reading of grid should be 1
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegParamProj{time}" uploaded' should have been shown
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegParamProj{time} gallery card
    Then the "BDDRegParams{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Script"
    And the table should have been reloaded by data sync
    Given the toolbox pane is shown
    Then "Type Name" input in Source pane in toolbox should have value "Script"
    When user enters "User" into "Type Name" input in Source pane in toolbox
    And user clicks on REFRESH button in Source pane in toolbox
    Then the "text of cell 1 of name" reading of grid should be "User"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDRegParamCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegParamCopy{time}" uploaded' should have been shown
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegParamCopy{time} gallery card
    Then the "BDDRegParams{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "User"
    And the table should have been reloaded by data sync
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegParamProj{time} gallery card
    Then the "BDDRegParams{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Script"
    And the table should have been reloaded by data sync
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:core#463
  Scenario: A project of a query is saved, and the query loses one of its parameters in its editor
    Given no project named "BDDRegParamLess{time}" is on the server
    And a query "BDDRegParamsLess{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: string pattern = "%"
      select name from entity_types where name = @typeName and name like @pattern
      """
    When user closes all views
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree
    Then the "BDDRegParamsLess{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Project"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamLess{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamLess{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamLess{time}" dialog
    Then the "Share BDDRegParamLess{time}" dialog should close
    When user closes all views
    When user collapses Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree
    And user picks "Edit..." from the context menu of Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user replaces the code of code editor with '--input: string typeName = "Project"'
    And user appends "select name from entity_types where name = @typeName" to code editor
    And user clicks on Save button
    And user closes the current view
    # opened again from the tree, the editor shows what was saved
    And user picks "Edit..." from the context menu of Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    And code editor should contain the text '--input: string typeName = "Project"'
    And code editor should not contain the text "pattern"
    When user closes the current view
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamLess{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDRegParamLess{time} gallery card should be visible

  @known-failure @realizes:core#463
  Scenario: The project opens after its query lost a parameter (core#463)
    When user double-clicks on BDDRegParamLess{time} gallery card
    Then the "BDDRegParamsLess{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Project"
    And "Data loading error" dialog should be absent

  @realizes:GROK-20174
  Scenario: A project whose query gained a parameter names the parameters, not "related entities were deleted"
    Given no project named "BDDRegParamMore{time}" is on the server
    And a query "BDDRegParamsMore{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    When user closes all views
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree
    Then the "BDDRegParamsMore{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Project"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamMore{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamMore{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamMore{time}" dialog
    Then the "Share BDDRegParamMore{time}" dialog should close
    When user closes all views
    When user collapses Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree
    And user picks "Edit..." from the context menu of Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user replaces the code of code editor with '--input: string pattern = "%"'
    And user appends '--input: string typeName = "Project"' to code editor
    And user appends "select name from entity_types where name = @typeName and name like @pattern" to code editor
    And user clicks on Save button
    And user closes the current view
    # opened again from the tree, the editor shows what was saved
    And user picks "Edit..." from the context menu of Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    And code editor should contain the text '--input: string pattern = "%"'
    When user closes the current view
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamMore{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegParamMore{time} gallery card
    Then "Data loading error" dialog should be visible
    And "Data loading error" dialog should contain text "BDDRegParamsMore{time}"
    And "Data loading error" dialog should contain text "required input parameters"
    And "Data loading error" dialog should not contain text "related entities"
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close
