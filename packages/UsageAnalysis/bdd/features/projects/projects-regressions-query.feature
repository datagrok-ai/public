@serial @realizes:views.projects
Feature: Projects regressions: projects of a parameterized query
  Regression guards for Jira bugs about projects built on a parameterized query: GROK-20006 (a
  changed parameter is not saved to the project or its copy), GROK-20627 (the project does not open
  after a parameter was removed from its query) and GROK-20174 (a misleading "related entities were
  deleted" error after a parameter was added). The reports use Northwind and Dbtests queries; here
  the feature makes its own query on System:Datagrok over public.entity_types (the platform's entity
  type names, the same on every server): one row per type name the parameter asks for. The query is
  run from Browse > Databases > Postgres > Datagrok, its parameter changed in Toolbox > Source,
  projects saved through the ribbon's Save dialog and reopened from the Dashboards gallery. Editing
  the query's SQL is done through the JS API (the edit prepares the scene; the claim is the reopen).

  What each scenario fails on: the project or its copy saved without data sync, with a creation
  script that lacks the changed value, or reopening with the old parameter value and row
  (GROK-20006); the reopen erroring or not opening the table after the query lost or gained a
  parameter (GROK-20627, GROK-20174). GROK-20174 asked only for an accurate error message; the
  product now opens such a project with the query's new parameter at its default, and that is the
  outcome the scenario holds it to: a refusal, even with an accurate message, fails it.

  GROK-20006's second case (parameters carried by a copied link, then saved) is not translated: no
  scenario here opens the project through a link.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after a save lets
  that one message through; every check after a reopen is strict. The query and every project are
  named with the run's time and removed when each scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-20006
  Scenario: A parameter changed in Toolbox > Source is saved to the project and to its copy
    Given no query named "BDDRegParams{time}" is on the server
    And no project named "BDDRegParamProj{time}" is on the server
    And no project named "BDDRegParamCopy{time}" is on the server
    And a query "BDDRegParams{time}" on the Datagrok connection is:
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
    And the only row of the table should read "Project" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamProj{time}" dialog
    Then the "Share BDDRegParamProj{time}" dialog should close
    And the "BDDRegParams{time}" table of the "BDDRegParamProj{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegParamProj{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParams{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    Given the toolbox pane is shown
    When user enters "Script" into "Type Name" text input inside toolbox
    And user clicks on REFRESH button inside toolbox
    Then the only row of the table should read "Script" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegParamProj{time}" uploaded' should have been shown
    And the "BDDRegParams{time}" table of the "BDDRegParamProj{time}" project should be saved with data sync
    And the creation script of the "BDDRegParams{time}" table of the "BDDRegParamProj{time}" project on the server should contain '("Script")'
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    Given user watches the task bar
    When user double-clicks on BDDRegParamProj{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParams{time}" view should be current
    And the table should have been reloaded by data sync
    And the only row of the table should read "Script" in the "name" column
    And "Type Name" text input inside toolbox should have value "Script"
    When user enters "User" into "Type Name" text input inside toolbox
    And user clicks on REFRESH button inside toolbox
    Then the only row of the table should read "User" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDRegParamCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegParamCopy{time}" uploaded' should have been shown
    And the "BDDRegParams{time}" table of the "BDDRegParamCopy{time}" project should be saved with data sync
    And the creation script of the "BDDRegParams{time}" table of the "BDDRegParamCopy{time}" project on the server should contain '("User")'
    And the creation script of the "BDDRegParams{time}" table of the "BDDRegParamProj{time}" project on the server should contain '("Script")'
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegParamCopy{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParams{time}" view should be current
    And the table should have been reloaded by data sync
    And the only row of the table should read "User" in the "name" column
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamProj{time}" into gallery search
    Given user watches the task bar
    When user double-clicks on BDDRegParamProj{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParams{time}" view should be current
    And the table should have been reloaded by data sync
    And the only row of the table should read "Script" in the "name" column
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:GROK-20627
  Scenario: A project opens after one of its query's parameters was removed
    Given no query named "BDDRegParamsLess{time}" is on the server
    And no project named "BDDRegParamLess{time}" is on the server
    And a query "BDDRegParamsLess{time}" on the Datagrok connection is:
      """
      --input: string typeName = "Project"
      --input: string pattern = "%"
      select name from entity_types where name = @typeName and name like @pattern
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree
    Then the "BDDRegParamsLess{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamLess{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamLess{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamLess{time}" dialog
    Then the "Share BDDRegParamLess{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And the query "BDDRegParamsLess{time}" on the server is changed to:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamLess{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegParamLess{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParamsLess{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:GROK-20174
  Scenario: A project opens after its query gained a parameter
    Given no query named "BDDRegParamsMore{time}" is on the server
    And no project named "BDDRegParamMore{time}" is on the server
    And a query "BDDRegParamsMore{time}" on the Datagrok connection is:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree
    Then the "BDDRegParamsMore{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegParamMore{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegParamMore{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegParamMore{time}" dialog
    Then the "Share BDDRegParamMore{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And the query "BDDRegParamsMore{time}" on the server is changed to:
      """
      --input: string pattern = "%"
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName and name like @pattern
      """
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegParamMore{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegParamMore{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegParamsMore{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    And no error or warning balloon should have been shown
    And no errors should have been logged
