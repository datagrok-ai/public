@realizes:views.projects @realizes:views.space
Feature: Projects regressions: a dashboard named like its query, moved into a space
  core#563: a query with a parameter on the Datagrok database's entity_types (made through the API:
  writing a query is the database features' subject) is run from the Browse tree, and its result is
  saved with Data sync through the ribbon's Save dialog under the name the dialog proposes, which is
  the query's own. The dashboard is moved into a space through Move to Space... on its Dashboards
  card, found in the space's gallery (a reload right after OK can come before the move lands, and the
  dashboard then opens as it was), the page is reloaded, and the dashboard is opened from its card. The table's creation script
  should still call the query; the move re-points the call at the dashboard's new namespace, and the
  open fails with a Data loading error. The scenario is a known failure until the dashboard opens.

  The query, the dashboard and the space carry the run's time and are removed when the feature
  starts and ends; the space is made through the API (the Create Space dialog is spaces-create's).

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And the browse panel is open
    And no project named "BDDSameName{time}" is on the server
    And no space named "BDDSameNameSpace{time}" is on the server
    And a space named "BDDSameNameSpace{time}" is on the server
    And a query "BDDSameName{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """

  # core#563: the move rewrites the creation script to call <space>:BDDSameName{time}, the dashboard itself
  @known-failure:67 @realizes:core#563
  Scenario: The dashboard named like its query opens after the move into a space
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDSameName{time} tree node inside browse tree and the query it runs completes
    Then the "BDDSameName{time}" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDSameName{time}" project table in "Save project" dialog should be checked
    And Name text input in "Save project" dialog should have value "BDDSameName{time}"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDSameName{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDSameName{time}" dialog
    Then the "Share BDDSameName{time}" dialog should close
    And the "BDDSameName{time}" table of the "BDDSameName{time}" project should be saved with data sync
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSameName{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Move to Space..." from the context menu of BDDSameName{time} gallery card
    Then "Move to space" dialog should be visible
    When user selects "BDDSameNameSpace{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    Given Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDSameNameSpace{time} tree node inside browse tree
    Then the "BDDSameNameSpace{time}" view should be current
    And "BDDSameName{time}" link in gallery should be visible
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSameName{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDSameName{time} gallery card and the project opens
    Then the "BDDSameName{time}" view should be current
    And "Data loading error" dialog should be absent
    And the table should have 1 rows
