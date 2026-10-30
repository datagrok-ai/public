@realizes:views.projects
Feature: Projects regressions: OK clicked as soon as the Save dialog opens
  core#945 (done): a query result saved by clicking OK in the ribbon's Save dialog right after it
  opened, before the dialog listed the table with its Data sync switch, was saved with no tables, and
  its reopen opened nothing. A query with a parameter on the Datagrok database's entity_types (made
  through the API: writing a query is the database features' subject) is run from the Browse tree in
  a fresh session of the page, as the ticket did, the Save dialog's OK is clicked as soon as the
  dialog is shown, and the dashboard is reopened from its Dashboards card with its table and rows.

  The query and the dashboard carry the run's time and are removed when the feature starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And no project named "BDDQuickSave{time}" is on the server
    And a query "BDDQuickSaveQ{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Script"
      select name from entity_types where name = @typeName
      """

  @realizes:core#945
  Scenario: The dashboard saved with an immediate OK reopens with its table
    When user reloads the page
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDQuickSaveQ{time} tree node inside browse tree and the query it runs completes
    Then the "BDDQuickSaveQ{time}" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDQuickSave{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDQuickSave{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDQuickSave{time}" dialog
    Then the "Share BDDQuickSave{time}" dialog should close
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDQuickSave{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDQuickSave{time} gallery card and the project opens
    Then the "BDDQuickSaveQ{time}" view should be current
    And the table should have 1 rows
    And the table should have been reloaded by data sync
