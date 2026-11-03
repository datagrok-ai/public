@realizes:views.projects
Feature: Projects regressions: a query result renamed from the status bar
  core#626: a query on the Datagrok database's entity_types (made through the API: writing a query is
  the database features' subject) is run from the Browse tree, and its result view is not saved as a
  project. The table name on the left of the status bar puts the table into the context panel, and
  Actions > Rename... there renames it. The table and its view should take the new name, as with
  Rename on the view tab. The scenario is a known failure until they do.

  The query carries the run's time and is removed when the feature starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And the browse panel is open
    And a query "BDDRenameQ{time}" on "System:Datagrok" is on the server:
      """
      select name from entity_types
      """

  # core#626: the rename from the status bar's table does not reach the table or its view
  @known-failure:38 @realizes:core#626
  Scenario: The table renamed through the status bar takes the new name with its view
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDRenameQ{time} tree node inside browse tree and the query it runs completes
    Then the "BDDRenameQ{time}" view should be current
    When user clicks on "BDDRenameQ{time}" text in status bar
    Then the context panel should show "BDDRenameQ{time}"
    Given Actions pane in context panel is expanded
    When user clicks on "Rename..." label in Actions pane in context panel
    Then "Rename table" dialog should be visible
    When user enters "BDDRenamedT{time}" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And the open table views should be exactly "BDDRenamedT{time}"
    And the "BDDRenamedT{time}" view should be current
