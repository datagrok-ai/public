@realizes:views.projects @realizes:views.functions
Feature: A dashboard on a query with a list parameter: the Source link after the save
  core#567: a query with a string parameter (typeName = "Project") and a list<string> parameter
  (otherNames, left empty) on the Datagrok database's entity_types (made through the API: writing a
  query is the database features' subject) is run from the Browse tree. Before the save, the copy
  icon of Toolbox > Source gives a /func/ link with typeName. The result is saved as a project
  through the ribbon's Save dialog under its own name. On a query whose parameters are all scalar the
  first save exposes them as URL parameters and the Source link of the dashboard carries them
  (projects-url-parameters.feature); here the list parameter turns that off for the scalar one too,
  and the link is the bare /p/ address. The scenario is a known failure until the link carries
  typeName.

  Names are letters and digits only and carry the run's time. The query and the project are removed
  when the feature starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And the browse panel is open
    And no project named "BDDUrlListProj{time}" is on the server
    And a query "BDDUrlListQuery{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: list<string> otherNames
      select name from entity_types where name = @typeName or name = any(@otherNames)
      """

  # core#567: one list<string> parameter leaves the first save without URL parameters
  @known-failure:54 @realizes:core#567
  Scenario: After the save, the Source link still carries the scalar parameter
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDUrlListQuery{time} tree node inside browse tree and the query it runs completes
    Then the "BDDUrlListQuery{time}" view should be current
    Given the toolbox pane is shown
    When user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "/func/"
    And tooltip should contain text "typeName=%22Project%22"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDUrlListProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDUrlListProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDUrlListProj{time}" dialog
    Then the "Share BDDUrlListProj{time}" dialog should close
    When user moves the pointer away from Source pane in toolbox
    And user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "/p/"
    And tooltip should contain text ".BDDUrlListProj{time}"
    And tooltip should contain text "typeName=Project"
