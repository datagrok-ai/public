@serial @realizes:views.projects
Feature: Projects regressions: the share link of the Save project dialog
  GROK-20930 (open): for a table of a parameterized query the Save project dialog shows a "Share
  link:" line built from the project name the dialog opened with; typing another name does not
  rebuild it, so the link offered points at a project that will not exist. The report used a
  Dbtests query on dev; here the feature makes its own query on System:Datagrok over
  public.entity_types, runs it from Browse > Databases > Postgres > Datagrok and opens the Save
  dialog from the ribbon.

  Reproduced on localhost, core 1.28.0 bc64f40e47 (a probe and two runs of this feature): the link
  keeps the query's name after a new name is typed, so the scenario is a known failure. The whole
  setup, the typed name included, is the Background; the scenario reads the link. Nothing is saved:
  the dialog is left open and closed by the shell reset; the query is named with the run's time and
  removed when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDDRegLinkQ{time}" is on the server
    And a query "BDDRegLinkQ{time}" on the Datagrok connection is:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegLinkQ{time} tree node inside browse tree
    Then the "BDDRegLinkQ{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And save dialog share link should contain text "BDDRegLinkQ{time}?typeName=Project"
    When user enters "BDDRegShareName{time}" into Name text input in "Save project" dialog
    And user presses Tab
    Then Name text input in "Save project" dialog should have value "BDDRegShareName{time}"

  @known-failure @realizes:GROK-20930
  Scenario: The share link follows the name typed into the Save dialog
    Then save dialog share link should contain text "BDDRegShareName{time}?typeName=Project"
