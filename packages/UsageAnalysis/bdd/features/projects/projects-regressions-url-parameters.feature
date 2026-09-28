@serial @realizes:views.projects
Feature: Projects regressions: the URL parameters icon of a just-saved dashboard
  GROK-20929 (open): Toolbox > Source of a saved dashboard built on a parameterized query shows,
  next to REFRESH, a sliders icon ("Choose which parameters the dashboard link carries"); in the tab
  where the project was just saved the icon is missing, and it shows only after the project is
  reopened. The report used a Dbtests query on dev; here the feature makes its own query on
  System:Datagrok over public.entity_types, runs it from Browse > Databases > Postgres > Datagrok,
  saves the project through the ribbon's Save dialog and reopens it from the Dashboards gallery.

  Reproduced on localhost, core 1.28.0 bc64f40e47 (a probe and two runs of this feature): no icon
  after the save, the icon after the reopen. The first scenario, the claim right after the save, is
  the known failure; its setup is the Background. The second one is the reopen, where the icon is
  there — the proof that the icon the first one misses exists at all.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after the save lets
  that one message through. The query and the project are named with the run's time and removed
  when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDDRegIconQ{time}" is on the server
    And no project named "BDDRegIcon{time}" is on the server
    And a query "BDDRegIconQ{time}" on the Datagrok connection is:
      """
      --input: string typeName = "Project"
      select name from entity_types where name = @typeName
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks Databases---Postgres---Datagrok---BDDRegIconQ{time} tree node inside browse tree
    Then the "BDDRegIconQ{time}" view should be current
    And the only row of the table should read "Project" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegIcon{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegIcon{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegIcon{time}" dialog
    Then the "Share BDDRegIcon{time}" dialog should close
    And 1 project named "BDDRegIcon{time}" should be on the server
    And no errors but the project preview's should have been logged
    Given the toolbox pane is shown
    Then REFRESH button inside toolbox should be visible

  @known-failure @realizes:GROK-20929
  Scenario: Right after the save, Toolbox > Source offers the URL parameters icon
    Then url parameters icon should be visible

  @realizes:GROK-20929
  Scenario: Reopened from the gallery, Toolbox > Source offers the URL parameters icon
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegIcon{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegIcon{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "BDDRegIconQ{time}" view should be current
    Given the toolbox pane is shown
    Then url parameters icon should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged
