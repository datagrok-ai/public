@journey @realizes:views.projects
Feature: Projects regressions: a project of a visual query with a parameter
  GROK-20535 (open): a visual query with a parameter, run with another value than the default and
  saved in a project, reopens with the default value. The report used Northwind; here the visual query
  is built in the UI on System:Datagrok's public.entity_types (Browse > Databases > Postgres >
  Datagrok > Schemas > public, New Visual Query...): its Where row takes "name" with the condition
  Project, exposed as a parameter, and its Order by row takes "name" (one row, so the order changes no
  claim). It is saved, run from Toolbox > Run query... with Script, the result saved as a project
  through the ribbon's Save dialog and reopened from the Dashboards gallery.

  The parameter is claimed where the UI shows it: its box checked in the builder, and the Run query
  dialog offering the Name field with the default Project. The reopen is the known failure: the
  reopened table holds the Project row, the default, where the project was saved with Script.

  The console is not claimed after the save: the Save dialog's preview logs "Unable to find element
  in cloned iframe" for some views (GROK-18606, won't fix). The query and the project are named with
  the run's time (letters and digits only: the Dashboards search misses "-" and "_") and removed when
  the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDDRegVQ{time}" is on the server
    And no project named "BDDRegVQProj{time}" is on the server

  Scenario: A visual query on entity_types filters name by a parameter
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New Visual Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user adds "name" to the "Where" row of the visual query
    And user enters "Project" into "name" where condition
    And user presses Enter
    And user checks "name" where parameter checkbox
    Then "name" where parameter checkbox should be checked
    When user adds "name" to the "Order by" row of the visual query
    Then the "Order by" row of the visual query should hold "name"
    And the visual query should preview 1 row
    When user enters "BDDRegVQ{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDRegVQ{time}" should be on the server

  Scenario: Run with another value, the result is saved as a project with Data sync
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then "BDDRegVQ{time}" dialog should be visible
    And Name text input in "BDDRegVQ{time}" dialog should have value "Project"
    When user enters "Script" into Name text input in "BDDRegVQ{time}" dialog
    And user clicks on OK button in "BDDRegVQ{time}" dialog and the query it runs completes
    Then the "BDDRegVQ{time}" dialog should close
    And the current view should be a TableView view
    And the table should have 1 row
    And the "text of cell 1 of name" reading of grid should be "Script"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDRegVQ{time}" project table in "Save project" dialog should be checked
    When user enters "BDDRegVQProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegVQProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegVQProj{time}" dialog
    Then the "Share BDDRegVQProj{time}" dialog should close
    And the "BDDRegVQ{time}" table of the "BDDRegVQProj{time}" project should be saved with data sync

  @known-failure:77 @realizes:GROK-20535
  Scenario: The reopened project holds the row of the value it was run with
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegVQProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegVQProj{time} gallery card and the project opens
    Then the current view should be a TableView view
    And the table should have been reloaded by data sync
    And the "text of cell 1 of name" reading of grid should be "Script"
