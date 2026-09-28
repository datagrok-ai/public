@serial @realizes:views.projects
Feature: Projects regressions: a project of a visual query with a parameter
  GROK-20535 (open): a visual query with a parameter, run with another value than the default and
  saved in a project, reopens with the default value. The report used Northwind; here the visual
  query is built in the UI on System:Datagrok's public.entity_types (Browse > Databases > Postgres
  > Datagrok > Schemas > public, New Visual Query...): its Where row takes "name" with the condition
  Project exposed as a parameter. It is saved, run from Toolbox > Run query... with Script, the
  result saved as a project through the ribbon's Save dialog and reopened from the Dashboards
  gallery.

  Reproduced on localhost, core 1.28.0 bc64f40e47 (a probe and two runs of this feature): the saved
  creation script carries the condition's default ("pattern":"Project"), and the reopened table holds
  the Project row. The whole setup, the reopen included, is the Background; the scenario reads the
  reopened row.

  After the condition and its parameter, one more property is set before Save: the Order by row
  takes "name" (one row, so the order changes no claim). The condition reaches the builder's query
  750 ms after its last keystroke (a debounced input), and the parameter checkbox with the builder's
  next run; the Order by gesture comes after both, and the builder's row count after it (1, the
  Project row) is read before Save, since the gesture's own duration is not bounded against the
  750 ms. The saved query is read back from the server holding the condition and its parameter.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after the save lets
  that one message through. The query and the project are named with the run's time and removed when
  the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And user clears the saved pivot table parameters
    And no query named "BDDRegVQ{time}" is on the server
    And no project named "BDDRegVQProj{time}" is on the server
    And the layout saved for the query "BDDRegVQ{time}" is deleted at the end
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New Visual Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity-types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user adds "name" to the "Where" row of the visual query
    And user enters "Project" into visual query name condition
    And user presses Enter
    And user checks visual query name parameter checkbox
    Then visual query name parameter checkbox should be checked
    When user adds "name" to the "Order by" row of the visual query
    Then the "Order by" row of the visual query should hold "name"
    And the visual query should have run to 1 row
    When user enters "BDDRegVQ{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDRegVQ{time}" should be on the server
    And the query "BDDRegVQ{time}" on the server should filter "name" by "Project" as a parameter
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then "BDDRegVQ{time}" dialog should be visible
    And Name text input in "BDDRegVQ{time}" dialog should have value "Project"
    When user enters "Script" into Name text input in "BDDRegVQ{time}" dialog
    And user clicks on OK button in "BDDRegVQ{time}" dialog
    Then the "BDDRegVQ{time}" dialog should close
    And the current view should be a TableView view
    And the only row of the table should read "Script" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegVQProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegVQProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegVQProj{time}" dialog
    Then the "Share BDDRegVQProj{time}" dialog should close
    And the "BDDRegVQ{time}" table of the "BDDRegVQProj{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegVQProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegVQProj{time} gallery card
    Then the task bar should have finished "Opening project"
    And the current view should be a TableView view
    And the table should have been reloaded by data sync

  @known-failure @realizes:GROK-20535
  Scenario: The reopened visual query project holds the row of the value it was run with
    Then the only row of the table should read "Script" in the "name" column
