@journey @serial @realizes:views.projects
Feature: Projects regressions: a project of a visual query with a pivot and a parameter
  GROK-20273 (done): a project of a visual query with a parameter and a pivot, run with another
  value than the default and saved, raised an error on reopening. The report used Northwind; here
  the visual query is built in the UI on System:Datagrok's public.entity_types (Browse > Databases
  > Postgres > Datagrok > Schemas > public, New Visual Query...): the Where row takes "name" with the
  condition Project exposed as a parameter, grouped by name and pivoted on is_package_entity with an
  aggregate of id. It is saved, run from Toolbox > Run query... with Script, the result saved as a
  project through the ribbon's Save dialog and reopened from the Dashboards gallery; the claim is the
  reopen without an error.

  The ticket's trigger is the non-default value reaching the reopen. While GROK-20535 (open) stands,
  it does not: the saved creation script carries the condition's default ("pattern":"Project"), so
  the reopen re-runs the pivot with Project, and "1 row, no error" holds for either value. The first
  scenario therefore guards only the error-free reopen of a pivoted visual query with a parameter;
  the second claims the reopened row is the Script one, and is a known failure linked to GROK-20535
  (reproduced on localhost, core 1.28.0 bc64f40e47, in projects-regressions-visual-query.feature).
  When GROK-20535 is fixed, the second scenario passes, the tag has to go, and from then on the
  whole path of GROK-20273 is under test. The feature is a journey: the known failure only reads
  the table the first scenario reopened and proved reloaded by data sync, so it never runs without
  that setup, and a failed setup fails the feature.

  The condition and its parameter are set first, then Group by, Aggregate and Pivot, then Save. The
  condition reaches the builder's query 750 ms after its last keystroke (a debounced input), and the
  parameter checkbox with the builder's next run; the builder's row count after the pivot (1, the
  Project row) is read before Save, since the gestures' own duration is not bounded against the
  750 ms. The saved query is read back from the server holding the condition and its parameter.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after the save lets
  that one message through; the checks after the reopen are strict. The query and the project are
  named with the run's time and removed when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And user clears the saved pivot table parameters

  @realizes:GROK-20273
  Scenario: A project of a visual query with a pivot and a parameter reopens without an error
    Given no query named "BDDRegVQPivot{time}" is on the server
    And no project named "BDDRegVQPivotProj{time}" is on the server
    And the layout saved for the query "BDDRegVQPivot{time}" is deleted at the end
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
    When user adds "name" to the "Group by" row of the visual query
    And user adds "id" to the "Aggregate" row of the visual query
    And user adds "is_package_entity" to the "Pivot" row of the visual query
    Then the "Group by" row of the visual query should hold "name"
    And the "Pivot" row of the visual query should hold "is_package_entity"
    And the visual query should have run to 1 row
    When user enters "BDDRegVQPivot{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDDRegVQPivot{time}" should be on the server
    And the query "BDDRegVQPivot{time}" on the server should filter "name" by "Project" as a parameter
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then "BDDRegVQPivot{time}" dialog should be visible
    And Name text input in "BDDRegVQPivot{time}" dialog should have value "Project"
    When user enters "Script" into Name text input in "BDDRegVQPivot{time}" dialog
    And user clicks on OK button in "BDDRegVQPivot{time}" dialog
    Then the "BDDRegVQPivot{time}" dialog should close
    And the current view should be a TableView view
    And the only row of the table should read "Script" in the "name" column
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegVQPivotProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegVQPivotProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegVQPivotProj{time}" dialog
    Then the "Share BDDRegVQPivotProj{time}" dialog should close
    And the "BDDRegVQPivot{time}" table of the "BDDRegVQPivotProj{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegVQPivotProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegVQPivotProj{time} gallery card
    Then the task bar should have finished "Opening project"
    And the current view should be a TableView view
    And the table should have been reloaded by data sync
    And the table should have 1 row
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @known-failure @realizes:GROK-20535
  Scenario: The reopened pivoted visual query holds the row of the value it was run with
    Then the only row of the table should read "Script" in the "name" column
