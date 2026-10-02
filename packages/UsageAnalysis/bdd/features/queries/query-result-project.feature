@journey @serial @realizes:views.queries
Feature: A query result saved as a project
  The result of a parameterized query, with a viewer added from the toolbox, saved as a project
  through the ribbon's Save dialog and opened again. Translated from the project part of the
  TestTrack Queries case browse-and-save-project (playwright-public/queries
  chembl-parameterized-and-project), on NorthwindTest's PostgresByStringChoices instead of CHEMBL —
  the save and reopen are the claim, and they do not need CHEMBL's slow search.

  The project is named with the run's time and removed, with its table and view, when the feature
  ends. It is @serial: the save uploads the table and the reopen reads it back, and under parallel
  features the stand answered the open slower than the step's budget.

  Not translated, and why: the case's "press + to add the result to the workspace" — the plus icon
  exists for a table query only, and a query run from the tree opens its own view already. The
  Share dialog that follows the save is cancelled: sharing is not what the case checks.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDD-Q-proj-{time}" is on the server

  Scenario: A query result with a viewer is saved as a project
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree
    And user selects "France" in "Ship Country" input in "PostgresByStringChoices" dialog
    And user clicks on OK button in "PostgresByStringChoices" dialog
    Then the current view should be a TableView view
    And the table should have 77 rows
    Given the toolbox pane is shown
    When user clicks on trellis plot icon in toolbox
    Then the open tableview should have 1 trellis plot viewer
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user enters "BDD-Q-proj-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Q-proj-{time}" should be on the server
    And "Share BDD-Q-proj-{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Q-proj-{time}" dialog
    Then the "Share BDD-Q-proj-{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The saved project opens with its table and its viewer
    When user closes all views
    And user opens the "BDD-Q-proj-{time}" project and waits for its table
    Then the table should have 77 rows
    And the open tableview should have 1 trellis plot viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown
