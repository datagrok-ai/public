@journey @serial @realizes:views.queries
Feature: A query's layout
  The Layout tab of a query edits the look of its result: viewers added there from the toolbox are
  saved with the query and come back when the saved query runs. Translated from the TestTrack
  Queries case query-layout and the layout half of query-postprocessing (playwright-public/queries
  query-layout), on the feature's own query over NorthwindTest's products (77 rows) — the case
  edits the shared PostgresAll query, which other suites read.

  The query is named with the run's time and removed when the feature ends; @serial for the same
  reason as the other features that save into NorthwindTest.

  The result of the saved query, with one more viewer, is then saved as a project through the
  ribbon's Save dialog, opened again with all three viewers, and the Source pane's REFRESH (the md's
  "Toolbox > File > Refresh") runs the query again without losing them (query-layout steps 9-13).
  The project is named with the run's time and removed when the feature ends.

  Not translated, and why: docking one viewer onto another inside the Layout tab — the dock zones
  carry no names.

  Background:
    Given user is logged in
    And the stand has a reachable "PostgresTest" connection
    And the browse panel is open
    And no query named "BDD-Q-layout-{time}" is on the server
    And no project named "BDD-Q-layout-proj-{time}" is on the server
    And the layouts named "BDD-Q-layout-{time}" are deleted when the feature ends

  Scenario: The Layout tab waits for a run, then takes viewers from the toolbox
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDD-Q-layout-{time}" into Name input
    And user replaces the code of code editor with "select * from products"
    And user clicks on Layout tab
    Then "Run query to get data and edit layout" text should be visible
    When user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 77
    Given the toolbox pane is shown
    When user moves the pointer away from play icon
    And user clicks on scatter plot icon in toolbox
    And user clicks on correlation plot icon in toolbox
    Then scatter plot viewer should be visible
    And correlation plot viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Save keeps the query with its layout, and the saved query runs with its viewers
    When user clicks on Save button
    Then 1 query named "BDD-Q-layout-{time}" should be on the server
    And the query "BDD-Q-layout-{time}" on the server should have a layout
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user refreshes the browse tree
    And user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 77 rows
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 correlation plot viewer

  Scenario: The result with one more viewer is saved as a project and opens with all three
    Given the toolbox pane is shown
    When user clicks on histogram icon in toolbox
    Then the open tableview should have 1 histogram viewer
    When user clicks on Save button
    Then "Save project" dialog should be visible
    When user enters "BDD-Q-layout-proj-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDD-Q-layout-proj-{time}" should be on the server
    And "Share BDD-Q-layout-proj-{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Q-layout-proj-{time}" dialog
    Then the "Share BDD-Q-layout-proj-{time}" dialog should close
    When user closes all views
    And user opens the "BDD-Q-layout-proj-{time}" project and waits for its table
    Then the table should have 77 rows
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 correlation plot viewer
    And the open tableview should have 1 histogram viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: REFRESH of the reopened result keeps its viewers
    Given the toolbox pane is shown
    When user clicks on REFRESH button in toolbox and the query it runs completes
    Then the table should have 77 rows
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 correlation plot viewer
    And the open tableview should have 1 histogram viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown
