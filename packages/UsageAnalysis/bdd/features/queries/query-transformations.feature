@journey @serial @realizes:views.queries
Feature: Transformations saved with a query
  A query's Transformations tab adds steps that run on its result every time it runs: a calculated
  column added there is in every result of the saved query, and is gone once the step is removed.
  Translated from the TestTrack Queries case transformations (playwright-public/queries
  postgres-query-transformations), on the feature's own query over NorthwindTest's products
  (77 rows) — the case edits the shared Products query, which takes a parameter and belongs to the
  Dbtests fixtures other suites read.

  The query is named with the run's time and removed when the feature ends; @serial for the same
  reason as the other features that save into NorthwindTest.

  Not translated, and why: the step list itself is read through its effect (the column in the
  result) — the steps carry no names a claim could read until the core names them. The CI spec's
  "hidden formula textarea" workaround is not needed: the formula is typed into the editor.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-Q-tr-{time}" is on the server

  Scenario: A column added in the Transformations tab is in the result of the saved query
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user hovers over Databases---Postgres---NorthwindTest tree node inside browse tree
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    And user picks "New Query..." from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDD-Q-tr-{time}" into Name input
    And user replaces the code of code editor with "select * from products"
    And user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 77
    When user moves the pointer away from play icon
    And user clicks on Transformations tab
    And user opens the "Add New Column" action of the transformations browser
    Then "Add New Column" dialog should be visible
    When user enters "doubled" into Name input in "Add New Column" dialog
    And user replaces the code of code editor in "Add New Column" dialog with "${productid} * 2"
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    When user clicks on Save button
    Then 1 query named "BDD-Q-tr-{time}" should be on the server
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then the current view should be a TableView view
    And the table should have 77 rows
    And the table should have a column "doubled"
    And the value of "doubled" column in row 1 should be "2"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Removing the step takes the column out of the saved query's result
    When user closes the current view
    Then the current view should be a DataQueryView view
    When user clicks on Transformations tab
    # the step list starts with the query's own run; the added column is the last step
    And user hovers over last "Remove step" icon
    And user clicks on last "Remove step" icon
    And user clicks on Save button
    And user clicks on "Run query..." action in toolbox
    Then the current view should be a TableView view
    And the table should have 77 rows
    And the table should not have a column "doubled"
    And no errors should have been logged
    And no error or warning balloon should have been shown
