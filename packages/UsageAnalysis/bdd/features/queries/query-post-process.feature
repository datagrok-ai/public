@journey @serial @realizes:views.queries
Feature: A query's post-process runs on its result
  The Post-Process tab of a query holds a script that runs on every result of the query: a line
  typed there and saved announces the row count whenever the saved query runs. Translated from the
  TestTrack Queries case query-postprocessing (playwright-public/queries query-postprocessing), on
  the feature's own query over NorthwindTest's products (77 rows).

  The query is named with the run's time and removed when the feature ends; @serial for the same
  reason as the other features that save into NorthwindTest.

  Not translated, and why: the case's layout half — two viewers on the Layout tab kept with the
  query — is query-layout.feature, since saving a layout fails for a non-admin today.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-Q-pp-{time}" is on the server

  Scenario: A line is typed into the Post-Process tab of a new query
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDD-Q-pp-{time}" into Name input
    And user replaces the code of code editor with "select * from products"
    And user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 77
    When user moves the pointer away from play icon
    And user clicks on Post-Process tab
    And user puts "grok.shell.info('PP' + result.rowCount);" on the first line of code editor
    Then code editor should contain the text "grok.shell.info('PP' + result.rowCount);"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # Candidate finding, ticket pending Olesia's manual walk: the Post-Process editor syncs its text
  # to the query 250 ms after the last key, and Save sends the synced copy
  # (data_query_view.dart:231), so a Save right after typing stores the old template. The core
  # fix (Save reads the editor) is in the separate core PR — remove the tag once it is deployed.
  @known-failure
  Scenario: Save right after typing keeps the line, and the saved query runs it
    When user clicks on Save button
    Then 1 query named "BDD-Q-pp-{time}" should be on the server
    And the query "BDD-Q-pp-{time}" on the server should have a post-process containing "grok.shell.info('PP' + result.rowCount);"
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    And user hovers over Databases---Postgres---NorthwindTest---BDD-Q-pp-{time} tree node inside browse tree
    And user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-pp-{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 77 rows
    And an info balloon containing "PP77" should have been shown
    And no errors should have been logged
