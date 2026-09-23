@serial @realizes:views.queries
Feature: A SQL query started from a table of the schema
  New SQL Query... on a table of a connection's schema opens the query editor already filled in:
  the table's name for the query and a select of the whole table for the SQL. Translated from the
  TestTrack Queries case new-sql-query (playwright-public/queries new-sql-query).

  The prefilled name is the table's, "products", which every run and every person would share; the
  feature claims it, then saves under its own run-named name and removes that query at its end.
  It is @serial, since it saves into the NorthwindTest connection other query features list.

  Not translated, and why: nothing of the case is left out; the save under another name is the
  only departure, for the reason above.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-Q-sql-{time}" is on the server

  Scenario: The editor opens on the table's select, runs it and saves it
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the tree says nothing when a connection has finished listing its children, and a Schemas row
    # opened before that comes back empty: the connection's Orders query is the sign it has
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree
    Then the current view should be a DataQueryView view
    And Name input should have value "products"
    And code editor should hold the code "select * from public.products"
    When user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 77
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then the current view should be a TableView view
    And the table should have 77 rows
    When user closes the current view
    Then the current view should be a DataQueryView view
    When user enters "BDD-Q-sql-{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDD-Q-sql-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown
