@realizes:views.queries
Feature: Get All and Get Top 100 on a table of the schema
  The two read commands of a table node: Get All opens the whole table, Get Top 100 its first
  hundred rows, each in a view of its own. Translated from the TestTrack Queries case
  get-all-get-top-100 (playwright-public/queries get-all-get-top-100), on NorthwindTest, whose
  orders table holds 830 rows. Read-only: nothing is saved, so nothing is cleaned up.

  Not translated, and why: the PostgresDart run — the commands are the same for every provider,
  and what differs is the provider's server side. Where the md walks Browse > Platform > Functions
  to reach a view, the feature takes the route directly (the walk itself is claimed in
  scripts-create).

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: Get All opens the whole table, Get Top 100 its first hundred rows
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the tree says nothing when a connection has finished listing its children, and a Schemas row
    # opened before that comes back empty: the connection's Orders query is the sign it has
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    Given user watches the task bar
    When user picks "Get All" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the task bar should have finished "orders"
    And table "orders" should have 830 rows
    Given the toolbox pane is hidden
    And the browse panel is open
    Given user watches the task bar
    When user picks "Get Top 100" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the task bar should have finished "orders"
    And table "orders (2)" should have 100 rows
    And table "orders" should have 830 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
