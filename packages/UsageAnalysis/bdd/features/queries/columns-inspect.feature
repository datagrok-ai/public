@realizes:views.queries
Feature: Every column of a schema in the context panel
  Clicking a column of a table in a connection's schema makes it the current object: the context
  panel shows it with its General, Actions, Inspect and Database meta panes. Translated from the
  TestTrack Queries case columns-inspect (playwright-public/queries columns-inspect) — every table
  of NorthwindTest's public schema, every column. Read-only.

  Not translated, and why: the PostgresDart run of the same walk — the panel is the same for every
  provider, and what differs is the provider's server side. The case's second part names
  "Postgres > Northwind", a connection only the public stand has; it runs on Postgres >
  NorthwindTest, the same database.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open

  Scenario: Every column of every table of the public schema is shown on click
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the tree says nothing when a connection has finished listing its children, and a Schemas row
    # opened before that comes back empty: the connection's Orders query is the sign it has
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks every column of every table under Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree
    Then every clicked column should have been shown in the context panel with "General, Actions, Inspect, Database meta"
    And no errors should have been logged
    And no error or warning balloon should have been shown
