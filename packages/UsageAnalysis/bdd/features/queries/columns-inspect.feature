@realizes:views.queries
Feature: Every column of a schema in the context panel
  Clicking a column of a table in a connection's schema makes it the current object: the context
  panel shows it with its General, Actions, Inspect and Database meta panes. Translated from the
  TestTrack Queries case columns-inspect (playwright-public/queries columns-inspect) — every table
  of NorthwindTest's public schema, every column, on PostgresDart and on Postgres.

  Read-only. PostgresDart is @full-stand — a minimal stack carries no PostgresDart connection.

  Not translated, and why: the case's second part names "Postgres > Northwind", a connection only
  the public stand has; it runs on Postgres > NorthwindTest, the same database.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open

  Scenario Outline: Every column of every table of <provider>'s public schema is shown on click
    Given Databases tree node inside browse tree is expanded
    When user collapses Databases---<other> tree node inside browse tree
    Given Databases---<provider> tree node inside browse tree is expanded
    And Databases---<provider>---NorthwindTest tree node inside browse tree is expanded
    # the tree says nothing when a connection has finished listing its children, and a Schemas row
    # opened before that comes back empty: the connection's Orders query is the sign it has
    Then Databases---<provider>---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---<provider>---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---<provider>---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks every column of every table under Databases---<provider>---NorthwindTest---Schemas---public tree node inside browse tree
    Then every clicked column should have been shown in the context panel with "General, Actions, Inspect, Database meta"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    @full-stand
    Examples:
      | provider     | other    |
      | PostgresDart | Postgres |

    Examples:
      | provider | other        |
      | Postgres | PostgresDart |
