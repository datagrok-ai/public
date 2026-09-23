@realizes:views.queries
Feature: Get All and Get Top 100 on a table of the schema
  The two read commands of a table node: Get All opens the whole table, Get Top 100 its first
  hundred rows, each in a view of its own. Translated from the TestTrack Queries case
  get-all-get-top-100 (playwright-public/queries get-all-get-top-100), for the Postgres and the
  PostgresDart providers on NorthwindTest, whose orders table holds 830 rows.

  Read-only: nothing is saved, so nothing is cleaned up. The PostgresDart row is @full-stand — a
  minimal stack carries no PostgresDart connection.

  Not translated, and why: nothing of the case is left out.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario Outline: <provider> opens the whole table and its first hundred rows
    Given Databases tree node inside browse tree is expanded
    # the tree keeps what other features opened: the other provider is folded so this one's
    # schema stays within the panel
    When user collapses Databases---<other> tree node inside browse tree
    Given Databases---<provider> tree node inside browse tree is expanded
    And Databases---<provider>---NorthwindTest tree node inside browse tree is expanded
    # the tree says nothing when a connection has finished listing its children, and a Schemas row
    # opened before that comes back empty: the connection's Orders query is the sign it has
    Then Databases---<provider>---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---<provider>---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---<provider>---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    Given user watches the task bar
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    When user hovers over Databases---<provider>---NorthwindTest---Schemas---public---orders tree node inside browse tree
    And user picks "Get All" from the context menu of Databases---<provider>---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the task bar should have finished "orders"
    And table "orders" should have 830 rows
    Given the toolbox pane is hidden
    And the browse panel is open
    Given user watches the task bar
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    When user hovers over Databases---<provider>---NorthwindTest---Schemas---public---orders tree node inside browse tree
    And user picks "Get Top 100" from the context menu of Databases---<provider>---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the task bar should have finished "orders"
    And table "orders (2)" should have 100 rows
    And table "orders" should have 830 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | provider | other        |
      | Postgres | PostgresDart |

    @full-stand
    Examples:
      | provider     | other    |
      | PostgresDart | Postgres |
