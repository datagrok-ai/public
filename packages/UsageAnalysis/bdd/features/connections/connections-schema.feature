@connections @journey
Feature: The schemas of a connection and the schema view
  A connection lists its schemas under a Schemas group; a schema's menu offers Browse, and Browse
  opens the schema view — one box per table, and a box's menu is the table's own menu. Translated
  from TestTrack Connections/schema.md (playwright-public connections/07-schema.test.ts) and the
  Browse-DB-03 case that browse-platform-and-databases.feature left for this file. The connection is
  CHEMBL (Postgres), which the Browse features already rely on; nothing is created on the server.

  The Schemas row is named by its path and a schema view's table box by its table, both since the
  core change that named them.

  Not translated, and why: the md's "Browse on the connection opens the schema" — Browse on a
  connection opens its queries gallery now (browse-context-panel-and-menus.feature claims the
  connection's menu); the schema view opens from Schemas > <schema> > Browse, which is what this
  feature does.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---CHEMBL tree node inside browse tree is expanded

  Scenario: The Schemas group lists the connection's schemas
    Given Databases---Postgres---CHEMBL---Schemas tree node inside browse tree is expanded
    Then the following elements should be visible:
      | Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree             |
      | Databases---Postgres---CHEMBL---Schemas---information-schema tree node inside browse tree |
    And no errors should have been logged

  Scenario: A schema's menu offers the schema view and the table actions
    When user opens the context menu of Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree
    Then the open menu should list "Browse"
    And the open menu should list "Open as table"
    And the open menu should list "New Table..."
    And the open menu should list "Import Table..."
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Browse opens the schema view with a box per table
    When user picks "Browse" from the context menu of Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree
    Then the "Schema: public" view should be current
    And "activities" schema table should be visible
    And "molecule_dictionary" schema table should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A table box's menu is the table's menu, and Get Top 100 reads the table
    When user opens the context menu of "activities" schema table
    Then the open menu should list "Get All"
    And the open menu should list "Get Top 100"
    And the open menu should list "New SQL Query..."
    And the open menu should list "New Visual Query..."
    When user picks "Get Top 100" from the open menu
    Then the "activities" view should be current
    And the table should have 100 rows
    And the table should have a column "activity_id"
    And no error or warning balloon should have been shown
