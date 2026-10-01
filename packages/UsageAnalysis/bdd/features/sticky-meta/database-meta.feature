@journey @serial @sticky-meta
Feature: Database meta of a database schema, of a table and of its column
  The Database meta pane of the context panel for the public schema of NorthwindTest, for its
  categories table and for that table's categoryid column: empty to begin with, filled and saved,
  found again after the page is reloaded, cleared and saved, and gone after another reload — and never
  shown on a table or a column it was not written for, a same-named column of another table included.
  Translated from TestTrack StickyMeta/04-database-meta.md (Tests 4.1, 4.2), its primary
  database-meta case (the schema-level part with its "test@#$!" value), and
  playwright-tests/e2e/stickymeta/04-database-meta.test.ts.

  GROK-19427 (string_list values of a column could not be saved) and GROK-19429 (Row Count could not
  be deleted) are fixed; the column's Values and Sample Values and the table's Row Count are written
  and cleared like the other fields. Every field is claimed empty before it is written, so a value
  left by a killed run cannot stand in for this run's save. The primary case's database-level meta of
  CHEMBL is not translated: the platform builds a Database meta pane for a connection only when its
  data source has no catalogs (df_properties.dart, the DataConnection branch), and Postgres has.

  A column of a database table is named in the panel as the database has it ("categoryid") while
  its entity is "Categoryid", which `the context panel should show` compares; the column scenarios
  name it by the panel's text and claim a field only a column's pane has first (MISSING.md).

  Clicking the schema node opens the schema's view, and a reload with that view in front logs
  NullError in DbSchemaView.saveStateMap (db_views.dart:192 — the view restored from its address has
  no connection yet); the schema scenario closes the views before each reload, so the claims are about
  the Database meta, and the error is reported on its own.

  NorthwindTest is the PostgresTest connection the DBTests package brings; a stand without it (the CI
  stack) skips the feature at its first step, with the reason.

  The values are the run's own and the last scenarios put every field back empty through the pane;
  the server-side restore the library does not have yet is in sticky-meta/MISSING.md. Serial: the
  schema, the table and the column are shared by every feature that reads NorthwindTest.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open

  Scenario: A schema's Database meta is saved, kept after a reload and cleared
    Given the stand has a reachable "PostgresTest" connection
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree
    Then context panel should contain text "public"
    Given "Database meta" pane in context panel is expanded
    Then Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    When user enters "test@#$! {time}" into Comment input in "Database meta" pane in context panel
    And user enters "test@#$! llm {time}" into "LLM Comment" input in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be enabled
    When user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    When user closes all views
    And user reloads the page
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree
    Then context panel should contain text "public"
    Given "Database meta" pane in context panel is expanded
    Then Comment input in "Database meta" pane in context panel should have value "test@#$! {time}"
    And "LLM Comment" input in "Database meta" pane in context panel should have value "test@#$! llm {time}"
    When user clears Comment input in "Database meta" pane in context panel
    And user clears "LLM Comment" input in "Database meta" pane in context panel
    And user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    When user closes all views
    And user reloads the page
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree
    Then context panel should contain text "public"
    Given "Database meta" pane in context panel is expanded
    Then Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A table's Database meta is saved (4.1)
    Given Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree
    Then the context panel should show "categories"
    Given "Database meta" pane in context panel is expanded
    Then the following elements should be visible:
      | Domains input in "Database meta" pane in context panel        |
      | "Row Count" input in "Database meta" pane in context panel    |
      | Comment input in "Database meta" pane in context panel        |
      | "LLM Comment" input in "Database meta" pane in context panel  |
      | Save button in "Database meta" pane in context panel          |
    And Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    And "Row Count" input in "Database meta" pane in context panel should have value ""
    When user enters "bdd-sm-db-{time} table" into Comment input in "Database meta" pane in context panel
    And user enters "bdd-sm-db-{time} table llm" into "LLM Comment" input in "Database meta" pane in context panel
    And user enters "8" into "Row Count" input in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be enabled
    When user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table's Database meta is there after a reload, and on no other table (4.1)
    When user reloads the page
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree
    Then the context panel should show "categories"
    Given "Database meta" pane in context panel is expanded
    Then Comment input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} table"
    And "LLM Comment" input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} table llm"
    And "Row Count" input in "Database meta" pane in context panel should have value "8"
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree
    Then the context panel should show "customers"
    And Comment input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} table"
    And "LLM Comment" input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} table llm"
    And "Row Count" input in "Database meta" pane in context panel should not have value "8"
    And no errors should have been logged

  Scenario: A column's Database meta is saved (4.2)
    Given Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree
    Then context panel should contain text "categoryid"
    Given "Database meta" pane in context panel is expanded
    Then the following elements should be visible:
      | "Is Unique" input in "Database meta" pane in context panel     |
      | Min input in "Database meta" pane in context panel             |
      | Max input in "Database meta" pane in context panel             |
      | Values input in "Database meta" pane in context panel          |
      | "Sample Values" input in "Database meta" pane in context panel |
      | "Unique Count" input in "Database meta" pane in context panel  |
      | Quality input in "Database meta" pane in context panel         |
      | Comment input in "Database meta" pane in context panel         |
      | "LLM Comment" input in "Database meta" pane in context panel   |
    And "Is Unique" input in "Database meta" pane in context panel should not be checked
    And Min input in "Database meta" pane in context panel should have value ""
    And Max input in "Database meta" pane in context panel should have value ""
    And Values input in "Database meta" pane in context panel should have value ""
    And "Sample Values" input in "Database meta" pane in context panel should have value ""
    And "Unique Count" input in "Database meta" pane in context panel should have value ""
    And Quality input in "Database meta" pane in context panel should have value ""
    And Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    When user checks "Is Unique" input in "Database meta" pane in context panel
    And user enters "1" into Min input in "Database meta" pane in context panel
    And user enters "8" into Max input in "Database meta" pane in context panel
    And user enters "1" into Values input in "Database meta" pane in context panel
    And user enters "2" into "Sample Values" input in "Database meta" pane in context panel
    And user enters "8" into "Unique Count" input in "Database meta" pane in context panel
    And user enters "good {time}" into Quality input in "Database meta" pane in context panel
    And user enters "bdd-sm-db-{time} column" into Comment input in "Database meta" pane in context panel
    And user enters "bdd-sm-db-{time} column llm" into "LLM Comment" input in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be enabled
    When user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The column's Database meta is there after a reload, and on no other column (4.2)
    When user reloads the page
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree
    Then context panel should contain text "categoryid"
    Given "Database meta" pane in context panel is expanded
    Then "Is Unique" input in "Database meta" pane in context panel should be checked
    And Min input in "Database meta" pane in context panel should have value "1"
    And Max input in "Database meta" pane in context panel should have value "8"
    And "Unique Count" input in "Database meta" pane in context panel should have value "8"
    And Quality input in "Database meta" pane in context panel should have value "good {time}"
    And Comment input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} column"
    And "LLM Comment" input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} column llm"
    And Values input in "Database meta" pane in context panel should have value "1"
    And "Sample Values" input in "Database meta" pane in context panel should have value "2"
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryname tree node inside browse tree
    Then context panel should contain text "categoryname"
    And Quality input in "Database meta" pane in context panel should not have value "good {time}"
    And Comment input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} column"
    Given Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---products---categoryid tree node inside browse tree
    Then context panel should contain text "categoryid"
    And Quality input in "Database meta" pane in context panel should not have value "good {time}"
    And Comment input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} column"
    And no errors should have been logged

  Scenario: Cleared and saved, the column's and the table's Database meta are gone after a reload (4.1, 4.2)
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree
    Then context panel should contain text "categoryid"
    And "Is Unique" input in "Database meta" pane in context panel should be checked
    And Comment input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} column"
    When user unchecks "Is Unique" input in "Database meta" pane in context panel
    And user clears Min input in "Database meta" pane in context panel
    And user clears Max input in "Database meta" pane in context panel
    And user clears Values input in "Database meta" pane in context panel
    And user clears "Sample Values" input in "Database meta" pane in context panel
    And user clears "Unique Count" input in "Database meta" pane in context panel
    And user clears Quality input in "Database meta" pane in context panel
    And user clears Comment input in "Database meta" pane in context panel
    And user clears "LLM Comment" input in "Database meta" pane in context panel
    And user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree
    Then the context panel should show "categories"
    And Comment input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} table"
    When user clears Comment input in "Database meta" pane in context panel
    And user clears "LLM Comment" input in "Database meta" pane in context panel
    And user clears "Row Count" input in "Database meta" pane in context panel
    And user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    When user reloads the page
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree
    Then the context panel should show "categories"
    Given "Database meta" pane in context panel is expanded
    Then Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    And "Row Count" input in "Database meta" pane in context panel should have value ""
    Given Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree
    Then context panel should contain text "categoryid"
    And "Is Unique" input in "Database meta" pane in context panel should not be checked
    And Min input in "Database meta" pane in context panel should have value ""
    And Max input in "Database meta" pane in context panel should have value ""
    And Values input in "Database meta" pane in context panel should have value ""
    And "Sample Values" input in "Database meta" pane in context panel should have value ""
    And "Unique Count" input in "Database meta" pane in context panel should have value ""
    And Quality input in "Database meta" pane in context panel should have value ""
    And Comment input in "Database meta" pane in context panel should have value ""
    And "LLM Comment" input in "Database meta" pane in context panel should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown
