@journey @serial @sticky-meta
Feature: Database meta of a database schema, of a table and of its column
  The Database meta pane of the context panel for the public schema of NorthwindTest, for its
  categories table and for that table's categoryid column: the fields each level offers, empty to begin
  with, filled and saved — the save claimed on the server — never shown on a table or a column it was
  not written for, a same-named column of another table included, and cleared and saved again.
  Translated from TestTrack StickyMeta/04-database-meta.md (Tests 4.1, 4.2), its primary
  database-meta case (the schema-level part with its "test@#$!" value), and
  playwright-tests/e2e/stickymeta/04-database-meta.test.ts.

  The md proves each save by reloading the page and reading the pane again; here the server is read
  back instead, and the round trip of every field through the store is the DBTests package's claim
  (src/db-annotations, Database Meta: DbSchemaInfo). GROK-19427 (string_list values of a column could
  not be saved) and GROK-19429 (Row Count could not be deleted) are fixed; the column's Values and
  Sample Values and the table's Row Count are written and cleared like the other fields. The primary
  case's database-level meta of CHEMBL is not translated: the platform builds a Database meta pane for
  a connection only when its data source has no catalogs (df_properties.dart, the DataConnection
  branch), and Postgres has.

  A column of a database table is named in the panel as the database has it ("categoryid") while
  its entity is "Categoryid", which `the context panel should show` compares; the column scenarios
  name it by the panel's text and claim a field only a column's pane has first.

  NorthwindTest is the PostgresTest connection the DBTests package brings; a stand without it (the CI
  stack) skips the feature at its first step, with the reason. Everything Datagrok keeps as Database
  meta for that connection is cleared through the API when the feature starts and when it ends, so a
  run killed after a save does not fail the next one. Serial: the schema, the table and the column are
  shared by every feature that reads NorthwindTest.

  Background:
    Given user is logged in
    And the stand has a reachable "PostgresTest" connection
    And the Database meta of the "PostgresTest" connection is cleared now and at feature end
    And the browse panel is open
    And the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded

  Scenario: A schema's Database meta is saved and cleared
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
    And the Database meta of "public" in the "PostgresTest" connection should be:
      | Comment     | test@#$! {time}     |
      | LLM Comment | test@#$! llm {time} |
    When user clears Comment input in "Database meta" pane in context panel
    And user clears "LLM Comment" input in "Database meta" pane in context panel
    And user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    And the Database meta of "public" in the "PostgresTest" connection should be:
      | Comment     |  |
      | LLM Comment |  |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A table's Database meta is saved, and on no other table (4.1)
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
    And the Database meta of "public.categories" in the "PostgresTest" connection should be:
      | Comment     | bdd-sm-db-{time} table     |
      | LLM Comment | bdd-sm-db-{time} table llm |
      | Row Count   | 8                          |
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree
    Then the context panel should show "customers"
    And Comment input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} table"
    And "LLM Comment" input in "Database meta" pane in context panel should not have value "bdd-sm-db-{time} table llm"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A column's Database meta is saved, and on no other column (4.2)
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
    And Comment input in "Database meta" pane in context panel should have value ""
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
    And the Database meta of "public.categories.categoryid" in the "PostgresTest" connection should be:
      | Is Unique     | true                        |
      | Min           | 1                           |
      | Max           | 8                           |
      | Values        | 1                           |
      | Sample Values | 2                           |
      | Unique Count  | 8                           |
      | Quality       | good {time}                 |
      | Comment       | bdd-sm-db-{time} column     |
      | LLM Comment   | bdd-sm-db-{time} column llm |
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

  Scenario: The column's and the table's Database meta are shown again, cleared and saved (4.1, 4.2)
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
    And the Database meta of "public.categories.categoryid" in the "PostgresTest" connection should be:
      | Min          |  |
      | Max          |  |
      | Values       |  |
      | Unique Count |  |
      | Quality      |  |
      | Comment      |  |
      | LLM Comment  |  |
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree
    Then the context panel should show "categories"
    And Comment input in "Database meta" pane in context panel should have value "bdd-sm-db-{time} table"
    When user clears Comment input in "Database meta" pane in context panel
    And user clears "LLM Comment" input in "Database meta" pane in context panel
    And user clears "Row Count" input in "Database meta" pane in context panel
    And user clicks on Save button in "Database meta" pane in context panel
    Then Save button in "Database meta" pane in context panel should be disabled
    And the Database meta of "public.categories" in the "PostgresTest" connection should be:
      | Comment     |  |
      | LLM Comment |  |
      | Row Count   |  |
    And no errors should have been logged
    And no error or warning balloon should have been shown
