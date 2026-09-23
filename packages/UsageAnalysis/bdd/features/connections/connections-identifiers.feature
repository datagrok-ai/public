@connections @needs-credentials @journey
Feature: Identifiers configured on a connection
  Configure Identifiers... asks for the primary schema and opens the identifiers view; an identifier
  names a semantic type for a table's column matched by a pattern, and after SAVE the column of that
  table comes with the semantic type; removing the configuration takes it away. Translated from
  TestTrack Connections/identifiers.md and identifiers-ui.md (playwright-public
  connections/02-identifiers.test.ts).

  The connection is the feature's own — BDD-Conn-Ident-{time}, db.datagrok.ai:54322/northwind,
  saved through the API and given the Northwind login and password in its Edit dialog from
  DG_PG_LOGIN / DG_PG_PASSWORD (typed from the environment, never printed) — so the configuration
  lives and dies with it; it is deleted at feature end and checked gone. {time} rather than {run}:
  the connection's name drops the dashes of its friendly name.

  Not translated, and why: "the values are highlighted in blue" — the semantic type the column
  carries is what the highlight is drawn from, and it is claimed directly. The md's reload is kept:
  a table opened in the same session as the SAVE carries no semantic type (checked on dev
  2026-09-22), so the identifiers reach a session that started after them.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    And Databases tree node inside browse tree is expanded
    And a "Postgres" connection named "BDD-Conn-Ident-{time}" is on the server
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree
    And user picks "Edit..." from the open menu
    And user enters the DG_PG_LOGIN secret into Login input in "Edit Connection" dialog
    And user enters the DG_PG_PASSWORD secret into Password input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close

  Scenario: The primary schema opens the identifiers view
    When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree
    And user picks "Configure Identifiers..." from the open menu
    Then "Select primary schema for Identifiers Configuration" dialog should be visible
    And Schema input in "Select primary schema for Identifiers Configuration" dialog should offer "public, pg_catalog, information_schema"
    When user selects "public" in Schema input in "Select primary schema for Identifiers Configuration" dialog
    And user clicks on OK button in "Select primary schema for Identifiers Configuration" dialog
    Then the "Select primary schema for Identifiers Configuration" dialog should close
    And "Add a new identifier" icon should be visible
    And no errors should have been logged

  Scenario: An identifier is added and saved
    When user clicks on "Add a new identifier" icon
    Then "Add Identifier" dialog should be visible
    When user enters "CUSTOMER_ID" into "Semantic Type" input in "Add Identifier" dialog
    And user selects "customers" in Table input in "Add Identifier" dialog
    And user selects "customerid" in Column input in "Add Identifier" dialog
    And user enters "[A-Z]{5}" into "Match Regexp" input in "Add Identifier" dialog
    And user clicks on Add button in "Add Identifier" dialog
    Then the "Add Identifier" dialog should close
    When user clicks on Save button
    Then an info balloon should have been shown
    And no error or warning balloon should have been shown

  Scenario: The table's column comes with the identifier's semantic type
    When user reloads the page
    Given the browse panel is open
    And Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree is expanded
    And Databases---Postgres---BDD-Conn-Ident-{time}---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public tree node inside browse tree is expanded
    When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public---customers tree node inside browse tree
    And user picks "Get All" from the open menu
    Then the "customers" view should be current
    And the table should have 91 rows
    And "customerid" column should have semantic type "CUSTOMER_ID"
    And no error or warning balloon should have been shown
