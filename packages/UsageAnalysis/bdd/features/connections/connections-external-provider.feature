@connections @full-stand @serial @needs-credentials @journey
Feature: Writing to an external Postgres through queries
  Queries on an external Postgres run DDL and DML and read back what they wrote. Translated from
  TestTrack Connections/external-provider.md (playwright-public connections/09-external-provider,
  whose CRUD part was hard-skipped: "no writable external Postgres on the CI stack"). On dev the
  It needs a login that may write in that database, which the suite does not have today: the stored
  credentials of the Dbtests PostgreSQLDBTests connection answer "permission denied for schema
  public", and the DG_PG_EXT_LOGIN / DG_PG_EXT_PASSWORD of the dev harness are refused by
  db.datagrok.ai:54327 ("password authentication failed", checked 2026-09-22) — the same wall the CI
  spec recorded when it skipped this case. The feature is therefore @needs-credentials: give the run
  a login with DDL rights on that database and it exercises the four statements.

  The table is bdd_tmp_{time}, its own per run, dropped before and after the feature; the queries are
  saved as BDD-Conn-Ext-…-{time} and deleted at feature end — both checked gone. The feature runs
  serially because every feature that saves queries lists the same connections; a full stand only.
  The old spec checked only that no error balloon followed each statement, which proves nothing: a
  query's error goes to the editor's Messages panel, not to a balloon. Here a SELECT reads back what
  the INSERT and the UPDATE wrote, which fails when any statement before it did nothing.

  Creating the connection itself from the dialog (the md's first step) needs the test database's
  login and password: it is connections-credentials.feature (@needs-credentials), so a run without
  secrets can leave it out.

  Not translated, and why: nothing of the md beyond that split.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And no table "bdd_tmp_{time}" is in the database of the "PostgreSQLDBTests" connection
    And no query named "BDD-Conn-Ext-Create-{time}" is on the server
    And no query named "BDD-Conn-Ext-Insert-{time}" is on the server
    And no query named "BDD-Conn-Ext-Update-{time}" is on the server
    And no query named "BDD-Conn-Ext-Select-{time}" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded

  Scenario: A CREATE TABLE query is saved and run
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree
    And user picks "New Query..." from the open menu
    Then the current view should be a DataQueryView view
    When user enters "BDD-Conn-Ext-Create-{time}" into Name input
    And user replaces the code of code editor with "create table bdd_tmp_{time} (id int, name varchar(50))"
    And user clicks on Save button
    Then 1 query named "BDD-Conn-Ext-Create-{time}" should be on the server
    When user clicks on play icon
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An INSERT query writes a row the SELECT query reads back
    When user closes the current view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree
    And user picks "New Query..." from the open menu
    And user enters "BDD-Conn-Ext-Insert-{time}" into Name input
    And user replaces the code of code editor with "insert into bdd_tmp_{time} (id, name) values (1, 'test')"
    And user clicks on Save button
    Then 1 query named "BDD-Conn-Ext-Insert-{time}" should be on the server
    When user clicks on play icon
    Then no error or warning balloon should have been shown
    When user closes the current view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree
    And user picks "New Query..." from the open menu
    And user enters "BDD-Conn-Ext-Select-{time}" into Name input
    And user replaces the code of code editor with "select id, name from bdd_tmp_{time}"
    And user clicks on Save button
    And user clicks on play icon
    Then the "rows" reading of grid should be 1
    And the "text of cell 1 of name" reading of grid should be "test"
    And no errors should have been logged

  Scenario: An UPDATE query changes the row
    When user closes the current view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree
    And user picks "New Query..." from the open menu
    And user enters "BDD-Conn-Ext-Update-{time}" into Name input
    And user replaces the code of code editor with "update bdd_tmp_{time} set name = 'bdd' where id = 1"
    And user clicks on Save button
    Then 1 query named "BDD-Conn-Ext-Update-{time}" should be on the server
    When user clicks on play icon
    Then no error or warning balloon should have been shown
    When user closes the current view
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree
    And user picks "New Query..." from the open menu
    And user replaces the code of code editor with "select id, name from bdd_tmp_{time}"
    And user clicks on play icon
    Then the "rows" reading of grid should be 1
    And the "text of cell 1 of name" reading of grid should be "bdd"
    And no errors should have been logged
