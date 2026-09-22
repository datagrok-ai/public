@connections @full-stand @serial @journey
Feature: The catalogs of an MS SQL connection
  An MS SQL connection lists its catalogs under a Catalogs group, a catalog its schemas and tables;
  a catalog opens on the context panel with a Database meta pane whose comment is saved with SAVE.
  Translated from TestTrack Connections/catalogs.md and catalogs-ui.md (playwright-public
  connections/10-catalogs.test.ts, skipped there: "no MS SQL connection for the playwright user" —
  the DBTests MS SQL NorthwindTest connection is visible on dev now). MS SQL is on a full stand only.

  The comment lands on a shared connection's catalog, so the feature is @serial and writes to tempdb,
  the catalog nothing else describes; the comment is claimed on the server (what SAVE wrote, not
  what the pane still shows), cleared by the last scenario and again at feature end. The Catalogs
  row carries the connection's own tree name — until the core names it, it is reached through its
  wrapper ("MS-SQL-MSSQLTest-Catalogs" tree group); switch to
  Databases---MS-SQL---NorthwindTest---Catalogs tree node once the core names the row.

  Not translated, and why: the md's icon claim (two databases for catalogs, tables for schemas) —
  the icons carry no name or label to read. The old spec's comment check typed and read the inputs
  back without SAVE, which proved nothing. Selecting another node and coming back, as the md does,
  cannot be the proof here: a click in the tree within 2 s of a property edit leaves the context
  panel on the edited object (the platform's propertyEdited guard), so the server is read instead.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---MS-SQL tree node inside browse tree is expanded
    And Databases---MS-SQL---NorthwindTest tree node inside browse tree is expanded

  Scenario: The Catalogs group lists the connection's databases, and a catalog its tables
    Given "MS-SQL-MSSQLTest-Catalogs" tree group inside browse tree is expanded
    Then the following elements should be visible:
      | Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree |
      | Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree    |
      | Databases---MS-SQL---NorthwindTest---Catalogs---master tree node inside browse tree    |
    Given Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree is expanded
    And Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo tree node inside browse tree is expanded
    Then Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo---orders tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: A catalog's menu opens it as a table
    When user right-clicks on Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree
    Then the open menu should list "Browse"
    And the open menu should list "Open as table"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: A comment on a catalog is saved on the server
    Given the "tempdb" catalog of the "MSSQLTest" connection has no comment
    When user clicks on Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree
    Then the context panel should show "tempdb"
    Given "Database meta" section in context panel is expanded
    When user enters "BDD comment {run}" into Comment input in context panel
    And user clicks on SAVE button in context panel
    Then the "tempdb" catalog of the "MSSQLTest" connection should have the comment "BDD comment {run}"
    And no error or warning balloon should have been shown

  # Candidate finding, ticket pending Olesia's manual walk: clearing the Comment input and pressing
  # Save leaves the old comment on the server (probed on dev 2026-09-22; typing a new comment does
  # save). The feature-end cleanup clears it through the API either way.
  @known-failure
  Scenario: Clearing the comment in the pane clears it on the server
    When user clears Comment input in context panel
    And user clicks on SAVE button in context panel
    Then the "tempdb" catalog of the "MSSQLTest" connection should have the comment ""
    And no error or warning balloon should have been shown
    And no errors should have been logged
