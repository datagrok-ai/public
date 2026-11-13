@skip @journey @serial @realizes:powerpack.db-explorer
Feature: Identifiers of a database with catalogs
  TODO: development only. The feature is set aside and does not run in CI: its database exists on no
  shared stand (a Docker container made by catalogs-db.mjs, see the README), it saves an identifiers
  configuration in the one file every configuration of the stand is kept in, and it reloads the page.
  Take @skip off while working on catalogs in Configure Identifiers or the db-explorer library, run it
  on a local stand with the fixture, and put the tag back before committing.

  Configure Identifiers... on a connection whose database has catalogs (Databricks, Snowflake, SQL Server)
  asks for the catalog before the schema, reads the schemas of that catalog alone, and saves a
  configuration that names the catalog: the editor's previews and the cards of the identifiers read
  catalog.schema.table.

  The stand needs the "Catalogs Demo" connection, a SQL Server that holds the catalogs archive (the
  connection's own), research and sales, with two hundred more beside them; a stand without it skips the
  feature. research and archive both hold chem.compounds with the same ids and other names (archive's
  are "ARCHIVED ..."), and only research holds chem.batches and the schema bio: a card or a list read
  from the connection's own catalog shows it at once. catalogs-db.mjs in the feature's folder makes the
  database from catalogs-fixture.sql; its login owns those 203 catalogs and sees them, master and tempdb.

  The dialog opens on the connection's own catalog and lists the others behind it; the claims on the
  list wait for it. No schema is picked for the user among several: the first by name is a system one
  on SQL Server. A configuration of the connection's own catalog is saved without the catalog's name,
  as one made before catalogs were asked for, and reads schema.table; another catalog is named. The
  saved configuration opened again from the Browse tree and removed there are the same for a database
  without catalogs and are claimed in identifiers.feature. The configuration the feature saves is
  removed when it ends.

  Background:
    Given user is logged in
    And the stand has a reachable "CatalogsDemo" connection
    And no identifiers configuration of the "CatalogsDemo" connection is on the server
    And the browse panel is open

  Scenario: The dialog asks for the catalog first and opens on the connection's own
    Given Databases tree node inside browse tree is expanded
    And Databases---MS-SQL tree node inside browse tree is expanded
    When user picks "Configure Identifiers..." from the context menu of Databases---MS-SQL---Catalogs-Demo tree node inside browse tree
    Then "Select catalog and schema for Identifiers Configuration" dialog should be visible
    And Catalog input in "Select catalog and schema for Identifiers Configuration" dialog should have value "archive"
    And Schema input in "Select catalog and schema for Identifiers Configuration" dialog should offer the choice "chem"
    And Schema input in "Select catalog and schema for Identifiers Configuration" dialog should not offer the choice "bio"
    And "Select catalog and schema for Identifiers Configuration" dialog should contain text "Select the schema that holds the identifier tables."
    And OK button in "Select catalog and schema for Identifiers Configuration" dialog should be disabled

  Scenario: Among two hundred catalogs one is found by a part of its name
    Then "Select catalog and schema for Identifiers Configuration" dialog should contain text "205 catalogs; archive is the connection's own."
    When user clicks on catalog search icon
    And user finds the catalog "tenant_150" by typing "tenant_15"
    Then Catalog input in "Select catalog and schema for Identifiers Configuration" dialog should have value "tenant_150"
    And Schema input in "Select catalog and schema for Identifiers Configuration" dialog should have value "dbo"
    And OK button in "Select catalog and schema for Identifiers Configuration" dialog should be enabled

  Scenario: The editor of the connection's own catalog names it and says nothing of links
    When user selects "archive" in Catalog input in "Select catalog and schema for Identifiers Configuration" dialog
    Then Schema input in "Select catalog and schema for Identifiers Configuration" dialog should offer the choice "core"
    When user selects "chem" in Schema input in "Select catalog and schema for Identifiers Configuration" dialog
    And user clicks on OK button in "Select catalog and schema for Identifiers Configuration" dialog
    Then the "CatalogsDemo.archive.chem Identifiers" view should be current
    And identifiers editor should contain text "Catalogs Demo: archive.chem"
    And identifiers editor should not contain text "add them under Explicit References"

  Scenario: Another catalog brings its own schemas, and its editor says where links come from
    Given the browse panel is open
    When user picks "Configure Identifiers..." from the context menu of Databases---MS-SQL---Catalogs-Demo tree node inside browse tree
    And user selects "research" in Catalog input in "Select catalog and schema for Identifiers Configuration" dialog
    Then Schema input in "Select catalog and schema for Identifiers Configuration" dialog should offer the choice "bio"
    When user selects "chem" in Schema input in "Select catalog and schema for Identifiers Configuration" dialog
    And user clicks on OK button in "Select catalog and schema for Identifiers Configuration" dialog
    Then the "CatalogsDemo.research.chem Identifiers" view should be current
    And identifiers editor should contain text "Catalogs Demo: research.chem"
    And identifiers editor should contain text "For research, add them under Explicit References."
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A configuration made for another catalog brings its identifiers and leaves the editor where it is
    When user clicks on Import JSON button in identifiers editor
    Then "Import Configuration" dialog should be visible
    When user enters '{"connectionName":"CatalogsDemo","nqName":"CatalogsDemo","dataSourceName":"MS SQL","catalogName":"sales","schemaName":"dbo","entryPoints":{"BDD_IMPORTED_ID":{"table":"batches","column":"batch_id"}}}' into configuration text
    And user clicks on Import button in "Import Configuration" dialog
    Then an info balloon containing "imported into CatalogsDemo.research.chem" should have been shown
    And identifiers editor should contain text "BDD_IMPORTED_ID"
    And identifiers editor should contain text "Catalogs Demo: research.chem"

  Scenario: An identifier is previewed from the picked catalog and saved with it
    When user clicks on add identifier icon
    Then "Add Identifier" dialog should be visible
    And Table input in "Add Identifier" dialog should offer the choice "batches"
    When user enters "BDD_COMPOUND_ID" into Semantic Type input in "Add Identifier" dialog
    And user selects "compounds" in Table input in "Add Identifier" dialog
    And user selects "compound_id" in Column input in "Add Identifier" dialog
    Then identifier preview should contain text "DGC-"
    And identifier preview should not contain text "ARCHIVED"
    When user enters "DGC-\d{6}" into Match Regexp input in "Add Identifier" dialog
    And user clicks on ADD button in "Add Identifier" dialog
    Then identifiers editor should contain text "BDD_COMPOUND_ID"
    When user clicks on Save button in identifiers editor
    Then an info balloon containing "Identifier configuration saved" should have been shown
    And the identifiers configurations of the "CatalogsDemo" connection on the server should cover "research.chem"
    And the saved identifiers of "research.chem" in the "CatalogsDemo" connection should be "BDD_COMPOUND_ID, BDD_IMPORTED_ID"
    And no errors should have been logged

  Scenario: The configuration is listed under the connection by catalog and schema
    Given the browse panel is open
    And Databases---MS-SQL---Catalogs-Demo tree node inside browse tree is expanded
    Then Databases---MS-SQL---Catalogs-Demo---Identifiers---research.chem tree node inside browse tree should be visible

  Scenario: After a reload the identifier opens the row of the picked catalog
    When user reloads the page
    Given the package autostarts have completed
    And the context panel is open
    And user opens a table "bdd-compounds" with:
      | compound   | note   |
      | DGC-000001 | first  |
      | DGC-000005 | second |
    Then "compound" column should have semantic type "BDD_COMPOUND_ID"
    When user clicks on the "cell 1 of compound" area of grid
    Then context panel should contain text "Aspirin"
    And context panel should not contain text "ARCHIVED"
    And no errors should have been logged
    And no error or warning balloon should have been shown
