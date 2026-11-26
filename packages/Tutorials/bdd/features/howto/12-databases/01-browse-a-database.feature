@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). Browse > Databases lists the connections by database type; the
  Biologics connection (PostgresDart, the stand's own Postgres, made by the Biologics package) holds its saved queries
  and a Schemas group. The clip opens the biologics schema, clicks the drugs table (200 rows: id, identifier, name,
  smiles), which shows a preview of its first rows and its details in the context panel, expands the table to its
  columns, clicks one, and opens the whole table with Get All. Read-only: nothing is made on the server.

  Scenario: Browse a database
    Find a database in Browse, look through its tables and columns, and open a table.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the browse panel is closed
    And the context panel is open
    When the video says "Databases are in Browse, under Databases"
    And user clicks on browse tab
    And user expands Databases tree node inside browse tree
    # caption: Connections are grouped by database type
    And user expands Databases---PostgresDart tree node inside browse tree
    # caption: Point at a connection to see its server and database
    And user hovers over Databases---PostgresDart---Biologics tree node inside browse tree
    # caption: A connection holds its saved queries and its schemas
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree
    # caption: A schema lists its tables
    And user expands Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree
    # caption: Click a table to preview its first rows
    And user clicks on Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree
    Then the context panel should show "drugs"
    # caption: The context panel lists its columns
    When user expands Columns section in context panel
    # caption: Or expand the table in the tree
    And user expands Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree
    Then Databases---PostgresDart---Biologics---Schemas---biologics---drugs---smiles tree node inside browse tree should be visible
    # caption: Get All opens the whole table
    When user picks "Get All" from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree
    Then the "drugs" view should be current
    And the table should have 200 rows
    And the table should have the columns "id, identifier, name, smiles"
    When the video says "Now it is a table like any other: sort, filter, chart" at the bottom
    And user hovers over "Rows: 200" text in status bar
    Then no errors should have been logged
