@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). Browse in a schema's right-click menu opens the schema view: one box per
  table with its columns (the biologics schema of the Biologics connection has 15 tables). A click on a box shows the
  table in the context panel, and a box's right-click menu is the table's own; Get Top 100 on drugs (200 rows) opens
  its first 100 rows as a table. Read-only: nothing is made on the server.

  Scenario: See a whole schema at once
    Lay out every table of a database schema on one page, then open one from there.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the context panel is open
    And the browse panel is open
    When the video says "One page with every table of a schema"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree
    # caption: Browse on a schema opens the schema view
    And user picks "Browse" from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree
    Then the "Schema: biologics" view should be current
    And "drugs" schema table should be visible
    # caption: Each box is a table and its columns
    When user hovers over "expression_batches" schema table
    # caption: Click a table's name for its details
    And user clicks on "sequences" text inside "sequences" schema table
    Then the context panel should show "sequences"
    # caption: Click a column for the column's details
    When user clicks on "heavy_chain" text inside "sequences" schema table
    Then the context panel should show "heavy_chain"
    # caption: Right-click a box for the table's commands
    When user picks "Get Top 100" from the context menu of "drugs" schema table
    Then the "drugs" view should be current
    And the table should have 100 rows
    When the video says "Get Top 100: the first hundred rows of the table" at the bottom
    And user hovers over "Rows: 100" text in status bar
    Then no errors should have been logged
