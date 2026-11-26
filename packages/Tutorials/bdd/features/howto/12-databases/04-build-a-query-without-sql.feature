@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). New Visual Query... on a table opens the query builder: rows for Data,
  Where, Group by, Aggregate (Pivot and Having show once there is an aggregation) and Order by, with the result under
  them, updated as the query changes. On biologics.expression_batches (1200 batches: expression_system CHO, E.coli or
  HEK293, chain HC or LC, yield_mg) the clip groups by expression_system (3 rows), aggregates yield_mg (the picker's
  default aggregation is avg) and pivots by chain (an HC and an LC column). Nothing is saved: the builder's view is
  closed by the shell reset. The steps for the builder are in bindings/howto-databases.ts.

  Scenario: Build a query without SQL
    Pick columns to group, aggregate and pivot a database table, and watch the result as you build.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "A query built by picking columns, no SQL"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree
    # caption: New Visual Query... on a table opens the builder
    And user picks "New Visual Query..." from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics---expression_batches tree node inside browse tree
    Then the current view should be a DataQueryView view
    When the video says "The result under the builder follows every change" at the bottom
    And user hovers over "Grid" grid
    # caption: Group by: one row per expression system
    And user adds "expression_system" to the "Group by" row of the visual query
    Then the "Group by" row of the visual query should hold "expression_system"
    And the visual query should preview 3 rows
    # caption: Aggregate: the average yield of each
    When user adds "yield_mg" to the "Aggregate" row of the visual query
    Then the "Aggregate" row of the visual query should hold "avg(yield_mg)"
    # caption: Pivot: a column for each chain
    When user adds "chain" to the "Pivot" row of the visual query
    Then the "Pivot" row of the visual query should hold "chain"
    And the visual query should preview 3 rows
    When the video says "Average yield by expression system and chain" at the bottom
    And user hovers over "Grid" grid
    Then no errors should have been logged
