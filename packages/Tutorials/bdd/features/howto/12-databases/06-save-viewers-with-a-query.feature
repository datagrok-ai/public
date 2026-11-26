@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). The Layout tab of the query editor sets how the query's result looks:
  run the query there, add viewers from the toolbox, and SAVE stores them with the query, so every run of it opens
  with them. The clip writes a select of expression_system, chain and yield_mg from biologics.expression_batches
  (1200 rows) on the Biologics connection, adds a box plot and a histogram in the Layout tab, saves the query as
  "Expression yields", and runs it from Browse: the result opens with both viewers.

  The query and its layout are the clip's own: deleted when the feature ends, and what an earlier run left is
  deleted at its start.

  Scenario: Save viewers with a query
    Arrange viewers in a query's Layout tab, so every run of the query opens with them.

    Given user is logged in
    And simple mode is off
    And no query named "Expression yields" is on the server
    And the user's own layouts named "Expression yields, ExpressionYields" are removed now and at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "A query can keep its viewers: its result opens with them"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    And user picks "New Query..." from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "Expression yields" into Name input
    And user replaces the code of code editor with "select expression_system, chain, yield_mg from biologics.expression_batches"
    # caption: The Layout tab sets how the result looks
    And user clicks on Layout tab
    Then "Run query to get data and edit layout" text should be visible
    # caption: Run the query to fill it
    When user runs the query with the play icon
    Then grid should be visible
    # caption: Add viewers from the toolbox
    When user clicks on box plot icon in toolbox
    And user clicks on histogram icon in toolbox
    Then box plot viewer should be visible
    And histogram viewer should be visible
    # caption: SAVE keeps them with the query
    When user clicks on Save button
    Then 1 query named "Expression yields" should be on the server
    And the query "Expression yields" on the server should have a layout
    # caption: Close the editor
    When user closes "Expression yields" view tab
    Then "Expression yields" view tab should be absent
    # the closed tab's tooltip stays on the page (the tab is removed under the pointer, its mouseleave never comes);
    # Browse comes back by itself when the last view closes
    Given the tooltip is hidden
    # caption: Run the saved query from Browse
    When user expands Databases---PostgresDart---Biologics tree node inside browse tree
    And user picks "Run" from the context menu of Databases---PostgresDart---Biologics---Expression-yields tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 1200 rows
    And the open tableview should have 1 box plot viewer
    And the open tableview should have 1 histogram viewer
    When the video says "The result opens with the viewers saved with the query" at the top
    And user hovers over box plot viewer
    Then no errors should have been logged
