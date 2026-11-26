@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). New Query... on a connection opens the query editor; the clip names the
  query, types a select over biologics.assay_types (23 assay types in 6 categories) that counts the types of each
  category, runs it with the play icon (the result shows under the editor), saves it, and finds it listed under the
  Biologics connection in Browse. The connection is the stand's own Postgres (no outside database).

  The saved query, "Assay types per category", is the clip's own: it is deleted when the feature ends, and one an
  earlier run left is deleted at its start.

  Scenario: Write a SQL query
    Write SQL against a database connection, run it, and save it as a query.

    Given user is logged in
    And simple mode is off
    And no query named "Assay types per category" is on the server
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "Write your own SQL on a connection"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    # caption: New Query... opens the SQL editor for this connection
    And user picks "New Query..." from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree
    Then the current view should be a DataQueryView view
    # caption: Give the query a name
    When user enters "Assay types per category" into Name input
    # caption: Type the SQL
    And user replaces the code of code editor with "select category, count(*) from biologics.assay_types group by category"
    # caption: Run it: the result shows under the editor
    And user runs the query with the play icon
    Then grid should be visible
    And the "rows" reading of grid should be 6
    # caption: Save keeps it on the server
    When user clicks on Save button
    Then 1 query named "Assay types per category" should be on the server
    When user clicks on browse tab
    # caption: It is listed under its connection, ready to run
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    Then Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree should be visible
    When user hovers over Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree
    Then no errors should have been logged
