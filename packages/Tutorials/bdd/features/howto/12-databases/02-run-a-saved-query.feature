@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). A saved query sits under its connection in Browse; Run asks for its
  parameters. The query is the Biologics package's "Antibody Profile by Organism and Target" on the stand's own
  Postgres: two choice parameters whose lists come from the database (target_organisms, targets), and a saved
  layout (antibody_profile.layout: filters, grid, form, pivot table, pie chart) that its result opens with. The clip
  runs it for E. coli and HER2, then picks PD-1 in the toolbox's Source pane and REFRESH runs it again. Nothing is
  saved: the result views are closed by the shell reset.

  Scenario: Run a saved query
    Run a saved query, pick its parameters, and run it again with other values.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "Saved queries are listed under their connection"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    # caption: Run asks for the query's parameters
    And user picks "Run" from the context menu of Databases---PostgresDart---Biologics---Antibody-Profile-by-Organism-and-Target tree node inside browse tree
    Then "Antibody Profile by Organism and Target" dialog should be visible
    # caption: The choices come from the database
    When user selects "E. coli" in Organism input in "Antibody Profile by Organism and Target" dialog
    And user selects "HER2" in Target input in "Antibody Profile by Organism and Target" dialog
    And user clicks on OK button in "Antibody Profile by Organism and Target" dialog
    Then the current view should be a TableView view
    And every value of "organism_name" column should contain "E. coli"
    And every value of "target_name" column should contain "HER2"
    And the open tableview should have 1 pie chart viewer
    When the video says "The result opens with the layout saved with the query" at the top
    And user hovers over pie chart viewer
    Then toolbox should be visible
    # caption: The toolbox keeps the parameters: pick another target
    When user selects "PD-1" in Target input in toolbox
    # caption: REFRESH runs the query again
    And user clicks on REFRESH button in toolbox and the query it runs completes
    Then every value of "target_name" column should contain "PD-1"
    And every value of "organism_name" column should contain "E. coli"
    When the video says "Same view, new rows from the database" at the top
    And user hovers over pie chart viewer
    Then no errors should have been logged
