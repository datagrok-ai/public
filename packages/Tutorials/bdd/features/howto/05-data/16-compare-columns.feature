@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The tennis players demo table (Files > Demo >
  tennis_players.csv, from DBpedia) has the country a player represents and their birthPlace. Data > Compare
  Columns (PowerPack) compares two columns row by row: with tennis_players as the table and country and
  birthPlace as the columns it counts 685 matches and 129 mismatches; Select selects the 129 rows. Some differ for
  real (row 2, Martina Navratilova: United States, born in the Czech Republic), some are DBpedia's mess (row 5,
  country "United Kingdom,United Kingdom"). The dialog has no OK: it is closed from its title bar. Off camera the
  image column (wikimedia URLs the grid would fetch) and the long abstract are removed and the grid shows name,
  country and birthPlace. Nothing is made on the server.

  Scenario: Compare two columns
    Check row by row where two columns agree, and select the rows where they do not.

    Given user is logged in
    And simple mode is off
    And user opens tennis-players dataset
    And user removes "image" column
    And user removes "abstract" column
    And the grid shows only the columns "name, country, birthPlace"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: country: who they play for; birthPlace: where they were born
    When user hovers over the "cell 2 of birthPlace" area of grid
    # caption: Compare Columns: where two columns agree, row by row
    And user picks "Data > Compare Columns..." from the top menu
    Then Compare Columns dialog should be visible
    When user selects "tennis_players" in "Tables" input in Compare Columns dialog
    # caption: The two columns to compare
    And user selects "country" in "Columns" input in Compare Columns dialog
    And user selects "birthPlace" in 4th input in Compare Columns dialog
    Then Compare Columns dialog should contain text "Matched: 685"
    And Compare Columns dialog should contain text "Mismatched: 129"
    # caption: Select: the rows where they differ
    When user clicks on "Select" link in Compare Columns dialog
    Then 129 rows should be selected
    When user clicks on "Close" icon in Compare Columns dialog
    Then Compare Columns dialog should be hidden
    And the value of "birthPlace" column in row 2 should be "Czech Republic"
    And the value of "country" column in row 5 should be "United Kingdom,United Kingdom"
    # caption: Some players were born in another country
    When user hovers over the "cell 2 of birthPlace" area of grid
    # caption: Others are messy data: a country written twice
    And user hovers over the "cell 5 of country" area of grid
    Then no errors should have been logged
