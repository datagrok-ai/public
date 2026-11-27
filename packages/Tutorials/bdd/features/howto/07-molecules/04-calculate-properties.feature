@guide @tutorial @viewport:1280x720
Feature: How-to · Molecules
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Properties as columns, then put to use: Chem > Calculate >
  Chemical Properties... opens on the molecule column with MW ticked; with the Chemical Properties calculator ticked,
  and Log P and PSA added, OK adds an MW, a LogP and a PSA column for every molecule, computed in the browser (Chem's
  Crux engine, RDKit's definitions). A line says what the Rule of Five asks of them; in the table's Search box
  (Ctrl+F, in the toolbox) "MW > 500" with Enter keeps the 211 molecules heavier than 500 and the status bar counts
  them; a double-click on the MW header then sorts them by it, heaviest first (the first double-click sorts
  descending). smiles (1000
  molecules) is opened off camera with the grid showing molregno and canonical_smiles; the new columns, the sort and
  the filter live in the page's table only and nothing is made on the server.

  Scenario: Calculate properties for every molecule
    Add molecular weight, LogP and polar surface area as columns, then filter and sort the molecules by them.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens smiles dataset
    And the grid shows only the columns "molregno, canonical_smiles"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Chem > Calculate > Chemical Properties
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And Molecules input in "Chemical Properties" dialog should contain text "canonical_smiles"
    And MW input in "Chemical Properties" dialog should be checked
    # caption: Tick the calculator and the properties you need
    When user checks "Chemical Properties" calculator in "Chemical Properties" dialog
    And user checks "Log P" input in "Chemical Properties" dialog
    And user checks PSA input in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then "Chemical Properties" dialog should be hidden
    And 3 new columns should have been added
    And a new column "MW" should have been added
    And a new column "LogP" should have been added
    And a new column "PSA" should have been added
    And every value of "MW" column should lie between 1 and 5000
    When the video says "A new column per property, computed in your browser" at the bottom
    And the video says "MW is the weight; LogP, oil versus water; PSA, the polar surface" at the bottom
    And the video says "Most oral drugs stay under MW 500 and LogP 5" at the bottom
    # caption: The Search box (Ctrl+F) takes a condition
    And user presses Control+f in grid overlay
    Then table search should be focused
    # caption: Enter keeps only the rows that match
    When user types "MW > 500" into table search
    And user presses Enter in table search
    Then 211 rows should pass the filter
    And the filter should pass exactly the rows where "MW" is between 500.001 and 5000
    When the video says "211 of the 1000 molecules are heavier than 500" at the bottom
    And user hovers over "Filtered: 211" text in status bar
    # caption: Double-click a header to sort by it: heaviest first
    And user double-clicks on the "header MW" area of grid
    Then the "sort column" reading of grid should be "MW"
    And the "sort direction" reading of grid should be "descending"
    When the video says "The heaviest molecules now come first" at the bottom
    Then 211 rows should pass the filter
    And the table should have 1000 rows
    And no errors should have been logged
