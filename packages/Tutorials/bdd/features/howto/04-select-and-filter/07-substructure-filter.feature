@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The substructure filter: a molecule column's card in the
  filter panel opens a sketcher on a click; a structure typed as SMILES (or drawn) and confirmed with OK keeps
  the molecules that contain it, and the card's settings switch the search to Not contains. The search runs in
  the browser (Chem's RDKit). spgi (100 molecules) is opened off camera with its default filter panel, and the
  OpenChemLib sketcher is pinned for the feature and put back after it; nothing is made on the server.

  Scenario: Filter molecules by a substructure
    Draw or paste a structure in the filter panel to keep only the molecules that contain it.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens the filter panel
    Then all rows should pass the filter
    # caption: Sketch opens the molecule editor
    When user clicks on "Sketch" text in "Structure" filter card
    Then sketcher dialog should be visible
    # caption: Draw a structure, or type or paste its SMILES
    When user clicks on molecule input of sketcher dialog
    And user types "C1CCNC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 33 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "C1CCNC1"
    When the video says "33 of 100 molecules contain the ring; it is highlighted in each"
    # caption: The search type can flip it to Not contains
    And user opens the settings of the "Structure" filter card
    And user picks search type "Not contains" in the "Structure" filter card
    Then the "search type of Structure" reading of filter panel should be "Not contains"
    And 67 rows should pass the filter
    And no errors should have been logged
