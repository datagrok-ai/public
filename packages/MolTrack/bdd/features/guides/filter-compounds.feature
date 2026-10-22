@guide
Feature: MolTrack EMS: filter compounds by structure and properties
  A guide: the answer to "how do I find the compounds that contain a substructure, or look like a
  molecule?". The Filters pane left of every list (the funnel opens it) builds the conditions row by row, as the legacy
  search did: a structure column offers "has substructure", "is substructure of", "is similar to" (with
  the similarity threshold beside the molecule) and "is exactly"; the other columns their usual
  operators. Every complete row narrows the list at once, and the address keeps the filter.

  Scenario: Phenols, then the compounds similar to bisphenol A
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And user opens the "Compounds" view of MolTrack EMS
    When user remembers the text of compound count
    # caption: Open the filters, and add a condition
    And user clicks on "Show the filters" button
    And user clicks on add condition in filters pane
    And user selects "Parent Structure" in condition property
    And user selects "has substructure" in condition operator
    # caption: Draw the substructure: a phenol
    And user clicks on molecule canvas in condition value
    And user types "c1ccccc1O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    # caption: Only the phenols are left
    Then compound count should not have the remembered text
    And compound list should contain text "Triclosan"
    # caption: Now the compounds similar to bisphenol A, at 0.4 or more
    When user selects "is similar to" in condition operator
    And user clicks on molecule canvas in condition value
    And user types "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    And user enters "0.4" into condition threshold
    Then compound list should contain text "4-(1,1,3,3-Tetramethylbutyl)phenol"
    # caption: Remove the condition: every compound is back
    When user clicks on condition remove
    Then compound count should have the remembered text
