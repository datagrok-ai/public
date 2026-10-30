@guide
Feature: MolTrack: filter compounds by structure and properties
  A guide: the answer to "how do I find the compounds that contain one substructure or another, but not some
  element?". The Filters pane left of every list builds the conditions row by row; Advanced nests groups, each
  with its own AND or OR. A structure column offers "has substructure", "is substructure of", "is similar to"
  and "is exactly", the other columns their usual operators. Every complete row narrows the list at once, and
  the address keeps the filter.

  Scenario: Phenols or anilines, without chlorine
    Given user is logged in
    And simple mode is off
    And the MolTrack demo data is on the server
    And user opens the "Compounds" view of MolTrack
    When user remembers the text of compound count
    # caption: The filters are on the left; Advanced nests groups of conditions
    And user clicks on advanced mode in filters pane
    # caption: A group inside: either of two substructures
    And user clicks on add group in filters pane
    And user clicks on last add condition in filters pane
    And user selects "Parent Structure" in first condition property
    And user selects "has substructure" in first condition operator
    # caption: The first: a phenol, drawn in the sketcher
    And user clicks on molecule canvas in first condition value
    And user types "c1ccccc1O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    # caption: The second, OR the first: an aniline
    And user clicks on last add condition in filters pane
    And user selects "Parent Structure" in second condition property
    And user selects "has substructure" in second condition operator
    And user clicks on molecule canvas in second condition value
    And user types "c1ccccc1N" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then compound count should not have the remembered text
    And compound list should contain text "Triclosan"
    # caption: AND, outside the group: no chlorine in the formula
    When user clicks on first add condition in filters pane
    And user selects "Formula" in third condition property
    And user selects "does not contain" in third condition operator
    And user types "Cl" into third condition value
    And user presses Tab in third condition value
    # caption: Phenols and anilines without chlorine: triclosan is gone
    Then compound list should contain text "Bisphenol A"
    And compound list should not contain text "Triclosan"
    # caption: Remove the chlorine condition: triclosan is back
    When user clicks on third condition remove
    Then compound list should contain text "Triclosan"
