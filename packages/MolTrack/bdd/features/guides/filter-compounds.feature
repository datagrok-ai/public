@guide
Feature: MolTrack EMS: filter compounds by structure and properties
  A guide: the answer to "how do I find the compounds that contain one substructure or another, but not some
  element?". The Filters pane left of every list builds the conditions row by row; Advanced nests groups, each
  with its own AND or OR. A structure column offers "has substructure", "is substructure of", "is similar to"
  and "is exactly", the other columns their usual operators. Every complete row narrows the list at once, and
  the address keeps the filter.

  Scenario: Phenols or anilines, without chlorine
    # caption: Sign in
    Given user is logged in
    # caption: Simple mode is off
    And simple mode is off
    # caption: The demo data is loaded
    And the MolTrack EMS demo data is on the server
    # caption: Open the Compounds view of MolTrack
    And user opens the "Compounds" view of MolTrack EMS
    # caption: Note how many compounds there are
    When user remembers the text of compound count
    # caption: Switch the filters to Advanced, to nest conditions
    And user clicks on advanced mode in filters pane
    # caption: Add a group: either of two substructures
    And user clicks on add group in filters pane
    # caption: Add its first condition
    And user clicks on last add condition in filters pane
    # caption: On the parent structure
    And user selects "Parent Structure" in first condition property
    # caption: Has substructure
    And user selects "has substructure" in first condition operator
    # caption: Open the sketcher
    And user clicks on molecule canvas in first condition value
    # caption: Type a phenol
    And user types "c1ccccc1O" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept it
    And user clicks on OK button in sketcher dialog
    # caption: Add the second condition to the group
    And user clicks on last add condition in filters pane
    # caption: Again on the parent structure
    And user selects "Parent Structure" in second condition property
    # caption: Has substructure
    And user selects "has substructure" in second condition operator
    # caption: Open the sketcher
    And user clicks on molecule canvas in second condition value
    # caption: Type an aniline
    And user types "c1ccccc1N" into molecule input of sketcher dialog
    # caption: Press Enter to draw it
    And user presses Enter in molecule input of sketcher dialog
    # caption: Accept it: phenols or anilines
    And user clicks on OK button in sketcher dialog
    # caption: The list narrows
    Then compound count should not have the remembered text
    # caption: Triclosan is among them
    And compound list should contain text "Triclosan"
    # caption: Add a condition outside the group
    When user clicks on first add condition in filters pane
    # caption: On the formula
    And user selects "Formula" in third condition property
    # caption: Does not contain
    And user selects "does not contain" in third condition operator
    # caption: Type Cl, for chlorine
    And user types "Cl" into third condition value
    # caption: Press Tab to apply it
    And user presses Tab in third condition value
    # caption: Bisphenol A stays
    Then compound list should contain text "Bisphenol A"
    # caption: Triclosan, with chlorine, is gone
    And compound list should not contain text "Triclosan"
    # caption: Remove the chlorine condition
    When user clicks on third condition remove
    # caption: Triclosan is back
    Then compound list should contain text "Triclosan"
