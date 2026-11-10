@guide @demo
Feature: Bioinformatics · Monomer libraries
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). Bio > Manage > Monomers opens the monomers of a library
  as a table (structure, symbol in its color, name, the cap of each attachment point) beside the monomer
  editor, which follows the current row: the structure in the sketcher, a chip per attachment point
  with its cap, and a click on a chip lists the common caps. The Library input shows another library's
  monomers, and Bio > Manage > Monomer Libraries lists the library files the platform loads together,
  with the symbols two of them define differently. Browsing only: nothing is saved, switched or picked,
  so no library and no setting on the server changes (the editor's Save is the monomer-editor
  feature's, the library switches the libraries feature's).

  Scenario: Browse the monomer libraries
    Every monomer of a library in a table, beside an editor that shows its structure, attachment points and caps.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And user opens HELM_sample dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "Monomers: the building blocks of every sequence"
    And user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the monomer sketcher of the Manage Monomers view should be ready
    And "R1" column should have semantic type "Molecule"
    Given the grid shows only the columns "Monomer, Symbol, Name, R1, R2"
    And the "R1" column of the grid is 190 px wide
    And the "R2" column of the grid is 190 px wide
    When user clicks on the "cell 4 of Name" area of grid
    Then "Name" input in monomer form should have value "Glutamic acid"
    And "R3" R-group chip in monomer form should contain the text "R3 OH"
    When user clicks on "R3" R-group chip in monomer form
    Then "OH" menu item should be checked
    When user hovers over "Azide" menu item
    And user presses Escape
    And user selects "polytool-lib" in "Library" input in monomer form
    Then the table should have 12 rows
    And "R1" column should have semantic type "Molecule"
    Given the grid shows only the columns "Monomer, Symbol, Name, R1, R2"
    And the "R1" column of the grid is 190 px wide
    And the "R2" column of the grid is 190 px wide
    When user clicks on the "cell 4 of Name" area of grid
    Then "Name" input in monomer form should have value "Azido Glycine"
    When user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the "Manage Monomer Libraries" view should be current
    And "HELMCoreLibrary.json" checkbox should be checked
    When the video says "Library files load together; clashing symbols are resolved here"
    Then no error or warning balloon should have been shown
    And no errors should have been logged
