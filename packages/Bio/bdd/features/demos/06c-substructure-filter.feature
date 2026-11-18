@guide @demo
Feature: 2 · Every sequence is a molecule · Substructure filter
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). The first 63 peptides of the BILN_complex sample, made molecules off camera
  (Bio | Transform | To Atomic Level, its defaults) and sorted by name from z to a, so half of the first screen
  (tyrocidine, the trivalent scaffold, tirzepatide) has no disulfide and drops out. The molecule column's
  substructure filter card, sketched on camera — a disulfide, C-S-S-C, typed as SMILES into the sketcher — keeps
  exactly the peptides whose molecule holds one (44 of the 63: Cys–Cys, D-Pen–D-Pen and Mpa–Cys bridges alike),
  whatever their notation wrote; their BILN beside them draws each disulfide as a bracket between its two Cys. The
  search runs on the molecules Bio built, in the browser (Chem's engine, Crux on this stand). Set up off camera:
  the conversion, the sort, the columns (the molecules 560 px wide, in 120 px rows), no current row, monomer names
  up to eight letters (D-Pen) and the filter panel with only the molecule column's card.

  Scenario: One chemical motif across every peptide
    A disulfide sketched once finds every peptide that holds one, whatever its notation.

    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "OpenChemLib"
    And the "Chem" package is installed
    And the Bio package is initialized
    And user opens BILN_complex dataset keeping the first 63 rows
    And the "biln" column is converted to atomic level
    And the grid is sorted by the "name" column in descending order
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "name, biln, molfile(biln)"
    And the "biln" column of the grid is 500 px wide
    And the "molfile(biln)" column of the grid is 560 px wide
    And the monomers of the "biln" column are drawn at 14 px
    And the monomer names of the "biln" column are cut at 8 characters
    And the rows of the grid are 120 px tall
    And no row of the table is current
    And user opens an empty filter panel
    And user adds a card for "molfile(biln)" to the filter panel
    When user clicks on "Sketch" text in "molfile(biln)" filter card
    Then sketcher dialog should be visible
    When user types "CSSC" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 44 rows should pass the filter
    When the video says "44 of 63 peptides carry a disulfide bridge"
    Then no error or warning balloon should have been shown
