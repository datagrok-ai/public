@guide @demo @viewport:1440x810
Feature: Bioinformatics · To atomic level
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio | Transform | To Atomic Level turns a whole column of sequences into
  molecules on Chem's macromolecule assembler: first HELM peptides of non-natural monomers with Non-linear off
  (plain chains drawn straight) and Highlight monomers on (each atom colored by the monomer it comes from), then
  disulfide-bridged peptide drugs in BILN in Non-linear mode (drawn as 2D depictions, rings and all). The pointer
  on a monomer of a sequence lights that monomer's atoms in the molecule of its row. Set up off camera: the
  molecule column left of the sequences, so the monomer tooltips never cover it, and no current row, so its color
  does not tint the molecules.

  Scenario: Sequences become molecules, monomer by monomer
    A column of peptides becomes molecules in one step, linear or cyclic; each monomer's atoms light up under the pointer.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_sample dataset keeping the first 5 rows
    And the toolbox pane is hidden
    And the context panel is closed
    And the "HELM" column of the grid is 560 px wide
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    And user unchecks "Non-linear" checkbox in "To Atomic Level" dialog
    And user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column "molfile(HELM)" should have been added
    And "molfile(HELM)" column should have no missing values
    Given the grid shows only the columns "molfile(HELM), HELM"
    And the "molfile(HELM)" column of the grid is 720 px wide
    And the "HELM" column of the grid is 560 px wide
    And the rows of the grid are 230 px tall
    And no row of the table is current
    When the video says "Every atom colored by the monomer it comes from"
    And user hovers over monomer 1 of the "cell 2 of HELM" area of grid
    And user hovers over monomer 7 of the "cell 2 of HELM" area of grid
    And user hovers over monomer 14 of the "cell 2 of HELM" area of grid
    And the video says "Non-linear mode: cycles and bridges drawn as rings"
    Given user opens BILN_complex dataset keeping the first 6 rows
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "name, biln"
    And the monomers of the "biln" column are drawn at 16 px
    And the "biln" column of the grid is 560 px wide
    And the rows of the grid are 60 px tall
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    Then "Non-linear" checkbox in "To Atomic Level" dialog should be checked
    When user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column "molfile(biln)" should have been added
    And "molfile(biln)" column should have no missing values
    Given the grid shows only the columns "name, molfile(biln), biln"
    And the "molfile(biln)" column of the grid is 560 px wide
    And the "biln" column of the grid is 560 px wide
    And the rows of the grid are 230 px tall
    And no row of the table is current
    When user hovers over monomer 1 of the "cell 1 of biln" area of grid
    And user hovers over monomer 4 of the "cell 2 of biln" area of grid
    Then no error or warning balloon should have been shown
