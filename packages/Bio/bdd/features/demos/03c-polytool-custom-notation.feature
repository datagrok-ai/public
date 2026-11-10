@guide @demo
Feature: Bioinformatics · Custom notation with PolyTool
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). SequenceTranslator's cyclized sample writes cyclic peptides in a
  custom notation, ring closures marked C(1)…C(1); Bio > PolyTool > Convert applies the notation's rules
  and writes each sequence as HELM, drawn as its ring, and as a molecule.

  Scenario: Custom notation to HELM and molecules
    PolyTool reads a custom peptide notation by its rules and writes each sequence as HELM and as a molecule.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens cyclized dataset
    Then "seqs" column should have semantic type "Macromolecule"
    Given the "seqs" column of the grid is 520 px wide
    And the monomers of the "seqs" column are drawn at 15 px
    And the rows of the grid are 150 px tall
    When user picks "Bio > PolyTool > Convert..." from the top menu
    Then "PolyTool Conversion" dialog should be visible
    When user clicks on OK button in "PolyTool Conversion" dialog
    Then the "PolyTool Conversion" dialog should close
    And a new column "transformed(seqs)" should have been added
    And a new column "molfile(seqs)" should have been added
    Given the "seqs" column of the grid is 420 px wide
    And the "transformed(seqs)" column of the grid is 300 px wide
    And the "molfile(seqs)" column of the grid is 340 px wide
    When the video says "Ring closures C(1)…C(1) become the HELM ring and the molecule's macrocycle"
    And user hovers over monomer 3 of the "cell 1 of transformed(seqs)" area of grid
    And user hovers over the "cell 4 of n" area of grid
    Then no errors should have been logged
