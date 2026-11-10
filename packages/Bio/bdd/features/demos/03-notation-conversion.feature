@guide @demo
Feature: Bioinformatics · Notation conversion
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio > Transform > Convert Sequence Notation rewrites a sequence
  column in another notation: the cyclic peptides of the HELM_cyclic sample, written in HELM and drawn as
  rings, become BILN, whose head-to-tail bond the renderer draws as a bracket, and separator text.

  Scenario: HELM to BILN and back to plain text
    One dialog converts between HELM, BILN, FASTA and separator notation, keeping non-natural monomers and ring bonds.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_cyclic dataset
    Then "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "helm"
    Given the grid shows only the columns "sequence_id, sequence"
    And the "sequence" column of the grid is 300 px wide
    And the rows of the grid are 130 px tall
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    When user selects "biln" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then a new column "biln(sequence)" should have been added
    And "biln(sequence)" column should have units "biln"
    Given the "biln(sequence)" column of the grid is 560 px wide
    And the monomers of the "biln(sequence)" column are drawn at 14 px
    When the video says "In BILN the head-to-tail bond is part of the line, drawn as a bracket"
    And user hovers over monomer 1 of the "cell 1 of biln(sequence)" area of grid
    And user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    And user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then a new column "separator(sequence)" should have been added
    Given the "separator(sequence)" column of the grid is 560 px wide
    And the monomers of the "separator(sequence)" column are drawn at 14 px
    Then no errors should have been logged
