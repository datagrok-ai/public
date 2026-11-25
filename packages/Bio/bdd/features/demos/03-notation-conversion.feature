@guide @demo
Feature: 1 · Any notation · Notation conversion
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio > Transform > Convert Sequence Notation rewrites a sequence
  column in another notation: the cyclic peptides of the HELM_cyclic sample, written in HELM and drawn as
  rings, become BILN, whose head-to-tail bond the renderer draws as a bracket, monomer names up to twelve
  letters in full. The second conversion, to separator text, is the same dialog run again off camera; the
  clip ends on the first four peptides side by side in the three notations, each a single line: the HELM
  ring, the BILN line with its bracket, the separator line without one (separator text has no way to write
  the ring bond). Set up off camera: the columns, their widths and fonts, the renderers' Max Monomer Length
  (seven letters in the last scene, so both lines fit beside the ring) and no row header there. Left out:
  separator text writes a dashed name (D-Cit) unbracketed and its renderer splits it in two, which a single
  line does not show.

  Scenario: HELM to BILN and back to plain text
    One dialog converts HELM, BILN, FASTA and separator, keeping non-natural monomers and rings.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_cyclic dataset
    Then "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "helm"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence_id, sequence"
    And the "sequence" column of the grid is 300 px wide
    And the rows of the grid are 130 px tall
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    When user selects "biln" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then a new column "biln(sequence)" should have been added
    And "biln(sequence)" column should have units "biln"
    Given the "biln(sequence)" column of the grid is 1080 px wide
    And the monomers of the "biln(sequence)" column are drawn at 14 px
    And the monomer names of the "biln(sequence)" column are cut at 12 characters
    When the video says "BILN writes the head-to-tail bond into the line"
    And user hovers over monomer 1 of the "cell 1 of biln(sequence)" area of grid
    And the video says "Separator: plain text, without the ring bond"
    Given user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    And user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then a new column "separator(sequence)" should have been added
    And "separator(sequence)" column should have units "separator"
    Given the tooltip is hidden
    And the grid shows only the columns "sequence, biln(sequence), separator(sequence)"
    And user sets "showRowHeader" property of grid to "false"
    And the "sequence" column of the grid is 210 px wide
    And the "biln(sequence)" column of the grid is 654 px wide
    And the "separator(sequence)" column of the grid is 654 px wide
    And the monomers of the "biln(sequence)" column are drawn at 12 px
    And the monomers of the "separator(sequence)" column are drawn at 12 px
    And the monomer names of the "biln(sequence)" column are cut at 7 characters
    And the monomer names of the "separator(sequence)" column are cut at 7 characters
    And the rows of the grid are 150 px tall
    And no row of the table is current
    Then no errors should have been logged
