@guide @demo
Feature: Bioinformatics · Nucleotides to atoms
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio | Transform | To Atomic Level on Bio's FASTA DNA and RNA samples: each
  letter becomes a whole nucleotide (sugar, base and phosphate, deoxyribose for DNA, ribose for RNA), the strand
  drawn straight (Non-linear off) and every nucleotide colored like its letter (Highlight monomers on); the
  pointer on a letter lights that nucleotide's atoms. The default window, not a zoomed one: a strand of 35
  nucleotides needs the width. Set up off camera: the molecule column left of the sequences, so the letter
  tooltips never cover it, the letters drawn larger, and no current row.

  Scenario: DNA and RNA, nucleotide by nucleotide
    Each letter of a DNA or RNA sequence becomes a whole nucleotide (sugar, base and phosphate), lit under the pointer.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_DNA dataset keeping the first 12 rows
    And the toolbox pane is hidden
    And the context panel is closed
    And the monomers of the "sequence" column are drawn at 14 px
    And the "sequence" column of the grid is 390 px wide
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    And user unchecks "Non-linear" checkbox in "To Atomic Level" dialog
    And user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column "molfile(sequence)" should have been added
    And "molfile(sequence)" column should have no missing values
    Given the grid shows only the columns "molfile(sequence), sequence"
    And the "sequence" column of the grid is 390 px wide
    And the "molfile(sequence)" column of the grid is 1110 px wide
    And the rows of the grid are 150 px tall
    And no row of the table is current
    When the video says "DNA: deoxyribose, base and phosphate per letter"
    And user hovers over monomer 4 of the "cell 2 of sequence" area of grid
    And user hovers over monomer 13 of the "cell 2 of sequence" area of grid
    And user hovers over monomer 25 of the "cell 2 of sequence" area of grid
    And the video says "RNA: ribose instead of deoxyribose, uracil for thymine"
    Given user opens FASTA_RNA dataset keeping the first 12 rows
    And the toolbox pane is hidden
    And the context panel is closed
    And the monomers of the "sequence" column are drawn at 14 px
    And the "sequence" column is converted to atomic level drawn straight with its monomers highlighted
    Then "molfile(sequence)" column should have no missing values
    Given the grid shows only the columns "molfile(sequence), sequence"
    And the "sequence" column of the grid is 390 px wide
    And the "molfile(sequence)" column of the grid is 1110 px wide
    And the rows of the grid are 150 px tall
    And no row of the table is current
    When user hovers over monomer 7 of the "cell 2 of sequence" area of grid
    And user hovers over monomer 16 of the "cell 2 of sequence" area of grid
    Then no error or warning balloon should have been shown
