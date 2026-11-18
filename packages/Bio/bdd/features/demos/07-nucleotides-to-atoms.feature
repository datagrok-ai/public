@guide @demo
Feature: 2 · Every sequence is a molecule · Nucleotides to atoms
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio | Transform | To Atomic Level on Bio's FASTA RNA and DNA samples, each
  strand cut off camera to its first 15 nucleotides (Bio's Extract Region), so a nucleotide is large: each letter
  becomes a whole nucleotide (sugar, base and phosphate: ribose for RNA, deoxyribose for DNA), the strand drawn
  straight (Non-linear off) and every nucleotide colored like its letter (Highlight monomers on); the pointer on an
  RNA letter lights that nucleotide's atoms. The DNA strands are converted off camera, under the closing note, and
  not hovered: a DNA letter's tooltip names its polymer type RNA. Set up off camera: the region, the molecule column
  left of the sequences, so the letter tooltips never cover it, the letters drawn larger, and no current row.

  Scenario: RNA and DNA, nucleotide by nucleotide
    Each letter of an RNA or DNA strand becomes a whole nucleotide.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_RNA dataset keeping the first 8 rows
    And the region 1 to 15 of the "sequence" column is extracted as "strand"
    And user removes "sequence" column
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the monomers of the "strand" column are drawn at 16 px
    And the "strand" column of the grid is 260 px wide
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    And user unchecks "Non-linear" checkbox in "To Atomic Level" dialog
    And user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column "molfile(strand)" should have been added
    And "molfile(strand)" column should have no missing values
    Given the grid shows only the columns "molfile(strand), strand"
    And the "strand" column of the grid is 230 px wide
    And the "molfile(strand)" column of the grid is 1220 px wide
    And the rows of the grid are 200 px tall
    And no row of the table is current
    When the video says "RNA: ribose, base and phosphate per letter"
    And user hovers over monomer 3 of the "cell 2 of strand" area of grid
    And user hovers over monomer 7 of the "cell 2 of strand" area of grid
    Given user opens FASTA_DNA dataset keeping the first 8 rows
    And the region 1 to 15 of the "sequence" column is extracted as "strand"
    And user removes "sequence" column
    And the toolbox pane is hidden
    And the context panel is closed
    And the monomers of the "strand" column are drawn at 16 px
    And the "strand" column is converted to atomic level drawn straight with its monomers highlighted
    Then "molfile(strand)" column should have no missing values
    Given the grid shows only the columns "molfile(strand), strand"
    And the "strand" column of the grid is 230 px wide
    And the "molfile(strand)" column of the grid is 1220 px wide
    And the rows of the grid are 200 px tall
    And no row of the table is current
    # caption: DNA: deoxyribose instead of ribose, thymine for uracil
    Then "strand" column should have tag "alphabet" equal to "DNA"
