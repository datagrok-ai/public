@guide @demo
Feature: 4 · A peptide series · Alignment
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). FASTA_PT is one peptide family: 100 variants, 31 to 38 residues
  long, that share two motifs (MTMFC near the start, EAM and a K near the end) at shifting offsets. Bio >
  Analyze > MSA aligns them in the browser (Kalign): the motifs fall into columns, and the aligned column gets
  a header with every position's conservation and WebLogo, whose tooltip shows the position's composition.
  A click on a peptide makes it the reference the other rows are drawn against: what they share with it
  fades, what differs stays in color; that note comes before the first click, so the fading is explained as it
  appears.

  Left out: the dialog's alignment parameters and its Clusters input (each cluster aligned on its own), and
  the UniProt sample (FASTA.csv), whose Cluster column groups unrelated fragments: aligned, they share too
  little for the fading to show anything.

  Scenario: Align a peptide family, then compare it with one peptide
    Kalign aligns a peptide family in the browser; one click shows what differs.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_PT dataset
    Then "sequence" column should have semantic type "Macromolecule"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "id, sequence"
    And the monomers of the "sequence" column are drawn at 16 px
    And the "sequence" column of the grid is 620 px wide
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    When user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(sequence)" should have been added
    And "msa(sequence)" column should have tag "aligned" equal to "SEQ.MSA"
    And every value of "msa(sequence)" column should have the same length
    Given the grid shows only the columns "id, msa(sequence)"
    And the monomers of the "msa(sequence)" column are drawn at 16 px
    And the "msa(sequence)" column of the grid is 820 px wide
    And the column headers of the grid are 160 px tall
    When the video says "Aligned: the shared motifs fall into columns"
    And user hovers over the "position 8 of msa(sequence) header" area of grid
    Then tooltip should contain text "Position: 8"
    When user hovers over the "position 34 of msa(sequence) header" area of grid
    Then tooltip should contain text "Position: 34"
    When the video says "What differs from the clicked peptide stays in color"
    And user clicks on the "cell 7 of msa(sequence)" area of grid
    And user hovers over the "cell 7 of id" area of grid
    And user clicks on the "cell 13 of msa(sequence)" area of grid
    And user hovers over the "cell 13 of id" area of grid
    Then no errors should have been logged
