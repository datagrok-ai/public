@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Bio's FASTA_PT sample is one peptide family: 100 variants, 31 to
  38 residues long, that share two motifs at shifting offsets. Bio > Analyze > MSA aligns them in the browser
  (Kalign) into a new column; the shared motifs fall into columns, and the aligned column's header shows each
  position's conservation and logo, whose tooltip gives the position's composition. A click on a peptide makes it the
  reference: what other rows share with it fades, what differs stays in color. Then Bio > Transform > Split to
  Monomers, on the aligned column, gives every position a column of its own (44 here), each residue drawn in its
  color, and the grid scrolls to them; position 8 is M in every peptide, the pointer on one shows the monomer, and a
  double-click on position 14's header sorts the peptides by the residue there. Off camera, before the first step: the columns chosen and sized, the headers made tall enough for the logo.
  Nothing is made on the server. Merged from the former clips "Align sequences" and "Split aligned sequences into
  positions" (which used Bio's MSA sample).

  Scenario: Align sequences and split them into positions
    MSA lines up a family of sequences in your browser; Split to Monomers then gives every position a column of its own.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_PT dataset
    Then "sequence" column should have semantic type "Macromolecule"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "id, sequence"
    And the "sequence" column of the grid stays 520 px wide
    And the column headers of the grid are 160 px tall
    When the video says "Related peptides, not lined up: a shared motif like EAM.K sits at different places"
    # caption: Bio > Analyze > MSA aligns them
    And user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    # caption: It runs in your browser
    When user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(sequence)" should have been added
    And "msa(sequence)" column should have tag "aligned" equal to "SEQ.MSA"
    And every value of "msa(sequence)" column should have the same length
    When the video says "Aligned: the shared motifs fall into columns" at the bottom
    And the video says "The header: how conserved each position is, and its letters" at the bottom
    And user hovers over the "position 8 of msa(sequence) header" area of grid
    Then tooltip should contain text "Position: 8"
    When user hovers over the "position 3 of msa(sequence) header" area of grid
    Then tooltip should contain text "Position: 3"
    When the video says "Click a peptide: what differs from it stays in color" at the bottom
    And user clicks on the "cell 7 of msa(sequence)" area of grid
    And user hovers over the "cell 7 of id" area of grid
    And user clicks on the "cell 13 of msa(sequence)" area of grid
    And user hovers over the "cell 13 of id" area of grid
    # caption: Bio > Transform > Split to Monomers
    And user picks "Bio > Transform > Split to Monomers..." from the top menu
    Then "Split to Monomers" dialog should be visible
    # caption: Pick the aligned column
    When user selects "msa(sequence)" in Sequence input in "Split to Monomers" dialog
    And user clicks on OK button in "Split to Monomers" dialog
    Then a new column "1" should have been added
    And the "cell type of 10" reading of grid should be "Monomer"
    And every value of "8" column should match "^M$"
    When the video says "Every position now has a column of its own" at the bottom
    # caption: Position 8 holds M in every peptide
    And user hovers over the "cell 3 of 8" area of grid
    Then tooltip should be visible
    # caption: Double-click a position's header to sort by it
    When user double-clicks on the "header 14" area of grid
    Then the "sort column" reading of grid should be "14"
    When the video says "Peptides with the same residue at position 14 now sit together" at the bottom
    Then no error or warning balloon should have been shown
    And no errors should have been logged
