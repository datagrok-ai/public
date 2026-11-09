@guide @demo
Feature: Bioinformatics · Alignment and differences
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio > Analyze > MSA aligns the peptides of the FASTA sample in the
  browser (Kalign); the aligned column gets a header with every position's conservation and WebLogo, whose
  tooltip shows the position's composition, and a click on a sequence makes it the reference the other
  rows are drawn against: what they share with it fades, what differs stays in color.

  Scenario: Align, then compare every sequence with the one you click
    Kalign aligns the peptides in the browser; the header shows each position's conservation, and the rows show how they differ from the sequence you click.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_sample dataset
    Then "Sequence" column should have semantic type "Macromolecule"
    Given the grid shows only the columns "Entry, Sequence, Activity"
    And the "Sequence" column of the grid is 330 px wide
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    When user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(Sequence)" should have been added
    And "msa(Sequence)" column should have tag "aligned" equal to "SEQ.MSA"
    Given the grid shows only the columns "Entry, msa(Sequence), Activity"
    And the "msa(Sequence)" column of the grid is 820 px wide
    And the monomers of the "msa(Sequence)" column are drawn at 16 px
    And the column headers of the grid are 160 px tall
    When the video says "Every position: its conservation and its composition"
    And user hovers over the "position 20 of msa(Sequence) header" area of grid
    And user clicks on the "cell 7 of msa(Sequence)" area of grid
    And user hovers over the "cell 7 of Entry" area of grid
    And the video says "Shared monomers fade: what differs from the clicked row stands out"
    And user clicks on the "cell 13 of msa(Sequence)" area of grid
    And user hovers over the "cell 13 of Entry" area of grid
    Then no errors should have been logged
