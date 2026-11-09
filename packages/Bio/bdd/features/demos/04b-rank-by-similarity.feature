@guide @demo
Feature: Bioinformatics · Similarity to a lead
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio > Calculate > Similarity scores the 99 peptides of
  FASTA_PT_activity against a lead typed in (c0_seq95, activity 4.93); sorted by the score, the lead comes
  first and its closest analog second, and a click on the lead draws every other row against it: c0_seq96
  differs at one position, and its activity is 1.19.

  Scenario: Rank the peptides by their similarity to a lead
    Every peptide scored against a lead and sorted by the score; one click shows where the closest analogs differ from it.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens FASTA_PT_activity dataset
    Then "sequence" column should have semantic type "Macromolecule"
    Given the grid shows only the columns "sequence_id, sequence, activity"
    And the "sequence" column of the grid is 300 px wide
    And the monomers of the "sequence" column are drawn at 16 px
    When user picks "Bio > Calculate > Similarity..." from the top menu
    Then Similarity dialog should be visible
    When user enters "DLKMMHQNTPGNMAQD" into Reference input in Similarity dialog
    And user clicks on OK button in Similarity dialog
    Then a new column "Similarity" should have been added
    And "Similarity" column should have its maximum in row 7
    When user double-clicks on the "header Similarity" area of grid
    Given user scrolls the grid to the "sequence_id" column
    When the video says "The lead first, then its closest analogs"
    And user clicks on the "cell 7 of sequence" area of grid
    And user hovers over the "cell 7 of sequence_id" area of grid
    And the video says "One substitution apart: activity 4.93 against 1.19"
    Then no errors should have been logged
