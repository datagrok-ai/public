@guide @demo @viewport:1440x810
Feature: Bioinformatics · Search
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity): Bio > Search > Subsequence Search puts a filter on the sequence
  column into the filter panel, and a motif typed there keeps the ten peptides that hold it; a click on
  a peptide makes the grid compare every other row with it, fading the residues they share; Bio >
  Search > Similarity Search docks a viewer listing the nearest neighbours of the current peptide in
  the whole table, drawn the same way, and it follows the current row; Bio > Search > Diversity Search
  docks the peptides most unlike each other, which come from all three families.

  Left out: the HELM filter (a query drawn in the HELM editor takes too long for a clip) and the viewers'
  distance settings.

  Scenario: Find peptides by motif and by similarity
    A motif typed into the filter keeps the peptides that hold it; the similarity viewer lists the nearest neighbours of the current one, the diversity viewer the most different ones.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 300 px wide
    When user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then the filter panel should have a filter on "sequence" column
    When user enters "NEMTK" into "Substructure" input in filters viewer
    Then 10 rows should pass the filter
    And the filter should pass exactly the rows where "sequence" contains "NEMTK"
    When the video says "Only the peptides that contain the motif remain"
    And user clicks on the "cell 42 of sequence" area of grid
    Then row 42 should be current
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then "Sequence Similarity Search" viewer should be visible
    And the "target row" reading of "Sequence Similarity Search" viewer should be 41
    Given user sets "limit" property of "Sequence Similarity Search" viewer to "6"
    Then the "neighbours" reading of "Sequence Similarity Search" viewer should be 7
    Given the "similar (sequence)" column of the "Sequence Similarity Search" viewer is 450 px wide
    When the video says "Its nearest neighbours: residues that differ stay in color"
    And user clicks on the "cell 55 of sequence" area of grid
    Then the "target row" reading of "Sequence Similarity Search" viewer should be 54
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 7
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then "Sequence Diversity Search" viewer should be visible
    And the subset of "Sequence Diversity Search" viewer should span 3 values of the "cluster" column
    When the video says "The most different peptides: all three families"
    Then no error or warning balloon should have been shown
    And no errors should have been logged
