@guide @demo
Feature: 4 · A peptide series · Search
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity) the grid shows the family of every peptide (the cluster column,
  colored off camera) and the sequences large (22 px residues, narrowed off camera when the first viewer
  docks beside the grid). Bio > Search > Subsequence Search puts a filter on the sequence column into the
  filter panel, and a motif typed there keeps the ten peptides that hold it, all of one family; a click on
  a peptide makes the grid compare every other row with it, fading the residues they share. Bio > Search >
  Similarity Search adds a viewer listing the nearest neighbours of the current peptide in the whole
  table, drawn the same way, and it follows the current row. With the motif cleared, Bio > Search >
  Diversity Search lists the peptides most unlike each other (three of them, set off camera), and a click
  on each makes it the current row of the grid: one from every family. The commands dock their viewers
  under the grid; off camera they move to its right, the diversity viewer under the similarity one.

  Left out: the HELM filter (a query drawn in the HELM editor takes too long for a clip) and the viewers'
  distance settings.

  Scenario: Find peptides by motif and by similarity
    A motif keeps one family; viewers list nearest neighbours and the most different peptides.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "cluster, sequence, activity"
    And user colors "cluster" column categorically:
      | 0 | #8e6cc0 |
      | 1 | #e9c46a |
      | 2 | #2a9d8f |
    And the "cluster" column of the grid is 60 px wide
    And the rows of the grid are 30 px tall
    And the monomers of the "sequence" column are drawn at 22 px
    And the "sequence" column of the grid is 500 px wide
    When user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then the filter panel should have a filter on "sequence" column
    When user enters "NEMTK" into "Substructure" input in filters viewer
    Then 10 rows should pass the filter
    And the filter should pass exactly the rows where "sequence" contains "NEMTK"
    And no rows where "cluster" is "0" should pass the filter
    And no rows where "cluster" is "2" should pass the filter
    When the video says "Only the peptides with the motif remain: one family"
    And user clicks on the "cell 42 of sequence" area of grid
    Then row 42 should be current
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then "Sequence Similarity Search" viewer should be visible
    And the "target row" reading of "Sequence Similarity Search" viewer should be 41
    Given "Sequence Similarity Search" viewer is docked to the right of the grid, 50% of the view wide
    And the "sequence" column of the grid is 380 px wide
    And user sets "limit" property of "Sequence Similarity Search" viewer to "6"
    Then the "neighbours" reading of "Sequence Similarity Search" viewer should be 7
    Given the "similar (sequence)" column of the "Sequence Similarity Search" viewer is 450 px wide
    When the video says "Its nearest neighbours: residues that differ stay in color"
    And user clicks on the "cell 55 of sequence" area of grid
    Then the "target row" reading of "Sequence Similarity Search" viewer should be 54
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 7
    When user clears "Substructure" input in filters viewer
    Then all rows should pass the filter
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then "Sequence Diversity Search" viewer should be visible
    Given "Sequence Diversity Search" viewer is docked under "Sequence Similarity Search" viewer
    And user sets "limit" property of "Sequence Diversity Search" viewer to "3"
    Then the "subset size" reading of "Sequence Diversity Search" viewer should be 3
    And the subset of "Sequence Diversity Search" viewer should span 3 values of the "cluster" column
    When the video says "The three most different peptides: one per family"
    And user clicks on the "row header 1" area of third grid viewer
    And user clicks on the "row header 2" area of third grid viewer
    And user clicks on the "row header 3" area of third grid viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
