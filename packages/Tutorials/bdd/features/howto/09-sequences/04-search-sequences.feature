@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The three searches of Bio's Search menu on FASTA_PT_activity
  (99 peptides of 16 residues in three families, the cluster column). Bio > Search > Subsequence Search puts a filter
  on the sequence column into the filter panel; the motif NEMTK typed there keeps the ten peptides that hold it, all
  of one family, and the status bar counts them; clearing the box brings all 99 back. Bio > Search > Similarity Search
  adds a viewer listing the nearest neighbours of the current peptide, drawn with what differs from it in color, and
  it follows the current row; Bio > Search > Diversity Search lists the ten peptides most unlike each other, from all
  three families (picked from a random start, so the rows vary from run to run), and a click on each of the first
  two makes it the grid's current row. Both viewers stay where the commands dock them, under the grid. Off camera:
  the family column colored, the sequences drawn large. Nothing is made on the server. Merged from the former clips
  "Find sequences by a motif" and "Find similar and diverse sequences".

  Scenario: Search sequences by motif, similarity and diversity
    Find the peptides that hold a motif, the nearest neighbours of one peptide, and a set as different as possible.

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
    And the "sequence" column of the grid is 420 px wide
    When the video says "99 peptides; the first column, colored, is each one's family"
    # caption: Bio > Search > Subsequence Search: find a motif
    And user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then the filter panel should have a filter on "sequence" column
    # caption: Type a motif into the sequence filter
    When user enters "NEMTK" into "Substructure" input in filters viewer
    Then 10 rows should pass the filter
    And the filter should pass exactly the rows where "sequence" contains "NEMTK"
    And no rows where "cluster" is "0" should pass the filter
    And no rows where "cluster" is "2" should pass the filter
    When the video says "Ten peptides hold the motif, all of one family"
    And user hovers over "Filtered: 10" text in status bar
    # caption: Clear the motif: all 99 are back
    And user clears "Substructure" input in filters viewer
    Then all rows should pass the filter
    # caption: Click a peptide
    When user clicks on the "cell 5 of sequence" area of grid
    Then row 5 should be current
    # caption: Bio > Search > Similarity Search: its nearest neighbours
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then "Sequence Similarity Search" viewer should be visible
    And the "target row" reading of "Sequence Similarity Search" viewer should be 4
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 11
    When the video says "Most similar first; what differs stays in color"
    # caption: Click another peptide: the list follows
    And user clicks on the "cell 10 of sequence" area of grid
    Then the "target row" reading of "Sequence Similarity Search" viewer should be 9
    # caption: Bio > Search > Diversity Search: the most different peptides
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then "Sequence Diversity Search" viewer should be visible
    And the "subset size" reading of "Sequence Diversity Search" viewer should be 10
    And the peptides "Sequence Diversity Search" viewer lists should come from 3 different "cluster" values
    When the video says "Ten peptides unlike each other, from all three families"
    # caption: Click one to find it in the table
    And user clicks on the "row header 1" area of third grid viewer
    And user clicks on the "row header 2" area of third grid viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
