@guide @demo
Feature: 4 · A peptide series · Sequence space and activity cliffs
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity): Bio > Analyze > Sequence Space, with the distance set to Monomer
  chemical distance in the dialog (how alike the monomers are as molecules, from their fingerprints, not
  only whether they are the same), embeds the peptides (UMAP) and clusters the embedding, and the scatter
  plot it docks shows the three families apart, colored by cluster; Bio > Analyze > Activity Cliffs, with
  the same distance, embeds them again and links every pair at least 80% alike, the peptides colored by
  activity and sized by how steep a cliff they sit on (SALI); its cliffs button lists the pairs, the
  residue that differs marked, and a click on a pair zooms the plot to it and makes its first peptide the
  grid's current row; the clip ends with the pointer on the pair's activities in the grid, 3.84 over -2.81,
  and no tooltip over anything (the link's own tooltip opens at the pointer, its first words under the
  arrow, so it is left out). The note on the distance comes once the sequence space is drawn, and holds
  it on screen before the next menu. Off camera: the cliffs table, which the button docks along the
  whole bottom of the view, moves under the grid, beside the cliffs plot the pointer goes to next; the
  cliffs plot is set not to filter the table by its zoom (the scatter plot's default), which would filter
  the table to the peptides in view and empty the sequence space plot.

  Left out: the dialogs' advanced settings (method and its options) and the Show only cliffs switch.

  Scenario: Map the sequence space and find the activity cliffs
    Peptides mapped by monomer chemistry; near-identical pairs far apart in activity are cliffs.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 260 px wide
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    When user selects "Monomer chemical distance" in Similarity input in "Sequence Space" dialog
    Then Similarity input in "Sequence Space" dialog should have value "Monomer chemical distance"
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And scatter plot viewer should be visible
    And "Description" property of scatter plot viewer should contain "similarity: Monomer chemical distance"
    And scatter plot viewer should be painted
    Given the grid shows only the columns "sequence, activity"
    When the video says "Placed by monomer chemistry: alike peptides land together"
    And user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    And editor of Activities input in "Sequence Activity Cliffs" dialog should have text "activity"
    When user selects "Monomer chemical distance" in Similarity input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And second scatter plot viewer should be visible
    And the "cliffs" reading of second scatter plot viewer should be 21
    And second scatter plot viewer should be painted
    Given the grid shows only the columns "sequence, activity"
    And user sets "zoomAndFilter" property of second scatter plot viewer to "no action"
    When user clicks on "cliffs" button
    Then second grid viewer should be visible
    Given second grid viewer is docked under the grid, 35% of the view tall
    When user clicks on the first "cell" area of second grid viewer
    Then row 70 should be current
    And all rows should pass the filter
    And second scatter plot viewer should be painted
    When the video says "A cliff: one residue apart, far apart in activity"
    And user hovers over the "cell 71 of activity" area of grid
    Then no error or warning balloon should have been shown
    And no errors should have been logged
