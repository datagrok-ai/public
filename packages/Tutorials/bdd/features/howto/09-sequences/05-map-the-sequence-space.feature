@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). On FASTA_PT_activity (99 peptides of 16 residues in three
  families, with an activity), two maps of the same peptides. Bio > Analyze > Sequence Space, with the similarity set
  to Monomer chemical distance (how alike the monomers are as molecules, not only whether they are the same), embeds
  the peptides in two dimensions (UMAP, in the browser) and clusters the embedding; the scatter plot it docks shows
  the three families apart, colored by cluster. The pointer on a dot shows its peptide; a click makes it the grid's
  current row. Then Bio > Analyze > Activity Cliffs, with the same similarity, embeds them again and links every pair
  at least 80% alike, the peptides colored by activity and sized by how steep a cliff they sit on (SALI); its cliffs
  button lists the 21 pairs, the residue that differs marked, and a click on a pair zooms the plot to it and makes its
  first peptide the grid's current row: 3.84 over -2.81. Off camera: before the first step the grid shows only the
  sequence and the activity; after the cliffs plot appears, that plot is set not to filter the table by its zoom
  (nothing on screen changes). Nothing is made on the server. Merged from the former clips "Map the sequence space"
  and "Find activity cliffs".

  Scenario: Map the sequence space and find the activity cliffs
    Place similar sequences together on a map, then find the pairs that are nearly identical but far apart in activity.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 260 px wide
    # caption: Bio > Analyze > Sequence Space
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    # caption: Similarity: how two sequences are compared
    When user selects "Monomer chemical distance" in Similarity input in "Sequence Space" dialog
    Then Similarity input in "Sequence Space" dialog should have value "Monomer chemical distance"
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And scatter plot viewer should be visible
    And scatter plot viewer should be painted
    When the video says "Each dot is a peptide; alike peptides land together, colored by cluster"
    # caption: The pointer on a dot shows its peptide
    And user hovers over the "marker of row 5" area of scatter plot viewer
    # caption: Click a dot to find its peptide in the table
    And user clicks on the "marker of row 15" area of scatter plot viewer
    Then row 15 should be current
    When the video says "Activity cliffs: alike peptides whose activity differs a lot"
    # caption: Bio > Analyze > Activity Cliffs
    And user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    And editor of Activities input in "Sequence Activity Cliffs" dialog should have text "activity"
    # caption: The sequences, their activity, and how to compare them
    When user selects "Monomer chemical distance" in Similarity input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And second scatter plot viewer should be visible
    And the "cliffs" reading of second scatter plot viewer should be 21
    And second scatter plot viewer should be painted
    Given user sets "zoomAndFilter" property of second scatter plot viewer to "no action"
    When the video says "Colored by activity; alike peptides linked, bigger dots on steeper cliffs"
    # caption: This button lists the cliffs: 21 pairs
    And user clicks on "cliffs" button
    Then second grid viewer should be visible
    # caption: Click a pair: the plot zooms to it
    When user clicks on the first "cell" area of second grid viewer
    Then row 70 should be current
    And all rows should pass the filter
    When the video says "A cliff: one residue apart, far apart in activity"
    And user hovers over the "cell 71 of activity" area of grid
    Then no error or warning balloon should have been shown
    And no errors should have been logged
