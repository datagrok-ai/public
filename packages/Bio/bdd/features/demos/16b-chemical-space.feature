@guide @demo
Feature: 4 · A peptide series · Chemical space
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). The peptides of FASTA_PT_activity (99, of 16 residues, in three families)
  are a sequence space already (Bio > Analyze > Sequence Space, run off camera with its defaults) and
  molecules (To Atomic Level, off camera: a V3000 molfile per peptide; the grid keeps only the sequences, so
  the plots get the width). Chem > Analyze > Chemical Space embeds the molecules by their fingerprints (UMAP
  over Tanimoto), as for any small molecules; off camera its plot moves beside the sequence space, takes the
  sequence space's cluster colors (purple, teal and blue, so the orange of a selection stands out) and draws
  no structure labels. A selection box around one family in the sequence space (the markers of its first 20
  rows) lights the same peptides up in the chemical space, where they keep together too: the families hold by
  chemistry as well. The one note comes once the chemical space is drawn, over the selection it announces. Chemistry runs in the browser
  (Chem's engine on this stand is Crux).

  Left out: the Chem Space dialog's settings (method, similarity, clustering of the embedding).

  Scenario: The same peptides in chemical space
    Molecules made from the sequences, mapped by chemistry beside their sequence space.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the browse panel is closed
    And the context panel is closed
    And the grid shows only the columns "sequence"
    And the "sequence" column of the grid is 250 px wide
    And user picks "Bio > Analyze > Sequence Space..." from the top menu
    And user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And scatter plot viewer should be painted
    Given the "sequence" column is converted to atomic level drawn straight with its monomers highlighted
    And the grid shows only the columns "sequence"
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    And Column input in "Chem Space" dialog should contain text "molfile(sequence)"
    When user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And second scatter plot viewer should be visible
    Given second scatter plot viewer is docked to the right of scatter plot viewer
    And the grid shows only the columns "sequence"
    And the rows of the grid are 25 px tall
    And user colors "Cluster (DBSCAN)" column categorically:
      | 1 | #8e6cc0 |
      | 2 | #2a9d8f |
      | 3 | #4575b4 |
    And user sets "colorColumnName" property of second scatter plot viewer to "Cluster (DBSCAN)"
    And user sets "displayLabels" property of second scatter plot viewer to "Never"
    Then second scatter plot viewer should be painted
    When the video says "Selected by sequence, the same peptides light up in chemical space"
    And user drags a selection box around the markers of rows 1 to 20 of scatter plot viewer
    Then some rows should be selected
    And no rows where "cluster" is "1" should be selected
    And no rows where "cluster" is "2" should be selected
    And second scatter plot viewer should show a selection highlight
    And no error or warning balloon should have been shown
    And no errors should have been logged
