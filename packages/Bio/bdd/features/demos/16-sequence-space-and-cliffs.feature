@guide @demo @viewport:1920x1080
Feature: Bioinformatics · Sequence space and activity cliffs
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity): Bio > Analyze > Sequence Space embeds the peptides by their
  pairwise distance (UMAP over Hamming) and clusters the embedding, and the scatter plot it docks
  shows the three families apart, colored by cluster; Bio > Analyze > Activity Cliffs embeds them
  again and links every pair at least 80% alike, the peptides colored by activity and sized by how
  steep a cliff they sit on (SALI); its cliffs button lists the pairs, the residue that differs
  marked, and a click on a pair zooms the plot to it, where the pointer on the link shows both
  peptides and their activities. The cliffs plot is set off camera not to filter the table by its
  zoom (the scatter plot's default): the zoom on the pair would otherwise filter the table to the
  five peptides in view and empty the sequence space plot.

  Left out: the dialogs' advanced settings (method, metric and their options) and the Show only
  cliffs switch.

  Scenario: Map the sequence space and find the activity cliffs
    Peptides embedded by sequence similarity fall into their families; near-identical pairs far apart in activity show up as cliffs.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 260 px wide
    When the video says "Sequence space: alike peptides land close together"
    And user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And scatter plot viewer should be visible
    And scatter plot viewer should be painted
    Given the grid shows only the columns "sequence, activity"
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    And editor of Activities input in "Sequence Activity Cliffs" dialog should have text "activity"
    When user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And second scatter plot viewer should be visible
    And the "cliffs" reading of second scatter plot viewer should be 3
    And second scatter plot viewer should be painted
    Given the grid shows only the columns "sequence, activity"
    And user sets "zoomAndFilter" property of second scatter plot viewer to "no action"
    When user clicks on "cliffs" button
    Then second grid viewer should be visible
    When user clicks on the first "cell" area of second grid viewer
    Then row 70 should be current
    And all rows should pass the filter
    And second scatter plot viewer should be painted
    When the video says "A cliff: one residue apart, far apart in activity"
    And user hovers over the "view" area of second scatter plot viewer
    Then exactly one tooltip should be shown
    And no error or warning balloon should have been shown
    And no errors should have been logged
