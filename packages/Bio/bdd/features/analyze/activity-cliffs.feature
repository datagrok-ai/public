@journey @realizes:bio.analyze.activity-cliffs
Feature: Sequence Activity Cliffs
  Bio | Analyze | Activity Cliffs... embeds the sequences, scores the pairs that are close in
  sequence and far in activity, and docks a scatter plot with the cliffs drawn over it. The
  fixture carries an Activity column: the editor needs a numeric column to fill Activities.

  Background:
    Given user is logged in
    And user opens FASTA_sample dataset
    And the Bio package is initialized
    Then "Sequence" column should have semantic type "Macromolecule"

  Scenario: The editor opens with the sequence column, an activity and a cutoff
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    And editor of Column input in "Sequence Activity Cliffs" dialog should have text "Sequence"
    And Method input in "Sequence Activity Cliffs" dialog should have value "UMAP"
    And Similarity input in "Sequence Activity Cliffs" dialog should have value "Hamming"
    And "Similarity cutoff" input in "Sequence Activity Cliffs" dialog should have value "80"

  Scenario: Running on the Activity column docks a cliff scatter plot
    When user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And "Sequence Activity Cliffs" dialog should be hidden
    And a new column "Embed_X_1" should have been added
    And a new column "Embed_Y_1" should have been added
    And a new column matching "sali|SALI" should have been added
    And scatter plot viewer should be visible
    And title of scatter plot viewer should have text "Activity cliffs"
    And "X" property of scatter plot viewer should be "Embed_X_1"
    And "Y" property of scatter plot viewer should be "Embed_Y_1"
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A second analysis with other settings docks alongside the first
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    When user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user selects "t-SNE" in Method input in "Sequence Activity Cliffs" dialog
    And user selects "Levenshtein" in Similarity input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_2" should have been added
    And second scatter plot viewer should be visible
    And "X" property of second scatter plot viewer should be "Embed_X_2"
    And "Description" property of second scatter plot viewer should contain "method: t-SNE"
    And "Description" property of second scatter plot viewer should contain "similarity: Levenshtein"
    And second scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged
