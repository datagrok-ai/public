@journey @realizes:bio.analyze.sequence-space
Feature: Sequence Space
  Bio | Analyze | Sequence Space... reduces a sequence column to two embedding columns and docks a
  scatter plot over them; run with the defaults, then again with another method and metric — the
  second run must be a second result, computed with what was edited, clustered again.

  The last scenario translates TestTrack Bio/bio-lifecycle-macromolecule-column.md scenario 3
  (GROK-19928): the table with both results is saved through the ribbon's Save project dialog and
  reopened, and the embedding columns, both scatter plots and the sequence column's type and units
  come back. The table opened by the feature has no creation script, so the dialog saves it as a
  snapshot; the md's Data sync switch is not claimed. The HELM and MSA runs of the manual cases are
  tested in Bio src/tests/sequence-space-test.ts ('sequenceSpaceHelm', 'sequenceSpaceMsa').

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset
    And the Bio package is initialized
    Then "fasta" column should have semantic type "Macromolecule"
    And "fasta" column should have units "fasta"

  Scenario: The editor opens on the sequence column with the default engine
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    And editor of Column input in "Sequence Space" dialog should have text "fasta"
    And Method input in "Sequence Space" dialog should have value "UMAP"
    And Similarity input in "Sequence Space" dialog should have value "Hamming"
    And "Plot embeddings" checkbox in "Sequence Space" dialog should be checked
    And "Cluster embeddings" checkbox in "Sequence Space" dialog should be checked

  Scenario: Running with the defaults appends the embeddings and docks the scatter plot
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And "Sequence Space" dialog should be hidden
    And a new column "Embed_X_1" should have been added
    And a new column "Embed_Y_1" should have been added
    And a new column matching "^Cluster \(DBSCAN\)" should have been added
    And "Embed_X_1" column should have no missing values
    And scatter plot viewer should be visible
    And title of scatter plot viewer should have text "Sequence space"
    And "X" property of scatter plot viewer should be "Embed_X_1"
    And "Y" property of scatter plot viewer should be "Embed_Y_1"
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The editor reopens over the first result and takes another method and metric
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    When user selects "t-SNE" in Method input in "Sequence Space" dialog
    And user selects "Levenshtein" in Similarity input in "Sequence Space" dialog
    Then Method input in "Sequence Space" dialog should have value "t-SNE"
    And Similarity input in "Sequence Space" dialog should have value "Levenshtein"

  Scenario: The edited run is a second result computed with the edited settings
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_2" should have been added
    And a new column "Embed_Y_2" should have been added
    And a new column matching "^Cluster \(DBSCAN\)" should have been added
    And second scatter plot viewer should be visible
    And "X" property of second scatter plot viewer should be "Embed_X_2"
    And "Description" property of second scatter plot viewer should contain "method: t-SNE"
    And "Description" property of second scatter plot viewer should contain "similarity: Levenshtein"
    And second scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Both results come back from a saved project (GROK-19928)
    Given no project named "bdd-bio-seqspace-{run}" is on the server
    When user opens the Save project dialog from the ribbon
    And user types "bdd-bio-seqspace-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "bdd-bio-seqspace-{run}" should be on the server
    And the "filter_FASTA" table of the "bdd-bio-seqspace-{run}" project should be saved as a snapshot
    When user closes all views
    And user opens the "bdd-bio-seqspace-{run}" project and waits for its table
    Then "fasta" column should have semantic type "Macromolecule"
    And "fasta" column should have units "fasta"
    And the table should have a column "Embed_X_1"
    And the table should have a column "Embed_Y_2"
    And scatter plot viewer should be visible
    And "X" property of scatter plot viewer should be "Embed_X_1"
    And "Y" property of scatter plot viewer should be "Embed_Y_1"
    And "X" property of second scatter plot viewer should be "Embed_X_2"
    And "Description" property of second scatter plot viewer should contain "method: t-SNE"
    And scatter plot viewer should be painted
    And second scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged
