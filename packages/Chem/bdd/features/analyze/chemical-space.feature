@journey @realizes:chem.cp.chemical-space
Feature: Chemical Space over SMILES, V2000 and V3000 molecules
  Chem | Analyze | Chemical Space... opens the Chem Space dialog on the molecule column with UMAP,
  Tanimoto, Plot embeddings and Cluster embeddings on and Cluster MCS off. OK adds a pair of Embed_X / Embed_Y
  columns, a Cluster (DBSCAN) column and a scatter plot of the embedding. A second run with
  t-SNE adds a second pair of embedding columns. The same holds for molecules read from a V2000 SDF
  and from a V3000 SDF.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed

  Scenario: The dialog opens on the molecule column of smiles-50
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    And Column input in "Chem Space" dialog should contain text "canonical_smiles"
    And Method input in "Chem Space" dialog should have value "UMAP"
    And Method input in "Chem Space" dialog should offer "UMAP, t-SNE"
    And Similarity input in "Chem Space" dialog should have value "Tanimoto"
    And Similarity input in "Chem Space" dialog should offer "Tanimoto, Asymmetric, Cosine, Sokal"
    And "Plot embeddings" input in "Chem Space" dialog should be checked
    And "Cluster MCS" input in "Chem Space" dialog should not be checked
    And no errors should have been logged

  Scenario: UMAP on SMILES adds the embedding columns and plots them
    When user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And a new column matching "^Embed_Y_" should have been added
    And a new column matching "^Cluster " should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And the newest column matching "^Embed_Y_" should have no missing values
    And the current view should hold at least 2 viewers
    And scatter plot viewer should be visible
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: A second run with t-SNE adds a second pair of embedding columns
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    And user selects "t-SNE" in Method input in "Chem Space" dialog
    And user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And 2 new columns matching "^Embed_[XY]_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: UMAP on V2000 molecules from an SDF
    Given user opens mol1K.sdf dataset
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    When user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And a new column matching "^Embed_Y_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And scatter plot viewer should be visible
    And no errors should have been logged

  Scenario: UMAP on V3000 molecules from an SDF
    Given user opens ApprovedDrugs2015 dataset
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    When user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And a new column matching "^Embed_Y_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And scatter plot viewer should be visible
    And no errors should have been logged
