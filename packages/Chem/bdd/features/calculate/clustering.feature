@journey @realizes:chem.cp.calculate-clustering
Feature: BitBIRCH clustering, Cluster MCS and the similarity matrix
  On smiles, BitBIRCH Clustering with its defaults appends a cluster column. On spgi-100, Cluster MCS
  over the BitBIRCH clusters appends an MCS column. Over the first 50 molecules, Similarity Matrix
  opens a table of the molecules' similarities.

  The clusters, the matrix and the MCS of each cluster are tested in Chem src/tests/calculate.ts ('clustering': 'bitbirch', 'similarity matrix', 'cluster MCS').

  Not here: Butina Cluster, a server-side Python script — see the bdd library's CLAUDE.md,
  "What never becomes a feature".

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: BitBIRCH clustering adds a cluster column
    When user picks "Chem > Calculate > BitBIRCH Clustering..." from the top menu
    Then "BitBIRCH Clustering" dialog should be visible
    And Molecules input in "BitBIRCH Clustering" dialog should contain text "canonical_smiles"
    And Threshold input in "BitBIRCH Clustering" dialog should have value "0.55"
    And "Fingerprint type" input in "BitBIRCH Clustering" dialog should have value "Morgan"
    When user clicks on OK button in "BitBIRCH Clustering" dialog
    Then the top menu command should have completed
    And a new column "Cluster (BitBIRCH)" should have been added
    And no errors should have been logged

  Scenario: Similarity Matrix opens the matrix of the molecules
    Given user opens smiles dataset keeping the first 50 rows as "clustering_matrix_subset"
    When user picks "Chem > Calculate > Similarity Matrix..." from the top menu
    Then "Similarity Matrix" dialog should be visible
    And Table input in "Similarity Matrix" dialog should have value "clustering_matrix_subset"
    And Molecules input in "Similarity Matrix" dialog should contain text "canonical_smiles"
    And Symbols input in "Similarity Matrix" dialog should contain text "molregno"
    When user clicks on OK button in "Similarity Matrix" dialog
    Then the top menu command should have completed
    And table "canonical_smiles similarity matrix" should be open
    And no errors should have been logged

  Scenario: Cluster MCS over the BitBIRCH clusters adds an MCS column
    Given user opens spgi dataset
    When user picks "Chem > Calculate > BitBIRCH Clustering..." from the top menu
    And user clicks on OK button in "BitBIRCH Clustering" dialog
    Then the top menu command should have completed
    And a new column "Cluster (BitBIRCH)" should have been added
    When user picks "Chem > Calculate > Cluster MCS..." from the top menu
    And user clicks on OK button in "Cluster MCS" dialog
    Then the top menu command should have completed
    And a new column "Cluster MCS" should have been added
    And no errors should have been logged
