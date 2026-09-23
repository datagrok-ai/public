@journey @realizes:chem.cp.calculate-clustering
Feature: BitBIRCH clustering, Cluster MCS and the similarity matrix
  On smiles, BitBIRCH Clustering with its defaults appends a cluster column that groups the
  molecules: more than one cluster, fewer clusters than rows. Butina Cluster, the Python script,
  groups the same fifty molecules its own way, and the two partitions do not agree row for row. On
  spgi-100, Cluster MCS over the BitBIRCH clusters appends a Molecule column with a structure on every
  row, fewer structures than rows. Over the first 50 molecules, Similarity Matrix builds a table with
  one similarity column per molecule, 1 down the diagonal and something lower off it.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: BitBIRCH clustering groups the molecules
    When user picks "Chem > Calculate > BitBIRCH Clustering..." from the top menu
    Then "BitBIRCH Clustering" dialog should be visible
    And Molecules input in "BitBIRCH Clustering" dialog should contain text "canonical_smiles"
    And Threshold input in "BitBIRCH Clustering" dialog should have value "0.55"
    And "Fingerprint type" input in "BitBIRCH Clustering" dialog should have value "Morgan"
    When user clicks on OK button in "BitBIRCH Clustering" dialog
    Then the top menu command should have completed
    And a new column "Cluster (BitBIRCH)" should have been added
    And "Cluster (BitBIRCH)" column should have no missing values
    And "Cluster (BitBIRCH)" column should have at least 2 distinct values
    And "Cluster (BitBIRCH)" column should have fewer distinct values than the table has rows
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: The similarity matrix reads 1 down its diagonal
    Given user opens smiles dataset keeping the first 50 rows as "clustering_matrix_subset"
    When user picks "Chem > Calculate > Similarity Matrix..." from the top menu
    Then "Similarity Matrix" dialog should be visible
    And Table input in "Similarity Matrix" dialog should have value "clustering_matrix_subset"
    And Molecules input in "Similarity Matrix" dialog should contain text "canonical_smiles"
    And Symbols input in "Similarity Matrix" dialog should contain text "molregno"
    When user clicks on OK button in "Similarity Matrix" dialog
    Then the top menu command should have completed
    And table "canonical_smiles similarity matrix" should be open
    And the similarity columns of table "canonical_smiles similarity matrix" should be symmetric, read 1 on the diagonal and less somewhere off it
    And no errors should have been logged

  Scenario: Cluster MCS over the BitBIRCH clusters writes a structure on every row
    Given user opens spgi dataset
    When user picks "Chem > Calculate > BitBIRCH Clustering..." from the top menu
    And user clicks on OK button in "BitBIRCH Clustering" dialog
    Then the top menu command should have completed
    And a new column "Cluster (BitBIRCH)" should have been added
    And "Cluster (BitBIRCH)" column should have at least 2 distinct values
    And "Cluster (BitBIRCH)" column should have fewer distinct values than the table has rows
    When user picks "Chem > Calculate > Cluster MCS..." from the top menu
    And user clicks on OK button in "Cluster MCS" dialog
    Then the top menu command should have completed
    And a new column matching "MCS|mcs|Scaffold" should have been added
    And the newest column matching "MCS|mcs|Scaffold" should have no missing values
    And "Cluster MCS" column should have fewer distinct values than the table has rows
    And the table should have 100 rows
    And no errors should have been logged

  Scenario: Butina clustering groups the same molecules its own way
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Butina Cluster..." from the top menu
    Then "Butina Molecules Clustering" dialog should be visible
    And Molecules input in "Butina Molecules Clustering" dialog should contain text "canonical_smiles"
    And "Distance Cutoff" input in "Butina Molecules Clustering" dialog should have value "0.4"
    When user clicks on OK button in "Butina Molecules Clustering" dialog
    Then the top menu command should have completed
    And a new column "cluster (Butina)" should have been added
    And "cluster (Butina)" column should have no missing values
    And "cluster (Butina)" column should have at least 2 distinct values
    And "cluster (Butina)" column should have fewer distinct values than the table has rows
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: BitBIRCH and Butina split the same table differently
    When user picks "Chem > Calculate > BitBIRCH Clustering..." from the top menu
    And user clicks on OK button in "BitBIRCH Clustering" dialog
    Then the top menu command should have completed
    And a new column "Cluster (BitBIRCH)" should have been added
    And "Cluster (BitBIRCH)" column should have at least 2 distinct values
    And "Cluster (BitBIRCH)" column should have fewer distinct values than the table has rows
    And some value of "Cluster (BitBIRCH)" column should differ from "cluster (Butina)" column in the same row
    And no errors should have been logged
