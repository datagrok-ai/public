@eda
@journey
Feature: ML > Cluster runs DBSCAN on two columns and Markov clustering on several
  ML | Cluster | DBSCAN... opens the DBSCAN dialog on the current table; with the X and Y columns picked, OK adds a
  "Cluster (DBSCAN)" column with one label per dense group of points and another for the points in no group.
  ML | Cluster | MCL... opens the Markov clustering dialog; OK docks the MCL viewer and adds the cluster, the
  embedding, the cluster size and the connectivity of every row as columns.

  DBSCAN runs on a nine-row table written here: two tight groups of four points and one point far from both. MCL
  runs on two number columns of iris: on number columns only, the dialog used to hand the clustering no column, and
  it never finished (public#4473). Left out on purpose:
  - which row goes to which cluster: EDA package test 'DBSCAN', and Peptides' 'Model: MCL clusters';
  - the dialogs' defaults and Use WebGPU (headless Chromium has no adapter): MISSING.md.
  Nothing is put on the server.

  Background:
    Given user is logged in

  Scenario: DBSCAN labels the two dense groups apart and the far point as noise
    Given user opens a table "points" with:
      | x    | y    |
      | 0    | 0    |
      | 0.01 | 0    |
      | 0    | 0.01 |
      | 0.01 | 0.01 |
      | 5    | 5    |
      | 5.01 | 5    |
      | 5    | 5.01 |
      | 5.01 | 5.01 |
      | 20   | 20   |
    When user picks "ML > Cluster > DBSCAN..." from the top menu
    Then "DBSCAN" dialog should be visible
    When user selects "x" in "X Col" input in "DBSCAN" dialog
    And user selects "y" in "Y Col" input in "DBSCAN" dialog
    And user clicks on OK button in "DBSCAN" dialog
    Then a new column "Cluster (DBSCAN)" should have been added
    And the newest column matching "^Cluster \(DBSCAN\)" should have 3 distinct values
    And no errors should have been logged

  @realizes:public#4473
  Scenario: MCL on two number columns docks its viewer and adds the cluster columns
    Given user opens iris dataset
    When user picks "ML > Cluster > MCL..." from the top menu
    Then "Markov clustering" dialog should be visible
    When user clicks on Columns input in "Markov clustering" dialog
    And user toggles the "Petal.Length" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user clicks on OK button in "Markov clustering" dialog
    Then 5 new columns should have been added
    And a new column "Cluster (MCL)" should have been added
    And the open tableview should have 1 MCL viewer
    And no errors should have been logged
