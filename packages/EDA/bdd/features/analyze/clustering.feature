@eda
Feature: ML > Cluster runs DBSCAN on two columns and Markov clustering on several
  ML | Cluster | DBSCAN... opens the DBSCAN dialog on the current table: the table, the X and Y columns,
  Epsilon (0.02) and Minimum points (4). OK adds a "Cluster (DBSCAN)" column with one label per dense group
  of points and another for the points in no group. ML | Cluster | MCL... opens the Markov clustering
  dialog: the table, the columns, Similarity Threshold (80), Inflation Factor (2), Max Iterations (5), Min
  Cluster Size (5) and Use WebGPU. OK docks the MCL viewer and adds the cluster, the embedding, the cluster
  size and the connectivity of every row as columns. Neither menu item has a help page; the DBSCAN
  parameters are the ones the help explains under Postprocessing in `explore/dim-reduction.md`.

  DBSCAN runs on a nine-row table written here: two tight groups of four points and one point far from
  both. MCL runs on iris (150 rows), over Species with Petal.Length and over the four measurements; with
  number columns only it never finishes (the viewer stays at "Updating..."), so that scenario is
  `@known-failure` (MISSING.md, Findings 30 Oct 2026). Not a journey: each run adds columns the
  next one's claims would read. Left out on purpose:
  - how many clusters MCL finds and which row goes where: EDA and the ml library's tests;
  - DBSCAN as the postprocessing of Reduce Dimensionality (`reduce-dimensionality.feature`), of Chemical
    Space and Sequence Space (Chem, Bio), and MCL inside the Peptides SAR viewer: their own features;
  - Use WebGPU: headless Chromium has no WebGPU adapter, so the input is off and disabled (MISSING.md).
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
    And Epsilon input in "DBSCAN" dialog should have value "0.02"
    And "Minimum points" input in "DBSCAN" dialog should have value "4"
    When user selects "x" in "X Col" input in "DBSCAN" dialog
    And user selects "y" in "Y Col" input in "DBSCAN" dialog
    And user clicks on OK button in "DBSCAN" dialog
    Then "DBSCAN" dialog should be hidden
    And a new column "Cluster (DBSCAN)" should have been added
    And the newest column matching "^Cluster \(DBSCAN\)" should have 3 distinct values
    And the value of "Cluster (DBSCAN)" column in row 1 should be "1"
    And the value of "Cluster (DBSCAN)" column in row 4 should be "1"
    And the value of "Cluster (DBSCAN)" column in row 5 should be "2"
    And the value of "Cluster (DBSCAN)" column in row 8 should be "2"
    And the value of "Cluster (DBSCAN)" column in row 9 should be "-1"
    And the table should have the columns "x, y, Cluster (DBSCAN)"
    And no errors should have been logged

  Scenario: MCL on a text and a number column docks its viewer and adds the cluster columns
    Given user opens iris dataset
    When user picks "ML > Cluster > MCL..." from the top menu
    Then "Markov clustering" dialog should be visible
    And Columns input in "Markov clustering" dialog should contain text "(0)"
    And "Similarity Threshold" input in "Markov clustering" dialog should have value "80"
    And "Inflation Factor" input in "Markov clustering" dialog should have value "2"
    And "Max Iterations" input in "Markov clustering" dialog should have value "5"
    And "Min Cluster Size" input in "Markov clustering" dialog should have value "5"
    When user clicks on Columns input in "Markov clustering" dialog
    And user toggles the "Species" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Length" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then Columns input in "Markov clustering" dialog should contain text "(2)"
    When user clicks on OK button in "Markov clustering" dialog
    Then "Markov clustering" dialog should be hidden
    And 5 new columns should have been added
    And a new column "Cluster (MCL)" should have been added
    And a new column "EmbedX (MCL)" should have been added
    And a new column "EmbedY (MCL)" should have been added
    And a new column "Cluster size (MCL)" should have been added
    And a new column "Connectivity (MCL)" should have been added
    And the newest column matching "^Cluster \(MCL\)" should have no missing values
    And the newest column matching "^EmbedX \(MCL\)" should have no missing values
    And the open tableview should have 1 MCL viewer
    And the table should have 150 rows
    And no errors should have been logged

  @known-failure:97
  Scenario: MCL on four number columns adds the cluster columns
    Given user opens iris dataset
    When user picks "ML > Cluster > MCL..." from the top menu
    And user clicks on Columns input in "Markov clustering" dialog
    And user toggles the "Sepal.Length" column in the column list of "Select columns..." dialog
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Length" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then Columns input in "Markov clustering" dialog should contain text "(4)"
    When user clicks on OK button in "Markov clustering" dialog
    Then "Markov clustering" dialog should be hidden
    And the open tableview should have 1 MCL viewer
    And a new column "Cluster (MCL)" should have been added
    And no errors should have been logged
