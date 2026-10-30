@eda
Feature: ML > Reduce Dimensionality embeds the chosen columns in two dimensions and plots them
  ML | Reduce Dimensionality... opens the Reduce Dimensionality dialog on the current table (help
  `explore/dim-reduction.md`): no column picked and OK disabled, Method UMAP (or t-SNE), Plot Embeddings on,
  Postprocessing DBSCAN clustering. With columns picked, OK adds an Embed_X / Embed_Y pair to the table,
  docks a scatter plot of them titled after the method and the columns, and the DBSCAN postprocessing adds
  a Cluster (DBSCAN) column that colours the plot. Without Plot Embeddings only the columns are added.

  Every scenario runs on cars (30 rows) with three number columns. Not a journey: every run adds columns
  and a scatter plot that the next run's claims would read. Left out on purpose:
  - the embeddings themselves, the method and DBSCAN settings behind the gear icons, the per-column encoding
    and the aggregation of mixed column types: EDA src/tests/dim-reduction-tests.ts and the ml library;
  - Radial Coloring, which the help lists under Postprocessing: it comes from the NLP package, not EDA;
  - Use WebGPU: headless Chromium has no WebGPU adapter (MISSING.md);
  - Chemical Space and Sequence Space, the same editor for molecules and sequences, have their features in
    Chem and Bio.
  The OK of a dialog with no column is greyed out by a CSS class only, so it is not reported as disabled:
  that scenario is `@known-failure` (MISSING.md, Findings 30 Oct 2026).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: The dialog opens with no column picked, and OK on three columns adds the embedding, its plot and clusters
    When user picks "ML > Reduce Dimensionality..." from the top menu
    Then "Reduce Dimensionality" dialog should be visible
    And Table input in "Reduce Dimensionality" dialog should have value "cars"
    And Columns input in "Reduce Dimensionality" dialog should contain text "(0)"
    And Method input in "Reduce Dimensionality" dialog should have value "UMAP"
    And Method input in "Reduce Dimensionality" dialog should offer "UMAP, t-SNE"
    And "Plot Embeddings" input in "Reduce Dimensionality" dialog should be checked
    And Postprocessing input in "Reduce Dimensionality" dialog should have value "DBSCAN clustering"
    When user clicks on Columns input in "Reduce Dimensionality" dialog
    And user toggles the "wheel.base" column in the column list of "Select columns..." dialog
    And user toggles the "length" column in the column list of "Select columns..." dialog
    And user toggles the "width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then Columns input in "Reduce Dimensionality" dialog should contain text "(3)"
    And OK button in "Reduce Dimensionality" dialog should be enabled
    Given user watches the task bar
    When user clicks on OK button in "Reduce Dimensionality" dialog
    Then "Reduce Dimensionality" dialog should be hidden
    And the task bar should have finished "dimensionality reduction"
    And 3 new columns should have been added
    And a new column "Embed_X_1" should have been added
    And a new column "Embed_Y_1" should have been added
    And a new column "Cluster (DBSCAN)" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And the newest column matching "^Embed_Y_" should have no missing values
    And the newest column matching "^Cluster \(DBSCAN\)" should have no missing values
    And scatter plot viewer should be visible
    And "title" property of scatter plot viewer should be "UMAP (wheel.base, length, width)"
    And "xColumnName" property of scatter plot viewer should be "Embed_X_1"
    And "yColumnName" property of scatter plot viewer should be "Embed_Y_1"
    And "colorColumnName" property of scatter plot viewer should be "Cluster (DBSCAN)"
    And scatter plot viewer should be painted
    And the table should have 30 rows
    And no errors should have been logged

  @known-failure:65
  Scenario: With no column picked OK is disabled
    When user picks "ML > Reduce Dimensionality..." from the top menu
    Then Columns input in "Reduce Dimensionality" dialog should contain text "(0)"
    And OK button in "Reduce Dimensionality" dialog should be disabled
    And no errors should have been logged

  Scenario: t-SNE with no postprocessing adds only the embedding and plots it under the method's name
    When user picks "ML > Reduce Dimensionality..." from the top menu
    And user clicks on Columns input in "Reduce Dimensionality" dialog
    And user toggles the "wheel.base" column in the column list of "Select columns..." dialog
    And user toggles the "length" column in the column list of "Select columns..." dialog
    And user toggles the "width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user selects "t-SNE" in Method input in "Reduce Dimensionality" dialog
    And user selects "None" in Postprocessing input in "Reduce Dimensionality" dialog
    Given user watches the task bar
    When user clicks on OK button in "Reduce Dimensionality" dialog
    Then the task bar should have finished "dimensionality reduction"
    And the newest column matching "^Embed_X_" should have no missing values
    And 2 new columns should have been added
    And scatter plot viewer should be visible
    And "title" property of scatter plot viewer should be "t-SNE (wheel.base, length, width)"
    And the table should not have a column "Cluster (DBSCAN)"
    And no errors should have been logged

  Scenario: Without Plot Embeddings the columns are added and no scatter plot is docked
    When user picks "ML > Reduce Dimensionality..." from the top menu
    And user clicks on Columns input in "Reduce Dimensionality" dialog
    And user toggles the "wheel.base" column in the column list of "Select columns..." dialog
    And user toggles the "length" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user unchecks "Plot Embeddings" input in "Reduce Dimensionality" dialog
    Given user watches the task bar
    When user clicks on OK button in "Reduce Dimensionality" dialog
    Then the task bar should have finished "dimensionality reduction"
    And the newest column matching "^Embed_Y_" should have no missing values
    And a new column "Embed_X_1" should have been added
    And the open tableview should have 0 scatter plot viewers
    And no errors should have been logged
