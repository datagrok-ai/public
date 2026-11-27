@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). A walk down the table view's ML menu, one line per group,
  nothing picked: Analyze (PCA, PLS, Multivariate Analysis and Group Comparison: T-test, ANOVA, Control
  Comparisons; EDA, in the browser), Reduce Dimensionality (UMAP, t-SNE), Impute Missing Values, Pareto Front,
  Cluster (Cluster Data, DBSCAN, MCL, Hierarchical Clustering) and Models (Train Model, Apply Model: the Forge
  package, shown only where it is installed). Tools > Random Data and Notebooks (Jupyter) are left out: Tools did
  not open from a rest on it in time (the walk then flew back to the bar), and pointing at Notebooks moved the pointer
  off the menu, which closed. A group is entered by resting on it, so the pointer stays inside the dropdown. demog is
  opened off camera; nothing is made on the server.

  Scenario: What the ML menu does
    The ML menu holds the statistics and machine learning tools that run on the table you are looking at.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "ML: statistics and machine learning on the open table" at the top
    # caption: Analyze: PCA, PLS and multivariate analysis
    And user points at "ML > Analyze > PCA..." in the top menu
    And user points at "ML > Analyze > Multivariate Analysis..." in the top menu
    # caption: Group Comparison: t-test, ANOVA, each group against a control
    And user points at "ML > Analyze > Group Comparison" in the top menu
    And user points at "ML > Analyze > Group Comparison > ANOVA..." in the top menu
    # caption: Reduce Dimensionality: a 2D map of the rows (UMAP, t-SNE)
    And user points at "ML > Reduce Dimensionality..." in the top menu
    # caption: Impute Missing Values: fill the gaps from similar rows
    And user points at "ML > Impute Missing Values..." in the top menu
    # caption: Pareto Front: the best trade-offs between several goals
    And user points at "ML > Pareto Front..." in the top menu
    # caption: Cluster: k-means, DBSCAN, MCL and a hierarchical tree
    And user points at "ML > Cluster" in the top menu
    And user points at "ML > Cluster > Cluster Data..." in the top menu
    # caption: Models: train a model on this table, apply it to another
    And user points at "ML > Models" in the top menu
    And user points at "ML > Models > Train Model..." in the top menu
    Then "Train Model..." menu item should be visible
    And no errors should have been logged
