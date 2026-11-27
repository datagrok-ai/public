@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Reduce Dimensionality... runs UMAP (or t-SNE) in a web
  worker of the EDA package, no script. On winequality the Columns picker takes All and leaves out quality and type,
  so the map is built from the eleven lab measurements only; OK adds Embed_X_1 / Embed_Y_1 and, from the default
  DBSCAN postprocessing, a Cluster (DBSCAN) column, and docks a scatter plot of the embedding colored by the
  clusters. The plot's color is then set to type: red and white wines land in different parts of the map (UMAP is
  random: the shape differs from run to run, the split by type does not). The type column is colored red / gold off
  camera (text only in the grid) so the legend reads naturally. winequality is opened off camera; nothing is made on
  the server.

  Scenario: Map similar rows with UMAP
    UMAP and t-SNE place similar rows close together on a 2D map, judging them by many columns at once.

    Given user is logged in
    And simple mode is off
    And user opens winequality dataset
    And user colors "type" column categorically:
      | red   | #b2182b |
      | white | #d9a520 |
    And the grid paints the colors of "type" column on its text
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "A map that puts wines with similar measurements close together" at the bottom
    # caption: Reduce Dimensionality is in the ML menu
    And user picks "ML > Reduce Dimensionality..." from the top menu
    Then "Reduce Dimensionality" dialog should be visible
    # caption: Columns: what to compare the rows by
    When user clicks on Columns input in "Reduce Dimensionality" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    # caption: Leave out quality and type: the map should not see them
    And user types "quality" into Search input in "Select columns..." dialog
    And user toggles the "quality" column in the column list of "Select columns..." dialog
    And user types "type" into Search input in "Select columns..." dialog
    And user toggles the "type" column in the column list of "Select columns..." dialog
    Then "Select columns..." dialog should contain text "11 checked"
    When user clicks on OK button in "Select columns..." dialog
    # caption: Method: UMAP, or t-SNE
    And user hovers over Method input in "Reduce Dimensionality" dialog
    Then Method input in "Reduce Dimensionality" dialog should have value "UMAP"
    Given user watches the task bar
    # caption: OK computes the map in the browser and plots it
    When user clicks on OK button in "Reduce Dimensionality" dialog
    Then the task bar should have finished "dimensionality reduction"
    And a new column "Cluster (DBSCAN)" should have been added
    And "xColumnName" property of scatter plot viewer should be "Embed_X_1"
    And "colorColumnName" property of scatter plot viewer should be "Cluster (DBSCAN)"
    When the video says "One dot per wine: close dots have similar measurements"
    # caption: Color by type: red and white wines land in different parts of the map
    And user picks "type" in the "color" column selector of scatter plot viewer
    Then "Color" property of scatter plot viewer should be "type"
    When user moves the pointer away from scatter plot viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
