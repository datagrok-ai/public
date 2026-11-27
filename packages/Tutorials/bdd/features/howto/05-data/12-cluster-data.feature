@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Cluster > Cluster Data runs k-means in the browser
  (dml.kMeans in the Dart client, no script). On iris the dialog takes every numeric column without missing values
  as a feature, the row number "col 1" included, so it is unticked in Select columns; 2 clusters, Show Scatter Plot
  colors a scatter plot by the clusters while the dialog is open, OK adds a "clusters" column, and the plot's X set
  to Petal.Width shows setosa apart from the other two. Two clusters, not three: the dialog starts k-means from a
  round-robin assignment (no k-means++), and on iris, sorted by species, 3 clusters split setosa in two and merge
  versicolor with virginica for every Normalize choice (in the report of the 05-data clips). iris is opened off
  camera; nothing is made on the server.

  Scenario: Cluster data
    Group similar rows with k-means: pick the columns to compare and the number of clusters.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Cluster Data: k-means, right in the browser
    When user picks "ML > Cluster > Cluster Data..." from the top menu
    Then "Cluster Data" dialog should be visible
    # caption: Features: the columns that make rows alike
    When user clicks on Features input in "Cluster Data" dialog
    And user toggles the "col 1" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    # caption: Clusters: how many groups to find
    And user enters "2" into Clusters input in "Cluster Data" dialog
    # caption: Show Scatter Plot: see the clusters as you change the settings
    And user checks "Show Scatter Plot" input in "Cluster Data" dialog
    Then the open tableview should have 1 scatter plot viewer
    # caption: OK keeps them as a column
    When user clicks on OK button in "Cluster Data" dialog
    Then "Cluster Data" dialog should be hidden
    And the table should have a column "clusters"
    And "clusters" column should have at least 2 distinct values
    # caption: By petal size, setosa forms a group of its own
    When user picks "Petal.Width" in the "x" column selector of scatter plot viewer
    Then "x Column Name" property of scatter plot viewer should be "Petal.Width"
    And no errors should have been logged
