@journey @realizes:chem.cp.activity-cliffs
Feature: Activity Cliffs over molecules and their activity
  Chem | Analyze | Activity Cliffs... opens on the molecule column of activity-cliffs (29 molecules in
  smiles with an Activity column) with UMAP, Tanimoto and a similarity cutoff of 80. OK embeds the
  molecules, plots them, and reports the cliffs it found on the plot; the cliff link opens the
  Activity cliffs panel. Show only cliffs keeps the rows that take part in one
  and switching it off lets every row through again. The similarity cutoff decides how many pairs count as
  cliffs: 52 of them at 20, two at the default 80 and one at 95. A fresh table carries each run,
  since a second run plots itself beside the first.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens activity-cliffs dataset

  Scenario: The dialog opens on the molecule column with its defaults
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    Then "Activity Cliffs" dialog should be visible
    And Column input in "Activity Cliffs" dialog should contain text "smiles"
    And Activities input in "Activity Cliffs" dialog should contain text "Activity"
    And Method input in "Activity Cliffs" dialog should have value "UMAP"
    And Similarity input in "Activity Cliffs" dialog should have value "Tanimoto"
    And no errors should have been logged

  Scenario: The run embeds the molecules, plots them and counts the cliffs
    When user clicks on OK button in "Activity Cliffs" dialog
    Then a new column matching "^Embed_X_" should have been added
    And a new column matching "^Embed_Y_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And scatter plot viewer should be visible
    And the "cliffs" reading of scatter plot viewer should be at least 1
    And the "only cliffs" reading of scatter plot viewer should be "false"
    And the table should have 29 rows
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: The cliff link opens the panel that lists the cliffs
    When user clicks on "cliffs" button
    Then "Activity cliffs" dock panel should be visible
    And no errors should have been logged

  Scenario: Show only cliffs keeps the rows that take part in one
    When user switches on "Show only cliffs" input
    Then the "only cliffs" reading of scatter plot viewer should be "true"
    And fewer than 29 rows should pass the filter
    When user switches off "Show only cliffs" input
    Then the "only cliffs" reading of scatter plot viewer should be "false"
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A lower similarity cutoff finds many more cliffs
    Given user opens activity-cliffs dataset
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    And user enters "20" into "Similarity cutoff" input in "Activity Cliffs" dialog
    Then "Similarity cutoff" input in "Activity Cliffs" dialog should have value "20"
    When user clicks on OK button in "Activity Cliffs" dialog
    Then the top menu command should have completed
    And the "cliffs" reading of scatter plot viewer should be at least 10
    And no errors should have been logged

  Scenario: A stricter similarity cutoff finds fewer
    Given user opens activity-cliffs dataset
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    And user enters "95" into "Similarity cutoff" input in "Activity Cliffs" dialog
    And user clicks on OK button in "Activity Cliffs" dialog
    Then the top menu command should have completed
    And the "cliffs" reading of scatter plot viewer should be 1
    And no errors should have been logged
