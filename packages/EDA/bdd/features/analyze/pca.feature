@journey @eda @realizes:ml.menu.analyze.pca
Feature: Principal component analysis
  ML | Analyze | PCA... over cars: the sixteen numeric columns, three components, then the same with
  Center and Scale on. Translated from files/TestTrack/EDA/pca.md and the package's
  playwright/pca.test.ts.

  The Features picker is the platform's "Select columns..." dialog; All checks every column it
  offers, and the picker counts them itself. The second run excludes PC1 to PC3 from its inputs,
  keeping the original sixteen predictors; the platform suffixes the repeated output names.
  Centered unit-variance scores over thirty rows have magnitude below six.

  Background:
    Given user is logged in
    And user opens cars dataset
    Then table "cars" should have 30 rows

  Scenario: Three components over every column add PC1 to PC3
    When user picks "ML > Analyze > PCA..." from the top menu
    Then "PCA" dialog should be visible
    When user clicks on editor of Features input in "PCA" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on All label in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "16 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input in "PCA" dialog should contain text "(16)"
    When user enters "3" into Components input in "PCA" dialog
    And user unchecks Center input in "PCA" dialog
    And user unchecks Scale input in "PCA" dialog
    Then Center input in "PCA" dialog should be unchecked
    And Scale input in "PCA" dialog should be unchecked
    When user clicks on OK button in "PCA" dialog
    Then the top menu command should have completed
    And "PCA" dialog should be hidden
    And 3 new columns should have been added
    And the table should have a column "PC1"
    And the table should have a column "PC2"
    And the table should have a column "PC3"
    And "PC1" column should have no missing values
    And "PC2" column should have no missing values
    And "PC3" column should have no missing values
    And "PC1" column should have at least 2 distinct values
    And "PC2" column should have at least 2 distinct values
    And "PC3" column should have at least 2 distinct values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Center and Scale add a second set of components beside the first
    When user picks "ML > Analyze > PCA..." from the top menu
    Then "PCA" dialog should be visible
    When user clicks on editor of Features input in "PCA" dialog
    And user clicks on All label in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "19 checked"
    When user types "PC" into Search input in "Select columns..." dialog
    Then the "rows shown" reading of grid viewer in "Select columns..." dialog should be 3
    And the "text of cell 17 of __name" reading of grid viewer in "Select columns..." dialog should be "PC1"
    And the "text of cell 18 of __name" reading of grid viewer in "Select columns..." dialog should be "PC2"
    And the "text of cell 19 of __name" reading of grid viewer in "Select columns..." dialog should be "PC3"
    When user clicks on the "cell 17 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 18 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 19 of x" area of grid viewer in "Select columns..." dialog
    Then the "text of cell 17 of x" reading of grid viewer in "Select columns..." dialog should be "false"
    And the "text of cell 18 of x" reading of grid viewer in "Select columns..." dialog should be "false"
    And the "text of cell 19 of x" reading of grid viewer in "Select columns..." dialog should be "false"
    And "Select columns..." dialog should contain text "16 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input in "PCA" dialog should contain text "(16)"
    When user enters "3" into Components input in "PCA" dialog
    And user checks Center input in "PCA" dialog
    And user checks Scale input in "PCA" dialog
    Then Center input in "PCA" dialog should be checked
    And Scale input in "PCA" dialog should be checked
    When user clicks on OK button in "PCA" dialog
    Then the top menu command should have completed
    And "PCA" dialog should be hidden
    And 3 new columns should have been added
    And the table should have a column "PC1 (2)"
    And the table should have a column "PC2 (2)"
    And the table should have a column "PC3 (2)"
    And "PC1 (2)" column should have no missing values
    And "PC2 (2)" column should have no missing values
    And "PC3 (2)" column should have no missing values
    And "PC1 (2)" column should have at least 2 distinct values
    And "PC2 (2)" column should have at least 2 distinct values
    And "PC3 (2)" column should have at least 2 distinct values
    And every value of "PC1 (2)" column should lie between -6 and 6
    And every value of "PC2 (2)" column should lie between -6 and 6
    And every value of "PC3 (2)" column should lie between -6 and 6
    And no error or warning balloon should have been shown
    And no errors should have been logged
