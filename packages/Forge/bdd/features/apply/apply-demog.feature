@journey @serial
Feature: Applying on data with natural gaps: demog's HEIGHT and WEIGHT, skipped by default or imputed
  ML | Forge | Apply... over demog.csv opened from the demo files, whose HEIGHT (751 empty cells) and WEIGHT (one) are
  the features of a regressor of AGE the Background trains and saves through Forge:trainModel (Linear Regression, the
  method the Train view suggests for this selection; no data stored). Translated from the manual scenario 3.15 (steps
  2-4) of core/docs/features/platform/predictive-modeling/forge/manual-checks.md. The gaps are the demo file's own:
  the feature reads their counts from its table first (751, 1, 752 rows with a gap in either, the empty AGE of row
  1193), so a changed file fails there, with the numbers; the predictions are claimed row by row against those gaps,
  and the counts of the user's columns are read again after both applications. Every application is counted on the
  server for this feature's model, so a count is exact. An application's end is its task bar entry gone (followed as
  an element: the entry shows the progress instead of "Predicting AGE" while it runs), the new column and the balloon,
  which comes after the record. Impute is claimed by what it predicts: the feature imputes copies of HEIGHT and WEIGHT
  with Eda:knnImpute and the dialog's settings, predicts that copy with the same model through Forge:applyModel and
  claims the column equals it in every row (an Impute that filled nothing would still predict every row, from the
  empty cells read as near zero); that check records one more application, after the counts are claimed. The user's
  columns keep their gaps row by row: the first prediction column is still empty exactly in the gap rows afterwards.

  Not claimed here: 3.15 step 1 (the Train view on demog: M1c; the Background's run record is claimed instead: 5097
  rows used, 753 skipped) and step 5 (the delete: the sweep); the grid and the context panel's statistics of steps 3-4
  (the cells and the missing-value counts are read from the table's columns); the layout: "Skip rows and Impute on one
  line", the Missing values input below Columns, Neighbors and Distance under it; the tooltip of a row "ends with" its
  missing-values sentence (claimed as contained); the prediction's values beyond the imputed copy's (Linear
  Regression on an int target gives decimals, nothing rounds them); the manual's command line over the newest
  application of the stand (the newest application of this feature's model is claimed instead).

  Background:
    Given user is logged in
    And no Forge models named "forge-test-apply-demog-{time}" are on the server
    And no Forge training runs of the table "forge-test-demog-{time}" are on the server
    And user opens the demo file "demog.csv" as "forge-test-demog-{time}"
    And a Forge model "forge-test-apply-demog-{time}" is trained on the current table through the API with:
      | target   | AGE               |
      | features | HEIGHT, WEIGHT    |
      | method   | Linear Regression |

  Scenario: Rows with a natural gap in HEIGHT or WEIGHT get no prediction by default, and are counted (3.15, steps 2-3)
    Given user counts the task bar entries
    Then the newest training run of the table "forge-test-demog-{time}" should have:
      | source                            | api  |
      | row_count                         | 5097 |
      | options.missingValues.mode        | skip |
      | options.missingValues.skippedRows | 753  |
    And "HEIGHT" column should have 751 missing values
    And "WEIGHT" column should have 1 missing value
    And the rows with a missing value in any of "HEIGHT, WEIGHT" columns should be 752
    And the value of "AGE" column in row 1193 should be ""
    And the value of "HEIGHT" column in row 965 should be ""
    And the value of "WEIGHT" column in row 2222 should be ""
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-demog-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-demog-{time}"
    And summary of "Columns" group should have the text "2 of 2 matched"
    And summary of "Columns" group should be valid
    And "Columns" group should be collapsed
    And "Missing values" input in "Apply predictive model" dialog should be visible
    And "Missing values" input in "Apply predictive model" dialog should have "Skip rows" chosen
    And "Neighbors" input in "Apply predictive model" dialog should be hidden
    And "Distance" input in "Apply predictive model" dialog should be hidden
    And OK button in "Apply predictive model" dialog should be enabled
    When user hovers over label of "Missing values" input in "Apply predictive model" dialog
    Then tooltip should contain text "Rows with missing values in: HEIGHT: 751, WEIGHT: 1."
    When user expands "Columns" group
    Then editor of "HEIGHT" input in "Apply predictive model" dialog should have text "HEIGHT"
    And editor of "WEIGHT" input in "Apply predictive model" dialog should have text "WEIGHT"
    When user hovers over label of "HEIGHT" input in "Apply predictive model" dialog
    Then tooltip should contain text "'HEIGHT' has 751 missing values."
    When user hovers over label of "WEIGHT" input in "Apply predictive model" dialog
    Then tooltip should contain text "'WEIGHT' has 1 missing value."
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting AGE" should have ended
    And the task bar should have shown 1 entry "Predicting AGE"
    And a new column "AGE (predicted)" should have been added
    And an info balloon containing 'Added the column "AGE (predicted)" to forge-test-demog-{time}; 752 rows skipped (missing values).' should have been shown
    And the table should have the columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted)"
    And the "AGE (predicted)" column should be empty exactly in the rows where any of "HEIGHT, WEIGHT" columns is empty
    And the value of "AGE (predicted)" column in row 965 should be ""
    And the value of "AGE (predicted)" column in row 2222 should be ""
    And the value of "AGE (predicted)" column in row 1193 should match "^-?\d"
    And "AGE (predicted)" column should have 752 missing values
    And "HEIGHT" column should have 751 missing values
    And "WEIGHT" column should have 1 missing value
    And the model "forge-test-apply-demog-{time}" should have 1 application on the server
    And the newest application of the model "forge-test-apply-demog-{time}" should have:
      | table_name   | forge-test-demog-{time} |
      | row_count    | 5850                    |
      | skipped_rows | 752                     |
      | column_name  | AGE (predicted)         |
      | source       | ui                      |
      | status       | completed               |
      | duration_ms  | number                  |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Impute fills the gaps on a copy and predicts every row; the table keeps its gaps (3.15, step 4)
    Given user counts the task bar entries
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-demog-{time}"
    And "Missing values" input in "Apply predictive model" dialog should be visible
    And "Missing values" input in "Apply predictive model" dialog should have "Skip rows" chosen
    When user selects "Impute" in "Missing values" input in "Apply predictive model" dialog
    Then "Missing values" input in "Apply predictive model" dialog should have "Impute" chosen
    And "Neighbors" input in "Apply predictive model" dialog should be visible
    And "Neighbors" input in "Apply predictive model" dialog should have value "4"
    And "Distance" input in "Apply predictive model" dialog should be visible
    And "Distance" input in "Apply predictive model" dialog should have value "Euclidean"
    And OK button in "Apply predictive model" dialog should be enabled
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting AGE" should have ended
    And the task bar should have shown 1 entry "Predicting AGE"
    And a new column "AGE (predicted 2)" should have been added
    And an info balloon containing 'Added the column "AGE (predicted 2)" to forge-test-demog-{time}.' should have been shown
    And the table should have the columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted), AGE (predicted 2)"
    And "AGE (predicted 2)" column should have no missing values
    And the value of "AGE (predicted 2)" column in row 965 should match "^-?\d"
    And the value of "AGE (predicted 2)" column in row 2222 should match "^-?\d"
    And "AGE (predicted)" column should have 752 missing values
    And "HEIGHT" column should have 751 missing values
    And "WEIGHT" column should have 1 missing value
    And the value of "HEIGHT" column in row 965 should be ""
    And the value of "WEIGHT" column in row 2222 should be ""
    And the "AGE (predicted)" column should be empty exactly in the rows where any of "HEIGHT, WEIGHT" columns is empty
    And the model "forge-test-apply-demog-{time}" should have 2 applications on the server
    And the newest application of the model "forge-test-apply-demog-{time}" should have:
      | table_name   | forge-test-demog-{time} |
      | row_count    | 5850                    |
      | skipped_rows | 0                       |
      | column_name  | AGE (predicted 2)       |
      | source       | ui                      |
      | status       | completed               |
      | duration_ms  | number                  |
    And the "AGE (predicted 2)" column should hold the predictions of the Forge model "forge-test-apply-demog-{time}" on "HEIGHT, WEIGHT" imputed with 4 neighbors and "Euclidean" distance
    And no error or warning balloon should have been shown
    And no errors should have been logged
