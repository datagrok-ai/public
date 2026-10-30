@eda
Feature: ML > Impute Missing Values fills the missing values of a column from its nearest neighbours
  ML | Impute Missing Values... opens the k-NN Imputation dialog on the current table (help
  `explore/missing-values-imputation.md`): Impute lists the columns that have missing values, Using the
  columns that find the neighbours, with Distance, Neighbors (4) and In-place (on). RUN fills each missing
  value with the average of its neighbours' values, in place or, with In-place off, in a new column. A table
  with no missing value gets an info balloon instead of the dialog.

  Every scenario runs on a five-row table written here: a = 1, 2, 5, 9, 10 and b = 10, 20, (missing), 90,
  100, so the four neighbours of row 3 average to 55 and its nearest one (a = 2) holds 20. Not a journey:
  each RUN changes the table the next scenario would start from. Left out on purpose:
  - the distance and weight matrix (Euclidean / Manhattan, per-feature metrics and weights, one-hot text
    features, Keep empty): EDA src/tests/mis-vals-imputation-tests.ts;
  - Impute missing as an option of Train Model: `models/train-options.feature`.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "gaps" with:
      | a  | b   |
      | 1  | 10  |
      | 2  | 20  |
      | 5  |     |
      | 9  | 90  |
      | 10 | 100 |

  Scenario: The dialog opens on the column with missing values, and RUN fills it in place
    When user picks "ML > Impute Missing Values..." from the top menu
    Then "k-NN Imputation" dialog should be visible
    And Impute input in "k-NN Imputation" dialog should contain text "(1)"
    And Using input in "k-NN Imputation" dialog should contain text "(2)"
    And Neighbors input in "k-NN Imputation" dialog should have value "4"
    And "In-place" input in "k-NN Imputation" dialog should be checked
    And RUN button in "k-NN Imputation" dialog should be enabled
    When user clicks on RUN button in "k-NN Imputation" dialog
    Then "k-NN Imputation" dialog should be hidden
    And "b" column should have no missing values
    And the value of "b" column in row 3 should be "55"
    And the value of "b" column in row 2 should be "20"
    And the value of "a" column in row 3 should be "5"
    And the table should have the columns "a, b"
    And no errors should have been logged

  Scenario: With In-place off and one neighbour, RUN adds the filled copy and leaves the column as it was
    When user picks "ML > Impute Missing Values..." from the top menu
    And user enters "1" into Neighbors input in "k-NN Imputation" dialog
    And user unchecks "In-place" input in "k-NN Imputation" dialog
    And user clicks on RUN button in "k-NN Imputation" dialog
    Then a new column "b(copy)" should have been added
    And the value of "b(copy)" column in row 3 should be "20"
    And the value of "b(copy)" column in row 4 should be "90"
    And "b" column should have missing values
    And the table should have the columns "a, b, b(copy)"
    And no errors should have been logged

  Scenario: Nothing to impute disables RUN
    When user picks "ML > Impute Missing Values..." from the top menu
    And user clicks on Impute input in "k-NN Imputation" dialog
    And user toggles the "b" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then RUN button in "k-NN Imputation" dialog should be disabled
    When user clicks on CANCEL button in "k-NN Imputation" dialog
    Then "k-NN Imputation" dialog should be hidden
    And "b" column should have missing values
    And no errors should have been logged

  Scenario: A table with no missing values gets an info balloon and no dialog
    Given user opens a table "full" with:
      | a | b  |
      | 1 | 10 |
      | 2 | 20 |
    When user picks "ML > Impute Missing Values..." from the top menu
    Then an info balloon containing "No missing values" should have been shown
    And "k-NN Imputation" dialog should be absent
    And no errors should have been logged
