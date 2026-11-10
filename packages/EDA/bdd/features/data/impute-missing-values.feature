@eda
@journey
Feature: ML > Impute Missing Values fills the missing values of a column from its nearest neighbours
  ML | Impute Missing Values... opens the k-NN Imputation dialog on the current table (help
  `explore/missing-values-imputation.md`): Impute lists the columns that have missing values, and RUN fills each
  missing value from its neighbours, in place by default. A table with no missing value gets an info balloon instead
  of the dialog.

  Runs on a five-row table written here: b misses its third value. The imputed values, the distances and In-place
  off are EDA's 'Missing values imputation' tests (MISSING.md). Nothing is put on the server.

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
    When user clicks on RUN button in "k-NN Imputation" dialog
    Then "b" column should have no missing values
    And the table should have the columns "a, b"
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
