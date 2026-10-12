@journey @eda @realizes:ml.menu.models.train-model @realizes:ml.menu.models.apply-model @realizes:eda.model.linear-regression
Feature: The preprocessing options of Train Model, kept by the saved model
  The options the Train Model view offers under the features — Ignore missing and Impute missing for
  missing values, Predict probability for a two-class target, One-hot encoding for text features — shown
  only where they apply, and remembered by a saved model, so Apply Model repeats them on new data.
  Translated from files/TestTrack/Models/models-testdemog-lifecycle-smoke.md blocks 1–2 and
  models-one-hot-suffix-collision.md, and their playwright-public specs. The old smoke spec only warned
  where Impute missing or Predict probability was out of place, and ticked Predict probability only
  when it happened to be shown.

  demog-1000 has missing HEIGHTs (128) and none elsewhere, so SEX predicted from HEIGHT and WEIGHT needs
  a missing-value option before any engine is offered; a change of the target clears the options
  ticked. Predict probability turns the two classes into 0 and 1 for a regression engine and turns its
  output back into the class names at the default cutoff (0.5), so the model applied to new rows writes
  A and B, not numbers; the cutoff is not moved, since its re-render sets no busy state to wait on
  (MISSING.md). One-hot encoding expands featureA and featureB into one 0/1 column per category inside the
  model; the model's inputs stay featureA and featureB, which is what Apply Model maps (2/2) on a fresh
  table. The old case expected the inputs listed as featureA=Yes, featureA=No, ... (MISSING.md, to ask).

  Not translated: the per-category column names inside the trained model (the old spec only re-read
  the name it had saved; what the model holds is the EDA package tests' claim), and the KNN dialog
  Impute missing opens (EDA's imputation feature).

  Background:
    Given user is logged in
    And no predictive model named "BDD-Probability-{run}" is on the server
    And no predictive model named "BDD-OneHot-{run}" is on the server

  Scenario: Ignore missing is offered for missing values, and hides Impute missing once ticked
    Given user opens demog-1000 dataset
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "SEX" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    Then the "text of cell 6 of __name" reading of grid viewer in "Select columns..." dialog should be "HEIGHT"
    And the "text of cell 7 of __name" reading of grid viewer in "Select columns..." dialog should be "WEIGHT"
    When user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 7 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should contain text "Column 'HEIGHT' contains missing values."
    And model preview should be invalid
    And "Ignore missing" input should be visible
    And "Impute missing" input should be visible
    And "Predict probability" input should be visible
    And "Model Engine" input should be absent
    When user checks "Ignore missing" input
    Then model preview should be ready
    And "Impute missing" input should be hidden
    And "Predict probability" input should be visible
    And "Model Engine" input should be visible
    And "Accuracy" table row should be visible
    And no errors should have been logged

  Scenario: A change to a numeric target clears the options ticked, and Predict probability goes
    When user selects "AGE" in Predict input
    Then editor of Predict input should have text "AGE"
    And "Ignore missing" input should be unchecked
    And model preview should be invalid
    And "Model Engine" input should be absent
    When user checks "Ignore missing" input
    Then model preview should be ready
    And "R squared" table row should be visible
    And "Predict probability" input should be hidden
    And no errors should have been logged

  Scenario: Predict probability trains a regression on the two classes, with a cutoff and a ROC curve
    Given user opens a table "readings" with:
      | f1 | f2 | target |
      | 1  | 7  | A      |
      | 2  | 3  | A      |
      | 3  | 9  | B      |
      | 4  | 1  | A      |
      | 5  | 8  | B      |
      | 6  | 2  | A      |
      | 7  | 10 | B      |
      | 8  | 4  | B      |
      | 9  | 6  | A      |
      | 10 | 5  | B      |
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    Then the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "f1"
    And the "text of cell 2 of __name" reading of grid viewer in "Select columns..." dialog should be "f2"
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then model preview should be ready
    And "Predict probability" input should be visible
    And "Positive class cutoff" input should be absent
    And model preview should not contain text "ROC Curve"
    When user checks "Predict probability" input
    Then model preview should be ready
    And "Positive class cutoff" input should be visible
    And model preview should contain text "ROC Curve"
    And "Eda: Linear Regression" heading should be visible
    When user clicks on Save button
    And user enters "BDD-Probability-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-Probability-{run}" should be on the server
    And no errors should have been logged

  Scenario: The probability model applied to new rows writes the class names
    Given user opens a table "new readings" with:
      | f1 | f2 |
      | 1  | 7  |
      | 3  | 9  |
      | 4  | 1  |
      | 7  | 10 |
      | 8  | 4  |
      | 9  | 6  |
    When user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    When user selects the predictive model "BDD-Probability-{run}" in Model input in "Apply predictive model" dialog
    Then Inputs input in "Apply predictive model" dialog should contain text "(2/2)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And a new column "target" should have been added
    And "target" column should have no missing values
    # one row of each class: an inverted cutoff or swapped class names would swap them
    And the "target" cell of row 1 should be displayed as "A"
    And the "target" cell of row 4 should be displayed as "B"
    And no errors should have been logged

  Scenario: One-hot encoding trains on two Yes/No features, and the model applies to a fresh table
    Given user opens a table "answers" with:
      | featureA | featureB | target |
      | Yes      | Yes      | 5.00   |
      | No       | Yes      | 3.07   |
      | Yes      | No       | 2.14   |
      | No       | No       | 0.21   |
      | Yes      | Yes      | 5.28   |
      | No       | Yes      | 3.35   |
      | Yes      | No       | 2.42   |
      | No       | No       | 0.49   |
      | Yes      | Yes      | 5.56   |
      | No       | Yes      | 3.63   |
      | Yes      | No       | 2.70   |
      | No       | No       | 0.77   |
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    Then the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "featureA"
    And the "text of cell 2 of __name" reading of grid viewer in "Select columns..." dialog should be "featureB"
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then model preview should contain text "Columns 'featureA, featureB' are categorical."
    And model preview should be invalid
    And "Model Engine" input should be absent
    When user checks "One-hot encoding" input
    Then model preview should be ready
    And "Model Engine" input should be visible
    When user selects "Eda: Linear Regression" in "Model Engine" input
    Then model preview should be ready
    And "Eda: Linear Regression" heading should be visible
    And "R squared" table row should be visible
    When user clicks on Save button
    And user enters "BDD-OneHot-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-OneHot-{run}" should be on the server
    Given user opens a table "new answers" with:
      | featureA | featureB |
      | No       | No       |
      | Yes      | No       |
      | No       | Yes      |
      | Yes      | Yes      |
    When user picks "ML > Models > Apply Model..." from the top menu
    And user selects the predictive model "BDD-OneHot-{run}" in Model input in "Apply predictive model" dialog
    Then Inputs input in "Apply predictive model" dialog should contain text "(2/2)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And a new column "target" should have been added
    And "target" column should have no missing values
    # Yes adds about 2 for featureA and about 3 for featureB: swapped inputs would swap rows 2 and 3
    And the "target" cell of row 2 should be displayed as "2.42"
    And the "target" cell of row 3 should be displayed as "3.35"
    And no error or warning balloon should have been shown
    And no errors should have been logged
